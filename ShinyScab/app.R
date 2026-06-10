library(shiny)
library(shinyjs)
library(SimInf)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shinyWidgets)

source("../Code/run_scenario.R")
source("../Code/plot_scenario.R")

setup_nodes <- function(Nfarms) {
  x=runif(Nfarms)*2500
  y=runif(Nfarms)*2500
  
  nodes=cbind(x,y)
  nodes
}


parse_numeric_vector <- function(x) {
  out <- strsplit(x, ",")[[1]]
  out <- trimws(out)
  as.numeric(out)
}


ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  titlePanel("Sheep Scab Transmission Model"),
  
  radioButtons(
    "farm_mode",
    "Farm setup mode:",
    choices = c("Manual (click)", "Random (generate)"),
    selected = "Manual (click)"
  ),
  
  plotOutput("farmEditor", height = "400px", click = "plot_click"),
  actionButton("clearFarms", "Clear farms"),
  
  actionButton("generateFarms", "Generate random farms"),
  radioButtons(
    "click_mode",
    "Click action:",
    choices = c("Add farm", "Toggle vaccination", "Toggle infection"),
    selected = "Add farm"
  ),
  
  sidebarLayout(
    sidebarPanel(
      numericInput("Nfarms", "Number of farms", value=5, min=2, max=10),
      numericInput("Nvacc", "Number of farms vaccinated", value=2, min=0, max=10),
      textInput("Vacc_times", "Vaccination times (comma-separated)", value = "10, 25"),
      numericInput("efficacy", "Vaccine efficacy", value=0.8, min=0, max=1),
      numericInput("coupling", "Spatial coupling between farms", value=0.5, min=0, max=1),
      numericInput("cutoff", "Max. distance at which scab can travel between farms (without trade) in km", value=1, min=0, max=50),
      radioButtons(
        "intro_type",
        "Introduce infestation on:",
        choices = c("Unvaccinated farm", "Vaccinated farm"),
        selected = "Unvaccinated farm"
      ),
      actionButton("toggleAdvanced", "Show/Hide Advanced Parameters"),
      div(id = "advancedParams",
          numericInput("alpha", "Daily shedding rate of infectious animal", value=0.005619795, min=0, max=1),
          numericInput("beta", "Decay rate of environmental infectious pressure", value=0.068, min=0, max=1),
          numericInput("upsilon", "Transmission rate from environment to susceptible sheep", value=0.002519243 * 100, min=0, max=1),
          numericInput("qprop", "Proportion of acute infections that become carriers", value=0.5, min=0, max=1),
          numericInput("tau", "Recovery rate for carriers (1/days)", value=1/653, min=0, max=1),
          numericInput("gamma", "Recovery rate for infectious sheep (1/days)", value=1/77, min=0, max=1),
          numericInput("epar", "Scaling rate for the contribution of carriers to infectious pressure", value=1/3, min=0, max=1),
          numericInput("wane", "Rate at which vaccine-based immunity wanes", value=1/653, min=0, max=1),
          numericInput("nsim", "Number of simulations", 100),
          plotOutput("gammaPlot")
      ),
      actionButton("run", "Run model")
    ),
    
    mainPanel(
      plotOutput("farmNetworkPlot", height = "300px"),
      plotOutput("plot")
    )
  )
)

server <- function(input, output, session) {
  
  # Hide advanced parameters initially
  shinyjs::hide("advancedParams")
  
  observeEvent(input$generateFarms, {
    
    if (input$farm_mode != "Random (generate)") return()
    
    Nfarms <- input$Nfarms
    Nvacc  <- min(input$Nvacc, Nfarms)
    
    # generate coordinates
    nodes <- setup_nodes(Nfarms)
    
    # assign vaccination
    vacc_idx <- sample(seq_len(Nfarms), Nvacc)
    
    # build df
    df <- data.frame(
      x = nodes[,1],
      y = nodes[,2],
      status = "Unvaccinated",
      infected = FALSE
    )
    
    df$status[vacc_idx] <- "Vaccinated"
    
    infect_idx <- sample(seq_len(Nfarms), 1)
    df$infected[infect_idx] <- TRUE
    
    df$farm_id <- seq_len(nrow(df))
    
    farms$df <- df
  })
  
  # Toggle visibility of advanced parameters
  observeEvent(input$toggleAdvanced, {
    shinyjs::toggle("advancedParams")
  })
  
  
  farms <- reactiveValues(
    df = data.frame(
      x = numeric(),
      y = numeric(),
      status = character(),
      infected = logical()
    )
  )
  
  
  observeEvent(input$plot_click, {
    
    if (input$farm_mode != "Manual (click)") return()
    
    click <- input$plot_click
    df <- farms$df
    
    if (input$click_mode == "Add farm") {
      
      df <- rbind(df, data.frame(
        x = click$x,
        y = click$y,
        status = "Unvaccinated",
        infected = FALSE
      ))
      
    } else {
      
      # find nearest farm
      if (nrow(df) > 0) {
        d <- sqrt((df$x - click$x)^2 + (df$y - click$y)^2)
        idx <- which.min(d)
        
        if (input$click_mode == "Toggle vaccination") {
          df$status[idx] <- ifelse(
            df$status[idx] == "Vaccinated",
            "Unvaccinated",
            "Vaccinated"
          )
        }
        
        if (input$click_mode == "Toggle infection") {
          df$infected[idx] <- !df$infected[idx]
        }
      }
    }
    farms$df <- df
  })
  
  
  
  observe({
    
    # only update if farms exist
    if (nrow(farms$df) > 0 && input$farm_mode == "Manual (click)") {
      
      updateNumericInput(session, "Nfarms",
                         value = nrow(farms$df))
      
      updateNumericInput(session, "Nvacc",
                         value = sum(farms$df$status == "Vaccinated"))
    }
  })
  
  
  
  observeEvent(input$clearFarms, {
    farms$df <- farms$df[0, ]
  })
  
  
  sim <- eventReactive(input$run, {
    
    # base params (copied from your notebook)
    params <- list(
      phi        = 0,
      upsilon    = input$upsilon,
      gamma      = input$gamma,
      alpha      = input$alpha,
      qprop      = input$qprop,
      tau        = input$tau,
      epar       = input$epar,
      wane       = input$wane,
      beta_t1    = input$beta,
      beta_t2    = input$beta,
      beta_t3    = input$beta,
      beta_t4    = input$beta,
      end_t1     = 91,
      end_t2     = 182,
      end_t3     = 273,
      end_t4     = 365,
      coupling   = input$coupling,
      Nfarms     = input$Nfarms,
      Nvacc      = input$Nvacc,
      cutoff     = input$cutoff*1000,
      source     = input$source,
      efficacy   = input$efficacy
    )
    
    s <- farms$df
    
    if (nrow(s) < 2) {
      showNotification("Please add at least 2 farms", type = "error")
      return(NULL)
    }
    
    vacc_times <- parse_numeric_vector(input$Vacc_times)
    
    if (length(vacc_times) == 0 || any(is.na(vacc_times))) {
      showNotification("Invalid vaccination times", type = "error")
      vaccination <- NULL
    } else {
      vaccination <- NULL  # initialise; will overwrite below
    }
    
  
    nodes <- as.matrix(s[, c("x", "y")])
    vaccinated_farms <- which(s$status == "Vaccinated")
  
    
    if (length(vaccinated_farms) == 0) {
      vaccination <- NULL
    } else {
      
      vaccination <- expand.grid(
        time = vacc_times,
        node = vaccinated_farms
      )
      
      vaccination$event <- "intTrans"
      vaccination$dest <- 0
      vaccination$n <- 0
      vaccination$proportion <- 1
      vaccination$select <- 1
      vaccination$shift <- 3
    }
    
    u0 <- data.frame(
      S = rep(100, nrow(s)),
      I = ifelse(s$infected, 1, 0),
      C = 0,
      V = 0,
      D = 0,
      H = 0
    )
    
    u0$S[u0$I == 1] <- 99
    
    vaccination <- data.frame(
      event = "intTrans", time = vacc_times, node = vaccinated_farms,
      dest = 0, n = 0, proportion = input$efficacy, select = 3, shift = 1
    )
    
    result <- run_scenario(
      params = params,
      nsim = input$nsim,
      tspan = seq(0, 730, 1),
      seed = 897342,
      cutoff = input$cutoff*1000,
      u0 = u0,
      vaccination = vaccination,
      nodes = nodes
    )
    
    result$status <- s$status[result$node]
    
    list(
      result = result,
      setup = s
    )
  })
  
  output$farmEditor <- renderPlot({
    
    df <- farms$df
    
    if (nrow(df) == 0) {
      return(
        ggplot() +
          xlim(0, 2500) +
          ylim(0, 2500) +
          coord_equal() +
          theme_minimal() +
          labs(title = "Click to place farms")
      )
    }
    
    df$farm_id <- seq_len(nrow(df))
    
    edges <- NULL
    
    if (nrow(df) > 1) {
      
      coords <- as.matrix(df[, c("x", "y")])
      dist_mat <- as.matrix(dist(coords))
      
      edges <- as.data.frame(as.table(dist_mat))
      names(edges) <- c("from", "to", "dist")
      
      edges$from <- as.integer(edges$from)
      edges$to   <- as.integer(edges$to)
      
      edges <- edges[edges$from < edges$to, ]
      
      cutoff_m <- input$cutoff * 1000
      edges <- edges[edges$dist <= cutoff_m, ]
      
      edges$strength <- input$coupling / edges$dist
      
      if (nrow(edges) > 0) {
        edges$linewidth <- scales::rescale(edges$strength, to = c(0.3, 2))
      }
      
      edges$x    <- df$x[edges$from]
      edges$y    <- df$y[edges$from]
      edges$xend <- df$x[edges$to]
      edges$yend <- df$y[edges$to]
    }
    
    ggplot() +
      
      {
        if (!is.null(edges) && nrow(edges) > 0)
          geom_segment(
            data = edges,
            aes(x = x, y = y, xend = xend, yend = yend,
                linewidth = linewidth),
            color = "grey40",
            alpha = 0.7
          )
      } +
      
      geom_point(
        data = df,
        aes(x = x, y = y, color = status, shape = infected),
        size = 4
      ) +
      
      geom_text(
        data = df,
        aes(x = x, y = y, label = farm_id),
        vjust = -1,
        size = 3
      ) +
      
      scale_color_manual(values = c(
        "Vaccinated" = "#1b9e77",
        "Unvaccinated" = "#d95f02"
      )) +
      
      scale_shape_manual(values = c(16, 17)) +
      scale_linewidth_identity() +
      
      coord_equal(xlim = c(0, 2500), ylim = c(0, 2500)) +
      theme_minimal() +
      labs(
        title = "Click to place farms",
        color = "Status",
        shape = "Infected"
      )
  })
  
  
  output$plot <- renderPlot({
    req(sim())
    plot_scenario(sim()$result)
  })
}

shinyApp(ui, server)
