library(SimInf)
library(shiny)
library(shinyjs)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)
library(shinyWidgets)

source("run_scenario.R")
source("plot_scenario.R")


setup_nodes <- function(Nfarms) {
  x=runif(Nfarms)*10000
  y=runif(Nfarms)*10000
  
  nodes=cbind(x,y)
  nodes
}


parse_numeric_vector <- function(x) {
  out <- strsplit(x, ",")[[1]]
  out <- trimws(out)
  as.numeric(out)
}

make_vacc_times <- function(n) {
  if (n == 0) return("")
  paste(rep(1, n), collapse = ", ")
}

calc_infections <- function(df) {
  
  df %>%
    arrange(sim, node, time) %>%
    group_by(sim, node) %>%
    mutate(dI = I - lag(I, default = first(I))) %>%
    
    summarise(
      infections = sum(pmax(dI, 0), na.rm = TRUE),
      .groups = "drop"
    )
}



ui <- fluidPage(
  useShinyjs(),  # Initialize shinyjs
  titlePanel("Sheep Scab Transmission Model"),
  
  # Description paragraph
  p(HTML("This app simulates the transmission of sheep scab mites between and within a network of farms based on an adaptation of the model model published by Nixon et al (2021). 
  Within a farm, animals are classed as susceptible, vaccinated, infested or carriers. Mites are transmitted between sheep and farms via contamination of an environmental reservoir that accumulates due to shedding of infested sheep and carriers.
  Transmission between farms is scaled depending on the distance between two farms with an upper distance limit for transmission imposed, beyond which transmission could only occur via direct animal movements.
  Vaccination is assumed to protect animals from becoming infested at clinical levels, meaning that vaccinated animals that are exposed to mites bypass the infested class and move directly to the carrier class. <br> <br>
  
  The user has the option to specify the exact location of farms by choosing the farm setup mode and then either pointing and clicking on the grid or selecting random mode and using the 'Generate random farms' button. <br>
  The vaccination and initial infection status of farms can then be toggled by selecting the relevant option from the menu and clicking on the farms. <br>
  The total number of farms and number of vaccinated farms will be updated automatically when the user interacts with the map or they can be set manually before generating random farms. <br>
  Vaccination times must be specified for each vaccinated farm in the order specified by the numbers assigned to the farms on the grid. <br>
  Users can set the vaccine efficacy and the maximum distance at which scab can travel between farms in the absence of trade. <br> <br>
         
  The model assumes that all farms are the same size and that the simulations run for 2 years. (Manual choice of simulation time would be easy - manual entry of farm sizes would be easy too but input intensive for the user).
  There is a further button to toggle advanced parameters. These parameters were extracted from the original Nixon et al (2021) model and it is suggested that these only be changed by experienced users.")),
  
  radioButtons(
    "farm_mode",
    "Farm setup mode:",
    choices = c("Manual (click)", "Random (generate)"),
    selected = "Manual (click)"
  ),
  
  fluidRow(
    column(
      width = 12,
      
      plotOutput("farmEditor", height = "500px", click = "plot_click"),
      
      actionButton("clearFarms", "Clear farms"),
      actionButton("generateFarms", "Generate random farms"),
      
      radioButtons(
        "click_mode",
        "Click action:",
        choices = c("Add farm", "Toggle vaccination", "Toggle infection"),
        selected = "Add farm"
      )
    )
  ),
  
  sidebarLayout(
    sidebarPanel(width=3,
      numericInput("Nfarms", "Number of farms", value=5, min=1, max=10),
      numericInput("Nvacc", "Number of farms vaccinated", value=2, min=0, max=10),
      numericInput("Nsheep", "Number of sheep per farm (all farms are assumed to be equal size)", value=100, min=10, max=10000),
      textInput("Vacc_times", "Vaccination times (comma-separated, one per vaccinated farm, in order)", value = "10, 25"),
      numericInput("efficacy", "Vaccine efficacy", value=0.8, min=0, max=1),
      numericInput("coupling", "Spatial coupling between farms", value=0.5, min=0, max=1),
      numericInput("cutoff", "Max. distance at which scab can travel between farms (without trade) in km", value=5, min=0, max=50),
      
      actionButton("toggleAdvanced", "Toggle Adv. Parameters"),
      div(id = "advancedParams",
          numericInput("alpha", "Daily shedding rate of infectious animal", value=signif(0.005619795,3), min=0, max=1),
          numericInput("beta", "Decay rate of environmental infectious pressure", value=0.0680, min=0, max=1),
          numericInput("upsilon", "Transmission rate from environment to susceptible sheep", value=signif(0.002519243*100,3), min=0, max=1),
          numericInput("qprop", "Proportion of acute infections that become carriers", value=0.5, min=0, max=1),
          numericInput("tau", "Recovery rate for carriers (1/days)", value=signif(1/653,3), min=0, max=1),
          numericInput("gamma", "Recovery rate for infectious sheep (1/days)", value=signif(1/77,3), min=0, max=1),
          numericInput("epar", "Scaling rate for the contribution of carriers to infectious pressure", value=signif(1/3,3), min=0, max=1),
          numericInput("wane", "Rate at which vaccine-based immunity wanes", value=signif(1/653,3), min=0, max=1),
          numericInput("nsim", "Number of simulations", 100)
      ),
      actionButton("run", "Run model")
    ),
    
    mainPanel(width=9,
      #plotOutput("farmNetworkPlot", height = "500px"),
      plotOutput("infectionReductionByFarm"),
      uiOutput("captionReductionByFarm"),
      plotOutput("infectionReductionGrouped"),
      uiOutput("captionReductionGrouped"),
      plotOutput("plot"),
      uiOutput("captionTrajectories")
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
    
    df <- data.frame(
      x = nodes[,1],
      y = nodes[,2],
      status = "Unvaccinated",
      infected = FALSE
    )
    
    # assign vaccination
    if (Nvacc > 0) {
      vacc_idx <- sample(seq_len(Nfarms), Nvacc)
      df$status[vacc_idx] <- "Vaccinated"
    }
    
    # assign infection
    infect_idx <- sample(seq_len(Nfarms), 1)
    df$infected[infect_idx] <- TRUE
    
    
    Nvacc <- sum(df$status == "Vaccinated")
    
    current_vals <- parse_numeric_vector(input$Vacc_times)
    
    if (length(current_vals) != Nvacc) {
      updateTextInput(
        session,
        "Vacc_times",
        value = make_vacc_times(Nvacc)
      )
    }
    
    
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
      
      new_row <- data.frame(
        x = click$x * 1000,
        y = click$y * 1000,
        status = "Unvaccinated",
        infected = FALSE
      )
      
      df <- dplyr::bind_rows(df, new_row)
    } else {
      
      # find nearest farm
      if (nrow(df) > 0) {
        d <- sqrt((df$x - click$x*1000)^2 + (df$y - click$y*1000)^2)
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
    
    
    
    Nvacc <- sum(df$status == "Vaccinated")
    
    current_vals <- parse_numeric_vector(input$Vacc_times)
    
    if (length(current_vals) != Nvacc) {
      updateTextInput(
        session,
        "Vacc_times",
        value = make_vacc_times(Nvacc)
      )
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
    
    # Validate inputs
    valid <- TRUE
    messages <- c()
    
    if((input$Nfarms > 10) || (input$Nfarms < 2)) messages <- c(messages, "Please specify 1-10 farms.")
    if(input$Nvacc > input$Nfarms) messages <- c(messages, "Number of vaccinated farms must be less than or equal to total number of farms.")
    if((input$Nsheep > 10000) || (input$Nsheep < 10)) messages <- c(messages, "Please between 10 and 10000 sheep per farm.")
    if((input$efficacy > 1) || (input$efficacy < 0)) messages <- c(messages, "Vaccine efficacy must be between 0 and 1.")
    if((input$coupling < 0) || (input$efficacy > 1)) messages <- c(messages, "Spatial coupling parameter must be between 0 and 1.")
    if(input$alpha < 0) messages <- c(messages, "Daily shedding rate of infectious animals must not be negative.")
    if(input$beta < 0) messages <- c(messages, "Decay rate of environmental infectious pressure must not be negative.")
    if(input$upsilon < 0) messages <- c(messages, "The transmission rate from environment to susceptible sheep must not be negative.")
    if((input$qprop < 0) || (input$qprop > 1)) messages <- c(messages, "The proportion of acute infections that become carriers must be between 0 and 1.")
    if(input$tau < 0) messages <- c(messages, "The recovery rate for carriers must not be negative.")
    if(input$gamma < 0) messages <- c(messages, "The recovery rate for infectious sheep must not be negative.")
    if((input$epar < 0) || (input$epar > 1)) messages <- c(messages, "The scaling rate for the contribution of carriers to infectious pressure must be between 0 and 1.")
    if(input$wane < 0) messages <- c(messages, "The rate at which vaccine-based immunity wanes must not be negative.")
    if(input$nsim < 100) messages <- c(messages, "The number of simulations must be at least 100.")
    
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
      Nsheep     = input$Nsheep,
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
    
    if(length(vacc_times) != input$Nvacc) messages <- c(messages, "Please specify exactly one vaccination event per vaccinated farm")
  
    nodes <- as.matrix(s[, c("x", "y")])
    vaccinated_farms <- which(s$status == "Vaccinated")
  
    
    if (length(vaccinated_farms) == 0) {
      vaccination <- NULL
    } else {
      vaccination <- data.frame(
        event = "intTrans", time = vacc_times, node = vaccinated_farms,
        dest = 0, n = 0, proportion = input$efficacy, select = 3, shift = 1
      )
    }
    
    u0 <- data.frame(
      S = rep(input$Nsheep, nrow(s)),
      I = ifelse(s$infected, 1, 0),
      C = 0,
      V = 0,
      D = 0,
      H = 0
    )
    
    u0$S[u0$I == 1] <- input$Nsheep-1
    
    if(length(messages) > 0){
      # Show messages in the model output
      output$modelOutput <- renderText({
        paste(messages, collapse = "\n")
      })
      # Stop execution — do not run simulation
      return()
    }
    
    result_vacc <- run_scenario(
      params = params,
      nsim = input$nsim,
      tspan = seq(0, 730, 1),
      seed = 897342,
      cutoff = input$cutoff*1000,
      u0 = u0,
      vaccination = vaccination,
      nodes = nodes
    )
    
    result_vacc$status <- s$status[result_vacc$node]
    
    # --- WITHOUT vaccination ---
    result_novacc <- run_scenario(
      params = params,
      nsim = input$nsim,
      tspan = seq(0, 730, 1),
      seed = 897342,
      cutoff = input$cutoff*1000,
      u0 = u0,
      vaccination = NULL,
      nodes = nodes
    )
    
    result_novacc$status <- s$status[result_novacc$node]
    
    inf_vacc <- calc_infections(result_vacc)
    inf_novacc <- calc_infections(result_novacc)
    
    total_vacc <- inf_vacc %>%
      group_by(sim) %>%
      summarise(total = sum(infections, na.rm = TRUE), .groups = "drop") %>%
      mutate(bin = ntile(total, 10))   # 10 = deciles (adjustable)
    
    total_novacc <- inf_novacc %>%
      group_by(sim) %>%
      summarise(total = sum(infections, na.rm = TRUE), .groups = "drop") %>%
      mutate(bin = ntile(total, 10))
    
    
    inf_vacc <- inf_vacc %>%
      left_join(total_vacc %>% select(sim, bin), by = "sim")
    
    inf_novacc <- inf_novacc %>%
      left_join(total_novacc %>% select(sim, bin), by = "sim")
    
    vacc_summary <- inf_vacc %>%
      group_by(bin, node) %>%
      summarise(
        inf_vacc = mean(infections, na.rm = TRUE),
        .groups = "drop"
      )
    
    novacc_summary <- inf_novacc %>%
      group_by(bin, node) %>%
      summarise(
        inf_novacc = mean(infections, na.rm = TRUE),
        .groups = "drop"
      )
    
    inf_compare <- vacc_summary %>%
      left_join(novacc_summary, by = c("bin", "node")) %>%
      mutate(
        inf_vacc   = coalesce(inf_vacc, 0),
        inf_novacc = coalesce(inf_novacc, 0),
        
        reduction  = inf_novacc - inf_vacc,
        status = s$status[node]
      )
    
    
    
    list(
      result = result_vacc,   # used for trajectory plot
      setup = s,
      infections = inf_compare
    )
    
  })
  
  output$farmEditor <- renderPlot({
    
    df <- farms$df
    
    df$farm_id <- seq_len(nrow(df))
    
    
    ggplot(df, aes(x = x/1000, y = y/1000)) +   
      geom_blank() +
          xlim(0, 10000) +
          ylim(0, 10000) +
          coord_equal() +
          theme_minimal() +
          labs(title = "Click to place farms")
    
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
            aes(x = x/1000, y = y/1000, xend = xend/1000, yend = yend/1000,
                linewidth = linewidth),
            color = "grey40",
            alpha = 0.7
          )
      } +
      
      geom_point(
        data = df,
        aes(x = x/1000, y = y/1000, color = status, shape = infected),
        size = 4
      ) +
      
      geom_text(
        data = df,
        aes(x = x/1000, y = y/1000, label = farm_id),
        vjust = -1,
        size = 5
      ) +
      
      scale_color_manual(values = c(
        "Vaccinated" = "#1b9e77",
        "Unvaccinated" = "#d95f02"
      )) +
      
      scale_shape_manual(values = c(16, 17)) +
      scale_linewidth_identity() +
      
      coord_equal(xlim = c(0, 10), ylim = c(0, 10)) +
      theme_minimal(base_size = 18) +
      labs(
        title = "Click to place farms",
        color = "Status",
        shape = "Infected",
        x = "x location (km)",
        y = "y location (km)"
      )
  })
  
  output$infectionReductionByFarm <- renderPlot({
    
    req(sim())
    
    df <- sim()$infections
    
    
    df_summary <- df %>%
      group_by(node, status) %>%
      summarise(
        mean_reduction = mean(reduction, na.rm = TRUE),
        lower = quantile(reduction, 0.1, na.rm = TRUE),   # 10th percentile
        upper = quantile(reduction, 0.9, na.rm = TRUE),   # 90th percentile
        .groups = "drop"
      )
    
    
    ggplot(df_summary, aes(x = factor(node), y = mean_reduction, color = status)) +
      geom_point(size = 3) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
      labs(
        x = "Farm",
        y = "Reduction in infections",
        title = "Effect of vaccination by farm"
      ) +
      theme_minimal(base_size=18) +
      theme(legend.position = "bottom")
    
    
  })
  
  
  output$captionReductionByFarm <- renderUI({
    req(sim())
    
    p(
      "The estimated reduction in the number of clinical cases of infestation over the duration of the simulation due to the specified vaccination programme is shown for each farm individually. Colours represent vaccination status. Points represent the mean reduction and error bars show the central 80% of simulated outcomes."
    )
  })
  
  
  output$infectionReductionGrouped <- renderPlot({
    
    req(sim())
    
    df <- sim()$infections
    
    
    df_summary_group <- df %>%
      group_by(status) %>%
      summarise(
        mean_reduction = mean(reduction, na.rm = TRUE),
        lower = quantile(reduction, 0.1, na.rm = TRUE),
        upper = quantile(reduction, 0.9, na.rm = TRUE),
        .groups = "drop"
      )
    
    
    ggplot(df_summary_group, aes(x = status, y = mean_reduction, color = status)) +
      geom_point(size = 4) +
      geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2) +
      scale_color_manual(values = c(
        "Vaccinated" = "#1b9e77",
        "Unvaccinated" = "#d95f02"
      )) +
      labs(
        x = "",
        y = "Reduction in infections",
        title = "Effect of vaccination by farm type"
      ) +
      theme_minimal(base_size=18) +
      theme(legend.position = "bottom")
    
  })
  
  
  output$captionReductionGrouped <- renderUI({
    req(sim())
    p(
      "The estimated reduction in the number of clinical cases of infestation is shown aggregated across all farms and split by vaccination status. Points represent the mean reduction and error bars show the central 80% of simulated outcomes."
    )
  })
  
  output$plot <- renderPlot({
    req(sim())
    plot_scenario(sim()$result,
                  pi_lower=0.1,
                  pi_upper=0.9
  )})
  
  
  output$captionTrajectories <- renderUI({
    req(sim())
    p(
      "The simulated trajectories of infested sheep (I) and carriers (C) is shown for each farm over the course of the simulation. Each thin line represents a single model simulation, the thick lines represent the medians and the shaded areas represent upper and low 10% quantiles, such that the shaded bands contain the middle 80% of simulations."
      )
    
  })
  
}

shinyApp(ui, server)
