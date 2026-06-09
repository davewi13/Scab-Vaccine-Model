library(SimInf)

compartments <- c("S", "I", "C", "V")

data("nodes", package = "SimInf")

#take just the first two nodes 
nodes <- nodes[1:3,]
# These nodes are so far apart from each other that they would never infect each other
# Make up some nodes with a small distance between
distance_between <- 1000
nodes[1,1] <- nodes[2,1]-distance_between
nodes[1,2] <- nodes[2,2]
nodes[3,1] <- nodes[2,1]+distance_between
nodes[3,2] <- nodes[2,2]

u0 <- data.frame(S = c(100,95,100), I= c(0,5,0), C = c(rep(0, nrow(nodes))), V = c(rep(0, nrow(nodes))), D = c(rep(0, nrow(nodes))), H= c(rep(0, nrow(nodes))))

d_ik <- distance_matrix(x = nodes$x , y= nodes$y, cutoff = 1999)

set.seed(123)

vaccination <- data.frame(event = "intTrans", time = c(1), node = 3, 
                          dest = 0, n = 100, proportion = 0, select = 3, shift = 1)

model <- SISe_sp(u0 = u0, tspan = seq(1, 2000, 1), events = vaccination, phi = 0,
                 upsilon = (0.002519243*100), gamma = (1/77), alpha = 0.005619795, qprop = (1/2), tau = (1/653), epar = (1/3), wane = (0),
                 beta_t1 = 0.06832791, beta_t2 = 0.06832791, beta_t3 = 0.06832791, beta_t4 = 0.06832791, end_t1 = 91,
                 end_t2 = 182, end_t3 = 273, end_t4 = 365, distance = d_ik,
                 coupling = 0.5642857)

nsim <- 100

results <- lapply(1:nsim, function(i) {
  res <- run(model = model, threads = 1)
  df <- trajectory(res)
  df$sim <- i
  df$vaccination <- TRUE
  df
})

library(data.table)
library(dplyr)
library(tidyr)

df <- rbindlist(results)

df_long <- df %>%
  pivot_longer(cols = c(S, I, C, V),
               names_to = "compartment",
               values_to = "count")

df_long$node <- factor(df_long$node, levels=c(1,2,3), labels=c("Unvaccinated", "Source", "Vaccinated"))

library(ggplot2)

ggplot(df_long[df_long$vaccination == T,],
       aes(x = time, y = count,
           group = interaction(sim, compartment), colour=compartment)) +
  geom_line(alpha = 0.2) +
  stat_summary(aes(group = interaction(node, compartment)),
               fun = mean,
               geom = "line",
               linewidth = 1.2) +
  facet_wrap(~node) +
  theme_minimal() +
  scale_y_continuous(breaks=seq(0,100,10)) +
  ylab("Within-flock prevalence (%)") +
  xlab("Time (days)")
  




# Code from Claude AI prompted with the following "Write some R wrapper code 
# so that I can run sensitivity analyses on this model for the specified 
# parameters. My current code is below (inserted code above).

# ============================================================
# Sensitivity Analysis Wrapper for SISe_sp Model
# ============================================================

library(SimInf)
library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

# ============================================================
# 1. BASE MODEL SETUP (unchanged from your original code)
# ============================================================

setup_nodes <- function(distance_between = 1000) {
  data("nodes", package = "SimInf")
  nodes <- nodes[1:3, ]
  nodes[1, 1] <- nodes[2, 1] - distance_between
  nodes[1, 2] <- nodes[2, 2]
  nodes[3, 1] <- nodes[2, 1] + distance_between
  nodes[3, 2] <- nodes[2, 2]
  nodes
}

base_params <- list(
  phi        = 0,
  upsilon    = 0.002519243 * 100,
  gamma      = 1 / 77,
  alpha      = 0.005619795,
  qprop      = 1 / 2,
  tau        = 1 / 653,
  epar       = 1 / 3,
  wane       = 0,
  beta_t1    = 0.06832791,
  beta_t2    = 0.06832791,
  beta_t3    = 0.06832791,
  beta_t4    = 0.06832791,
  end_t1     = 91,
  end_t2     = 182,
  end_t3     = 273,
  end_t4     = 365,
  coupling   = 0.5642857
)

# ============================================================
# 2. CORE HELPER: build & run a single parameterisation
# ============================================================

#' Build and run the SISe_sp model for one parameter set
#'
#' @param params  Named list of parameters (merged over base_params)
#' @param nsim    Number of stochastic replicates
#' @param tspan   Time span vector
#' @param seed    RNG seed (set NULL to skip)
#' @return        data.table with columns: time, node, S, I, C, V, sim
run_scenario <- function(params      = list(),
                         nsim        = 100,
                         tspan       = seq(1, 730, 1),
                         seed        = 123,
                         cutoff      = 1999,
                         u0          = NULL,
                         vaccination = NULL) {
  
  p <- modifyList(base_params, params)   # overlay user params on defaults
  
  nodes <- setup_nodes()
  d_ik  <- distance_matrix(x = nodes$x, y = nodes$y, cutoff = cutoff)
  
  model <- SISe_sp(
    u0      = u0,
    tspan   = tspan,
    events  = vaccination,
    phi     = p$phi,
    upsilon = p$upsilon,
    gamma   = p$gamma,
    alpha   = p$alpha,
    qprop   = p$qprop,
    tau     = p$tau,
    epar    = p$epar,
    wane    = p$wane,
    beta_t1 = p$beta_t1, beta_t2 = p$beta_t2,
    beta_t3 = p$beta_t3, beta_t4 = p$beta_t4,
    end_t1  = p$end_t1,  end_t2  = p$end_t2,
    end_t3  = p$end_t3,  end_t4  = p$end_t4,
    distance = d_ik,
    coupling = p$coupling
  )
  
  if (!is.null(seed)) set.seed(seed)
  
  results <- lapply(seq_len(nsim), function(i) {
    df      <- trajectory(run(model, threads = 1))
    df$sim  <- i
    df
  })
  
  rbindlist(results)
}

# ============================================================
# 3. ONE-AT-A-TIME (OAT) SENSITIVITY ANALYSIS
# ============================================================

#' Run OAT sensitivity analysis
#'
#' @param param_ranges  Named list; each element is a numeric vector of values
#'                      to test for that parameter, e.g.
#'                      list(gamma = c(1/100, 1/77, 1/50),
#'                           coupling = c(0.1, 0.5, 0.9))
#' @param nsim          Replicates per scenario
#' @param tspan         Time span
#' @param seed          RNG seed
#' @return              data.table with a `scenario` and `param_value` column
run_oat_sensitivity <- function(param_ranges,
                                nsim  = 100,
                                tspan = seq(1, 730, 1),
                                seed  = 123) {
  
  all_results <- list()
  
  for (param_name in names(param_ranges)) {
    values <- param_ranges[[param_name]]
    message(sprintf("OAT: varying '%s' (%d values)", param_name, length(values)))
    
    for (val in values) {
      params_i           <- setNames(list(val), param_name)
      scenario_label     <- sprintf("%s = %s", param_name, signif(val, 4))
      message("  Running: ", scenario_label)
      
      dt                 <- run_scenario(params_i, nsim, tspan, seed)
      dt$param_name      <- param_name
      dt$param_value     <- val
      dt$param_value_chr <- as.character(signif(val, 4))
      dt$scenario        <- scenario_label
      all_results[[length(all_results) + 1]] <- dt
    }
  }
  
  rbindlist(all_results)
}

# ============================================================
# 4. GRID / FACTORIAL SENSITIVITY ANALYSIS
# ============================================================

#' Run a full factorial grid over a set of parameters
#'
#' @param param_grid  Named list of vectors — all combinations are tested, e.g.
#'                    list(gamma = c(1/100, 1/77), coupling = c(0.3, 0.56))
#' @param nsim        Replicates per scenario
#' @param tspan       Time span
#' @param seed        RNG seed
#' @return            data.table tagged with one column per parameter
run_grid_sensitivity <- function(param_grid,
                                 nsim  = 100,
                                 tspan = seq(1, 730, 1),
                                 seed  = 123) {
  
  grid        <- expand.grid(param_grid, stringsAsFactors = FALSE)
  all_results <- list()
  
  message(sprintf("Grid sensitivity: %d combinations × %d sims each", nrow(grid), nsim))
  
  for (i in seq_len(nrow(grid))) {
    params_i       <- as.list(grid[i, , drop = FALSE])
    scenario_label <- paste(
      mapply(function(k, v) sprintf("%s=%s", k, signif(v, 3)),
             names(params_i), params_i),
      collapse = ", "
    )
    message(sprintf("  [%d/%d] %s", i, nrow(grid), scenario_label))
    
    dt             <- run_scenario(params_i, nsim, tspan, seed)
    dt$scenario    <- scenario_label
    dt$grid_row    <- i
    for (nm in names(params_i)) dt[[nm]] <- params_i[[nm]]   # one col per param
    
    all_results[[i]] <- dt
  }
  
  rbindlist(all_results, fill = TRUE)
}

# ============================================================
# 5. SUMMARY HELPER
# ============================================================

#' Summarise raw trajectory data for plotting
#'
#' @param dt            data.table from run_*_sensitivity()
#' @param group_cols    Additional grouping columns (besides time, node, compartment)
#' @return              Summarised data.table with mean / 5th / 95th percentiles
summarise_sensitivity <- function(dt, group_cols = c("scenario", "param_name", "param_value")) {
  
  # keep only group cols that actually exist in dt
  group_cols <- intersect(group_cols, names(dt))
  
  dt_long <- pivot_longer(dt, cols = c(S, I, C, V),
                          names_to = "compartment", values_to = "count")
  setDT(dt_long)
  
  dt_long[, node := factor(node, levels = 1:3,
                           labels = c("Unvaccinated", "Source", "Vaccinated"))]
  
  grp <- c("time", "node", "compartment", group_cols)
  
  dt_long[, .(
    mean  = mean(count),
    lo95  = quantile(count, 0.05),
    hi95  = quantile(count, 0.95)
  ), by = grp]
}

# ============================================================
# 6. PLOTTING HELPERS
# ============================================================

#' Plot OAT sensitivity — faceted by param_name × node,
#' coloured by relative position (low / mid / high) within each parameter
#'
#' @param summary_dt    Output of summarise_sensitivity() on OAT results
#' @param compartments  Which compartments to show
#' @param ribbon        Show 90% CI ribbon?
plot_oat <- function(summary_dt,
                     compartments = c("S", "I", "C", "V"),
                     ribbon = TRUE) {
  
  dt <- copy(summary_dt[compartment %in% compartments])
  
  # Assign a relative rank label within each parameter
  dt[, value_rank := {
    uvals  <- sort(unique(param_value))
    n      <- length(uvals)
    labels <- if (n == 1) {
      "mid"
    } else if (n == 2) {
      c("low", "high")
    } else {
      # Evenly spread labels: low, ..., high
      c("low",
        if (n > 2) paste0("mid", if (n > 3) seq_len(n - 2) else ""),
        "high")
    }
    factor(labels[match(param_value, uvals)], levels = labels)
  }, by = param_name]
  
  # Build a palette that maps rank positions to colours
  # (reused identically across every parameter panel)
  all_ranks  <- sort(unique(dt$value_rank))
  n_ranks    <- length(all_ranks)
  rank_pal <- setNames(
    viridisLite::viridis(n_ranks, option = "plasma", end = 0.9),
    all_ranks
  )
  
  p <- ggplot(dt, aes(x = time, y = mean,
                      colour = value_rank,
                      linetype = compartment,
                      group = interaction(param_value_chr, compartment))) +
    geom_line(linewidth = 0.8)
  
  if (ribbon)
    p <- p + geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = value_rank),
                         alpha = 0.15, colour = NA)
  
  p +
    scale_colour_manual(
      values = rank_pal,
      labels = function(x) x,   # show "low", "mid", "high" etc.
      name   = "Relative value"
    ) +
    scale_fill_manual(
      values = rank_pal,
      labels = function(x) x,
      name   = "Relative value"
    ) +
    # Secondary annotation: show the actual numeric values per panel
    geom_text(
      data = dt[time == max(time),
                .(mean = mean(mean)), 
                by = .(node, param_name, param_value_chr, compartment, value_rank)],
      aes(x = max(dt$time), y = mean, label = param_value_chr,
          colour = value_rank),
      hjust = -0.05, size = 2.8, show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    facet_grid(param_name ~ node) +
    coord_cartesian(clip = "off") +
    labs(
      colour   = "Relative value",
      fill     = "Relative value",
      linetype = "Compartment",
      y        = "Within-flock prevalence (%)",
      x        = "Time (days)",
      title    = "One-at-a-time sensitivity analysis"
    ) +
    theme_minimal() +
    theme(plot.margin = margin(5, 60, 5, 5))  # right margin for labels
}

#' Plot grid sensitivity — faceted by node and scenario
#'
#' @param summary_dt    Output of summarise_sensitivity() on grid results
#' @param compartments  Which compartments to show
#' @param ribbon        Show 90 % CI ribbon?
plot_grid <- function(summary_dt,
                      compartments = c("S", "I", "C", "V"),
                      ribbon = TRUE) {
  
  dt <- summary_dt[compartment %in% compartments]
  
  p <- ggplot(dt, aes(x = time, y = mean,
                      colour = compartment,
                      group = compartment)) +
    geom_line(linewidth = 0.8)
  
  if (ribbon)
    p <- p + geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = compartment),
                         alpha = 0.15, colour = NA)
  
  p +
    facet_grid(scenario ~ node) +
    labs(colour = "Compartment", fill = "Compartment",
         y = "Within-flock prevalence (%)", x = "Time (days)",
         title = "Grid sensitivity analysis") +
    theme_minimal() +
    theme(strip.text.y = element_text(size = 7))
}

# ============================================================
# 7. EXAMPLE USAGE
# ============================================================

## --- One-at-a-time ---
oat_ranges <- list(
  upsilon  = c(0.001*100, 0.002519243*100, 0.006*100),
  alpha    = c(0.0005619795, 0.005619795, 0.012)
)

oat_raw     <- run_oat_sensitivity(oat_ranges, nsim = 100)
oat_summary <- summarise_sensitivity(oat_raw,
                                     group_cols = c("scenario", "param_name",
                                                    "param_value", "param_value_chr"))
plot_oat(oat_summary, compartments = c("I", "C"))

## --- One-at-a-time ---
oat_ranges <- list(
  gamma    = c(1/100, 1/77, 1/50),
  tau      = c(1/2000, 1/653, 1/200)
)

oat_raw     <- run_oat_sensitivity(oat_ranges, nsim = 100)
oat_summary <- summarise_sensitivity(oat_raw,
                                     group_cols = c("scenario", "param_name",
                                                    "param_value", "param_value_chr"))
plot_oat(oat_summary, compartments = c("I", "C"))

## --- One-at-a-time ---
oat_ranges <- list(
  qprop    = c(0.1, 0.5, 0.9),
  epar     = c(0.5, 1/3, 0.01)
)

oat_raw     <- run_oat_sensitivity(oat_ranges, nsim = 100)
oat_summary <- summarise_sensitivity(oat_raw,
                                     group_cols = c("scenario", "param_name",
                                                    "param_value", "param_value_chr"))
plot_oat(oat_summary, compartments = c("I", "C"))

## --- One-at-a-time ---
oat_ranges <- list(
  coupling = c(0.1, 0.4, 0.8)
)

oat_raw     <- run_oat_sensitivity(oat_ranges, nsim = 100)
oat_summary <- summarise_sensitivity(oat_raw,
                                     group_cols = c("scenario", "param_name",
                                                    "param_value", "param_value_chr"))
plot_oat(oat_summary, compartments = c("I", "C"))

## --- Grid (factorial) ---
grid_ranges <- list(
  gamma    = c(1/100, 1/77, 1/50),
  coupling = c(0.1, 0.5642857, 0.9)
)

grid_raw     <- run_grid_sensitivity(grid_ranges, nsim = 100)
grid_summary <- summarise_sensitivity(grid_raw,
                                      group_cols = c("scenario", "gamma", "coupling"))
plot_grid(grid_summary, compartments = c("I", "C"))
