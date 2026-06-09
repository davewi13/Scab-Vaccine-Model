#' Run OAT sensitivity analysis
#'
#' @param param_ranges  Named list; each element is a numeric vector of values
#'                      to test for that parameter, e.g.
#'                      list(gamma = c(1/100, 1/77, 1/50),
#'                           coupling = c(0.1, 0.5, 0.9))
#' @param nsim          Replicates per scenario
#' @param tspan         Time span
#' @param seed          RNG seed
#' @param cutoff        The distance (in metres) beyond which no transmission between farms occurs
#' @param u0            Data frame (columns S, I, C, V, D, H and one row per farm) giving initial number in each compartment
#' @param vaccination   Data frame of the form specified by SimInf scheduling vaccination events
#' @return              data.table with a `scenario` and `param_value` column
run_oat_sensitivity <- function(param_ranges,
                                nsim  = 100,
                                tspan = seq(1, 730, 1),
                                seed  = 123,
                                cutoff      = 1999,
                                u0          = NULL,
                                vaccination = NULL) {
  
  all_results <- list()
  
  for (param_name in names(param_ranges)) {
    values <- param_ranges[[param_name]]
    message(sprintf("OAT: varying '%s' (%d values)", param_name, length(values)))
    
    for (val in values) {
      params_i           <- setNames(list(val), param_name)
      scenario_label     <- sprintf("%s = %s", param_name, signif(val, 3))
      message("  Running: ", scenario_label)
      
      dt                 <- run_scenario(params_i, nsim, tspan, seed, cutoff, u0, vaccination)
      dt$param_name      <- param_name
      dt$param_value     <- val
      dt$param_value_chr <- as.character(signif(val, 3))
      dt$scenario        <- scenario_label
      all_results[[length(all_results) + 1]] <- dt
    }
  }
  
  rbindlist(all_results)
}