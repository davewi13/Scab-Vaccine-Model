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