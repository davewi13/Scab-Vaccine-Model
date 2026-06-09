#' Plot compartment trajectories with median and prediction interval ribbon
#'
#' @param data        Wide data frame with columns node, time, sim and compartment
#'                    columns (S, I, C, V), or an already-pivoted long data frame
#'                    with columns node, time, sim, compartment, count
#' @param compartments Character vector of compartment names to plot
#'                    (default: c("C", "I"))
#' @param alpha_lines  Transparency of individual simulation lines (default: 0.1)
#' @param linewidth_individual Line width for individual simulation lines (default: 0.3)
#' @param linewidth_median     Line width for median trajectory line (default: 1)
#' @param alpha_ribbon Transparency of the prediction interval ribbon (default: 0.3)
#' @param pi_lower     Lower quantile for prediction interval (default: 0.025)
#' @param pi_upper     Upper quantile for prediction interval (default: 0.975)
#' @return             A ggplot object facetted by node, with individual simulation
#'                     lines, a shaded 95% prediction ribbon, and a median trajectory
#'                     line, coloured by compartment
plot_scenario <- function(data,
                          compartments = c("C", "I"),
                          alpha_lines  = 0.1,
                          linewidth_individual = 0.3,
                          linewidth_median     = 1,
                          alpha_ribbon = 0.3,
                          pi_lower     = 0.025,
                          pi_upper     = 0.975) {
  
  # Pivot to long if needed (accepts either wide or already-long data)
  if (any(c("S", "I", "C", "V") %in% names(data))) {
    data <- pivot_longer(data,
                         cols      = any_of(c("S", "I", "C", "V")),
                         names_to  = "compartment",
                         values_to = "count")
  }
  
  # Filter to requested compartments
  data_filt <- data %>% filter(compartment %in% compartments)
  
  # Summary statistics
  q_summary <- data_filt %>%
    group_by(node, time, compartment) %>%
    summarise(
      median = median(count),
      lower  = quantile(count, pi_lower),
      upper  = quantile(count, pi_upper),
      .groups = "drop"
    )
  
  ggplot() +
    geom_line(
      data = data_filt,
      aes(x = time, y = count, colour = compartment,
          group = interaction(sim, compartment)),
      alpha     = alpha_lines,
      linewidth = linewidth_individual
    ) +
    geom_ribbon(
      data = q_summary,
      aes(x = time, ymin = lower, ymax = upper, fill = compartment),
      alpha = alpha_ribbon
    ) +
    geom_line(
      data = q_summary,
      aes(x = time, y = median, colour = compartment),
      linewidth = linewidth_median
    ) +
    facet_wrap(~ node) +
    labs(
      x      = "Time (days)",
      y      = "Count",
      colour = "Compartment",
      fill   = "Compartment"
    ) +
    theme_bw() +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey90")
    )
}