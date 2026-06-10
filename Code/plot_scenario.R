plot_scenario <- function(data,
                          compartments = c("C", "I"),
                          alpha_lines  = 0.1,
                          linewidth_individual = 0.3,
                          linewidth_median     = 1,
                          alpha_ribbon = 0.3,
                          pi_lower     = 0.025,
                          pi_upper     = 0.975) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Pivot to long if needed
  if (any(c("S", "I", "C", "V") %in% names(data))) {
    data <- pivot_longer(data,
                         cols      = any_of(c("S", "I", "C", "V")),
                         names_to  = "compartment",
                         values_to = "count")
  }
  
  # Filter to requested compartments
  data_filt <- data %>%
    filter(compartment %in% compartments)
  
  # ✅ Create readable labels (requires 'status' column in data)
  if ("status" %in% names(data_filt)) {
    data_filt <- data_filt %>%
      mutate(node_label = paste0("Farm ", node, " (", status, ")"))
  } else {
    data_filt <- data_filt %>%
      mutate(node_label = paste0("Farm ", node))
  }
  
  # ✅ Summary WITH node_label included (this was the key fix)
  q_summary <- data_filt %>%
    group_by(node, node_label, time, compartment) %>%
    summarise(
      median = median(count),
      lower  = quantile(count, pi_lower),
      upper  = quantile(count, pi_upper),
      .groups = "drop"
    )
  
  # ✅ Plot
  ggplot() +
    geom_line(
      data = data_filt,
      aes(x = time, y = count, colour = compartment,
          group = interaction(sim, compartment, node)),
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
    facet_wrap(~ node_label) +
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