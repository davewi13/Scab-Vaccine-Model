#' Plot OAT sensitivity — faceted by param_name × node,
#' coloured by relative position (low / mid / high) within each parameter
#'
#' @param summary_dt    Output of summarise_sensitivity() on OAT results
#' @param compartments  Which compartments to show
#' @param ribbon        Show 90% CI ribbon?
#' @param base_size     Base font size passed to theme_minimal() (default: 11).
#'                      All text elements scale relative to this value.
plot_oat <- function(summary_dt,
                     compartments = c("S", "I", "C", "V"),
                     ribbon       = TRUE,
                     base_size    = 11) {
  
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
      c("low",
        if (n > 2) paste0("mid", if (n > 3) seq_len(n - 2) else ""),
        "high")
    }
    factor(labels[match(param_value, uvals)], levels = labels)
  }, by = param_name]
  
  # Build palette
  all_ranks <- sort(unique(dt$value_rank))
  n_ranks   <- length(all_ranks)
  rank_pal  <- setNames(
    viridisLite::viridis(n_ranks, option = "plasma", end = 0.9),
    all_ranks
  )
  
  p <- ggplot(dt, aes(x = time, y = mean,
                      colour   = value_rank,
                      linetype = compartment,
                      group    = interaction(param_value_chr, compartment))) +
    geom_line(linewidth = 0.8)
  
  if (ribbon)
    p <- p + geom_ribbon(aes(ymin = lo95, ymax = hi95, fill = value_rank),
                         alpha = 0.15, colour = NA)
  
  p +
    scale_colour_manual(values = rank_pal, name = "Relative value") +
    scale_fill_manual(values = rank_pal,   name = "Relative value") +
    facet_grid(param_name ~ node) +
    coord_cartesian(clip = "off") +
    labs(
      colour   = "Relative value",
      fill     = "Relative value",
      linetype = "Compartment",
      y        = "Within-flock prevalence (%)",
      x        = "Time (days)"
    ) +
    theme_minimal(base_size = base_size) +
    theme(
      plot.margin      = margin(5, 60, 5, 5),
      strip.text       = element_text(size = rel(0.9)),
      axis.title       = element_text(size = rel(1.0)),
      axis.text        = element_text(size = rel(0.85)),
      legend.title     = element_text(size = rel(0.9)),
      legend.text      = element_text(size = rel(0.85))
    )
}