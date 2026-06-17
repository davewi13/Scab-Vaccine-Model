#' Build and run the SISe_sp model for one parameter set
#'
#' @param params      Named list of parameters (merged over base_params)
#' @param nsim        Number of stochastic replicates
#' @param tspan       Time span vector
#' @param seed        RNG seed (set NULL to skip)
#' @param cutoff      The distance (in metres) beyond which no transmission between farms occurs
#' @param u0          Data frame (columns S, I, C, V, D, H and one row per farm) giving initial number in each compartment
#' @param vaccination Data frame of the form specified by SimInf scheduling vaccination events
#' @return            data.table with columns: time, node, S, I, C, V, sim
run_scenario <- function(params      = list(),
                         nsim        = 100,
                         tspan       = seq(1, 730, 1),
                         seed        = 123,
                         cutoff      = 1999,
                         u0          = NULL,
                         vaccination = NULL,
                         nodes       = NULL) {
  
  p <- params   # overlay user params on defaults
  
  
  d_ik <- distance_matrix(
    x = nodes[,1],
    y = nodes[,2],
    cutoff = cutoff
  )
  
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