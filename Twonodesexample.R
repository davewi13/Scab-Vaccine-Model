library(SimInf)


compartments <- c("S", "I", "C", "V")

data("nodes", package = "SimInf")

#take just the first two nodes 
nodes <- nodes[1:2,]


u0 <- data.frame(S = c(95,100), I= c(5,0), C = c(rep(0, nrow(nodes))), V = c(rep(0, nrow(nodes))), D = c(rep(0, nrow(nodes))), H= c(rep(0, nrow(nodes))))

d_ik <- distance_matrix(x = nodes$x , y= nodes$y, cutoff = 2500)

set.seed(123)

vaccination <- data.frame(event = "intTrans", time = c(5,10), node = 1:2, 
                          dest = 0, n = 0, proportion = 0.8, select = 3, shift = 1)

model <- SISe_sp(u0 = u0, tspan = 1:365, events = vaccination, phi = 0,
                 upsilon = (0.2519243), gamma = (1/77), alpha = 0.005619795, qprop = (1/2), tau = (1/653), epar = (1/3), wane = (1/653), vaccshed=(1/6),
                 beta_t1 = 0.06832791, beta_t2 = 0.06832791, beta_t3 = 0.06832791, beta_t4 = 0.06832791, end_t1 = 91,
                 end_t2 = 182, end_t3 = 273, end_t4 = 365, distance = d_ik,
                 coupling = 0.5642857)

model_novac <- SISe_sp(u0 = u0, tspan = 1:365, events = NULL, phi = 0,
                 upsilon = (0.2519243), gamma = (1/77), alpha = 0.005619795, qprop = (1/2), tau = (1/653), epar = (1/3), wane = (1/653), vaccshed=(1/6), 
                 beta_t1 = 0.06832791, beta_t2 = 0.06832791, beta_t3 = 0.06832791, beta_t4 = 0.06832791, end_t1 = 91,
                 end_t2 = 182, end_t3 = 273, end_t4 = 365, distance = d_ik,
                 coupling = 0.5642857)

nsim <- 100




############################
############################
##  Prevalence #############
############################
############################

results_vac <- lapply(1:nsim, function(i) {
  res <- run(model = model, threads = 1)
  df <- prevalence(res, formula = I + C ~ S + I + C + V, type = "wnp")
  df$sim <- i
  df$vaccination <- TRUE
  df
})

results_novac <- lapply(1:nsim, function(i) {
  res <- run(model = model_novac, threads = 1)
  df <- prevalence(res, formula = I + C ~ S + I + C + V, type = "wnp")
  df$sim <- i
  df$vaccination <- FALSE
  df
})

withinnodeprev <- do.call(rbind, c(results_vac, results_novac))
withinnodeprev$node <- factor(withinnodeprev$node)

library(ggplot2)
ggplot(withinnodeprev, aes(x = time, y = prevalence, group = interaction(node, sim), colour = node)) +
  geom_line(alpha = 0.2) +
  stat_summary(aes(group=NULL),
               fun = mean,
               geom = "line",
               linewidth = 1.2) +
  facet_wrap(~vaccination) +
  theme_minimal()

plot(cases, main = "", xlim = c(0, 365), xlab = "Time", ylab = "Number of cases", 
     do.points = FALSE)




###############################
###############################
## Raw numbers ################
###############################
###############################

results_vac <- lapply(1:nsim, function(i) {
  res <- run(model = model, threads = 1)
  df <- trajectory(res)
  df$sim <- i
  df$vaccination <- TRUE
  df
})

results_novac <- lapply(1:nsim, function(i) {
  res <- run(model = model_novac, threads = 1)
  df <- trajectory(res)
  df$sim <- i
  df$vaccination <- FALSE
  df
})

library(dplyr)
library(tidyr)

df_all <- bind_rows(results_vac, results_novac)

df_long <- df_all %>%
  pivot_longer(cols = c(S, I, C, V),
               names_to = "compartment",
               values_to = "count")

df_long$node <- factor(df_long$node)

ggplot(df_long,
       aes(x = time, y = count,
           group = interaction(sim, node, compartment), colour=compartment)) +
  geom_line(alpha = 0.2) +
  stat_summary(aes(linetype = node,
                   group = interaction(node, compartment)),
               fun = mean,
               geom = "line",
               linewidth = 1.2) +
  facet_wrap(~vaccination) +
  theme_minimal()

## Create an 'epicurve' function to estimate the average number of new cases per
## day from n = 1000 realizations. To clear infection that was introduced in the
## previous trajectory, animals are first moved to the susceptible compartment.
## Then, one infected individual is introduced into a randomly sampled node from
## the population.  Note that we use the 'L' suffix to create an integer value
## rather than a numeric value.  Run the model and accumulate 'Icum'.  For
## efficiency, we use 'as.is = TRUE', the internal matrix format, to extract
## 'Icum' in every node at each time-point in 'tspan'.
epicurve <- function(model, n = 1000) {
  Icum <- numeric(length(model@tspan))
  for (i in seq_len(n)) {
    ## Move all infected individuals to the susceptible compartment. This is to remove
    ## infection from the node with one infected individual from the previous
    ## trajectory.
    model@u0["S", ] <- model@u0["S", ] + model@u0["I", ] + model@u0["C", ] + model@u0["V", ]
    model@u0["I", ] <- 0L
    model@u0["C", ] <- 0L
    model@u0["V", ] <- 0L
    
    ## Sample one node in the population where to introduce infection with one
    ## infected individual.
    j <- sample(seq_len(Nn(model)), 1)
    model@u0["I", j] <- 1L
    model@u0["S", j] <- model@u0["S", j] - 1L
    
    ## Run the model
    result <- run(model = model)
    
    ## Accumulate Icum.
    traj <- trajectory(model = result, compartments = "I", as.is = TRUE)
    Icum <- Icum + colSums(traj)
  }
  
  stepfun(model@tspan[-1], diff(c(0, Icum/n)))
}
