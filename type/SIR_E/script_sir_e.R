# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.r")
Rcpp::sourceCpp("mcmc_functions_sir_eg.cpp")


# starting options
num_initial_infectives <- 1
num_individuals <- 100
min_fs <- 60
max_fs <- 90

# epi parameters
beta <- 0.02
gamma <- 5
kappa <- 0
delta <- 5

simulate_data <- F
if(simulate_data) {
  simulated_data <- SimulateOutbreak(num_individuals,
                                     num_initial_infectives,
                                     D,
                                     beta,
                                     gamma,
                                     kappa,
                                     delta,
                                     min_fs,
                                     max_fs)
  paste("final size =",sum(simulated_data$epi_data$exposure_times != -1))
  saveRDS(simulated_data, "simulated_data.rds")
} else {
  simulated_data <- readRDS("simulated_data.rds")
}

truth <- list("beta" = beta,
              "gamma" = gamma,
              "kappa" = kappa,
              "delta" = delta,
              "infection_sum" = sum(simulated_data$infection_times
                                    [simulated_data$infection_times != -1]))





# Debug flags are
# beta, delta, gamma, kappa, augmented data
# prior parameters are beta shape, beta rate, gamma shape, gamma rate
# other parameters have uniform priors
MCMC_options <- list("initial_chain_state" = c(beta,
                                               delta,
                                               gamma,
                                               kappa),
                     "iterations" = 2000,
                     "prior_parameters" = list("beta_shape" = 1e-3,
                                               "beta_rate" = 1e-3,
                                               "gamma_shape" = 1e-3,
                                               "gamma_rate" = 1e-3,
                                               "delta_lower" = 0,
                                               "delta_higher" = 5,
                                               "kappa_lower" = 0,
                                               "kappa_higher" = 1),
                     "debug_flags" = list("beta" = 0,
                                          "gamma" = 0,
                                          "delta" = 1,
                                          "kappa" = 1,
                                          "aug_data" = 0),
                     "output_file" = "output.dat",
                     "debug_file" = "debug_file.dat",
                     "proposal_variance" = list("beta" = 0.01, 
                                                "gamma" = 0.5, 
                                                "delta" = 0.01, 
                                                "kappa" = 0.01),
                     "num_aug_updates" = 50,
                     "proposal_type" = 1)


MCMC(MCMC_options, 
     simulated_data$epi_data$infection_times, 
     simulated_data$epi_data$removal_times,
     simulated_data$epi_data$source,
     simulated_data$epi_data$x,
     simulated_data$epi_data$y)

res <- ProcessOutputMCMC("output.dat")
PlotTrace(res,truth,"Gamma proposals")









MCMC_options <- list("initial_chain_state" = c(beta,
                                               delta,
                                               gamma,
                                               kappa),
                     "iterations" = 250000,
                     "prior_parameters" = list("beta_shape" = 1e-3,
                                               "beta_rate" = 1e-3,
                                               "gamma_shape" = 1e-3,
                                               "gamma_rate" = 1e-3,
                                               "delta_lower" = 0,
                                               "delta_higher" = 5,
                                               "kappa_lower" = 0,
                                               "kappa_higher" = 1),
                     "debug_flags" = list("beta" = 0,
                                          "gamma" = 0,
                                          "delta" = 1,
                                          "kappa" = 1,
                                          "aug_data" = 0),
                     "output_file" = "output.dat",
                     "debug_file" = "debug_file.dat",
                     "proposal_variance" = list("beta" = 0.01, 
                                                "gamma" = 0.5, 
                                                "delta" = 0.01, 
                                                "kappa" = 0.01),
                     "num_aug_updates" = 50,
                     "proposal_type" = 0)


MCMC(MCMC_options, 
     simulated_data$infection_times, 
     simulated_data$removal_times,
     simulated_data$source,
     simulated_data$x,
     simulated_data$y)

res_uni <- ProcessOutputMCMC("output.dat")


MCMC_options$proposal_type <- 1
MCMC(MCMC_options, 
     simulated_data$infection_times, 
     simulated_data$removal_times,
     simulated_data$source,
     simulated_data$x,
     simulated_data$y)
res_gam <- ProcessOutputMCMC("output.dat")



PlotTrace(res_uni,truth,"Uniform proposals")
PlotTrace(res_gam,truth,"Gamma proposals")





num_iterations <- 250000



df <- data.frame("parameter" = factor(c(rep(c("beta","gamma","inf_sum"),each=num_iterations),
                                        rep(c("beta","gamma","inf_sum"),each=num_iterations))),
                 "value" = c(res_uni$beta,res_uni$gamma,res_uni$infection_sum,
                             res_gam$beta,res_gam$gamma,res_gam$infection_sum),
                 "iteration" = rep(1:num_iterations,6),
                 "type" = factor(rep(c("uniform","gamma"),each=num_iterations*3)))

library(ggplot2)
library(dplyr)
df %>% 
  filter(parameter=="beta") %>% 
  ggplot(aes(x=iteration,y=value)) + geom_line(aes(color=type))

df %>% 
  filter(parameter=="beta") %>% 
  ggplot(aes(x=value,fill=type)) + geom_histogram(aes(y=stat(count/sum(count)))) +
  ylab("density") + title("beta")

df %>% 
  filter(parameter=="beta") %>% 
  ggplot(aes(x=value,fill=type)) + geom_density(alpha=0.25, aes(y=stat(count/sum(count)))) +
  ylab("density") + title("beta")

df %>% 
  filter(parameter=="gamma") %>% 
  ggplot(aes(x=value,fill=type)) + geom_histogram(aes(y=stat(count/sum(count)))) +
  ylab("density") + title("gamma")

df %>% 
  filter(parameter=="gamma") %>% 
  ggplot(aes(x=value,fill=type)) + geom_density(alpha=0.25, aes(y=stat(count/sum(count)))) +
  ylab("density") + title("gamma")


