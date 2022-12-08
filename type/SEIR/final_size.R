# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("../../simulation_functions.r")
Rcpp::sourceCpp("mcmc_functions_dev.cpp")

## Simulation settings
num_individuals <- 50
num_initial_infectives <- 1
transmission_kernel <- "Spatial"
inf_period_distribution <- "Gamma"
lat_period_distribution <- "Gamma"
obs_model <- "Uniform"
mut_model <- "JC"

## Parameter values
alpha <- 0.5
beta <- 0.004
gamma <- 0.5
kappa <- 0.02
delta <- 5
zeta <- 5
lambda <- 4*10^(-6)/10000
fs_lim <- NA


# Set pars for simulation
transmission_kernel_pars <- c(beta, kappa)
inf_period_pars <- c(delta, gamma)
lat_period_pars <- c(zeta, alpha)
mutation_pars <- c(lambda)
sequence_length <- 8000


max_iter <- 100
fs_vec <- numeric(max_iter)
for(i in 1:max_iter) {
  simulated_data <- SimulateData(num_individuals,
                                 num_initial_infectives,
                                 transmission_kernel,
                                 transmisison_kernel_pars,
                                 inf_period_distribution,
                                 inf_period_pars,
                                 lat_period_distribution,
                                 lat_period_pars,
                                 obs_model,
                                 mut_model,
                                 mutation_pars,
                                 fs_lim)
  fs_vec[i] <- sum(simulated_data$epi_data$infection_times != -1)
}
hist(fs_vec)