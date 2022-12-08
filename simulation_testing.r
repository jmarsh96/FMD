# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions.r")


num_individuals <- 100
num_initial_infectives <- 1
transmission_kernel <- "Spatial"
beta <- 0.02
kappa <- 0.02
transmission_kernel_pars <- c(beta, kappa)
inf_period_distribution <- "Exp"
gamma <- 0.5
inf_period_pars <- c(gamma)
obs_model <- "Uniform"
mut_model <- "JC"
mutation_rate <- 4.4*10^(-6)
mutation_pars <- c(mutation_rate)

simulated_data <- SimulateData(num_individuals,
                               num_initial_infectives,
                               transmission_kernel,
                               transmisison_kernel_pars,
                               inf_period_distribution,
                               inf_period_pars,
                               lat_period_distribution = NA,
                               lat_period_pars = NA,
                               obs_model,
                               mut_model,
                               mutation_pars,
                               fs_lim = NA) 



