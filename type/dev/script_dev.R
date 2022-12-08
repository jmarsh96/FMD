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
mut_model <- "Kimura"

## Parameter values
alpha <- 5
beta <- 0.004
gamma <- 5
kappa <- 0.02
delta <- 25
zeta <- 25
transition_rate <- 0.00003
transversion_rate <- 0.000001
fs_lim <- c(4,100)


# Set pars for simulation
transmission_kernel_pars <- c(beta, kappa)
inf_period_pars <- c(delta, gamma)
lat_period_pars <- c(zeta, alpha)
mutation_pars <- c(transition_rate, transversion_rate)
sequence_length <- 8000


simulate_data <- F
if(simulate_data) {
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
  print(paste0("final size = ", sum(simulated_data$epi_data$exposure_times != -1)))
  saveRDS(simulated_data, "simulated_data.rds")
} else {
  if(file.exists("simulated_data.rds")) {
    simulated_data <- readRDS("simulated_data.rds")
  }
}




exp_times <- simulated_data$epi_data$exposure_times
inf_times <- simulated_data$epi_data$infection_times
exp_sum <- sum(exp_times[exp_times != -1])
inf_sum <- sum(inf_times[inf_times != -1])

truth <- list("alpha" = alpha,
              "beta" = beta,
              "gamma" = gamma,
              "kappa" = kappa,
              "delta" = delta,
              "zeta" = zeta,
              "transition_rate" = transition_rate,
              "transversion_rate" = transversion_rate,
              "exposure_sum" = exp_sum,
              "infection_sum" = inf_sum)
par_names <- names(truth)

MCMC_options <- list("initial_chain_state" = as.numeric(truth[1:8]),
                     "iterations" = 1000,
                     "prior_parameters" = list("alpha_shape" = 1e-3,
                                               "alpha_rate" = 1e-3,
                                               "beta_shape" = 1e-3,
                                               "beta_rate" = 1e-3,
                                               "gamma_shape" = 1e-3,
                                               "gamma_rate" = 1e-3,
                                               "kappa_lower" = 0,
                                               "kappa_higher" = 1,
                                               "delta_lower" = 0,
                                               "delta_higher" = 10,
                                               "zeta_lower" = 0,
                                               "zeta_higher" = 10,
                                               "transition_rate_lower" = 0,
                                               "transition_rate_higher" = 1,
                                               "transversion_rate_lower" = 0,
                                               "transversion_rate_higher" = 1),
                     "debug_flags" = list("alpha" = 0,
                                          "beta" = 0,
                                          "gamma" = 0,
                                          "kappa" = 1,
                                          "delta" = 1,
                                          "zeta" = 1,
                                          "transition_rate" = 0,
                                          "transversion_rate" = 0,
                                          "aug_data" = 0),
                     "output_file" = "output.dat",
                     "debug_file" = "debug_file.dat",
                     "proposal_variance" = list("kappa" = 0.05, 
                                                "delta" = 0.5,
                                                "zeta" = 0.5,
                                                "transition_rate" = 0.00001,
                                                "transversion_rate" = 0.000005),
                     "num_aug_updates" = 50)

set.seed(2)
MCMC_SEIR_K80(MCMC_options,
              simulated_data$gen_data$observed$N,
              simulated_data$epi_data$exposure_times,
              simulated_data$epi_data$infection_times, 
              simulated_data$epi_data$removal_times,
              simulated_data$epi_data$source,
              simulated_data$epi_data$x,
              simulated_data$epi_data$y,
              simulated_data$gen_data$observed$genetic_ids,
              simulated_data$gen_data$observed$sample_times,
              simulated_data$gen_data$observed$subtype_numbers,
              simulated_data$gen_data$observed$T_matrix,
              simulated_data$gen_data$observed$V_matrix)


res <- ProcessOutputMCMC_SEIR_K80("output.dat", par_names) 
PlotTrace(res, par_names, truth) 













#PlotOutbreak(simulated_data)


