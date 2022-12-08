# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions_dev.r")
Rcpp::sourceCpp("mcmc_functions_dev.cpp")


# starting options
num_initial_infectives <- 1
num_individuals <- 50
min_fs <- 13
max_fs <- 20

# epi parameters
beta <- 0.003
gamma <- 0.5
kappa <- 0.01
delta <- 5

## genetic parameters
mutation_rate <- 4*10^(-6)
sequence_length <- 8176
average_import_distance <- 10
average_variant_distance <- 40


simulate_data <- T
if(simulate_data) {
  simulated_data <- SimulateOutbreak(num_individuals,
                                     num_initial_infectives,
                                     sequence_length,
                                     beta,
                                     gamma,
                                     kappa,
                                     delta,
                                     mutation_rate,
                                     average_import_distance, 
                                     average_variant_distance,
                                     min_fs,
                                     max_fs)
  saveRDS(simulated_data, "simulated_data.rds")
} else {
  simulated_data <- readRDS("simulated_data.rds")
}
paste("final size =",sum(simulated_data$epi_data$infection_times != -1))

truth <- list("beta" = beta,
              "gamma" = gamma,
              "kappa" = kappa,
              "delta" = delta,
              "lambda" = mutation_rate,
              "infection_sum" = sum(simulated_data$epi_data$infection_times
                                    [simulated_data$epi_data$infection_times != -1]))




# Debug flags are
# beta, delta, gamma, kappa, augmented data
# prior parameters are beta shape, beta rate, gamma shape, gamma rate
# other parameters have uniform priors
MCMC_options <- list("initial_chain_state" = c(beta,
                                               delta,
                                               gamma,
                                               kappa,
                                               mutation_rate),
                     "iterations" = 2000,
                     "prior_parameters" = list("beta_shape" = 1e-3,
                                               "beta_rate" = 1e-3,
                                               "gamma_shape" = 1e-3,
                                               "gamma_rate" = 1e-3,
                                               "delta_lower" = 0,
                                               "delta_higher" = 5,
                                               "kappa_lower" = 0,
                                               "kappa_higher" = 1,
                                               "lambda_lower" = 0,
                                               "lambda_higher" = 1),
                     "debug_flags" = list("beta" = 0,
                                          "gamma" = 0,
                                          "delta" = 1,
                                          "kappa" = 1,
                                          "lambda" = 0,
                                          "aug_data" = 0),
                     "output_file" = "output.dat",
                     "debug_file" = "debug_file.dat",
                     "proposal_variance" = list("beta" = 0.01, 
                                                "gamma" = 0.5, 
                                                "delta" = 0.01, 
                                                "kappa" = 0.01,
                                                "lambda" = 0.000001),
                     "num_aug_updates" = 25)

set.seed(2)
MCMC_SIR_JC(MCMC_options,
            simulated_data$gen_data$observed$N,
            simulated_data$epi_data$infection_times, 
            simulated_data$epi_data$removal_times,
            simulated_data$epi_data$source,
            simulated_data$epi_data$x,
            simulated_data$epi_data$y,
            simulated_data$gen_data$observed$genetic_ids,
            simulated_data$gen_data$observed$sample_times,
            simulated_data$gen_data$observed$subtype_numbers,
            simulated_data$gen_data$observed$gen_matrix)


res <- ProcessOutputMCMC("output.dat")
PlotTrace(res,truth)



Rcpp::sourceCpp("../SIR_EG/mcmc_functions_sir_eg.cpp")
MCMC(MCMC_options,
     simulated_data$gen_data$observed$N,
     simulated_data$epi_data$infection_times, 
     simulated_data$epi_data$removal_times,
     simulated_data$epi_data$source,
     simulated_data$epi_data$x,
     simulated_data$epi_data$y,
     simulated_data$gen_data$observed$genetic_ids,
     simulated_data$gen_data$observed$sample_times,
     simulated_data$gen_data$observed$subtype_numbers,
     simulated_data$gen_data$observed$gen_matrix)
res_old <- ProcessOutputMCMC("output.dat")
PlotTrace(res_old,truth)






