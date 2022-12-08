# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions_sir_eg.r")
Rcpp::sourceCpp("mcmc_functions_sir_eg.cpp")


# starting options
num_initial_infectives <- 1
num_individuals <- 50
min_fs <- 4
max_fs <- 50

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
                     "iterations" = 1000,
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
                                                "lambda" = 0.000003),
                     "num_aug_updates" = 50)

set.seed(1)
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


res <- ProcessOutputMCMC("output.dat")
PlotTrace(res,truth)





source_table <- CalculateSourceTable(res, simulated_data)
source_table
source_accuracy <- mean(source_table$inf_source==source_table$true_source)
source_accuracy


ever_infected <- which(simulated_data$epi_data$infection_times!=-1)
inf_mat <- res$infection_times[,ever_infected]
par(mfrow=c(3,3))
for(i in 1:length(ever_infected)) {
  person <- ever_infected[i]
  inf_time <- round(simulated_data$epi_data$infection_times[person],3)
  rem_time <- round(simulated_data$epi_data$removal_times[person],3)
  if(simulated_data$epi_data$infection_times[person] > 0) {
    plot(inf_mat[,i],type="l",ylab=paste0("Infection time trace for ",person),
         main=paste0("inf time = ",inf_time,
                     ", rem time = ",rem_time))
    abline(h=simulated_data$epi_data$infection_times[person],col=2)
  }
}
plot(res$source[,25],type="l")
plot(res$source[,30],type="l")






Rcpp::sourceCpp("../SIR_E/mcmc_functions_sir_e.cpp")
MCMC_options_SIR_E <- list("initial_chain_state" = c(beta,
                                               delta,
                                               gamma,
                                               kappa),
                     "iterations" = 5000,
                     "prior_parameters" = list("beta_shape" = 1e-3,
                                               "beta_rate" = 1e-3,
                                               "gamma_shape" = 1e-3,
                                               "gamma_rate" = 1e-3,
                                               "delta_lower" = 0,
                                               "delta_higher" = 5,
                                               "kappa_lower" = 0,
                                               "kappa_higher" = 1),
                     "debug_flags" = list("beta" = 1,
                                          "gamma" = 1,
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
MCMC_SIR_E(MCMC_options_SIR_E,
           simulated_data$epi_data$infection_times, 
           simulated_data$epi_data$removal_times,
           simulated_data$epi_data$source,
           simulated_data$epi_data$x,
           simulated_data$epi_data$y)

res_SIR_E <- ProcessOutputMCMC_SIR_E("output.dat")
PlotTrace_SIR_EG(res_SIR_E, truth)

ever_infected <- which(simulated_data$epi_data$infection_times!=-1)
inf_mat <- res_SIR_E$infection_times[,ever_infected]
par(mfrow=c(3,3))
for(i in 1:length(ever_infected)) {
  person <- ever_infected[i]
  inf_time <- round(simulated_data$epi_data$infection_times[person],3)
  rem_time <- round(simulated_data$epi_data$removal_times[person],3)
  if(simulated_data$epi_data$infection_times[person] > 0) {
    plot(inf_mat[,i],type="l",ylab=paste0("Infection time trace for ",person),
         main=paste0("inf time = ",inf_time,
                     ", rem time = ",rem_time))
    abline(h=simulated_data$epi_data$infection_times[person],col=2)
  }
}
plot(res_SIR_E$source[,25],type="l")
plot(res_SIR_E$source[,30],type="l")