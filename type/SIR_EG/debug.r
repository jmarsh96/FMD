# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions_sir_eg.r")
Rcpp::sourceCpp("mcmc_functions_sir_eg.cpp")


# starting options
num_initial_infectives <- 1
num_individuals <- 50
min_fs <- 3
max_fs <- 50

# epi parameters
beta <- 0.002
gamma <- 0.5
kappa <- 0.02
delta <- 5

## genetic parameters
mutation_rate <- 4*10^(-6)#/1000
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
  paste("final size =",sum(simulated_data$epi_data$infection_times != -1))
  saveRDS(simulated_data, "simulated_data.rds")
} else {
  simulated_data <- readRDS("simulated_data.rds")
}

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
                     "iterations" = 5000,
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



par(mfrow=c(2,1))
target <- 2
plot(res$infection_times[,target], type="l")
abline(h=simulated_data$epi_data$infection_times[target],col=2)
plot(res_SIR_E$infection_times[,target], type="l")
abline(h=simulated_data$epi_data$infection_times[target],col=2)



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
           simulated_data$epi_data$symptom_times,
           simulated_data$epi_data$source,
           simulated_data$epi_data$x,
           simulated_data$epi_data$y)

res_SIR_E <- ProcessOutputMCMC_SIR_E("output.dat")
PlotTrace_SIR_EG(res_SIR_E,truth)


plot_individual <- F
if(plot_individual) {
  ever_infected <- which(simulated_data$epi_data$infection_times != -1)
  dev.new()
  par(mfrow=c(4,2))
  plot(res$loglik,type="l",ylab="loglik",main="Epi + Genetic data")
  plot(res$infection_sum,type="l",ylab="Infection time sum")
  abline(h=truth$infection_sum,col=2)
  for(person in ever_infected) {
    if(simulated_data$epi_data$infection_times[person] != 0) {
      plot(res$infection_times[,person],type="l",ylab="Infection time",main="Infection time trace for ID = 25")
      abline(h=simulated_data$epi_data$infection_times[person],col=2)
      plot(res$source[,person],type="l",ylab=paste0("source for ",person))
    }
  }
  
  
  
  
  dev.new()
  par(mfrow=c(4,2))
  plot(res_SIR_E$loglik,type="l",ylab="loglik",main="Epi data only")
  plot(res_SIR_E$infection_sum,type="l",ylab="Infection time sum")
  abline(h=truth$infection_sum,col=2)
  for(person in ever_infected) {
    if(simulated_data$epi_data$infection_times[person] != 0) {
      plot(res_SIR_E$infection_times[,person],type="l",ylab="Infection time",main="Infection time trace for ID = 25")
      abline(h=simulated_data$epi_data$infection_times[person],col=2)
      plot(res_SIR_E$source[,person],type="l",ylab=paste0("source for ",person))
    }
  }
} else {
  dev.new()
  par(mfrow=c(2,2))
  plot(res$loglik,type="l",ylab="loglik",main="Epi + Genetic data")
  plot(res_SIR_E$loglik,type="l",ylab="loglik",main="Epi data only")
  plot(res$infection_sum,type="l",ylab="Infection time sum")
  abline(h=truth$infection_sum,col=2)
  plot(res_SIR_E$infection_sum,type="l",ylab="Infection time sum")
  abline(h=truth$infection_sum,col=2)
}



plot_inf_trace_to_file <- T
if(plot_inf_trace_to_file) {
  pdf("inf_trace_plots.pdf")
  ever_infected <- which(simulated_data$epi_data$infection_times != -1)
  par(mfrow=c(4,2))
  for(person in ever_infected) {
    if(simulated_data$epi_data$infection_times[person] != 0) {
      plot(res$infection_times[,person],type="l",ylab="Infection time",main=paste0("EPI+GEN - Infection time trace for ID = ",person))
      abline(h=simulated_data$epi_data$infection_times[person],col=2)
      plot(res_SIR_E$infection_times[,person],type="l",ylab="Infection time",main=paste0("EPI - Infection time trace for ID = ",person))
      abline(h=simulated_data$epi_data$infection_times[person],col=2)
    }
  }
  dev.off()
}



