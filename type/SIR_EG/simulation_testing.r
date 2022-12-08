# Here we wish to investigate the performance of the forward simulation for the variant analysis
rm(list=ls())
setwd(dirname(rstudioapi::getSourceEditorContext()$path))
source("simulation_functions_sir_eg.r")


num_individuals <- 100
num_initial_infectives <- 1
transmission_kernel <- "Spatial"
beta <- 0.02
kappa <- 0.02
transmission_kernel_pars <- c(beta, kappa)
inf_period_distribution <- "Exp"
gamma <- 0.5
inf_period_pars <- c(gamma)

x <- SimulateOutbreak(num_individuals,
                      num_initial_infectives,
                      transmission_kernel,
                      transmisison_kernel_pars,
                      inf_period_distribution,
                      inf_period_pars,
                      lat_period_distribution = NA,
                      lat_period_pars = NA) 
sum(!is.na(x$infection_times))

## Gamma inf period
x <- SimulateOutbreak(num_individuals,
                      num_initial_infectives,
                      transmission_kernel,
                      transmisison_kernel_pars,
                      inf_period_distribution = "Gamma",
                      inf_period_pars = c(0.5, 0.5),
                      lat_period_distribution = NA,
                      lat_period_pars = NA) 
sum(!is.na(x$infection_times))

x <- SimulateOutbreak(num_individuals,
                      num_initial_infectives,
                      transmission_kernel = "Spatial",
                      transmisison_kernel_pars = 0.02,
                      inf_period_distribution = "Weibull",
                      inf_period_pars = c(0.5, 0.5),
                      lat_period_distribution = "Weibull",
                      lat_period_pars = c(2,2)) 
sum(!is.na(x$infection_times))





