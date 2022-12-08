SimulateSIR <- function(num_individuals,
                        num_initial_infectives,
                        D,
                        beta,
                        gamma,
                        kappa,
                        delta) {
  # individual coordinates
  x_coords <- runif(num_individuals)
  y_coords <- runif(num_individuals)
  
  # indicator vectors, ith element = 1 iff they belong in that class
  S <- rep(1,num_individuals)
  I <- numeric(num_individuals)
  R <- numeric(num_individuals)
  
  source <- rep(-1,num_individuals)
  
  # matrix of transition times, 
  # first column is infection time
  # second column is removal time
  times <- matrix(NA,nrow=num_individuals,ncol=2)
  
  ## Update vectors with initial conditions
  initial_infective <- sample(1:num_individuals,size=num_initial_infectives)
  S[initial_infective] <- 0
  I[initial_infective] <- 1
  times[initial_infective,1] <- 0
  rem_time <- rgamma(1,shape=delta,rate=gamma)
  times[initial_infective,2] <- rem_time
  source[initial_infective] <- -1
  
  # record time
  t <- 0
  while(sum(I) > 0 && sum(S) > 0) {
    # Generate times until next event
    S_idx <- which(S==1)
    I_idx <- which(I==1)
    

    ## For each susceptible, genetic the time to infection
    next_infection <- NA
    next_infection_time <- Inf
    next_infection_source <- NA
    #browser()
    currently_infective <- which(times[,1] <= t & times[,2] > t)
    for(i in S_idx) {
      if(length(currently_infective) > 0) {
        #browser()
        infectious_pressure <- numeric(length(currently_infective))
        for(j in 1:length(currently_infective)) {
          cur_infective <- currently_infective[j]
          euclidian_dist <- sqrt((x_coords[cur_infective]-x_coords[i])^2+
                                   (y_coords[cur_infective]-y_coords[i])^2)
          infectious_pressure[j] <- beta*exp(-kappa*euclidian_dist)
        }
        
        cur_infection_time <- t + rexp(1, rate=sum(infectious_pressure))
        if(cur_infection_time < next_infection_time) {
          #browser()
          next_infection <- i
          next_infection_time <- cur_infection_time
          if(length(currently_infective)==1) {
            next_infection_source <- currently_infective
          } else {
            next_infection_source <- sample(currently_infective, 
                                            size=1,
                                            prob=infectious_pressure/
                                              sum(infectious_pressure))
          }
        }
      }
    }
    
    
    #browser()
    # look at the next removal time (I->R transition)
    next_removal <- I_idx[which.min(times[I_idx,2])]
    if(length(next_removal)==0) {
      next_removal_time <- Inf
    } else {
      next_removal_time <- times[next_removal,2]
    }
    #browser()
    
    event_times <- c(next_infection_time,
                     next_removal_time)
    
    next_event <- which.min(event_times)
    
    if(next_event==1) {
      # Infection has occurred (S -> E)
      S[next_infection] <- 0
      I[next_infection] <- 1
      times[next_infection,1] <- next_infection_time
      times[next_infection,2] <- next_infection_time + rgamma(1,shape=delta,rate=gamma)
      source[next_infection] <- next_infection_source
      t <- next_infection_time
    } else if(next_event==2) {
      # An exposed individual is now infective (E->I)
      I[next_removal] <- 0
      R[next_removal] <- 1
      t <- next_removal_time
    } 
  }
  out <- data.frame("infection_times" = times[,1],
                    "removal_times" = times[,2],
                    "source" = source,
                    "x" = x_coords,
                    "y" = y_coords)
  return(out)
}

SimulateOutbreak <- function(num_individuals,
                             num_initial_infectives,
                             D,
                             beta,
                             gamma,
                             kappa,
                             delta,
                             fs_min=0,
                             fs_max=0)
{
  if(fs_min == 0) {
    epi_data <- SimulateSIR(num_individuals,
                            num_initial_infectives,
                            D,
                            beta,
                            gamma,
                            kappa,
                            delta)
  } else {
    data_generated <- F
    while(!data_generated) {
      epi_data <- SimulateSIR(num_individuals,
                              num_initial_infectives,
                              D,
                              beta,
                              gamma,
                              kappa,
                              delta)
      final_size <- sum(!is.na(epi_data$infection_times))
      if(final_size >= fs_min && final_size <= fs_max) data_generated <- T
    }
  }
  
  
  epi_data$infection_times[is.na(epi_data$infection_times)] <- -1
  epi_data$removal_times[is.na(epi_data$removal_times)] <- -1
  return(epi_data)
}

ProcessOutputMCMC_SIR_E <- function(filename) {
  num_parameters <- 4
  num_parameters <- num_parameters + 1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  N <- (ncol(res)-num_parameters)/2
  browser()
  out <- list("beta" = res$V1,
              "delta" = res$V2,
              "gamma" = res$V3,
              "kappa" = res$V4,
              "loglik" = res$V5,
              "infection_sum" = apply(res[,(num_parameters+1):(num_parameters+N)],1,
                                     function(x) sum(x[x!=-1])),
              "source" = res[,((num_parameters+N)+1):ncol(res)])
  return(out)
}

PlotTrace <- function(res, truth = NA, plot_title = "") {
  #dev.new()
  par(mfrow=c(3,3))
  plot(res$beta, type="l", ylab="beta")
  abline(h=truth$beta,col=2)
  plot(res$gamma, type="l", ylab="gamma")
  abline(h=truth$gamma,col=2)
  plot(res$kappa, type="l", ylab="kappa")
  abline(h=truth$kappa,col=2)
  plot(res$delta, type="l", ylab="delta")
  abline(h=truth$delta,col=2)
  plot(res$loglik, type="l", ylab="loglik")
  plot(res$beta/res$gamma,type="l", ylab="beta/gamma")
  abline(h=truth$beta/truth$gamma,col=2)
  plot(res$infection_sum, type="l", ylab="infection_sum", main=plot_title)
  abline(h=truth$infection_sum,col=2)
  par(mfrow=c(1,1))
}


