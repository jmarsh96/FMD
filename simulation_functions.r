## Simulates an outbreak
#
# Parameters
# transmission_kernel - Spatial, Homogenous
# tranmission_kernel_pars - Spatial(beta, kappa), homogenous(beta)
# inf_period_distribution - (R_i-I_i) Exponential, Gamma, Weibull
# inf_period_pars - Exp(rate), Gamma(shape,rate), Weibull(shape,scale)
# lat_period_distribution - (I_i-E_i) Exponential, Gamma, Weibull
# lat_period_pars - Exp(rate), Gamma(shape,rate), Weibull(shape,scale)
# kernel - h_ij in the distance - Spatial kernel, homogenous
# kernel_pars - Spatial(kappa)
SimulateOutbreak <- function(num_individuals,
                             num_initial_infectives,
                             transmission_kernel,
                             transmisison_kernel_pars,
                             inf_period_distribution,
                             inf_period_pars,
                             lat_period_distribution = NA,
                             lat_period_pars = NA) {
  
  
  total_pop <- num_individuals+num_initial_infectives
  
  ## Settings for the infectious period distribution
  if(inf_period_distribution == "Exp") {
    if(length(inf_period_pars) == 1) {
      inf_rate = inf_period_pars
      rand_inf_period <- function() rexp(1, inf_rate)
    } else {
      stop("Incorrect number of infectious period parameters")
    }
  } else if(inf_period_distribution == "Gamma") {
    if(length(inf_period_pars) == 2) {
      inf_shape = inf_period_pars[1]
      inf_rate = inf_period_pars[2]
      rand_inf_period <- function() rgamma(1, shape = inf_shape, rate = inf_rate)
    } else {
      stop("Incorrect number of infectious period parameters")
    }
  } else if(inf_period_distribution == "Weibull") {
    if(length(inf_period_pars) == 2) {
      inf_shape = inf_period_pars[1]
      inf_scale = inf_period_pars[2]
      rand_inf_period <- function() rweibull(1, inf_shape, inf_scale)
    } else {
      stop("Incorrect number of infectious period parameters")
    }
  } else {
    stop("Infectious period distribution not known")
  }
  
  ## Settings for the transmission kernel
  if(transmission_kernel == "Homogenous") {
    if(length(transmission_kernel_pars) == 1) {
      trans_kernel <- function(i,j) beta
    } else {
      stop("Incorrect number of parameters for transmission kernel")
    }
  } else if(transmission_kernel == "Spatial") {
    if(length(transmission_kernel_pars) == 2) {
      trans_kernel <- function(i,j,x_coords,y_coods) {
        beta <- transmission_kernel_pars[1]
        kappa <- transmission_kernel_pars[2]
        euclidian_distance <- sqrt((x_coords[i]-x_coords[j])^2+
                                     (y_coords[i]-y_coords[j])^2)
        return(beta*exp(-kappa*euclidian_distance))
      }
    } else {
      stop("Incorrect number of parameters for transmission kernel")
    }
  }
  
  
  #browser()
  if(is.na(lat_period_distribution)) {
    ## Simulate from the SIR model
    # individual coordinates
    
    x_coords <- runif(total_pop)
    y_coords <- runif(total_pop)
    
    # indicator vectors, ith element = 1 iff they belong in that class
    S <- rep(1,total_pop)
    I <- numeric(total_pop)
    R <- numeric(total_pop)
    
    source <- rep(-1,total_pop)
    
    # matrix of transition times, 
    # first column is infection time
    # second column is removal time
    times <- matrix(NA,nrow=total_pop,ncol=2)
    
    ## Update vectors with initial conditions
    initial_infective <- sample(1:total_pop,size=num_initial_infectives)
    S[initial_infective] <- 0
    I[initial_infective] <- 1
    times[initial_infective,1] <- 0
    rem_time <- rand_inf_period()
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
            infectious_pressure[j] <- trans_kernel(i, cur_infective,
                                                   x_coords, y_coords)
          }
          
          cur_infection_time <- t + rexp(1, rate=sum(infectious_pressure))
          if(cur_infection_time < next_infection_time) {
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
        times[next_infection,2] <- next_infection_time + rand_inf_period()
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
  } else {
    ## Simulate from the SEIR model
    
    ## Settings for latent period
    if(lat_period_distribution == "Exp") {
      if(length(lat_period_pars) == 1) {
        lat_rate = lat_period_pars
        rand_lat_period <- function() rexp(1, lat_rate)
      } else {
        stop("Incorrect number of latent period parameters")
      }
    } else if(lat_period_distribution == "Gamma") {
      if(length(lat_period_pars) == 2) {
        lat_shape = lat_period_pars[1]
        lat_rate = lat_period_pars[2]
        rand_lat_period <- function() rgamma(1, shape = lat_shape, rate = lat_rate)
      } else {
        stop("Incorrect number of latent period parameters")
      }
    } else if(lat_period_distribution == "Weibull") {
      if(length(lat_period_pars) == 2) {
        lat_shape = lat_period_pars[1]
        lat_scale = lat_period_pars[2]
        rand_lat_period <- function() rweibull(1, lat_shape, lat_scale)
      } else {
        stop("Incorrect number of latent period parameters")
      }
    } else {
      stop("Latent period distribution not known")
    }
    
    ## Simulate data
    x_coords <- runif(total_pop)
    y_coords <- runif(total_pop)
    
    # indicator vectors, ith element = 1 iff they belong in that class
    S <- rep(1,total_pop)
    E <- numeric(total_pop)
    I <- numeric(total_pop)
    R <- numeric(total_pop)
    
    source <- rep(-1,total_pop)
    
    # matrix of transition times, first column in exposure time
    # second column is infection time
    # third column is removal time
    # fourth is lesion time
    times <- matrix(NA,nrow=total_pop,ncol=3)
    
    ## Update vectors with initial conditions
    initial_infective <- sample(1:total_pop,size=num_initial_infectives)
    S[initial_infective] <- 0
    E[initial_infective] <- 1
    times[initial_infective,1] <- 0
    inf_time <- rand_lat_period()
    #rem_time <- inf_time + rgamma(1,delta,gamma)
    rem_time <- inf_time + rand_inf_period()
    times[initial_infective,2] <- inf_time
    times[initial_infective,3] <- rem_time
    source[initial_infective] <- -1
    
    # record time
    t <- 0
    while((sum(I) > 0 || sum(E) > 0) && sum(S) > 0) {
      # Generate times until next event
      #browser()
      S_idx <- which(S==1)
      E_idx <- which(E==1)
      I_idx <- which(I==1)
      
      ## For each susceptible, genetic the time to infection
      next_exposure <- NA
      next_exposure_time <- Inf
      next_exposure_source <- NA
      #browser()
      currently_infective <- which(times[,2] <= t & times[,3] > t)
      for(i in S_idx) {
        if(length(currently_infective) > 0) {
          #browser()
          infectious_pressure <- numeric(length(currently_infective))
          for(j in 1:length(currently_infective)) {
            cur_infective <- currently_infective[j]
            infectious_pressure[j] <- trans_kernel(i, cur_infective,
                                                   x_coords, y_coords)
          }
          
          cur_exposure_time <- t + rexp(1, rate=sum(infectious_pressure))
          if(cur_exposure_time < next_exposure_time) {
            #browser()
            next_exposure <- i
            next_exposure_time <- cur_exposure_time
            if(length(currently_infective)==1) {
              next_exposure_source <- currently_infective
            } else {
              next_exposure_source <- sample(currently_infective, size=1,
                                             prob=infectious_pressure/
                                               sum(infectious_pressure))
            }
          }
        }
      }
      
      #browser()
      ## find the time to the next E->I transition
      next_infective <- E_idx[which.min(times[E_idx,2])]
      if(length(next_infective)==0) {
        next_infective_time <- Inf
      } else {
        next_infective_time <- times[next_infective,2]
      }
      
      
      # look at the next removal time (I->R transition)
      next_removal <- I_idx[which.min(times[I_idx,3])]
      if(length(next_removal)==0) {
        next_removal_time <- Inf
      } else {
        next_removal_time <- times[next_removal,3]
      }
      #browser()
      
      event_times <- c(next_exposure_time, 
                       next_infective_time,
                       next_removal_time)
      
      next_event <- which.min(event_times)
      
      if(next_event==1) {
        # Exposure has occured (S -> E)
        #browser()
        S[next_exposure] <- 0
        E[next_exposure] <- 1
        next_inf_time <- next_exposure_time + rand_lat_period()
        #next_rem_time <- next_inf_time + rgamma(1,delta,gamma)
        next_rem_time <- next_inf_time + rand_inf_period()
        times[next_exposure,1] <- next_exposure_time
        times[next_exposure,2] <- next_inf_time
        times[next_exposure,3] <- next_rem_time
        source[next_exposure] <- next_exposure_source
        t <- next_exposure_time
      } else if(next_event==2) {
        # An exposed individual is now infective (E->I)
        E[next_infective] <- 0
        I[next_infective] <- 1
        t <- next_infective_time
      } else {
        # Removal has occured (I -> R)
        I[next_removal] <- 0
        R[next_removal] <- 1
        t <- next_removal_time
      }
      #browser()
    }
    
    
    out <- data.frame("exposure_times" = times[,1],
                      "infection_times" = times[,2],
                      "removal_times" = times[,3],
                      "source" = source,
                      "x" = x_coords,
                      "y" = y_coords)
    return(out)
  }
}

## Simulates an outbreak with background transmission
#
# Parameters
# transmission_kernel - Spatial, Homogenous
# tranmission_kernel_pars - Spatial(beta, kappa), homogenous(beta)
# inf_period_distribution - (R_i-I_i) Exponential, Gamma, Weibull
# inf_period_pars - Exp(rate), Gamma(shape,rate), Weibull(shape,scale)
# lat_period_distribution - (I_i-E_i) Exponential, Gamma, Weibull
# lat_period_pars - Exp(rate), Gamma(shape,rate), Weibull(shape,scale)
# kernel - h_ij in the distance - Spatial kernel, homogenous
# kernel_pars - Spatial(kappa)
SimulateOutbreak_bkg <- function(num_individuals,
                                 T_end,
                             transmission_kernel,
                             transmisison_kernel_pars,
                             inf_period_distribution,
                             inf_period_pars,
                             lat_period_distribution = NA,
                             lat_period_pars = NA) {
  
  
  total_pop <- num_individuals
  
  ## Settings for the infectious period distribution
  if(inf_period_distribution == "Exp") {
    if(length(inf_period_pars) == 1) {
      inf_rate = inf_period_pars
      rand_inf_period <- function() rexp(1, inf_rate)
    } else {
      stop("Incorrect number of infectious period parameters")
    }
  } else if(inf_period_distribution == "Gamma") {
    if(length(inf_period_pars) == 2) {
      inf_shape = inf_period_pars[1]
      inf_rate = inf_period_pars[2]
      rand_inf_period <- function() rgamma(1, shape = inf_shape, rate = inf_rate)
    } else {
      stop("Incorrect number of infectious period parameters")
    }
  } else if(inf_period_distribution == "Weibull") {
    if(length(inf_period_pars) == 2) {
      inf_shape = inf_period_pars[1]
      inf_scale = inf_period_pars[2]
      rand_inf_period <- function() rweibull(1, inf_shape, inf_scale)
    } else {
      stop("Incorrect number of infectious period parameters")
    }
  } else {
    stop("Infectious period distribution not known")
  }
  
  ## Settings for the transmission kernel
  if(transmission_kernel == "Homogenous") {
    if(length(transmission_kernel_pars) == 1) {
      trans_kernel <- function(i,j) beta
    } else {
      stop("Incorrect number of parameters for transmission kernel")
    }
  } else if(transmission_kernel == "Spatial") {
    if(length(transmission_kernel_pars) == 2) {
      trans_kernel <- function(i,j,x_coords,y_coods) {
        beta <- transmission_kernel_pars[1]
        kappa <- transmission_kernel_pars[2]
        euclidian_distance <- sqrt((x_coords[i]-x_coords[j])^2+
                                     (y_coords[i]-y_coords[j])^2)
        return(beta*exp(-kappa*euclidian_distance))
      }
    } else {
      stop("Incorrect number of parameters for transmission kernel")
    }
  } else if(transmission_kernel == "Spatial_b") {
    ## spatial transmission kernel with background transmission
    if(length(transmission_kernel_pars) == 3) {
      trans_kernel <- function(i,j,x_coords,y_coods) {
        beta <- transmission_kernel_pars[1]
        beta0 <- transmission_kernel_pars[2]
        kappa <- transmission_kernel_pars[3]
        euclidian_distance <- sqrt((x_coords[i]-x_coords[j])^2+
                                     (y_coords[i]-y_coords[j])^2)
        return(beta*exp(-kappa*euclidian_distance))
      }
    }
  }
  
  
  if(is.na(lat_period_distribution)) {
    ## Simulate from the SIR model
    # individual coordinates
    
    x_coords <- runif(total_pop)
    y_coords <- runif(total_pop)
    
    # indicator vectors, ith element = 1 iff they belong in that class
    S <- rep(1,total_pop)
    I <- numeric(total_pop)
    R <- numeric(total_pop)
    
    source <- rep(-2,total_pop)
    
    # matrix of transition times, 
    # first column is infection time
    # second column is removal time
    times <- matrix(NA,nrow=total_pop,ncol=2)
    
    # record time
    t <- 0
    while(sum(S) > 0 && t < T_end) {
      # Generate times until next event
      S_idx <- which(S==1)
      I_idx <- which(I==1)
      
      ## For each susceptible, genetic the time to infection
      next_infection <- NA
      next_infection_time <- Inf
      next_infection_source <- NA
      #browser()
      currently_infective <- c(which(times[,1] <= t & times[,2] > t), -1)
      
      for(i in S_idx) {
        infectious_pressure <- numeric(length(currently_infective))
        if(length(currently_infective > 1)) {
          for(j in 1:length(currently_infective)) {
            cur_infective <- currently_infective[j]
            if(cur_infective == -1) {
              infectious_pressure[j] <- beta0
            } else {
              infectious_pressure[j] <- trans_kernel(i, cur_infective,
                                                     x_coords, y_coords)
            }

              
          }
        }
        cur_infection_time <- t + rexp(1, rate=sum(infectious_pressure))
        if(cur_infection_time < next_infection_time) {
          next_infection <- i
          next_infection_time <- cur_infection_time
          next_infection_source <- sample(currently_infective, 
                                          size=1,
                                          prob=infectious_pressure/
                                            sum(infectious_pressure))
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
      #browser()
      if(next_event==1) {
        # Infection has occurred (S -> E)
        S[next_infection] <- 0
        I[next_infection] <- 1
        times[next_infection,1] <- next_infection_time
        times[next_infection,2] <- next_infection_time + rand_inf_period()
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
  } else {
    ## Simulate from the SEIR model
    
    ## Settings for latent period
    if(lat_period_distribution == "Exp") {
      if(length(lat_period_pars) == 1) {
        lat_rate = lat_period_pars
        rand_lat_period <- function() rexp(1, lat_rate)
      } else {
        stop("Incorrect number of latent period parameters")
      }
    } else if(lat_period_distribution == "Gamma") {
      if(length(lat_period_pars) == 2) {
        lat_shape = lat_period_pars[1]
        lat_rate = lat_period_pars[2]
        rand_lat_period <- function() rgamma(1, shape = lat_shape, rate = lat_rate)
      } else {
        stop("Incorrect number of latent period parameters")
      }
    } else if(lat_period_distribution == "Weibull") {
      if(length(lat_period_pars) == 2) {
        lat_shape = lat_period_pars[1]
        lat_scale = lat_period_pars[2]
        rand_lat_period <- function() rweibull(1, lat_shape, lat_scale)
      } else {
        stop("Incorrect number of latent period parameters")
      }
    } else {
      stop("Latent period distribution not known")
    }
    
    ## Simulate data
    x_coords <- runif(total_pop)
    y_coords <- runif(total_pop)
    
    # indicator vectors, ith element = 1 iff they belong in that class
    S <- rep(1,total_pop)
    E <- numeric(total_pop)
    I <- numeric(total_pop)
    R <- numeric(total_pop)
    
    source <- rep(-1,total_pop)
    
    # matrix of transition times, first column in exposure time
    # second column is infection time
    # third column is removal time
    # fourth is lesion time
    times <- matrix(NA,nrow=total_pop,ncol=3)
    
    ## Update vectors with initial conditions
    initial_infective <- sample(1:total_pop,size=num_initial_infectives)
    S[initial_infective] <- 0
    E[initial_infective] <- 1
    times[initial_infective,1] <- 0
    inf_time <- rand_lat_period()
    #rem_time <- inf_time + rgamma(1,delta,gamma)
    rem_time <- inf_time + rand_inf_period()
    times[initial_infective,2] <- inf_time
    times[initial_infective,3] <- rem_time
    source[initial_infective] <- -1
    
    # record time
    t <- 0
    while((sum(I) > 0 || sum(E) > 0) && sum(S) > 0) {
      # Generate times until next event
      #browser()
      S_idx <- which(S==1)
      E_idx <- which(E==1)
      I_idx <- which(I==1)
      
      ## For each susceptible, genetic the time to infection
      next_exposure <- NA
      next_exposure_time <- Inf
      next_exposure_source <- NA
      #browser()
      currently_infective <- which(times[,2] <= t & times[,3] > t)
      for(i in S_idx) {
        if(length(currently_infective) > 0) {
          #browser()
          infectious_pressure <- numeric(length(currently_infective))
          for(j in 1:length(currently_infective)) {
            cur_infective <- currently_infective[j]
            infectious_pressure[j] <- trans_kernel(i, cur_infective,
                                                   x_coords, y_coords)
          }
          
          cur_exposure_time <- t + rexp(1, rate=sum(infectious_pressure))
          if(cur_exposure_time < next_exposure_time) {
            #browser()
            next_exposure <- i
            next_exposure_time <- cur_exposure_time
            if(length(currently_infective)==1) {
              next_exposure_source <- currently_infective
            } else {
              next_exposure_source <- sample(currently_infective, size=1,
                                             prob=infectious_pressure/
                                               sum(infectious_pressure))
            }
          }
        }
      }
      
      #browser()
      ## find the time to the next E->I transition
      next_infective <- E_idx[which.min(times[E_idx,2])]
      if(length(next_infective)==0) {
        next_infective_time <- Inf
      } else {
        next_infective_time <- times[next_infective,2]
      }
      
      
      # look at the next removal time (I->R transition)
      next_removal <- I_idx[which.min(times[I_idx,3])]
      if(length(next_removal)==0) {
        next_removal_time <- Inf
      } else {
        next_removal_time <- times[next_removal,3]
      }
      #browser()
      
      event_times <- c(next_exposure_time, 
                       next_infective_time,
                       next_removal_time)
      
      next_event <- which.min(event_times)
      
      if(next_event==1) {
        # Exposure has occured (S -> E)
        #browser()
        S[next_exposure] <- 0
        E[next_exposure] <- 1
        next_inf_time <- next_exposure_time + rand_lat_period()
        #next_rem_time <- next_inf_time + rgamma(1,delta,gamma)
        next_rem_time <- next_inf_time + rand_inf_period()
        times[next_exposure,1] <- next_exposure_time
        times[next_exposure,2] <- next_inf_time
        times[next_exposure,3] <- next_rem_time
        source[next_exposure] <- next_exposure_source
        t <- next_exposure_time
      } else if(next_event==2) {
        # An exposed individual is now infective (E->I)
        E[next_infective] <- 0
        I[next_infective] <- 1
        t <- next_infective_time
      } else {
        # Removal has occured (I -> R)
        I[next_removal] <- 0
        R[next_removal] <- 1
        t <- next_removal_time
      }
      #browser()
    }
    
    
    out <- data.frame("exposure_times" = times[,1],
                      "infection_times" = times[,2],
                      "removal_times" = times[,3],
                      "source" = source,
                      "x" = x_coords,
                      "y" = y_coords)
    return(out)
  }
}


SimulateObservations <- function(outbreak, observational_model) {
  if(observational_model == "Uniform") {
    symptom_times <- rep(NA,length(outbreak$infection_times))
    ever_infected <- which(!is.na(outbreak$infection_times))
    for(i in ever_infected) {
      inf_time <- outbreak$infection_times[i]
      rem_time <- outbreak$removal_times[i]
      #upper <- inf_time + (rem_time-inf_time)/3
      symp_time <- runif(1, min=inf_time, max=rem_time)
      #symp_time <- max(inf_time, rem_time-4)
      symptom_times[i] <- symp_time
    }
    outbreak$symptom_times <- symptom_times
    return(outbreak)
  } else {
    stop("Unknown observational model")
  }
}

K80_prob <- function(alpha, beta, t) {
  prob_t <- 0.25 + 0.25*exp(-4*beta*t) - 0.5*exp(-2*(alpha+beta)*t)
  prob_v <- (0.25 - 0.25*exp(-4*beta*t))*2
  out <- list("t" = prob_t,
              "v" = prob_v,
              "n" = 1-prob_t-prob_v)
  return(out)
}

JC_prob_mutation <- function(rate,t) {
  prob <- 0.75*(1-exp(-4*rate*t))
  return(prob)
}

CalculateDistanceBetweenNodes <- function(i,j,distance_matrix,gen_source) {
  if(i==j) return(0)
  path1 <- ReturnPathToRoot(i,gen_source)
  path2 <- ReturnPathToRoot(j,gen_source)
  
  ## check if they are in the same pathway
  same_pathway <- F
  if(length(path1) > length(path2)) {
    ## look for j in path1
    j_loc <- which(path1==j)
    if(length(j_loc)>0) {
      same_pathway <- T
    }
  } else {
    ## look for i in path 2
    i_loc = which(path2==i)
    if(length(i_loc)>0) {
      same_pathway <- T
    }
  }
  #browser()
  if(same_pathway) {
    ## they are in the same pathway so just add nodes
    if(length(path1)>length(path2)) {
      j_loc = which(path1 == j)
      nodes_to_sum <- path1[j_loc:length(path1)]
    } else {
      i_loc = which(path2 == i)
      nodes_to_sum <- path2[i_loc:length(path2)]
    }
    
  } else {
    ## try and find a common node
    common_node <- ReturnCommonNode(path1,path2)
    if(!is.null(common_node)) {
      ## common nodes found
      common_node_path1_loc <- which(common_node == path1)
      common_node_path2_loc <- which(common_node == path2)
      nodes_to_sum <- c(rev(path1[(common_node_path1_loc+1):length(path1)]),path2[common_node_path2_loc:length(path2)])
    } else {
      ## completely separate pathway, sum all the way to the roots and then look at the distance between roots
      nodes_to_sum <- c(rev(path1),path2)
    }
  }
  dist <- 0
  for(i in 1:(length(nodes_to_sum)-1)) {
    node_i <- nodes_to_sum[i]
    node_j <- nodes_to_sum[i+1]
    dist <- dist + distance_matrix[node_i, node_j]
  }
  return(dist)
}

ReturnPathToRoot <- function(node, gen_source) {
  #browser()
  path <- c(node)
  source <- gen_source[node] + 1
  while(source != 0) {
    path <- c(source, path)
    source <- gen_source[source] + 1
  }
  return(path)
}

ReturnCommonNode <- function(path1, path2) {
  ## here the paths are backwards
  if(length(path2) > length(path1)) {
    temp <- path1
    path1 <- path2
    path2 <- temp
  }
  common_node <- NULL
  ## path 1 should be longer
  for(i in length(path2):1) {
    ## check all the nodes in the smaller path and determine if it is in the bigger path
    current_node = path2[i]
    if(current_node %in% path1) {
      common_node <- current_node
      break
    }
  }
  return(common_node)
}

SimulateGeneticData <- function(outbreak_data,
                                N,
                                mutation_model,
                                mutation_pars) {
  
  ## do something with these later
  average_import_distance <- 10
  average_variant_distance <- 40
  #browser()
  ## Settings
  if(is.na(mutation_model)) {
    return(NA)
  }
  if(mutation_model == "JC") {
    if(length(mutation_pars) == 1) {
      rand_mutations <- function(t) {
        mutation_rate <- mutation_pars[1]
        mut_prob <- JC_prob_mutation(mutation_rate, t)
        num_mutations <- rbinom(1, prob = mut_prob, size=sequence_length)
        return(num_mutations)
      }
      
    } else {
      stop("Incorrect number of parameters for the mutation model")
    }
  } else if(mutation_model == "Kimura") {
    if(length(mutation_pars) == 2) {
      rand_mutations <- function(t) {
        transition_rate <- mutation_pars[1]
        transversion_rate <- mutation_pars[2]
        mut_prob <- K80_prob(transition_rate, transversion_rate, t)
        total_mutations <- rbinom(1, prob=mut_prob$t + mut_prob$v, size = sequence_length)
        total_transitions <- rbinom(1, prob=mut_prob$t/(mut_prob$t + mut_prob$v), size = total_mutations)
        out <- list("t" = total_transitions,
                    "v" = total_mutations-total_transitions)
        return(out)
      }
    } else {
      stop("Incorrect number of parameters for the mutation model")
    }
  } else {
    stop("Mutation model not known")
  }
  
  
  ever_infected <- which(!is.na(outbreak_data$infection_times))
  if("exposure_times" %in% colnames(outbreak_data)) {
    branching_times <- outbreak_data$exposure_times
  } else {
    branching_times <- outbreak_data$infection_times
  }
  genetic_ids <- c(ever_infected,ever_infected)
  sample_times <- c(outbreak_data$symptom_times[ever_infected],
                    branching_times[ever_infected])
  
  ## Consider implementing multiple subtypes at a later date
  subtype_numbers <- rep(1,length(genetic_ids))
  num_variants <- max(subtype_numbers)
  
  ## now add everyone at the time of exposure with each variant
  observed_length <- length(genetic_ids)/2
  #ever_infected <- which(epi_data$true_coltimes != -1 & epi_data$hcw_ind == 0)
  for(i in 1:num_variants) {
    for(person in ever_infected) {
      #browser()
      gen_loc <- which(person == genetic_ids & i == subtype_numbers & 
                         sample_times == branching_times[person])
      if(length(gen_loc)==0) {
        genetic_ids <- c(genetic_ids, person)
        sample_times <- c(sample_times, branching_times[person])
        subtype_numbers <- c(subtype_numbers, i)
      }
    }
  }
  data_list <- list("branching_times" = branching_times,
                    "t_i" = branching_times,
                    "source" = outbreak_data$source,
                    "genetic_ids" = genetic_ids,
                    "sample_times" = sample_times,
                    "subtype_numbers" = subtype_numbers)
  
  gen_source <- ReturnGenSourceVector_R(data_list)
  #gen_source[gen_source >= 0] <- gen_source[gen_source >= 0] + 1
  
  if(mutation_model == "JC") {
    full_matrix <- matrix(NA,
                          nrow=length(genetic_ids),
                          ncol=length(genetic_ids))
    
    
    #browser()
    ## Fill in direct links
    for(i in 1:length(gen_source)) {
      current_source <- gen_source[i] + 1
      full_matrix[i,i] <- 0
      if(current_source == 0) {
        ## first introduction of the pathogen in this chain, compare to master distance
        
      } else {
        time_diff <- sample_times[i] - sample_times[current_source]
        if(time_diff != 0) {
          num_mutations <- rand_mutations(time_diff)
        } else {
          num_mutations <- 0
        }
        full_matrix[current_source,i] <- num_mutations
        full_matrix[i,current_source] <- num_mutations
        
      }
    }
    
    ## Fill in imported sequences
    imports <- which(gen_source == -1)
    if(length(imports)>1) {
      for(i in 1:(length(imports)-1)) {
        for(j in (i+1):length(imports)) {
          variant_i <- variant_numbers[imports[i]]
          variant_j <- variant_numbers[imports[j]]
          if(variant_i == variant_j) {
            draw <- rpois(1,average_import_distance)
          } else {
            draw <- rpois(1,average_variant_distance)
          }
          
          #full_distance_matrix[imports[j],imports[i]] <- draw
          #full_distance_matrix[imports[i],imports[j]] <- draw
        }
      }
    }
    
    for(i in 1:nrow(full_matrix)) 
    {
      for(j in 1:ncol(full_matrix)) 
      {
        if(is.na(full_matrix[i,j])) 
        {
          distance <- CalculateDistanceBetweenNodes(i,j,full_matrix, 
                                                    gen_source)
          #print(c(i,j))
          if(is.null(distance)) return(NULL)
          full_matrix[i,j] <- distance
        }
      }
    }
    #browser()
    observed <- list("genetic_ids" = genetic_ids[1:observed_length],
                     "sample_times" = sample_times[1:observed_length],
                     "gen_matrix" = full_matrix[1:observed_length, 1:observed_length],
                     "subtype_numbers" = subtype_numbers[1:observed_length],
                     "N" =  N)
    truth <- list("genetic_ids" = genetic_ids,
                  "sample_times" = sample_times,
                  "gen_matrix" = full_matrix,
                  "subtype_numbers" = subtype_numbers,
                  "N" =  N)
    
    out <- list("observed" = observed,
                "truth" = truth)
    
    return(out)
  } else {
    full_T_matrix <- matrix(NA,
                            nrow=length(genetic_ids),
                            ncol=length(genetic_ids))
    full_V_matrix <- matrix(NA,
                            nrow=length(genetic_ids),
                            ncol=length(genetic_ids))
    
    
    #browser()
    ## Fill in direct links
    for(i in 1:length(gen_source)) {
      current_source <- gen_source[i] + 1
      full_T_matrix[i,i] <- 0
      full_V_matrix[i,i] <- 0
      if(current_source == 0) {
        ## first introduction of the pathogen in this chain, compare to master distance
        
      } else {
        time_diff <- sample_times[i] - sample_times[current_source]
        if(time_diff != 0) {
          num_mutations <- rand_mutations(time_diff)
        } else {
          num_mutations <- 0
        }
        full_T_matrix[current_source,i] <- num_mutations$t
        full_T_matrix[i,current_source] <- num_mutations$t
        full_V_matrix[current_source,i] <- num_mutations$v
        full_V_matrix[i,current_source] <- num_mutations$v
        
      }
    }
    
    ## Fill in imported sequences
    imports <- which(gen_source == -1)
    if(length(imports)>1) {
      for(i in 1:(length(imports)-1)) {
        for(j in (i+1):length(imports)) {
          variant_i <- variant_numbers[imports[i]]
          variant_j <- variant_numbers[imports[j]]
          if(variant_i == variant_j) {
            draw <- rpois(1,average_import_distance)
          } else {
            draw <- rpois(1,average_variant_distance)
          }
          
          #full_distance_matrix[imports[j],imports[i]] <- draw
          #full_distance_matrix[imports[i],imports[j]] <- draw
        }
      }
    }
    
    for(i in 1:nrow(full_T_matrix)) 
    {
      for(j in 1:ncol(full_T_matrix)) 
      {
        if(is.na(full_T_matrix[i,j])) 
        {
          distance_T <- CalculateDistanceBetweenNodes(i,j,full_T_matrix, 
                                                    gen_source)
          distance_V <- CalculateDistanceBetweenNodes(i,j,full_V_matrix, 
                                                      gen_source)
          
          #print(c(i,j))
          if(is.null(distance_T) || is.null(distance_V)) return(NULL)
          full_T_matrix[i,j] <- distance_T
          full_V_matrix[i,j] <- distance_V
        }
      }
    }
    #browser()
    observed <- list("genetic_ids" = genetic_ids[1:observed_length],
                     "sample_times" = sample_times[1:observed_length],
                     "T_matrix" = full_T_matrix[1:observed_length, 1:observed_length],
                     "V_matrix" = full_V_matrix[1:observed_length, 1:observed_length],
                     "subtype_numbers" = subtype_numbers[1:observed_length],
                     "N" =  N)
    truth <- list("genetic_ids" = genetic_ids,
                  "sample_times" = sample_times,
                  "T_matrix" = full_T_matrix,
                  "V_matrix" = full_V_matrix,
                  "subtype_numbers" = subtype_numbers,
                  "N" =  N)
    
    out <- list("observed" = observed,
                "truth" = truth)
    
    return(out)
  }
  
  
}



## Simulate data from the model
SimulateData <- function(num_individuals,
                         num_initial_infectives,
                         transmission_kernel,
                         transmisison_kernel_pars,
                         inf_period_distribution,
                         inf_period_pars,
                         lat_period_distribution,
                         lat_period_pars,
                         observation_model,
                         mutation_model,
                         mutation_pars,
                         fs_lim = NA) {
  
  ## First simulate epi data
  #browser()
  counter <- 0
  repeat {
    counter <- counter + 1
    outbreak_data <- SimulateOutbreak(num_individuals,
                                      num_initial_infectives,
                                      transmission_kernel,
                                      transmisison_kernel_pars,
                                      inf_period_distribution,
                                      inf_period_pars,
                                      lat_period_distribution,
                                      lat_period_pars)
    final_size <- sum(!is.na(outbreak_data$infection_times))
    #browser()
    if(any(is.na(fs_lim)) || (final_size >= fs_lim[1] && final_size <= fs_lim[2])) {
      break
    }
    
    if(counter > 200) {
      stop("final size limit too restrictive with parameter choices")
    }
  }
  #browser()
  ## Simulate observations
  outbreak_data <- SimulateObservations(outbreak_data, observation_model)
  #browser()
  if(is.na(mutation_model)) {
    gen_data <- NA
  } else {
    ## Simulate genetic data
    gen_data <- SimulateGeneticData(outbreak_data,
                                    sequence_length,
                                    mutation_model,
                                    mutation_pars) 
  }
  
  
  
  # Final formatting for MCMC
  if("exposure_times" %in% colnames(outbreak_data)) {
    outbreak_data$exposure_times[is.na(outbreak_data$exposure_times)] <- -1
  }
  outbreak_data$infection_times[is.na(outbreak_data$infection_times)] <- -1
  outbreak_data$removal_times[is.na(outbreak_data$removal_times)] <- -1
  outbreak_data$symptom_times[is.na(outbreak_data$symptom_times)] <- -1
  
  out <- list("epi_data" = outbreak_data,
              "gen_data" = gen_data)
  return(out)
} 



## Simulate data from the model (with background)
SimulateData_bkg <- function(num_individuals,
                             T_end,
                         transmission_kernel,
                         transmisison_kernel_pars,
                         inf_period_distribution,
                         inf_period_pars,
                         lat_period_distribution,
                         lat_period_pars,
                         observation_model,
                         mutation_model,
                         mutation_pars,
                         fs_lim = NA) {
  
  ## First simulate epi data
  #browser()
  counter <- 0
  repeat {
    counter <- counter + 1
    outbreak_data <- SimulateOutbreak_bkg(num_individuals,
                                          T_end,
                                      transmission_kernel,
                                      transmisison_kernel_pars,
                                      inf_period_distribution,
                                      inf_period_pars,
                                      lat_period_distribution,
                                      lat_period_pars)
    final_size <- sum(!is.na(outbreak_data$infection_times))
    #browser()
    if(any(is.na(fs_lim)) || (final_size >= fs_lim[1] && final_size <= fs_lim[2])) {
      break
    }
    
    if(counter > 200) {
      stop("final size limit too restrictive with parameter choices")
    }
  }
  #browser()
  ## Simulate observations
  outbreak_data <- SimulateObservations(outbreak_data, observation_model)
  #browser()
  if(is.na(mutation_model)) {
    gen_data <- NA
  } else {
    ## Simulate genetic data
    gen_data <- SimulateGeneticData(outbreak_data,
                                    sequence_length,
                                    mutation_model,
                                    mutation_pars) 
  }
  
  
  
  # Final formatting for MCMC
  if("exposure_times" %in% colnames(outbreak_data)) {
    outbreak_data$exposure_times[is.na(outbreak_data$exposure_times)] <- -1
  }
  outbreak_data$infection_times[is.na(outbreak_data$infection_times)] <- -1
  outbreak_data$removal_times[is.na(outbreak_data$removal_times)] <- -1
  outbreak_data$symptom_times[is.na(outbreak_data$symptom_times)] <- -1
  
  out <- list("epi_data" = outbreak_data,
              "gen_data" = gen_data)
  return(out)
} 

ProcessOutputMCMC_SEIR_K80 <- function(filename, parameters) {
  num_parameters <- length(parameters)
  
  ## work out how to do this better in future with different types
  num_parameters <- num_parameters -1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  if(!is.na(lat_period_distribution)) {
    N <- (ncol(res)-num_parameters)/3  
  } else {
    N <- (ncol(res)-num_parameters)/2
  }
  #browser()
  out <- lapply(1:num_parameters, function(i) res[,i])
  names(out) <- c(parameters[1:8],"loglik")
  if(!is.na(lat_period_distribution)) {
    out$exposure_sum <- apply(res[,(num_parameters+1):(N+num_parameters)], 1,
                              function(x) sum(x[x!=-1]))
    out$infection_sum <- apply(res[,(N+num_parameters+1):(2*N+num_parameters)], 1,
                               function(x) sum(x[x!=-1]))
    out$exposure_times <- res[,(num_parameters+1):(N+num_parameters)]
    out$infection_times <- res[,(N+num_parameters+1):(2*N+num_parameters)]
    out$source <- res[,(2*N+num_parameters+1):(3*N+num_parameters)]
  } else {
    ## implement functionality for this later 
    out$infection_sum <- apply(res[,(num_parameters+1):(num_parameters+N)],1,
                               function(x) sum(x[x!=-1]))
    out$infection_times <- res[,(num_parameters+1):(num_parameters+N)]
    out$source <- res[,((num_parameters+N)+1):ncol(res)]
  } 
  return(out)
}

ProcessOutputMCMC_SIR_B <- function(filename, parameters) {
  num_parameters <- length(parameters)
  
  ## work out how to do this better in future with different types
  num_parameters <- num_parameters - 1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  if(!is.na(lat_period_distribution)) {
    N <- (ncol(res)-num_parameters-1)/3  
  } else {
    N <- (ncol(res)-num_parameters-1)/2
  }
  #browser()
  out <- lapply(1:(num_parameters+1), function(i) res[,i])
  names(out) <- c(parameters[1:5],"loglik")
  if(!is.na(lat_period_distribution)) {
    out$exposure_sum <- apply(res[,(num_parameters+1):(N+num_parameters)], 1,
                              function(x) sum(x[x!=-1]))
    out$infection_sum <- apply(res[,(N+num_parameters+1):(2*N+num_parameters)], 1,
                               function(x) sum(x[x!=-1]))
    out$exposure_times <- res[,(num_parameters+1):(N+num_parameters)]
    out$infection_times <- res[,(N+num_parameters+1):(2*N+num_parameters)]
    out$source <- res[,(2*N+num_parameters+1):(3*N+num_parameters)]
  } else {
    ## implement functionality for this later 
    out$infection_sum <- apply(res[,(num_parameters+2):(num_parameters+1+N)],1,
                               function(x) sum(x[x!=-1]))
    out$infection_times <- res[,(num_parameters+2):(num_parameters+1+N)]
    out$source <- res[,((num_parameters+N)+1+1):ncol(res)]
  } 
  return(out)
}

ProcessOutputMCMC <- function(filename, parameters) {
  num_parameters <- length(parameters)
  
  ## work out how to do this better in future with different types
  num_parameters <- num_parameters -1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  if(!is.na(lat_period_distribution)) {
    N <- (ncol(res)-num_parameters)/3  
  } else {
    N <- (ncol(res)-num_parameters)/2
  }
  #browser()
  out <- lapply(1:num_parameters, function(i) res[,i])
  names(out) <- c(parameters[1:7],"loglik")
  if(!is.na(lat_period_distribution)) {
    out$exposure_sum <- apply(res[,(num_parameters+1):(N+num_parameters)], 1,
                              function(x) sum(x[x!=-1]))
    out$infection_sum <- apply(res[,(N+num_parameters+1):(2*N+num_parameters)], 1,
                               function(x) sum(x[x!=-1]))
    out$exposure_times <- res[,(num_parameters+1):(N+num_parameters)]
    out$infection_times <- res[,(N+num_parameters+1):(2*N+num_parameters)]
    out$source <- res[,(2*N+num_parameters+1):(3*N+num_parameters)]
  } else {
    ## implement functionality for this later 
    out$infection_sum <- apply(res[,(num_parameters+1):(num_parameters+N)],1,
                               function(x) sum(x[x!=-1]))
    out$infection_times <- res[,(num_parameters+1):(num_parameters+N)]
    out$source <- res[,((num_parameters+N)+1):ncol(res)]
  } 
  return(out)
}


ProcessOutputMCMC_SEIR <- function(filename, parameters) {
  num_parameters <- length(parameters)
  
  ## work out how to do this better in future with different types
  num_parameters <- num_parameters -1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  if(!is.na(lat_period_distribution)) {
    N <- (ncol(res)-num_parameters)/3  
  } else {
    N <- (ncol(res)-num_parameters)/2
  }
  #browser()
  out <- lapply(1:num_parameters, function(i) res[,i])
  names(out) <- c(parameters[1:6],"loglik")
  if(!is.na(lat_period_distribution)) {
    out$exposure_sum <- apply(res[,(num_parameters+1):(N+num_parameters)], 1,
                              function(x) sum(x[x!=-1]))
    out$infection_sum <- apply(res[,(N+num_parameters+1):(2*N+num_parameters)], 1,
                               function(x) sum(x[x!=-1]))
    out$exposure_times <- res[,(num_parameters+1):(N+num_parameters)]
    out$infection_times <- res[,(N+num_parameters+1):(2*N+num_parameters)]
    out$source <- res[,(2*N+num_parameters+1):(3*N+num_parameters)]
  } else {
    ## implement functionality for this later 
    out$infection_sum <- apply(res[,(num_parameters+1):(num_parameters+N)],1,
                               function(x) sum(x[x!=-1]))
    out$infection_times <- res[,(num_parameters+1):(num_parameters+N)]
    out$source <- res[,((num_parameters+N)+1):ncol(res)]
  } 
  return(out)
}

ProcessOutputMCMC_SIR <- function(filename, parameters) {
  num_parameters <- length(parameters)
  
  ## work out how to do this better in future with different types
  num_parameters <- num_parameters  # correction for the loglik
  res <- read.table(filename)
  #browser()
  if(!is.na(lat_period_distribution)) {
    N <- (ncol(res)-num_parameters)/3  
  } else {
    N <- (ncol(res)-num_parameters)/2
  }
  #browser()
  out <- lapply(1:num_parameters, function(i) res[,i])
  names(out) <- c(parameters[1:5],"loglik")
  if(!is.na(lat_period_distribution)) {
    out$exposure_sum <- apply(res[,(num_parameters+1):(N+num_parameters)], 1,
                              function(x) sum(x[x!=-1]))
    out$infection_sum <- apply(res[,(N+num_parameters+1):(2*N+num_parameters)], 1,
                               function(x) sum(x[x!=-1]))
    out$exposure_times <- res[,(num_parameters+1):(N+num_parameters)]
    out$infection_times <- res[,(N+num_parameters+1):(2*N+num_parameters)]
    out$source <- res[,(2*N+num_parameters+1):(3*N+num_parameters)]
  } else {
    ## implement functionality for this later 
    out$infection_sum <- apply(res[,(num_parameters+1):(num_parameters+N)],1,
                               function(x) sum(x[x!=-1]))
    out$infection_times <- res[,(num_parameters+1):(num_parameters+N)]
    out$source <- res[,((num_parameters+N)+1):ncol(res)]
  } 
  return(out)
}

PlotTrace <- function(res, par_names, truth) {
  par(mfrow=c(3,4))
  for(par in par_names) {
    res_loc <- which(names(res) == par)
    truth_loc <- which(names(truth) == par)
    plot(res[[res_loc]], type="l", ylab=par)
    abline(h=truth[truth_loc],col="red")
  }
  
  plot(res$loglik, type="l", ylab="loglik")
  

  
  par(mfrow=c(1,1))
}

PlotOutbreak <- function(simulated_data) {
  ever_infected <- which(simulated_data$epi_data$exposure_times != -1)
  ever_infected <- ever_infected[order(simulated_data$epi_data$exposure_times[ever_infected])]
  exp_times <- simulated_data$epi_data$exposure_times
  inf_times <- simulated_data$epi_data$infection_times
  rem_times <- simulated_data$epi_data$removal_times
  source <- simulated_data$epi_data$source
  
  plot(0,type="n",ylim=c(0,length(ever_infected)),
       xlim=c(0,max(simulated_data$epi_data$removal_times[ever_infected])),
       yaxt="n", ylab="Person")
  for(i in 1:length(ever_infected)) {
    person <- ever_infected[i]
    segments(exp_times[person],i,inf_times[person],i,col="blue")
    segments(inf_times[person],i,rem_times[person],i,col="green")
    cur_source <- source[person]
    if(cur_source != -1) {
      source_loc <- which(cur_source == ever_infected)
      arrows(exp_times[person],source_loc,exp_times[person],i,col="red")
    }
  }
  
  gen_ids <- simulated_data$gen_data$observed$genetic_ids
  sample_times <- simulated_data$gen_data$observed$sample_times
  for(i in 1:length(gen_ids)) {
    cur_id <- gen_ids[i]
    loc <- which(ever_infected == cur_id)
    cur_time <- sample_times[i]
    text(cur_time, loc, i-1)
  }
  axis(2, at=1:length(ever_infected), labels=ever_infected-1)
}

## Return a list where each entry contains the marginal posterior distribution for the ith
## ever_infected
CalculateSourceList <- function(res, simulated_data) {
  ever_infected <- which(simulated_data$epi_data$exposure_times != -1)
  source_list <- vector('list',length(ever_infected))
  for(i in 1:length(ever_infected)) {
    post_source_distribution <- table(res$source[,ever_infected[i]])/sum(table(res$source[,ever_infected[i]]))
    source_list[[i]] <- post_source_distribution
  } 
  return(source_list)
}

ReturnSourceTable <- function(res, simulated_data) {
  ever_infected <- which(simulated_data$epi_data$exposure_times != -1)
  source_data <- data.frame("id" = ever_infected,
                            "inferred_source" = NA,
                            "posterior_probability" = NA,
                            "true_source" = simulated_data$epi_data$source[ever_infected],
                            "match" = NA)
  source_list <- CalculateSourceList(res, simulated_data)
  for(i in 1:length(ever_infected)) {
    post_source_dist <- source_list[[i]]
    max_loc <- which.max(post_source_dist)
    post_prob <- post_source_dist[max_loc]
    inferred_source <- as.numeric(names(post_source_dist)[max_loc])
    if(inferred_source != -1) inferred_source <- inferred_source + 1
    source_data$inferred_source[i] <- inferred_source
    source_data$posterior_probability[i] <- post_prob
    source_data$match[i] <- ifelse(source_data$inferred_source[i] == source_data$true_source[i], 1, 0)
  }
  return(source_data)
}

CalculateSourceInformation <- function(res, simulated_data) {
  source_table <- ReturnSourceTable(res, simulated_data)
  accuracy <- mean(source_table$match)
  out <- list("source_table" = source_table,
              "accuracy" = accuracy)
  return(out)
}
