SimulateSIR <- function(num_individuals,
                        num_initial_infectives,
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

# Simulates symptom times for when genetic sequences are taken
SimulateObservations <- function(epi_data) {
  symptom_times <- rep(NA,length(epi_data$infection_times))
  ever_infected <- which(!is.na(epi_data$infection_times))
  for(i in ever_infected) {
    inf_time <- epi_data$infection_times[i]
    rem_time <- epi_data$removal_times[i]
    #upper <- inf_time + (rem_time-inf_time)/3
    symp_time <- runif(1, min=inf_time, max=rem_time)
    #symp_time <- max(inf_time, rem_time-4)
    symptom_times[i] <- symp_time
  }
  epi_data$symptom_times <- symptom_times
  return(epi_data)
}

# Simulates genetic data
SimulateGeneticData <- function(epi_data, 
                                lambda, 
                                N,
                                average_import_distance, 
                                average_variant_distance) {

  ever_infected <- which(!is.na(epi_data$infection_times))
  genetic_ids <- c(ever_infected,ever_infected)
  sample_times <- c(epi_data$symptom_times[ever_infected],
                    epi_data$infection_times[ever_infected])
  subtype_numbers <- rep(1,length(genetic_ids))
  
  
  num_variants <- max(subtype_numbers)
  
  ## now add everyone at the time of exposure with each variant
  observed_length <- length(genetic_ids)/2
  #ever_infected <- which(epi_data$true_coltimes != -1 & epi_data$hcw_ind == 0)
  for(i in 1:num_variants) {
    for(person in ever_infected) {
      gen_loc <- which(person == genetic_ids & i == subtype_numbers & 
                         sample_times == epi_data$infection_times[person])
      if(length(gen_loc)==0) {
        genetic_ids <- c(genetic_ids, person)
        sample_times <- c(sample_times, epi_data$true_coltimes[person])
        subtype_numbers <- c(subtype_numbers, i)
      }
    }
  }
  #browser()
  data_list <- list("t_i" = epi_data$infection_times,
                    "source" = epi_data$source,
                    "genetic_ids" = genetic_ids,
                    "sample_times" = sample_times,
                    "subtype_numbers" = subtype_numbers)
  #browser()
  gen_source <- ReturnGenSourceVector_R(data_list)
  #gen_source[gen_source >= 0] <- gen_source[gen_source >= 0] + 1
  
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
        #browser()
        mutation_prob <- JC_prob_mutation(lambda, time_diff)
        num_mutations <- rbinom(1, size=N, prob=mutation_prob)
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

SimulateOutbreak <- function(num_individuals,
                             num_initial_infectives,
                             sequence_length,
                             beta,
                             gamma,
                             kappa,
                             delta,
                             lambda,
                             average_import_distance, 
                             average_variant_distance,
                             fs_min=0,
                             fs_max=0)
{
  #browser()
  if(fs_max == 0) fs_max <- num_individuals
  if(fs_min == 0) {
    epi_data <- SimulateSIR(num_individuals,
                            num_initial_infectives,
                            beta,
                            gamma,
                            kappa,
                            delta)
  } else {
    data_generated <- F
    while(!data_generated) {
      epi_data <- SimulateSIR(num_individuals,
                              num_initial_infectives,
                              beta,
                              gamma,
                              kappa,
                              delta)
      final_size <- sum(!is.na(epi_data$infection_times))
      if(final_size >= fs_min && final_size <= fs_max) data_generated <- T
    }
  }
  
  ## Simulate symptoms
  epi_data <- SimulateObservations(epi_data)
  
  ## Simulate genetic data
  #browser()
  gen_data <- SimulateGeneticData(epi_data, 
                                  lambda, 
                                  sequence_length,
                                  average_import_distance, 
                                  average_variant_distance)
  
  
  epi_data$symptom_times[is.na(epi_data$symptom_times)] <- -1
  epi_data$infection_times[is.na(epi_data$infection_times)] <- -1
  epi_data$removal_times[is.na(epi_data$removal_times)] <- -1
  
  out <- list("epi_data" = epi_data,
              "gen_data" = gen_data)
  return(out)
}

ProcessOutputMCMC <- function(filename) {
  num_parameters <- 5
  num_parameters <- num_parameters + 1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  N <- (ncol(res)-num_parameters)/2
  #browser()
  out <- list("beta" = res$V1,
              "delta" = res$V2,
              "gamma" = res$V3,
              "kappa" = res$V4,
              "lambda" = res$V5,
              "loglik" = res$V6,
              "infection_sum" = apply(res[,(num_parameters+1):(num_parameters+N)],1,
                                      function(x) sum(x[x!=-1])),
              "infection_times" = res[,(num_parameters+1):(num_parameters+N)],
              "source" = res[,((num_parameters+N)+1):ncol(res)])
  return(out)
}

CalculateSourceTable <- function(res, simulated_data) {
  df <- data.frame("id" = integer(),
                   "inf_source" = integer(),
                   "true_source" = integer(),
                   "prob" = numeric())
  ever_infected <- which(simulated_data$epi_data$infection_times != -1)
  for(i in 1:length(ever_infected)) {
    person <- ever_infected[i]
    post_source_dist <- table(res$source[,person])/length(res$source[,person])
    true_source <- simulated_data$epi_data$source[person]
    max_loc <- which.max(post_source_dist)
    post_mode <- as.numeric(names(post_source_dist)[max_loc])
    if(post_mode != -1) post_mode <- post_mode + 1
    post_prob <- post_source_dist[max_loc]
    cur_df <- data.frame("id" = person,
                         "inf_source" = post_mode,
                         "true_source" = true_source,
                         "prob" = post_prob)
    df <- rbind(df,cur_df)
  }
  rownames(df) <- 1:nrow(df)
  return(df)
}

CalculateSourceAccuracy <- function(res, simulated_data)
{
  ever_infected <- which(simulated_data$epi_data$infection_times != -1)
  correct_sources <- 0
  for(i in 1:length(ever_infected)) {
    person <- ever_infected[i]
    post_source_dist <- table(res$source[,person])/length(res$source[,person])
    true_source <- simulated_data$epi_data$source[person] - 1
    post_mode <- as.numeric(names(post_source_dist)[which.max(post_source_dist)])
    if(post_mode == true_source) {
      correct_sources <- correct_sources + 1
    }
  }
  return(correct_sources/length(ever_infected))
}

PlotTrace <- function(res, truth = NA) {
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
  plot(res$lambda, type="l", ylab="lambda")
  abline(h=truth$lambda,col=2)
  plot(res$beta/res$gamma, type="l", ylab="beta/gamma")
  abline(h=truth$beta/truth$gamma,col=2)
  plot(res$loglik, type="l", ylab="loglik")
  plot(res$infection_sum, type="l", ylab="infection_sum")
  abline(h=truth$infection_sum,col=2)
  par(mfrow=c(1,1))
}

ProcessOutputMCMC_SIR_E <- function(filename) {
  num_parameters <- 4
  num_parameters <- num_parameters + 1 # correction for the loglik
  res <- read.table(filename)
  #browser()
  N <- (ncol(res)-num_parameters)/2
  #browser()
  out <- list("beta" = res$V1,
              "delta" = res$V2,
              "gamma" = res$V3,
              "kappa" = res$V4,
              "loglik" = res$V5,
              "infection_sum" = apply(res[,(num_parameters+1):(num_parameters+N)],1,
                                      function(x) sum(x[x!=-1])),
              "infection_times" = res[,(num_parameters+1):(num_parameters+N)],
              "source" = res[,((num_parameters+N)+1):ncol(res)])
  return(out)
}

PlotTrace_SIR_EG <- function(res, truth = NA, plot_title = "") {
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
