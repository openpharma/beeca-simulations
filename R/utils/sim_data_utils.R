
library(dplyr)
library(clustermq)


## Simulate covariates and stratification variables 
sim_covars <- function(df, scenario, N=NULL) {
  # simulate covariates
  if (is.null(N)) 
    N <- scenario$N
  for (i in scenario$covars$name) {
    this_covar <- scenario$covars[scenario$covars$name == i,]
    if (this_covar$dist == "rnorm") {
      rdist_func <- get(this_covar$dist)
      rdist_mean <- this_covar$arg1
      rdist_sd <- this_covar$arg2
      df[[i]] <- rdist_func(N, rdist_mean, rdist_sd)
    } else if (this_covar$dist == "rbinom") {
      rdist_func <- get(this_covar$dist)
      rdist_size <- this_covar$arg1
      rdist_prob <- this_covar$arg2
      df[[i]] <- rdist_func(N, rdist_size, rdist_prob)
    } else {
      stop(sprintf("unknown dist func %s", this_covar$dist))
    }
  }
  
  # derive strata variables if needed
  if (!is.null(scenario$make_strata)) {
    for (i in scenario$make_strata$name) {
      this_strata <- scenario$make_strata[scenario$make_strata$name == i,]
      
      # if creating binary strata, label as 0/1
      if (length(this_strata$breaks) == 1) {
        df[[i]] <- cut(df[[this_strata$covar]], c(-Inf, this_strata$breaks, Inf), labels = F)-1
      } else {
        # TODO: if > 2 strata, create additional dummy variables
        stop("Attempting to create >2 strata levels. Not yet supported.")
      }
      
    }
    
  }
  
  return(df)
}


## Simulate treatment assignment 
sim_trt <- function(df, scenario) {
  N <- scenario$N
  if (is.null(scenario$stratified_randomization) || (scenario$stratified_randomization==FALSE)) {
    
    df[["trt"]] <- rbinom(N, 1, scenario$trt_assignment)
    
  } else {
    
    df[["trt"]] <- NA
    
    # define all strata if multiple strata vars
    if (length(scenario$strata) > 1) {
      all_strata <- apply(as.data.frame(df)[, scenario$strata], 
                          1, paste, collapse = "-" )
    } else {
      all_strata <- df[[scenario$strata]]
    }
    df[["all_strata"]] <- all_strata
    
    for (stratum in unique(all_strata)) {
      n_stratum <- sum(all_strata == stratum)
      
      trt_stratum <- rbinom(n_stratum, 1, scenario$trt_assignment)
      
      # assign stratum trt to dataframe
      df[["trt"]][which(all_strata == stratum)] <- trt_stratum
    }
    
  }
  
  return(df)
}


# Simulate complete trial dataset: covariates, treatment assignment and outcome
sim_dat <- function(scenario, base_data=NULL) {
  
  N <- scenario$N
  
  # if no base data supplied, simulate new data from scratch
  if (is.null(base_data)) {
    
    df <- list()
    
    # simulate covariates
    df <- sim_covars(df, scenario)
    
    # simulate treatment assignment
    df <- sim_trt(df, scenario)
    
    # convert to dataframe
    df <- as.data.frame(df)
    
  } else {
    df <- base_data
  }
  
  ## if no trt in base data, simulate treatment assignment
  if (!"trt" %in% colnames(df)) {
    df <- sim_trt(df, scenario)
  }
  
  # if no outcome in base data, simulate y conditional on baseline df
  if (!"y" %in% colnames(df)) {
    # simulate outcome
    df[["y"]] <- rbinom(N, 1, plogis(scenario$response_func(df)))
  }
  
  # otherwise, if x, trt and y already in base_data, no new data simulated
  
  return(df)
  
}


get_true_diff <- function(scenario, 
                          estimand = c("PATE", "CPATE", "SATE"),
                          base_data=NULL){
  
  if (estimand == "PATE") {
    N <- 1e6
    
    df <- list()
    
    # simulate population-level covariates
    df <- sim_covars(df, scenario, N=N)
    
    df <- as.data.frame(df)
    
    p_1 <- mean(plogis(scenario$response_func(df %>% dplyr::mutate(trt=1)))) 
    p_0 <- mean(plogis(scenario$response_func(df %>% dplyr::mutate(trt=0))))
    
  } else if (estimand == "CPATE") {
    # use base_data x and predict counterfactuals
    p_1 <- mean(plogis(scenario$response_func(base_data %>% dplyr::mutate(trt=1)))) 
    p_0 <- mean(plogis(scenario$response_func(base_data %>% dplyr::mutate(trt=0))))
    
  } else if (estimand == "SATE") {
    # TODO
    stop("SATE not implemented")
  }
  
  data.frame(scenario = scenario$scenario,
             true_diff = p_1 - p_0)
}



## Main simulation function.
## Simulate a single dataset based on scenario specification,
## fit working model and estimate variance with each method
one_result <- function(scenario, base_data=NULL){
  
  dat <- sim_dat(scenario, base_data)
  
  dat$trt <- factor(dat$trt)
  fit <- glm(as.formula(scenario$working_model), family = "binomial", data=dat)
  
  data.frame(point_estimate = beeca_wrapper(fit)$marginal_est,
             
             # Ye method without stratified randomization adjustment
             ye_beeca_var = beeca_wrapper(fit, "Ye")$marginal_var,
             
             # Ye method with stratified randomization adjustment
             ye_beeca_strata_var = beeca_wrapper(fit, "Ye", strata=scenario$strata)$marginal_var,

             # RobinCar implementation
             #ye_robincar_var = robincar_wrapper(fit),
             
             # Ge method
             ge_var = beeca_wrapper(fit, "Ge")$marginal_var,
             
             # Ye method original implementation (based on paper - deviates from RobinCar development)
             ye_beeca_mod_var = beeca_wrapper(fit, "Ye", mod=T)$marginal_var,
             
             # Ge method with HC0 sandwich variance 
             ge_hc0_var = beeca_wrapper(fit, "Ge", type="HC0")$marginal_var,
             
             
             true_diff_pate = get_true_diff(scenario, "PATE")$true_diff,
             true_diff_cpate = get_true_diff(scenario, "CPATE", dat)$true_diff,
             #true_diff_sate = get_true_diff(scenario, "SATE", dat)$true_diff,
             
             scenario = scenario[["scenario"]])
  
}



## Simulation wrapper for clustermq
run_sim_wrapper <- function(scenario, base_data=NULL){
  
  if(!exists('worker_setup_complete', mode='logical')) {
    cat("Calling setup function:\n")
    
    source("R/utils/sim_variance_estimators.R")
    source("R/utils/sim_data_utils.R")
    
    worker_setup_complete <<- TRUE
  } else {
    cat("Calling GC:\n")
    print(gc())
  }
  
  one_result(scenario, base_data)
  
}
