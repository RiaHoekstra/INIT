#' 
#' Function to test for homogeneity between idiographic network models
#'
#' @param data Data set.
#' @param vars String indicating the variables included in the model. 
#' @param beepvar String indicating the beep. 
#' @param dayvar String indicating the day.
#' @param idvar String indicating the subject id variable name.
#' @param estimator The estimator to be used. "ML" for maximum likelihood estimation 
#' or "FIML", for full-information maximum likelihood estimation. By default set to ML
#' @param network_type Type of network to estimate. Currently implemented are "saturated" or "pruned". 
#' By default set to "saturated". When "pruned" is specified, non significant edges are pruned from the model with an alpha rate of 0.05.
#' @param homogeneity_test Type of homogeneity test. "homogeneity_overall", refers to 
#' test both contemporaneous and temporal network structures for homogeneity, 
#' "homogeneity_contemporaneous" to test contemporaneous network structures for homogeneity, 
#' "homogeneity_temporal" to test temporal network structures for homogeneity.
#' By default set to "homogeneity_overall".
#' @param save_models If TRUE, all models will be saved. By default set to FALSE to save memory.
#' 
#' @importFrom dplyr %>%
#' @importFrom psychonetrics runmodel
#' 
INIT <- function(
  data, 
  vars,  
  dayvar, 
  beepvar,
  idvar, 
  estimator = "ML",
  network_type = "saturated",
  homogeneity_test = "overall",
  save_models = FALSE
  ){
  
  # ------ Set defaults -----
  
  # Add id if missing:
  if (missing(idvar)) {
    
    idvar <- "id"
    data[[idvar]] <- 1
    
  } else if (!is.character(idvar) || length(idvar) != 1 || !idvar %in% names(data)){
    
    stop("'idvar' must be a string indicating the column name of the data.")
    
  } else {
    
    id <- data[,idvar]
    id <- as.numeric(as.factor(id)) # get ID's
    
  } 
  
  # Add day if missing:
  if (missing(dayvar)){
    
    dayvar <- "day"
    data[[dayvar]] <- 1
    
  } else if (!is.character(dayvar) || length(dayvar) != 1 || !dayvar %in% names(data)){
    
    stop("'dayvar' must be a string indicating a column name of the data.")
    
  } 
  
  # Add beep if missing:
  if (missing(beepvar)){
    
    beepvar <- "beep"
    data[[beepvar]] <- ave(seq_len(nrow(data)), 
                           data[[idvar]], 
                           data[[dayvar]], 
                           FUN = seq_along)
    
  } else if (!is.character(beepvar) || length(beepvar) != 1 || !beepvar %in% names(data)){
    
    stop("'beepvar' must be a string indicating a column name of the data.")
    
  } 
  
  # Add estimator if missing:
  if(missing(estimator)){
    estimator <- "ML"
  }
  
  # Add network_type if missing:
  if(missing(network_type)){
    network_type <- "saturated"
  }
  
  # Add homogeneity_test if missing:
  if(missing(homogeneity_test)){
    homogeneity_test <- "homogeneity_overall"
  }
  
  # ------ Input Checks -----
  
  # Check if data is a data frame
  if (!(is.data.frame(data) || is.matrix(data))){
    stop("'data' argument must be a data frame")
  }
  
  # If data is a matrix coerce to data frame:
  if (is.matrix(data)){
    data <- as.data.frame(data)
  }
  
  # Check dayvar:
  if (!is.numeric(data[[dayvar]])){
    stop("'dayvar' is not numeric")
  }
  
  # Check beepvar:
  if (!is.numeric(data[[beepvar]])){
    stop("'beepvar' is not numeric")
  }
  
  # Check argument inputs:
  if (!(estimator %in% c("ML", "FIML"))) {
    stop("estimator invalid; needs to be 'ML' or 'FIML'")
  }
  
  if (!(network_type %in% c("saturated", "pruned"))) {
    stop("network_type invalid; needs to be 'saturated' or 'pruned'")
  }
  
  if (!(homogeneity_test %in% c("homogeneity_overall", "homogeneity_contemporaneous", "homogeneity_temporal"))) {
    stop("homogeneity_test invalid; needs to be 'homogeneity_overall', 'homogeneity_contemporaneous', 'homogeneity_temporal'")
  }
  
  # ------ Collect info -----
  n_person <- length(unique(id)) # number of individuals
  n_vars <- length(vars) # number of variables
  
  # ------ Specify network model using psychonetrics -----
  
  if(estimator == "ML"){
    
    mod <- psychonetrics::gvar(data, 
                               vars = vars,
                               beepvar = beepvar,
                               dayvar = dayvar,
                               groups = idvar, 
                               estimator = "ML") 
    
  } else { # if estimator == "FIML"
    
    mod <- psychonetrics::gvar(data, 
                               vars = vars,
                               beepvar = beepvar,
                               dayvar = dayvar,
                               groups = idvar, 
                               estimator = "FIML")
  } # end: estimator
  
  # ------ Estimate network model using psychonetrics -----
  
  mod <- mod %>% psychonetrics::runmodel()
  
  # ------ Test for homogeneity -----
  
  if(network_type == "saturated"){
    
    if(homogeneity_test == "homogeneity_overall"){
      
      # Constrain  model: 
      mod_constrained <- mod %>% 
        psychonetrics::groupequal("omega_zeta") %>% 
        psychonetrics::groupequal("beta") %>% 
        psychonetrics::runmodel()
      
      # Get model fit indices: 
      fit_indicies <- psychonetrics::compare(different = mod, 
                                             equal = mod_constrained)
      
      # AIC
      mod_AIC <- fit_indicies$AIC[1]
      mod_constrained_AIC <- fit_indicies$AIC[2]
      delta_AIC <- fit_indicies$AIC[1] - fit_indicies$AIC[2]
     
      # Results 
      res <- list(different_AIC = mod_AIC, 
                  equal_AIC = mod_constrained_AIC, 
                  delta_AIC = delta_AIC)
      
    } else if(homogeneity_test == "homogeneity_contemp") {
      
      # Constrain  model: 
      mod_constrained <- mod %>% 
        psychonetrics::groupequal("omega_zeta") %>% 
        psychonetrics::runmodel()
      
      # Get model fit indices: 
      fit_indicies <- psychonetrics::compare(different = mod, 
                                             equal = mod_constrained)
      
      # AIC
      mod_AIC <- fit_indicies$AIC[1]
      mod_constrained_AIC <- fit_indicies$AIC[2]
      delta_AIC <- fit_indicies$AIC[1] - fit_indicies$AIC[2]
      
      # Results 
      res <- list(different_AIC = mod_AIC, 
                  equal_AIC = mod_constrained_AIC, 
                  delta_AIC = delta_AIC)
      
    } else { # if homogeneity_test == "homogeneity_temp"
      
      # Constrain  model: 
      mod_constrained <- mod %>% 
        psychonetrics::groupequal("beta") %>% 
        psychonetrics::runmodel()
      
      # Get model fit indices: 
      fit_indicies <- psychonetrics::compare(different = mod, 
                                             equal = mod_constrained)
      
      # AIC
      mod_AIC <- fit_indicies$AIC[1]
      mod_constrained_AIC <- fit_indicies$AIC[2]
      delta_AIC <- fit_indicies$AIC[1] - fit_indicies$AIC[2]
      
      # Results 
      res <- list(different_AIC = mod_AIC, 
                  equal_AIC = mod_constrained_AIC, 
                  delta_AIC = delta_AIC)
      
    } # End: homogeneity_test
    
  } else { # if network_type == "pruned"
    
    # Prune network model:
    mod_pruned <- mod %>% 
      psychonetrics::prune(alpha = 0.05, recursive = FALSE) 
    
    # Create union model: 
    mod_union <- mod_pruned %>% 
      psychonetrics::unionmodel() %>% 
      psychonetrics::runmodel()
    
    if(homogeneity_test == "homogeneity_overall"){
      
      # Constrain union model: 
      mod_constrained <- mod_union %>% 
        psychonetrics::groupequal("omega_zeta") %>% 
        psychonetrics::groupequal("beta") %>% 
        runmodel()
      
      # Get fit indices: 
      fit_indicies <- psychonetrics::compare(different = mod_union, 
                                             equal = mod_constrained)
      
      # BIC
      mod_BIC <- fit_indicies$BIC[1]
      mod_constrained_BIC <- fit_indicies$BIC[2]
      delta_BIC <- fit_indicies$BIC[1] - fit_indicies$BIC[2]
      
      # Results 
      res <- list(different_BIC = mod_BIC, 
                  equal_BIC = mod_constrained_BIC, 
                  delta_BIC = delta_BIC)
      
    } else if(homogeneity_test == "homogeneity_contemp"){
      
      # Constrain contemporaneous union model: 
      mod_constrained <- mod_union %>% 
        psychonetrics::groupequal("omega_zeta") %>% 
        runmodel()
      
      # Get fit indices: 
      fit_indicies <- psychonetrics::compare(different = mod_union, 
                                             equal = mod_constrained)
      
      # BIC
      mod_BIC <- fit_indicies$BIC[1]
      mod_constrained_BIC <- fit_indicies$BIC[2]
      delta_BIC <- fit_indicies$BIC[1] - fit_indicies$BIC[2]
      
      # Results 
      res <- list(different_BIC = mod_BIC, 
                  equal_BIC = mod_constrained_BIC, 
                  delta_BIC = delta_BIC)
      
    } else { # if homogeneity_test == "homogeneity_temp"
      
      # Constrain temporal union model: 
      mod_constrained <- mod_union %>% 
        psychonetrics::groupequal("beta") %>% 
        runmodel()
      
      # Get fit indices: 
      fit_indicies <- psychonetrics::compare(different = mod_union, 
                                             equal = mod_constrained)
      
      # BIC
      mod_BIC <- fit_indicies$BIC[1]
      mod_constrained_BIC <- fit_indicies$BIC[2]
      delta_BIC <- fit_indicies$BIC[1] - fit_indicies$BIC[2]
      
      # Results 
      res <- list(different_BIC = mod_BIC, 
                  equal_BIC = mod_constrained_BIC, 
                  delta_BIC = delta_BIC)
      
    } # End: homogeneity_test 
    
  } # End: network_type 
    
  
  # ------ Output -----
  
  # Add input to results:
  res[['input']] <- list(
    vars = vars, 
    estimator = estimator,
    network_type = network_type,
    homogeneity_test = homogeneity_test
  )
  
  return(res)
  
}


