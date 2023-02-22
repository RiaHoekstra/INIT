#' 
#' Function to test homogeneity in idiograhic network models
#'
#' @param data dataset 
#' @param vars string indicating the variables included in the model 
#' @param beepvar string indicating the beep  
#' @param dayvar string indicating the day
#' @param idvar string indicating the subject id variable name
#' @param estimator type of estimator, i.e., "ML" or "FIML". By default set to ML
#' @param test type of homogeneity test, i.e., "homogeneity", "homogeneity_contemp", "homogeneity_temp". By default set to homogeneity
#' @param networkType type of network, i.e., saturated or pruned network model. By default set to saturated
#' @param saveModels if TRUE, all models will be saved. By default set to FALSE to save memory
#' @param progressBar if TRUE progress bar is included 
#' 
#' @return
#' @export
#'
#' @examples
#' 
INIT <- function(
  data, 
  vars,  
  beepvar = NULL,   
  dayvar = NULL, 
  idvar, 
  estimator = "ML",
  test = "homogeneity", # different test options 
  networkType = "saturated", # type of network estimated 
  saveModels = FALSE,
  progressbar = TRUE
  ){

  # ------ To do add more input Checks -----
  
  # check modelType string
  if (!(networkType %in% c("pruned", "saturated"))) {
    stop("networkType invalid; needs to be 'pruned' or 'saturated'")
  }
  
  # check estimator string
  if (!(estimator %in% c("ML", "FIML"))) {
    stop("estimator is invalid; needs to be 'ML' or 'FIML'")
  }
  
  
  # ------ Collect basic info -----
  id <- data[,idvar]
  id <- as.numeric(as.factor(idvar)) # get ID's
  n_person <- length(unique(idvar)) # number of individuals
  n_vars <- length(vars) # number of variables
  
  
  # Progress bar
  if(progressbar == TRUE) pb <- txtProgressBar(min = 0, max=it, initial = 0, char="-", style = 3)
  # if (progressbar==TRUE) pb <- txtProgressBar(max=it, style = 3)
  

  # to do add  dayvar & beepvar
  
  # ------ Specify network model -----
  
  if(estimator == "FIML"){
    
    mod <- psychonetrics::gvar(data, 
                               vars = vars,
                               #beepvar = beepvar,
                               #dayvar = dayvar,
                               groups = idvar, 
                               estimator = "FIML") 
  } else {
    
    mod <- psychonetrics::gvar(data, 
                               vars = vars,
                               #beepvar = beepvar,
                               #dayvar = dayvar,
                               groups = idvar, 
                               estimator = "ML")
  } # end if: estimator is specified

  
  # ------ Estimate network model -----
  
  mod <- mod %>% psychonetrics::runmodel()
  
  if(networkType == "pruned"){
    
    # Prune network model:
    mod_pruned <- mod %>% 
      psychonetrics::prune(alpha = 0.05, recursive = FALSE) 
    
    # Create union model: 
    mod_union <- mod_p %>% 
      psychonetrics::unionmodel() %>% 
      psychonetrics::runmodel()
    
    if(test == "homogeneity_contemp"){
      
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
      
    } else if(test == "homogeneity_temp"){
      
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
      
    } else {
      
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
      
    } # End: test 
    
  } else {
    
    if(test == "homogeneity_contemp"){
      
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
      
    } else if(test == "homogeneity_temp") {
      
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
      
    } else {
      
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
      
    } # End: test
    
  } # End: networkType 
    
  
  return(list(homogeneity = homogeneity, 
              homogeneity_contemp = homogeneity_contemp,
              homogeneity_union = homogeneity_union,
              homogeneity_union_contemp = homogeneity_union_contemp)
         )
  
}


