INIT <- function(
  data, 
  vars,  
  dayvar, 
  beepvar,
  idvar, 
  network_type = c("saturated", "pruned"),
  homogeneity_test = c("overall", "contemporaneous", "temporal"),
  save_models = FALSE,
  alpha = 0.05,
  ...
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
  
  # # # Add estimator if missing:
  # if(missing(estimator)){
  #   estimator <- "ML"
  # }
  # 
  # # Add network_type if missing:
  # if(missing(network_type)){
  #   network_type <- "saturated"
  # }
  # 
  # # Add homogeneity_test if missing:
  # if(missing(homogeneity_test)){
  #   homogeneity_test <- "overall"
  # }
  # 
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
  # estimator <- match.arg(estimator)
  # if (!(estimator %in% c("ML", "FIML"))) {
  #   stop("estimator invalid; needs to be 'ML' or 'FIML'")
  # }
  
  network_type <- match.arg(network_type)
  # if (!(network_type %in% c("saturated", "pruned"))) {
  #   stop("network_type invalid; needs to be 'saturated' or 'pruned'")
  # }
  # 
  
  homogeneity_test <- match.arg(homogeneity_test)
  # if (!(homogeneity_test %in% c("overall", "contemporaneous", "temporal"))) {
  #   stop("homogeneity_test invalid; needs to be 'overall', 'contemporaneous', 'temporal'")
  # }
  
  # ------ Collect info -----
  n_person <- length(unique(id)) # number of individuals
  n_vars <- length(vars) # number of variables
  
  # ------ Specify network model using psychonetrics -----
  mod <- psychonetrics::gvar(data, 
                             vars = vars,
                             beepvar = beepvar,
                             dayvar = dayvar,
                             groups = idvar, 
                             ...) 
  # if(estimator == "ML"){
  #   
  #   mod <- psychonetrics::gvar(data, 
  #                              vars = vars,
  #                              beepvar = beepvar,
  #                              dayvar = dayvar,
  #                              groups = idvar, 
  #                              estimator = "ML") 
  #   
  # } else { # if estimator == "FIML"
  #   
  #   mod <- psychonetrics::gvar(data, 
  #                              vars = vars,
  #                              beepvar = beepvar,
  #                              dayvar = dayvar,
  #                              groups = idvar, 
  #                              estimator = "FIML")
  # } # end: estimator
  
  # ------ Estimate network model using psychonetrics -----
  
  mod <- mod %>% psychonetrics::runmodel()
  
  # ------ Test for homogeneity -----
  
  # FIXME: SE: Hieronder zit veel code repititie. Dat kan problematisch zijn als je wilt debuggen. Ik heb dit nu 
  # een beetje omgegooid zodat het minder code repititie gebruikt.
  
  # Prune the network if needed first:
  if (network_type == "pruned"){
    
    # Prune network model:
    mod <- mod %>% psychonetrics::prune(alpha = alpha, recursive = FALSE)  %>% 
      psychonetrics::unionmodel() %>% 
      psychonetrics::runmodel()
  
  }

  # Define homogeneity model:
  if (homogeneity_test == "overall"){
    
    mod_constrained <- mod %>% 
      psychonetrics::groupequal("omega_zeta") %>% 
      psychonetrics::groupequal("beta") %>% runmodel
    
  } else if (homogeneity_test == "contemporaneous"){
    
    mod_constrained <- mod %>% 
      psychonetrics::groupequal("omega_zeta") %>% runmodel
    
  } else {
    
    mod_constrained <- mod %>% 
      psychonetrics::groupequal("beta") %>% runmodel
    
  }
  
  # Get model fit indices: 
  fit_indicies <- psychonetrics::compare(different = mod, 
                                         equal = mod_constrained)
  
  # AIC:
  mod_AIC <- fit_indicies$AIC[fit_indicies$model == "different"]
  mod_constrained_AIC <- fit_indicies$AIC[fit_indicies$model == "equal"] 
  delta_AIC <- mod_AIC - mod_constrained_AIC
  
  # Results:
  res <- list(different_network_AIC = mod_AIC, 
              equal_network_AIC = mod_constrained_AIC, 
              delta_AIC = delta_AIC)
  
  # Save model: 
  if(save_models){
    
    if(mod_AIC < mod_constrained_AIC){
      
      # get matrix from psychonetrics:
      mod_contemp <- mod %>% 
        psychonetrics::getmatrix("omega_zeta")
      mod_temp <- mod %>% 
        psychonetrics::getmatrix("beta")
      
      # Change names to correspond to idvar:
      names(mod_contemp) <- unique(data[[idvar]])
      names(mod_temp) <- unique(data[[idvar]])
      
      res[["model"]] <- list(
        contemporaneous = mod_contemp,
        temporal = mod_temp)
      
    } else # if mod_AIC$AIC > mod_constrained_AIC$AIC
      
      # get matrix from psychonetrics:
      mod_constrained_contemp <- mod_constrained %>% 
        psychonetrics::getmatrix("omega_zeta")
    
    mod_constrained_temp <- mod_constrained %>% 
      psychonetrics::getmatrix("beta")
    
    # Change names to correspond to idvar:
    names(mod_constrained_contemp) <- unique(data[[idvar]])
    names(mod_constrained_temp) <- unique(data[[idvar]])
    
    res[["model"]] <- list(
      contemporaneous = mod_constrained_contemp,
      temporal = mod_constrained_temp)
    
  } # End: save_model
  
  
  # Some checks:
  if (mod_AIC < -1e+10 |  mod_constrained_AIC < -1e+10 | abs(delta_AIC) > 1e08){
    warning("POSSIBLE OPTIMIZATION ERRORS DETECTED: Check the results manually by setting save_models = FALSE.")
  }
  
  # 
  # if(network_type == "saturated"){
  #   
  #   if(homogeneity_test == "overall"){
  #     
  #     # Constrain  model: 
  #     mod_constrained <- mod %>% 
  #       psychonetrics::groupequal("omega_zeta") %>% 
  #       psychonetrics::groupequal("beta") %>% 
  #       psychonetrics::runmodel()
  #     
  #     # Get model fit indices: 
  #     fit_indicies <- psychonetrics::compare(different = mod, 
  #                                            equal = mod_constrained)
  #     
  #     # AIC:
  #     mod_AIC <- fit_indicies %>% 
  #       filter(model == "different")
  #     
  #     mod_constrained_AIC <- fit_indicies %>% 
  #       filter(model == "equal")
  #     
  #     delta_AIC <- mod_AIC$AIC - mod_constrained_AIC$AIC
  #     
  #     # Results:
  #     res <- list(different_network_AIC = mod_AIC$AIC, 
  #                 equal_network_AIC = mod_constrained_AIC$AIC, 
  #                 delta_AIC = delta_AIC)
  #     
  #     # Save model: 
  #     if(save_models == TRUE){
  #       
  #       if(mod_AIC$AIC < mod_constrained_AIC$AIC){
  #         
  #         # get matrix from psychonetrics:
  #         mod_contemp <- mod %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #         mod_temp <- mod %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_contemp) <- unique(data[[idvar]])
  #         names(mod_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_contemp,
  #           temporal = mod_temp)
  #         
  #       } else # if mod_AIC$AIC > mod_constrained_AIC$AIC
  #         
  #         # get matrix from psychonetrics:
  #         mod_constrained_contemp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #       
  #       mod_constrained_temp <- mod_constrained %>% 
  #         psychonetrics::getmatrix("beta")
  #       
  #       # Change names to correspond to idvar:
  #       names(mod_constrained_contemp) <- unique(data[[idvar]])
  #       names(mod_constrained_temp) <- unique(data[[idvar]])
  #       
  #       res[["model"]] <- list(
  #         contemporaneous = mod_constrained_contemp,
  #         temporal = mod_constrained_temp)
  #       
  #     } # End: save_model
  #     
  #   } else if(homogeneity_test == "contemporaneous") {
  #     
  #     # Constrain  model: 
  #     mod_constrained <- mod %>% 
  #       psychonetrics::groupequal("omega_zeta") %>% 
  #       psychonetrics::runmodel()
  #     
  #     # Get model fit indices: 
  #     fit_indicies <- psychonetrics::compare(different = mod, 
  #                                            equal = mod_constrained)
  #     
  #     # AIC:
  #     mod_AIC <- fit_indicies %>% 
  #       filter(model == "different")
  #     
  #     mod_constrained_AIC <- fit_indicies %>% 
  #       filter(model == "equal")
  #     
  #     delta_AIC <- mod_AIC$AIC - mod_constrained_AIC$AIC
  #     
  #     # Results:
  #     res <- list(different_network_AIC = mod_AIC$AIC, 
  #                 equal_network_AIC = mod_constrained_AIC$AIC, 
  #                 delta_AIC = delta_AIC)
  #     
  #     # Save model: 
  #     if(save_models == TRUE){
  #       
  #       if(mod_AIC$AIC < mod_constrained_AIC$AIC){
  #         
  #         # get matrix from psychonetrics:
  #         mod_contemp <- mod %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #         mod_temp <- mod %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_contemp) <- unique(data[[idvar]])
  #         names(mod_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_contemp,
  #           temporal = mod_temp)
  #         
  #       } else # if mod_AIC$AIC > mod_constrained_AIC$AIC
  #         
  #         # get matrix from psychonetrics:
  #         mod_constrained_contemp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #       
  #       mod_constrained_temp <- mod_constrained %>% 
  #         psychonetrics::getmatrix("beta")
  #       
  #       # Change names to correspond to idvar:
  #       names(mod_constrained_contemp) <- unique(data[[idvar]])
  #       names(mod_constrained_temp) <- unique(data[[idvar]])
  #       
  #       res[["model"]] <- list(
  #         contemporaneous = mod_constrained_contemp,
  #         temporal = mod_constrained_temp)
  #       
  #     } # End: save_model
  #     
  #   } else { # if homogeneity_test == "homogeneity_temp"
  #     
  #     # Constrain  model: 
  #     mod_constrained <- mod %>% 
  #       psychonetrics::groupequal("beta") %>% 
  #       psychonetrics::runmodel()
  #     
  #     # Get model fit indices: 
  #     fit_indicies <- psychonetrics::compare(different = mod, 
  #                                            equal = mod_constrained)
  #     
  #     # AIC:
  #     mod_AIC <- fit_indicies %>% 
  #       filter(model == "different")
  #     
  #     mod_constrained_AIC <- fit_indicies %>% 
  #       filter(model == "equal")
  #     
  #     delta_AIC <- mod_AIC$AIC - mod_constrained_AIC$AIC
  #     
  #     # Results:
  #     res <- list(different_network_AIC = mod_AIC$AIC, 
  #                 equal_network_AIC = mod_constrained_AIC$AIC, 
  #                 delta_AIC = delta_AIC)
  #     
  #     # Save model: 
  #     if(save_models == TRUE){
  #       
  #       if(mod_AIC$AIC < mod_constrained_AIC$AIC){
  #         
  #         # get matrix from psychonetrics:
  #         mod_contemp <- mod %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #         mod_temp <- mod %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_contemp) <- unique(data[[idvar]])
  #         names(mod_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_contemp,
  #           temporal = mod_temp)
  #         
  #       } else # if mod_AIC$AIC > mod_constrained_AIC$AIC
  #         
  #         # get matrix from psychonetrics:
  #         mod_constrained_contemp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #       
  #        mod_constrained_temp <- mod_constrained %>% 
  #          psychonetrics::getmatrix("beta")
  #       
  #        # Change names to correspond to idvar:
  #        names(mod_constrained_contemp) <- unique(data[[idvar]])
  #        names(mod_constrained_temp) <- unique(data[[idvar]])
  #       
  #        res[["model"]] <- list(
  #          contemporaneous = mod_constrained_contemp,
  #          temporal = mod_constrained_temp)
  #       
  #     } # End: save_model
  #     
  #   } # End: homogeneity_test
  #   
  # } else { # if network_type == "pruned"
  #   
  #   # Prune network model:
  #   mod_pruned <- mod %>% 
  #     psychonetrics::prune(alpha = 0.05, recursive = FALSE) 
  #   
  #   # Create union model: 
  #   mod_union <- mod_pruned %>% 
  #     psychonetrics::unionmodel() %>% 
  #     psychonetrics::runmodel()
  #   
  #   if(homogeneity_test == "overall"){
  #     
  #     # Constrain union model: 
  #     mod_constrained <- mod_union %>% 
  #       psychonetrics::groupequal("omega_zeta") %>% 
  #       psychonetrics::groupequal("beta") %>% 
  #       runmodel()
  #     
  #     # Get fit indices: 
  #     fit_indicies <- psychonetrics::compare(different = mod_pruned, 
  #                                            equal = mod_constrained)
  #     
  #     # BIC
  #     mod_BIC <- fit_indicies %>% 
  #       filter(model == "different")
  #   
  #     mod_constrained_BIC <- fit_indicies %>% 
  #       filter(model == "equal")
  #     
  #     delta_BIC <- mod_BIC$BIC - mod_constrained_BIC$BIC
  #     
  #     # Results
  #     res <- list(different_network_BIC = mod_BIC$BIC, 
  #                 equal_network_BIC = mod_constrained_BIC$BIC, 
  #                 delta_BIC = delta_BIC)
  #     
  #     # Save model: 
  #     if(save_models == TRUE){
  #       
  #       if(mod_BIC$BIC < mod_constrained_BIC$BIC){
  #         
  #         # get matrix from psychonetrics:
  #         mod_contemp <- mod %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #         mod_temp <- mod %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_contemp) <- unique(data[[idvar]])
  #         names(mod_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_contemp,
  #           temporal = mod_temp)
  #         
  #       } else # if mod_BIC$BIC > mod_constrained_BIC$BIC
  #         
  #         # get matrix from psychonetrics:
  #         mod_constrained_contemp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #       
  #         mod_constrained_temp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_constrained_contemp) <- unique(data[[idvar]])
  #         names(mod_constrained_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_constrained_contemp,
  #           temporal = mod_constrained_temp)
  #         
  #     } # End: save_model
  #     
  #   } else if(homogeneity_test == "homogeneity_contemp"){
  #     
  #     # Constrain contemporaneous union model: 
  #     mod_constrained <- mod_union %>% 
  #       psychonetrics::groupequal("omega_zeta") %>% 
  #       runmodel()
  #     
  #     # Get fit indices: 
  #     fit_indicies <- psychonetrics::compare(different = mod_pruned, 
  #                                            equal = mod_constrained)
  #   
  #     # BIC:
  #     mod_BIC <- fit_indicies %>% 
  #       filter(model == "different")
  #     
  #     mod_constrained_BIC <- fit_indicies %>% 
  #       filter(model == "equal")
  #     
  #     delta_BIC <- mod_BIC$BIC - mod_constrained_BIC$BIC
  #     
  #     # Results:
  #     res <- list(different_network_BIC = mod_BIC$BIC, 
  #                 equal_network_BIC = mod_constrained_BIC$BIC, 
  #                 delta_BIC = delta_BIC)
  #     
  #     # Save model: 
  #     if(save_models == TRUE){
  #       
  #       if(mod_BIC$BIC < mod_constrained_BIC$BIC){
  #         
  #         # get matrix from psychonetrics:
  #         mod_contemp <- mod %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #         mod_temp <- mod %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_contemp) <- unique(data[[idvar]])
  #         names(mod_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_contemp,
  #           temporal = mod_temp)
  #         
  #       } else # if mod_BIC$BIC > mod_constrained_BIC$BIC
  #         
  #         # get matrix from psychonetrics:
  #         mod_constrained_contemp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #       
  #       mod_constrained_temp <- mod_constrained %>% 
  #         psychonetrics::getmatrix("beta")
  #       
  #       # Change names to correspond to idvar:
  #       names(mod_constrained_contemp) <- unique(data[[idvar]])
  #       names(mod_constrained_temp) <- unique(data[[idvar]])
  #       
  #       res[["model"]] <- list(
  #         contemporaneous = mod_constrained_contemp,
  #         temporal = mod_constrained_temp)
  #       
  #     } # End: save_model
  #     
  #   } else { # if homogeneity_test == "homogeneity_temp"
  #     
  #     # Constrain temporal union model: 
  #     mod_constrained <- mod_union %>% 
  #       psychonetrics::groupequal("beta") %>% 
  #       runmodel()
  #     
  #     # Get fit indices: 
  #     fit_indicies <- psychonetrics::compare(different = mod_pruned, 
  #                                            equal = mod_constrained)
  #     
  #     # BIC:
  #     mod_BIC <- fit_indicies %>% 
  #       filter(model == "different")
  #     
  #     mod_constrained_BIC <- fit_indicies %>% 
  #       filter(model == "equal")
  #     
  #     delta_BIC <- mod_BIC$BIC - mod_constrained_BIC$BIC
  #     
  #     # Results:
  #     res <- list(different_network_BIC = mod_BIC$BIC, 
  #                 equal_network_BIC = mod_constrained_BIC$BIC, 
  #                 delta_BIC = delta_BIC)
  #     
  #     # Save model: 
  #     if(save_models == TRUE){
  #       
  #       if(mod_BIC$BIC < mod_constrained_BIC$BIC){
  #         
  #         # get matrix from psychonetrics:
  #         mod_contemp <- mod %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #         mod_temp <- mod %>% 
  #           psychonetrics::getmatrix("beta")
  #         
  #         # Change names to correspond to idvar:
  #         names(mod_contemp) <- unique(data[[idvar]])
  #         names(mod_temp) <- unique(data[[idvar]])
  #         
  #         res[["model"]] <- list(
  #           contemporaneous = mod_contemp,
  #           temporal = mod_temp)
  #         
  #       } else # if mod_BIC$BIC > mod_constrained_BIC$BIC
  #         
  #         # get matrix from psychonetrics:
  #         mod_constrained_contemp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("omega_zeta")
  #       
  #         mod_constrained_temp <- mod_constrained %>% 
  #           psychonetrics::getmatrix("beta")
  #       
  #         # Change names to correspond to idvar:
  #         names(mod_constrained_contemp) <- unique(data[[idvar]])
  #         names(mod_constrained_temp) <- unique(data[[idvar]])
  #       
  #         res[["model"]] <- list(
  #           contemporaneous = mod_constrained_contemp,
  #           temporal = mod_constrained_temp)
  #       
  #     } # End: save_model
  #     
  #   } # End: homogeneity_test 
  #   
  # } # End: network_type 
  
  # ------ Output -----  
  
  # Add input to results:
  res[['input']] <- list(
    vars = vars, 
    network_type = network_type,
    homogeneity_test = homogeneity_test,
    psychonetrics_args = list(...))
  
  # Assign a class for S3 methods:
  class(res) <- c("INIT","list")
  
  return(res)
  
}


