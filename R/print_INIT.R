#' 
#' Print function for INIT objects
#'
#' @param x A INIT object
#' @param ... 
#'
#' @return
#' @export
#'
print.INIT <- function(x,...){
  
  name <- deparse(substitute(x))[[1]]
  
  cat("\n ----- INIT results -----")
  
  cat("\n")
  
  if(x$input$network_type == "saturated"){
    cat("\n The model", ifelse(x$delta_AIC < 0,
                               "without equality constraints", 
                           "with equality constraints"),
        "is preferred.", "\n")
    
    cat("\n")

    cat(" INIT summary \n",
        "  - AIC homogeneity model:", x$homogeneity_model_AIC,"\n",
        "  - AIC heterogeneity model:", x$heterogeneity_model_AIC,"\n",
        "  - delta AIC:", x$delta_AIC,"\n")
    
  } else {
    
    cat(        "\n The model", ifelse(x$delta_BIC < 0,
                                      "without equality constraints",
                                      "with equality constraints"), 
                "is preferred.", "\n")
    
    cat("\n")
    
    cat(" INIT summary \n",
        "  - BIC homogeneity model:", x$homogeneity_model_BIC,"\n",
        "  - BIC heterogeneity model:", x$heterogeneity_model_BIC,"\n",
        "  - delta BIC:", x$delta_BIC,"\n")
  }

  cat("\n")
  
  cat(" input summary \n",
      "  - number of variables in the network:", x$input$vars, "\n",
      "  - estimator used:", x$input$estimator, "\n",
      "  - estimated network type:",x$input$network_type,"\n",
      "  - homogeneity test:",x$input$homogeneity_test,"\n")

  if(x$input$save_models){
    cat(paste0("\n Estimated networks are stored in ", name,"$model"))
  }
  
  # Warning?
  # if (x$homogeneity_model_AIC < -1e+10 |  x$heterogeneity_model_AIC < -1e+10 | abs(x$delta_AIC) > 1e08){
  #   warning("POSSIBLE OPTIMIZATION ERRORS DETECTED: Check the results manually by setting save_models = FALSE.")
  # }
  
}

summary.INIT <- function(object,...) print(object)