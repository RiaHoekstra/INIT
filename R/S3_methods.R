print.INIT <- function(x,...){
  
  cat("INIT test result:\n",
      "  - Network type:",x$input$network_type,"\n",
      "  - Homogeneity test:",x$input$homogeneity_test,"\n",
      "  - AIC homogeneity model:",x$equal_network_AIC,"\n",
      "  - AIC heterogeneity model:",x$different_network_AIC,"\n",
      "  - The",ifelse(x$delta_AIC<0,"HETEROGENEITY","HOMOGENEITY"),"model is preferred!","\n")
  
  # Warning?
  if (x$equal_network_AIC < -1e+10 |  x$different_network_AIC < -1e+10 | abs(x$delta_AIC) > 1e08){
    warning("POSSIBLE OPTIMIZATION ERRORS DETECTED: Check the results manually by setting save_models = FALSE.")
  }
  
}

summary.INIT <- function(object,...) print(object)