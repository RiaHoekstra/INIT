# INIT

The `R` package `INIT` (Individual Network Invariance Test) provides an easy way to test for (in)equality between idiographic network structures. INIT extends standard model comparison techniques from Structural Equation Modeling (SEM) to network psychometrics and can be used to test if the same network structure holds across different individuals or if the same network structure holds for one individual over time. By inspecting model selection criteria such as the AIC and BIC, INIT provides an easy way to compare idiographic network structures.

# Installation 

You can install the current developmental version in R using the remotes-package:

```
remotes::install_github("RiaHoekstra/INIT")
```

# Example 

## Compare networks using default arguments
This example shows how to compare idiographic networks using INIT's default arguments. 

```
library(INIT)

# Load data of two individuals: 
load("INIT_data.rda") 

# Compare networks with default arguments:
INIT(data = data, idvar = colnames(data)[1], vars = colnames(data)[2:7])
```

## Changing arguments
This example shows how default arguments such as 'estimator', 'network_type', and 'homogeneity_test' can be specified. 

```
library(INIT)
# Load data of two individuals: 
load("INIT_data.rda") 

# Compare networks with default arguments:
INIT(data = data, idvar = colnames(data)[1],vars = colnames(data)[2:7], estimator = "FIML", network_type = "pruned", homogeneity_test = "temporal", save_models = TRUE)
```

## Inspect estimated networks
INIT makes use of the `psychonetrics` package to estimate idiographic network structures. This example shows how the estimated networks can be saved and inspected. 

```
library(INIT)
# Load data of two individuals: 
load("INIT_data.rda") 

# Save INIT results in an object:
res <- INIT(data = data, idvar = colnames(data)[1],vars = colnames(data)[2:7], save_models = TRUE)

# Inspect summary of results: 
res 

# Inspect estimated network structures: 
res$network
```

