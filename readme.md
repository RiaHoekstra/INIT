# INIT

The `R` package `INIT` (Individual Network Invariance Test), provides an easy way to test for (in)equality between idiographic network structures. INIT extends standard model comparison techniques from Structural Equation Modeling (SEM) to network psychometrics. Similar to testing for Measurement Invariance (MI) using multi-group SEM methods to determine if the same concept is measured across individuals, INIT can be used to test if the same network structure applies across different individuals, or if the same network structure holds for one individual over time. INIT can be used to compare saturated (i.e., fully connected networks) or pruned idiographic network structures by inspecting model selection criteria such as the AIC and BIC, respectively.

# Installation 

You can install the current developmental version in R use the remotes-pakage:

```
remotes::install_github("RiaHoekstra/INIT")
```

# Example

## Comparing two individuals

This basic example shows how to apply INIT to compare the idiographic network structure of two individuals.

```
# Generate two homogenous VAR models:
set.seed(42)

library("graphicalVAR")

# Simulate model:
Mod <- randomGVARmodel(6,probKappaEdge = 0.8,probBetaEdge = 0.8)

# Simulate dataset 1:
Data1 <- graphicalVARsim(100,Mod$beta,Mod$kappa)

# Simulate dataset 2:
Data2 <- graphicalVARsim(100,Mod$beta,Mod$kappa)

# Combine:
Data1 <- as.data.frame(Data1)
Data1$id <- 1

Data2 <- as.data.frame(Data2)
Data2$id <- 2

vars <- paste0("V",1:4)

Data <- rbind(Data1,Data2)

# Run INIT:
init_res <- INIT(
  data = Data,
  vars = vars,
  idvar = "id",
  estimator = "FIML",
  network_type = "saturated",
  homogeneity_test = "homogeneity_overall"
  )

# Print results:
init_res
```

`network_type` can be adjusted to estimate "saturated" or "pruned" network structures. `homogeneity_test` can be adjusted to test for homogeneity between both the temporal and contemporaneous network structures, or to test for either the homogeneity between temporal or contemporaneous network structures. 
