# INIT

The `R` package `INIT` (Individual Network Invariance Test), provides an easy way to test for (in)equality between idiographic network structures. INIT extends standard model comparison techniques from Structural Equation Modeling (SEM) to network psychometrics. Similar to testing for Measurement Invariance (MI) using multi-group SEM methods to determine if the same concept is measured across individuals, INIT can be used to test if the same network structure applies across different individuals, or if the same network structure holds for one individual over time. INIT can be used to compare saturated (i.e., fully connected networks) or pruned idiographic network structures by inspecting model selection criteria such as the AIC and BIC, respectively.

# Installation 

You can install the current developmental version in R use the remotes-pakage:

```
remotes::install_github("RiaHoekstra/INIT")
```

# Example

## Comparing two individuals

This basic example shows how to apply INIT to compare the idiographic network structure of two individuals. To this end we will use the example dataset that is loaded with the package. This dataset contains simulated data for 6 variables, where n = 2 and t = 400. 

```
INIT(data = data, vars = colnames(data)[1:6], idvar = colnames(data)[7], estimator = "FIML", network_type = "saturated", homogeneity_test = "homogeneity_overall")
```

`network_type` can be adjusted to estimate "saturated" or "pruned" network structures. `homogeneity_test` can be adjusted to test for homogeneity between both the temporal and contemporaneous network structures, or to test for either the homogeneity between temporal or contemporaneous network structures. 
