set.seed(5)

library("bootnet")
library("graphicalVAR")

rewireFun <- function(x, p, directed){
  if (missing(directed)){
    directed <- !all(x == t(x))
  }
  
  if (directed){
    ind <- diag(1,ncol(x)) != 1
  } else {
    ind <- upper.tri(x)
  }
  
  # Current edges:
  curEdges <- which(x!=0 & ind, arr.ind = TRUE)
  
  # Select to rewire:
  toRewire <- which(runif(nrow(curEdges)) < p)
  
  for (i in seq_along(toRewire)){
    curZeros <- which(x==0 & ind, arr.ind=TRUE)
    dest <- sample(seq_len(nrow(curZeros)),1)
    
    x[curZeros[dest,1],curZeros[dest,2]] <- x[curEdges[toRewire[i],1],curEdges[toRewire[i],2]]
    x[curEdges[toRewire[i],1],curEdges[toRewire[i],2]] <- 0
    
    if (!directed){
      x[curZeros[dest,2],curZeros[dest,1]] <- x[curEdges[toRewire[i],2],curEdges[toRewire[i],1]]
      x[curEdges[toRewire[i],2],curEdges[toRewire[i],1]] <- 0
    }
  }
  
  return(x)
}

pcor2prec <- function(x){
  n <- ncol(x)
  out <- solve(cov2cor(solve(diag(n) - x)))
  out[abs(x) < sqrt(.Machine$double.eps) & !diag(n)==1] <- 0
  out
}

# ------------------------------------------------------------------------------
# Simulate time-series for two individuals under null hypothesis, 
# i.e., underlying network structures are similar
# ------------------------------------------------------------------------------

nNode <- 6 # number of nodes
nTime <- 300 # number of time-series per individual

# GVAR parameters: 

# Generate contemporaneous effects: 
pcor1 <- genGGM(nNode, p = 0, propPositive = 0.8)
pcor2 <- rewireFun(pcor1, p = 0, FALSE)

kappa1 <- pcor2prec(pcor1)
kappa2 <- pcor2prec(pcor2)

# Generate temporal effects
trueBeta <- diag(1, nNode)
for (i in 1:nNode){
  trueBeta[(i+(2-1)) %%nNode+1, i] <- sample(c(-1, 1), 1, 0.5)
}

beta1 <- trueBeta * rnorm(nNode^2, 0.3, 0.1)
beta2 <- rewireFun(beta1, p = 0)

# Add between subject effects? 

# Simulate data:
data1 <- as.data.frame(graphicalVARsim(nTime, beta1, kappa1))
data1$id <- 1
data2 <- as.data.frame(graphicalVARsim(nTime, beta2, kappa2))
data2$id <- 2

# Combine data for multi group:
data <- rbind(data1, data2)


# ------------------------------------------------------------------------------
# Compare networks using INIT
# ------------------------------------------------------------------------------

# Idiographic networks can be compared on (1) overall homogeneity, (2) contemporal network 
# homogeneity, (3) temporal network homogeneity for (1) saturated networks or 
# (2) pruned networks. 
# Below all options are illustrated. 

INIT(data = data, 
     idvar = "id", 
     estimator = "FIML", 
     network = "saturated")


data = data
idvar = "id"
estimator = "FIML"
vars = colnames(data)[1:6]
networkType = "pruned"


# To do: 
# if model is already estimated? 

