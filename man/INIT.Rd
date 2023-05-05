\name{INIT}
\alias{INIT}
\title{Function to test for homogeneity between idiographic network models}
\usage{
INIT(
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
)
}
\arguments{
\item{data}{Data set.}

\item{vars}{String indicating the variables included in the model.}

\item{dayvar}{String indicating the day.}

\item{beepvar}{String indicating the beep.}

\item{idvar}{String indicating the subject id variable name.}

\item{network_type}{Type of network to estimate. Currently implemented are "saturated" or "pruned".
By default set to "saturated". When "pruned" is specified, non significant edges are pruned from the model with an alpha rate of 0.05.}

\item{homogeneity_test}{Type of homogeneity test. "overall", refers to
test both contemporaneous and temporal network structures for homogeneity,
"contemporaneous" to test contemporaneous network structures for homogeneity,
"temporal" to test temporal network structures for homogeneity.
By default set to "overall".}

\item{save_models}{If TRUE, model output for contemporaneous and temporal network models will be saved and returned as output.
By default set to FALSE to save memory.}

\item{alpha}{The alpha level used in pruning}

\item{...}{any argument sent to psychonetrics::gvar}
}
\description{
Function to test for homogeneity between idiographic network models
}
\author{
Ria Hoekstra <h.a.hoekstra@uva.nl>
}
\examples{

# Generate two homogenous VAR models:
library("graphicalVAR")

# Simulate model:
set.seed(42)

Mod <- randomGVARmodel(6,probKappaEdge = 0.8,probBetaEdge = 0.8)

# Simulate data 1:
Data1 <- graphicalVARsim(100,Mod$beta,Mod$kappa)

# Simulate data 2:
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
  Data,
  vars,
  idvar = "id",
  standardize = "z"
)

# Print results:
init_res

}
