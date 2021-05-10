### This is the code to calculate the delta_F (increase or decrease) in feasibility in
#### a 2-species community by addding a non-resident species
rm(list=ls())
### load libraries
library('pcalg')
source("functions.r")
if(!require(deSolve)) {install.packages("deSolve"); library(deSolve)}
if(!require(matrixcalc)) {install.packages("matrixcalc"); library(matrixcalc)}

### This is the 5-species interaction matrix inferred from Gould et al. PNAS (2018) using the methodology described in Deng et al. 2021
## Note that using the inferred pairwise interaction matrix in Gould et al. yield qualitatively similar results
alpha5 <- as.matrix(read.table(paste("Gould_Matrix_2.txt")))
S5 <- 5
q <- combn(1:5,3)  ### the 10 combinations of 3-species communities

deltaF <- matrix(0,10,3) ## initialize
S <- 3  # number of residents plus non-resident
counter <- 0
for(j in 1:10){ 
  index <- q[,j]  # obtain combinations
  alpha <- alpha5[index,index] # extract the 3-species interaction matrix

  for(i in 1:S){ ## Each species in the community becomes a non-resident at a time
    
    ## calculate the feasibility of the 3 species
    C3 <- diag(c(-1,-1,-1), 3) # constraints to use only the positive orthant for r's
    e <- try(om <- Omega_overlap(alpha,C3),silent = T) # calculates Omega^(1/S) using triangulation
    if (class(e)!="try-error") {f3 <- om} # if trinagulation is possible
    if (class(e)=="try-error") { # if triangulation is not possible then use Monte Carlo as approximation
      om <- 0
      for(x in 1:30){ ## to reduce numerical instabilities by convergence
          om <- om + (Omega(alpha)/2) # this is to penalize for the positive orthant
      }
      f3 <- om / (30)
    }
    #f3 <- Omega_overlap(alpha,C3) # calculates Omega^(1/S)
    
    counter <- counter + 1
    alpham <- alpha5[index,index]
    alpham <- alpham[-i,-i]  # makes the 2-species commuunity

    C2 <- diag(c(-1,-1), 2)  # constraints to use only the positive orthant for r's
    e <- try(om <- Omega_overlap(alpham,C2),silent = T) # calculates Omega^(1/S) using triangulation
    if (class(e)!="try-error") {f2 <- om} # if trinagulation is possible
    if (class(e)=="try-error") {  # if triangulation is not possible then use Monte Carlo as approximation
      om <- 0
      for(x in 1:30){ ## to reduce numerical instabilities by convergence
        om <- om + (Omega(alpham)/2) # this is to penalize for the positive orthant
      }
      f2 <- om / (30)
    }
    #f2 <- Omega_overlap(alpham,C2)

    deltaF[counter] <- f3 - f2  # calculate delta_F for each case from the combinations q
    print(c(index,i,deltaF[counter])) # index = 3-species community; i = non-resident, delta_F 
  }
}