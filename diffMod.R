##################################
## Allen Roberts
## 5 June 2017
## Dynamic TaSP model for DO ART
##################################

rm(list = ls())

library(tidyverse)
library(plyr)
library(deSolve)

theme_set(theme_bw())

## Index
age <- seq(1, 12) ## First dimension
age_names <- c("0-4", "5-9", "10-14", "15-19", "20-24", "25-29", "30-34", "35-39", "40-44", "45-49", "50-54", "55-59")
sex <- seq(1, 2) ## Second dimension
sex_names <- c("male", "female")
hiv <- seq(1, 5) ## Third dimension - could consider renaming this to "disease" or something like that
hiv_names <- c("negative", "stage1", "stage2", "stage3", "stage4")
hiv_states <- c("stage1", "stage2", "stage3", "stage4")

pop <- array(data = rep(0, length(age) * length(sex) * length(hiv)),
             dim = c(length(age), length(sex), length(hiv)),
             dimnames = list(age_names, sex_names, hiv_names))

diff <- pop


## Initial parameters
## Function that makes a list of disease parameters with default values
disease_params <- function(Beta = 0.3
                           , alpha = 4 ## rate of beta decline with prevalence
                           , q = 2 ## Effect of behavior change in response to mortality
                           , progRt = (1/10)*4 ## rate of of progression through each of the I classes, for 10 years total
                           , birthRt = rep(.03*2, length(age)) ## age-specific birth rate, 3% of people give birth per year
                           , deathRt = 1/(60*2) ## 60 year natural life expectancy
                           , ageRt = 1/5 ## In a given year, 1/5 of a 5-year age group will progress to the next age group
                           , state = pop ## state variables
                           , deriv = diff ## change in state variables
                           
)
  return(as.list(environment()))

tseqMonth <- seq(1975, 2020, by = 1/12)

initInf <-  exp(-7)
initSusc <- 5000

pop[,,"negative"] <- initSusc/(length(age) * length(sex))
pop[,,hiv_states] <- initInf/(length(age) * length(sex) * length(hiv_states))

SImod <- function(tt, nn, parms) with(c(parms, as.list(tt)), {
  
  ## browser()
  ## Derived state variables
  N <- sum(state)
  S <- sum(state[,,"negative"])
  I <- sum(state[,, hiv_states])
  
  ## Derived parameters
  mortResponse <- exp(-q*progRt*sum(state[,,"stage4"])/N) ## Factor by which observed mortality reduces beta
  transmissionCoef <- mortResponse * Beta * exp(-alpha * I/N) ## Infectious contact rate
  
  ## Derivatives ## Maybe can make this faster by pre-allocating arrays for births, deaths, aging, etc?
  ## Births
  deriv["0-4",,"negative"] <- deriv["0-4",,"negative"] + sum(birthRt * state[, "female", ])/length(sex) ## Divide between male and female births
  
  ## Background mortality
  deriv <- deriv - deathRt * state
  
  ## Aging
  deriv <- deriv - ageRt * state ## People leaving the compartment due to age
  deriv[2:length(age),,] <- deriv[2:length(age),,] + ageRt * state[1:(length(age) - 1),,] ## People entering the compartment due to age
  
  ## Disease Progression
  deriv[,,hiv_states] <- deriv[,,hiv_states] - progRt * state[,,hiv_states] ## Progression out of the HIV compartments
  deriv[,,c("stage2", "stage3", "stage4")] <- deriv[,,c("stage2", "stage3", "stage4")] + progRt * state[,,c("stage1", "stage2", "stage3")] ## Progression into the HIV compartments from previous compartments
  
  ## Transmission matrices
  prev_mat <- apply(state[,,hiv_states], 1:2, sum)/apply(state, 1:2, sum) ## Age- and sex-specific prevalence
  susc_mat <- apply(state[,,"negative"], 1:2, sum) ## Age- and sex-specific number of susceptibles
  
  ## Transmission
  deriv[,,"stage1"] <- deriv[,,"stage1"] + transmissionCoef*susc_mat*prev_mat[, rev(colnames(susc_mat))] ## Note that I'm reversing the order of the columns of the prevalence matrix, which ensures that the prevalence in females is multipled by the number of susceptibles in males (and vice versa). This implies perfect assortativity by age.
  deriv[,,"negative"] <- deriv[,,"negative"] - transmissionCoef*susc_mat*prev_mat[, rev(colnames(susc_mat))]
  
  return(deriv)
})


for(ii in tseqMonth) {
  diff <- SImod(tseqMonth, as.vector(pop), disease_params())
  pop <- pop + diff
  diff <- array(data = rep(0, length(age) * length(sex) * length(hiv)),
                dim = c(length(age), length(sex), length(hiv)),
                dimnames = list(age_names, sex_names, hiv_names))
}



reconstruct <- function(vec) {
  array(data = vec,
        dim = c(length(age), length(sex), length(hiv)),
        dimnames = list(age_names, sex_names, hiv_names))
  
}

SImod(tseqMonth, as.vector(pop), disease_params())
reconstruct(unlist(SImod(tseqMonth, as.vector(pop), disease_params())))



# init <- c(S=initPop, I1 = initPrev, I2 = emptyArray, I3 = emptyArray, I4 = emptyArray, CI = emptyArray, CD =  emptyArray) ## modeling proportion of population
# Is <- paste0('I',1:4) ## for easy indexing
