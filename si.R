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
index <- expand.grid(hiv = seq(0, 4), male = c(0, 1))

## Initial parameters
## Function that makes a list of disease parameters with default values
disease_params <- function(Beta = 0.3
                           , alpha = 4 ## rate of beta decline with prevalence
                           , q = 2 ## Effect of behavior change in response to mortality
                           , progRt = (1/10)*4 ## rate of of progression through each of the I classes, for 10 years total
                           , birthRt = .03 ## birth rate, 3% of people give birth per year
                           , deathRt = 1/60 ## 60 year natural life expectancy
                           , ind = index

)
  return(as.list(environment()))

disease_params()
tseqMonth <- seq(1975, 2020, by = 1/12)

initInf <- exp(-7)
initSusc <- 5000

initial <- rep(0, nrow(index))
initial[index$hiv == 0] <- initSusc
initial[index$hiv == 1] <- initInf

SImod <- function(tt, nn, parms) with(c(parms, as.list(tt)), {
  
  ## browser()
  N <- sum(nn)
  I <- sum(nn[ind$hiv > 0]) ## Currently length 1: all HIV-infected people
  S <- nn[ind$hiv == 0] ## Currently length 2: males and females
  
  mortResponse <- exp(-q*progRt*sum(nn[ind$hiv == 4])/N)
  
  transmissionCoef <- mortResponse * Beta * exp(-alpha * sum(I)/N) ## Infectious contact rate
  
  ## state variable derivatives (ODE system)
  deriv <- rep(NA, length(nn))
  
  deriv[ind$hiv == 0] <- birthRt*N - deathRt*S - transmissionCoef*S*I/N ## Instantaneous rate of change: Susceptibles
  deriv[ind$hiv == 1] <-	transmissionCoef*S*I/N - progRt*nn[ind$hiv == 1] - deathRt*nn[ind$hiv == 1] ## Instantaneous rate of change: Infection class I1
  deriv[ind$hiv == 2] <-	progRt*nn[ind$hiv == 1] - progRt*nn[ind$hiv == 2] - deathRt*nn[ind$hiv == 2] ## Instantaneous rate of change:  Infection class I2
  deriv[ind$hiv == 3] <-	progRt*nn[ind$hiv == 2] - progRt*nn[ind$hiv == 3] - deathRt*nn[ind$hiv == 3] ## Instantaneous rate of change:  Infection class I2
  deriv[ind$hiv == 4] <-	progRt*nn[ind$hiv == 3] - progRt*nn[ind$hiv == 4] - deathRt*nn[ind$hiv == 4] ## Instantaneous rate of change:  Infection class I2
  # deriv[6] <-	transmissionCoef*S*I/N ## Instantaneous rate of change: Cumulative incidence
  # deriv[7] <-	progRt*I4 ## Instantaneous rate of change: Cumulative mortality
  
  return(list(deriv))
})

SImod(tseqMonth, initial, disease_params())


# init <- c(S=initPop, I1 = initPrev, I2 = emptyArray, I3 = emptyArray, I4 = emptyArray, CI = emptyArray, CD =  emptyArray) ## modeling proportion of population
# Is <- paste0('I',1:4) ## for easy indexing

## Function to run the deterministic model simulation, based on the ODE system defined in SImod().
simEpidemic <- function(tseq = tseqMonth, init = initial, modFunction=SImod, parms = disease_params()) {
  simDat <- as.data.frame(lsoda(init, tseq, modFunction, parms=parms))
  simDat$I <- rowSums(simDat[, 1 + (which(index$hiv > 0))])
  simDat$N <- rowSums(simDat[, 2:ncol(simDat)])
  simDat$P <- with(simDat, I/N)
  return(simDat)
}

simEpidemic()

## Prevalence estimates
kzn_prev <- read_csv("data/kzn_hiv_prev_survey.csv")
kzn_prev$source <- "HSRC 15-49"
zaf_prev <- read_csv("data/zaf_hiv_prev_ests.csv")
zaf_prev$source <- "Spectrum"
zaf_prev$numSamp <- 300 ## Make up "sample size" for the UNAIDS estimates for now
# fake_prev <- data.frame(time = seq(1980, 1985, by = 1), source = "Fake", mean = 0.001, numSamp = 300)
prev <- rbind.fill(kzn_prev, zaf_prev)
prev$numPos <- round(prev$mean * prev$numSamp, 0)
prev <- prev[prev$time < 2005, ] ## For now fit pre-ART data
prev <- prev[prev$source == "Spectrum", ]
rm(kzn_prev, zaf_prev)

## Since we need to be able to easily separate fitted and fixed parameters,
## let's make a function that takes fitted and fixed parameters together and
## puts them back in a parameter list (similar to the output of
## disease_params()). We also want it to be able to take logged parameter values
## and unlog them automatically back into the parameter list, so that we can fit
## on a log scale, but run our model with the unlogged values.
subsParms <- function(fit.params, fixed.params=disease_params())
  within(fixed.params, {
    loggedParms <- names(fit.params)[grepl('log_', names(fit.params))]
    unloggedParms <- names(fit.params)[!grepl('log_', names(fit.params))]        
    for(nm in unloggedParms) assign(nm, as.numeric(fit.params[nm]))
    for(nm in loggedParms) assign(gsub('log_','',nm), exp(as.numeric(fit.params[nm])))
    rm(nm, loggedParms, unloggedParms)
  })

## Likelihood
nllikelihood <- function(parms = disease_params(), obsDat=prev) {
  simDat <- simEpidemic(init = initial, parms=parms)
  ## What are the rows from our simulation at which we have observed data?
  matchedTimes <- simDat$time %in% obsDat$time
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P[matchedTimes], log = T)
  return(sum(nlls))
}

guess.params <- c(log_Beta = log(5), log_alpha = log(8), log_q = log(2))
subsParms(guess.params, disease_params())

## Make likelihood a function of fixed and fitted parameters.
objFXN <- function(fit.params ## parameters to fit
                   , fixed.params = disease_params() ## fixed paramters
                   , obsDat=prev) {
  parms <- subsParms(fit.params, fixed.params)
  nllikelihood(parms, obsDat = obsDat) ## then call likelihood
}
objFXN(guess.params, disease_params())

## Optimization
trace <- 3

## SANN: This is stochastic, be CAREFUL sometimes it gets stuck at local minima
## for unreasonble parameters. If you see this happen, run it again!
init.pars <- c(log_alpha = log(4), log_Beta = log(.9), log_q = log(2))
optim.vals <- optim(par = init.pars
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = prev
                    , control = list(trace = trace, maxit = 150)
                    , method = "SANN")
exp(optim.vals$par)

## Normally we use SANN first and then follow with Nelder-Mead since SANN is stochastic and will
## make sure to help you be sure that you aren't at local minima. We feed the last parameters of
## SANN in as the first values of Nelder-Mead
optim.vals <- optim(par = optim.vals$par
                    , objFXN
                    , fixed.params = disease_params()
                    , obsDat = prev
                    , control = list(trace = trace, maxit = 800, reltol = 10^-7)
                    , method = "Nelder-Mead"
                    , hessian = T)
optim.vals
MLEfits <- optim.vals$par
exp(MLEfits)

## Run simulation
out <- simEpidemic(init = initial, parms = disease_params("Beta" = exp(MLEfits['log_Beta']), "alpha" = exp(MLEfits['log_alpha']), "q" = exp(MLEfits['log_q'])))

# disease_params <- function(Beta = 0.45
#                            , alpha = 1 ## rate of beta decline with prevalence
#                            , q = 80 ## Effect of behavior change in response to mortality
#                            , progRt = (1/10)*4 ## rate of of progression through each of the I classes, for 10 years total
#                            , birthRt = .03 ## birth rate, 3% of people give birth per year
#                            , deathRt = 1/60 ## 60 year natural life expectancy
#                            
# )
#   return(as.list(environment()))
# 
# disease_params()
# 
# initPrev <- exp(3)
# 
# out <- simEpidemic(init = init)

out$S <- rowSums(out[, 1 + which(index$hiv == 0)])
long <- (out
  %>% gather(., state, n, -time)
)

## Plots
ggplot(data = prev, aes(x = time, y = mean)) +
  geom_point(aes(shape = source)) +
  ## geom_errorbar(aes(ymin = lower, ymax = upper)) +
  geom_line(data = out, aes(x = time, y = P)) +
  labs(x = "Year", y = "Prevalence") +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  ggtitle("HIV Prevalence")


ggplot(data = long[long$state %in% c("S", "I", "N"), ], aes(x = time, y = n)) +
  geom_line(aes(colour = state)) +
  labs(x = "Year", y = "Number") +
  ggtitle("Number of people in S, I, and N compartments")





