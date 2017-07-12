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

zeroes <- array(data = rep(0, length(age) * length(sex) * length(hiv)),
             dim = c(length(age), length(sex), length(hiv)),
             dimnames = list(age_names, sex_names, hiv_names))

pop <- zeroes

## Function to reconstruct multi-dimensional array from single dimension output vector
reconstruct <- function(vec) {
  array(data = vec,
        dim = c(length(age), length(sex), length(hiv)),
        dimnames = list(age_names, sex_names, hiv_names))
  
}

## Function to create data.frame from single dimension output vector. Note that this function depends on the order by which the array is initiated in the first place
createOutputFrame <- function(n) {
  
  df <- expand.grid(age = age_names, sex = sex_names, hiv = hiv_names)
  df$n <- n
  return(df)
}

## Initial parameters
## Function that makes a list of disease parameters with default values
disease_params <- function(Beta = 0.7
                           , alpha = 7 ## rate of beta decline with prevalence
                           , q = 2 ## Effect of behavior change in response to mortality
                           , progRt = (1/10)*4 ## rate of of progression through each of the I classes, for 10 years total
                           , birthRt = rep(0.05, length(age)) ## age-specific birth rate, 3% of people give birth per year
                           , deathRt = 1/(60) ## 60 year natural life expectancy
                           , ageRt = 1/5 ## In a given year, 1/5 of a 5-year age group will progress to the next age group
                           , births = zeroes
                           , deaths = zeroes
                           , age_in = zeroes
                           , age_out = zeroes
                           , prog_in = zeroes
                           , prog_out = zeroes
                           , trans = zeroes
                

)
  return(as.list(environment()))

tseqMonth <- seq(1975, 2020, by = 1/12)

initInf <- exp(-7)
initSusc <- 5000

pop[,,"negative"] <- initSusc/(length(age) * length(sex))
pop[,,"stage1"] <- initInf/(length(age) * length(sex))

SImod <- function(tt, yy, parms) with(as.list(parms), {
  
  ## browser()
  
  state <- reconstruct(yy)

  ## Derived state variables
  N <- sum(state)
  S <- sum(state[,,"negative"])
  I <- sum(state[,, hiv_states])
  
  ## Derived parameters
  mortResponse <- exp(-q*progRt*sum(state[,,"stage4"])/N) ## Factor by which observed mortality reduces beta
  transmissionCoef <- mortResponse * Beta * exp(-alpha * I/N) ## Infectious contact rate
  
  ## Derivatives ## Maybe can make this faster by pre-allocating arrays for births, deaths, aging, etc?
  ## Births
  births["0-4",,"negative"] <- sum(birthRt * state[, "female", ])/length(sex) ## Divide between male and female births
  
  ## Background mortality
  deaths <- - deathRt * state
    
  ## Aging
  age_out <- -ageRt * state
  age_in[2:length(age),,] <- ageRt * state[1:(length(age) - 1),,]
  
  ## Disease Progression
  prog_out[,,hiv_states] <- -progRt * state[,,hiv_states]
  prog_in[,, c("stage2", "stage3", "stage4")] <- progRt * state[,,c("stage1", "stage2", "stage3")]

  ## Transmission matrices
  prev_mat <- apply(state[,,hiv_states], 1:2, sum)/apply(state, 1:2, sum) ## Age- and sex-specific prevalence
  susc_mat <- apply(state[,,"negative"], 1:2, sum) ## Age- and sex-specific number of susceptibles
  
  ## Transmission
  trans[,,"stage1"] <- transmissionCoef*susc_mat*prev_mat[, rev(colnames(susc_mat))]
  trans[,,"negative"] <- -transmissionCoef*susc_mat*prev_mat[, rev(colnames(susc_mat))]
  ## Note that I'm reversing the order of the columns of the prevalence matrix, which ensures that the prevalence in females is multipled by the number of susceptibles in males (and vice versa). This implies perfect assortativity by age.

  deriv <- births + deaths + age_in + age_out + prog_in + prog_out + trans
  return(list(c(deriv)))
})


## SImod(yy = c(pop), tt = tseqMonth, disease_params())
# reconstruct(unlist(SImod(tseqMonth, as.vector(pop), disease_params())))

## Function to run the deterministic model simulation, based on the ODE system defined in SImod().
simEpidemic <- function(tseq = tseqMonth, init = c(pop), modFunction=SImod, parms = disease_params(), detail = "fit") {
  
  simDat <- as.data.frame(lsoda(y = init, times = tseq, modFunction, parms=parms))
  
  if(detail == "fit") { ## Calculate prevalence for fitting purposes
    
    out <- apply(simDat[, 2:ncol(simDat)], 1, function(x){
      
      state <- reconstruct(x)
      N <- sum(state)
      S <- sum(state[,,"negative"])
      I <- sum(state[,, hiv_states])
      
      return(c(N, S, I))
      
      
    })
    
    out <- as.data.frame(t(out))
    names(out) <- c("N", "S", "I")
    out$time <- tseqMonth
    out$P <- with(out, I/N)
    
  } else if(detail == "output") { ## Return full output data.frame to calculate statistics
    
    out_list <- alply(simDat, 1, function(x) {
      df <- createOutputFrame(c(t(x[2:ncol(x)])))
      df$time <- x[1, 1]
      return(df)
    })
    
    out <- do.call("rbind", out_list)
    
    
  }
  
  return(out)
}

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
  simDat <- simEpidemic(init = c(pop), parms=parms, detail = "fit")
  ## What are the rows from our simulation at which we have observed data?
  matchedTimes <- simDat$time %in% obsDat$time
  nlls <- -dbinom(obsDat$numPos, obsDat$numSamp, prob = simDat$P[matchedTimes], log = T)
  return(sum(nlls))
}

guess.params <- c(log_Beta = log(0.5), log_alpha = log(8), log_q = log(2))
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
init.pars <- c(log_alpha = log(7), log_Beta = log(0.7), log_q = log(2))
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
out <- simEpidemic(init = c(pop), parms = disease_params("Beta" = exp(MLEfits['log_Beta']), "alpha" = exp(MLEfits['log_alpha']), "q" = exp(MLEfits['log_q'])), detail = "output")

## Calculate statistics
## Sex-specific population
detach(package:plyr) ## Apparently this is necessary for dplyr to work correctly
sex_pop <- (out
            %>% group_by(., sex, time)
            %>% summarise(., sex_total = sum(n))
)

sex_pop_plot <- ggplot(data = sex_pop, aes(x = time, y = sex_total/10000)) +
  geom_line(aes(colour = sex)) +
  labs(x = "Year", y = "Total population (tens of thousands)")

## Age/sex-specific population
age_sex_pop <- (out
            %>% group_by(., sex, age, time)
            %>% summarise(., age_sex_total = sum(n))
)

age_sex_pop_plot <- ggplot(data = age_sex_pop, aes(x = time, y = age_sex_total/10000)) +
  geom_line(aes(colour = sex)) +
  labs(x = "Year", y = "Total population (tens of thousands)") +
  facet_wrap(~age)

## Sex-specific prevalence
sex_inf_pop <- (out
            %>% filter(., hiv != "negative")
            %>% group_by(., sex, time)
            %>% summarise(., sex_inf = sum(n))
)

sex_prev <- inner_join(sex_pop, sex_inf_pop)
sex_prev$prev <- with(sex_prev, sex_inf/sex_total)
  
sex_prev_plot <- ggplot(data = sex_prev, aes(x = time, y = prev)) +
  geom_line(aes(colour = sex)) +
  labs(x = "Year", y = "Prevalence") 

## Age/sex-specific prevalence
age_sex_inf_pop <- (out
                %>% filter(., hiv != "negative")
                %>% group_by(., age, sex, time)
                %>% summarise(., age_sex_inf = sum(n))
)

age_sex_prev <- inner_join(age_sex_pop, age_sex_inf_pop)
age_sex_prev$prev <- with(age_sex_prev, age_sex_inf/age_sex_total)

age_sex_prev_plot <- ggplot(data = age_sex_prev, aes(x = time, y = prev)) +
  geom_line(aes(colour = sex)) +
  labs(x = "Year", y = "Prevalence") +
  facet_wrap(~age)





