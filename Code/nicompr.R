#### Non-Informative COM-Poisson priors using RStan or NICOMPR ####
#   by:         Meyer, MJ, Graye, A, & Sellers, KF
#   modified:   03/23/23
# NICOMPR

#### load libraries ####
library(rstan)
library(shinystan)

#### compile CMP stan models ####

##### flat and conjugate prior model #####
# see file compois_conjugate.stan for details
cmpmodel  <- stan_model(file = 'Code/compois_conjugate.stan')

##### jeffreys' prior model #####
# see file compois_jeffreys.stan for details
cmpmodel_j  <- stan_model(file = 'Code/compois_jeffreys.stan')

#### diagnostic functions ####

##### WAIC #####
WAIC <- function(x, ...){
  UseMethod('WAIC')
}

WAIC.stanfit <- function(model){
  psamp       <- extract(model)
  logLikMat   <- psamp$log_lik
  lppd        <- sum(log(apply(exp(logLikMat), 2, mean)))
  pWAIC       <- sum(apply(logLikMat, 2, var))
  waic_out    <- -2*lppd + 2*pWAIC
}

WAIC.cmp <- function(model){
  model$diagnostics$WAIC
}

##### loo #####
loo.cmp <- function(model){
  model$diagnostics$loo$estimates[3,1]
}

#### model functions ####
paramCheck <- function(a,b,c){ 
  b/c > log(factorial(floor(a/c))) + (a/c - floor(a/c))*log(floor(a/c) + 1) 
}

Z   <- function(lambda, nu, J){
  log_Z_terms <- vector('numeric', length = J + 1)
  for(j in 0:J){
    log_Z_terms[j] <- j * log(lambda) - nu * lgamma(j + 1)
  }
  Z_out   <- sum(exp(log_Z_terms))

  return(Z_out)
}

cmp   <- function(X, priors = list(a = 0, b = 0, c = 0), 
                  type = c('conjugate', 'flat', 'Jeffreys'), 
                  chains = 4, iter = 2000, ...){
  if(length(type) > 1){
    type <- 'conjugate'
  }
  
  if(type != 'flat' & type != 'conjugate' & type != 'Jeffreys'){
    stop('prior type must be flat, conjugate, or Jeffreys')
  }
  
  if(type == 'conjugate' & is.null(priors)){
    warning('conjugate prior requires selection of hyper-parameters, defaults to flat prior when arugment priors is exlcuded')
  }
  
  if(type == 'flat'){
    a <- b <- c <- 0
  } else{
    if(is.null(priors)){
      a <- b <- c <- 1
    } else{
      a <- priors$a
      b <- priors$b
      c <- priors$c
    }
  }
  
  
  if(type == 'conjugate'){
    checkPriors <- paramCheck(a, b, c)
    if(!checkPriors){
      stop('hyper-parameters must staisfy paramCheck(a, b, c)')
    }
  }
  
  if(type == 'Jeffreys'){
    dl    <- list(n = length(X), # sample size
                  X = X, # data
                  S1 = sum(X), S2 = sum(lfactorial(X))) # sufficient statistics
    fit   <- sampling(cmpmodel_j, data = dl, chains = chains, iter = iter, ...)
  } else{
    dl    <- list(n = length(X), # sample size
                  X = X, # data
                  S1 = sum(X), S2 = sum(lfactorial(X)), # sufficient statistics
                  a = a, b = b, c = c)                  # hyper-parameters
    fit   <- sampling(cmpmodel, data = dl, chains = chains, iter = iter, ...)
  }
  
  ## diagnostics ##
  fit_loo     <- suppressWarnings(loo(fit))
  fit_waic    <- WAIC(fit)
  lambdaMat   <- matrix(unlist(extract(fit, pars = 'lambda')), nrow = chains,
                        byrow = TRUE)
  nuMat       <- matrix(unlist(extract(fit, pars = 'nu')), nrow = chains,
                        byrow = TRUE)
  fit_Rhats   <- list(lambda = Rhat(t(lambdaMat)), nu = Rhat(t(nuMat)))
  
  ## hmc details ##
  hmc         <- list(chains = chains, iter = iter, 
                      numDiverge = sum(get_divergent_iterations(fit)))
  
  out <- list(fit = fit, data = dl, priors = priors,
              diagnostics = list(loo = fit_loo, WAIC = fit_waic, Rhat = fit_Rhats),
              type = type, hmc = hmc)
  
  class(out) <- 'cmp'
  
  return(out)
}

print.cmp <- function(model, ...){
  print(model$fit, pars = c('lambda', 'nu'), ...)
}



coef.cmp <- function(model, ...){
  posts           <- extract(model$fit, pars = c('lambda', 'nu'), ...)
  lambda_med      <- median(posts$lambda)
  nu_med          <- median(posts$nu)
  out             <- matrix(c(lambda_med, nu_med), nrow = 1)
  colnames(out)   <- c('lambda', 'nu')
  rownames(out)   <- ''
    
  return(out)
}

summary.cmp <- function(model, alpha = 0.05, digits = 4, ...){
  posts           <- extract(model$fit, pars = c('lambda', 'nu'), ...)
  lambda_q        <- quantile(posts$lambda, probs = c(0.5, alpha/2, 1-alpha/2))
  nu_q            <- quantile(posts$nu, probs = c(0.5, alpha/2, 1-alpha/2))
  sums            <- rbind(lambda_q, nu_q)
  oldColnames     <- colnames(sums)
  tab             <- cbind(sums, unlist(model$diagnostics$Rhat))
  rownames(tab)   <- c('lambda', 'nu')
  colnames(tab)   <- c(oldColnames, 'Rhat')
  type            <- model$type
  n               <- model$data$n
  cat(paste("COM-Poisson using", type, "prior\n"))
  cat(paste("Model fit on", n, "total subjects\n\n"))
  
  cat("Model fit criteria\n")
  cat(paste("WAIC:", round(WAIC(model))), "\n")
  cat(paste("loo: ", round(loo(model)), "\n\n"))
  
  print(round(tab, digits = digits))
  cat(paste("\nPosterior samples based on", model$hmc$chains, "chains of", model$hmc$iter, "samples\n"))
  cat(paste("after a warm-up of size", model$hmc$iter/2, "(per chain, discarded)"))
  
  if(model$hmc$numDiverge > 0){
    cat(paste("\n\nCaution:", model$hmc$numDiverge, "divergent chains detected"))
  }
}

posts <- function(x, ...){
  UseMethod('posts')
}

posts.cmp <- function(model, ...){
  post           <- extract(model$fit, ...)
  
  return(post)
}

credint <- function(x, ...){
  UseMethod('credint')
}

credint.cmp <- function(model, alpha = 0.05, ...){
  posts           <- extract(model$fit, pars = c('lambda', 'nu'), ...)
  lambda_q        <- quantile(posts$lambda, probs = c(alpha/2, 1-alpha/2))
  nu_q            <- quantile(posts$nu, probs = c(alpha/2, 1-alpha/2))
  
  out <- list(lambda = lambda_q, nu = nu_q)
  
  return(out)
  
}
