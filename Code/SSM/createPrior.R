# Priors

density <- function(par){
  
  d1 <- sum(dunif(par[1:9], min = lower[1:9], max = upper[1:9], log = T))
  
  d2 <- sum(dgamma(par[10:11], shape = 2, scale = 0.1, log = T))
  
  d3 <- sum(dgamma(par[12:14], shape = 2, scale = 0.1, log = T))
  
  return(d1 + d2 +d3)
  
  
}


sampler <- function(n = 1){
  
  d1 <- runif(n, lower[1], upper[1])
  d2 <- runif(n, lower[2], upper[2])
  d3 <- runif(n, lower[3], upper[3])
  d4 <- runif(n, lower[4], upper[4])
  d5 <- runif(n, lower[5], upper[5])
  d6 <- runif(n, lower[6], upper[6])
  d7 <- runif(n, lower[7], upper[7])
  d8 <- runif(n, lower[8], upper[8])
  d9 <- runif(n, lower[9], upper[9])
  
  d10 <- rgamma(n, shape = 2, scale = 0.1)
  d11 <- rgamma(n, shape = 2, scale = 0.1)
  
  d12 <- rgamma(n, shape = 2, scale = 0.1)
  d13 <- rgamma(n, shape = 2, scale = 0.1)
  d14 <- rgamma(n, shape = 2, scale = 0.1)
  
  
  return(cbind(d1,d2,d3,d4,d5,d6,d7,d8,d9,
               d10,d11,d12,d13,d14))
  
}



prior <- createPrior(density = density, sampler = sampler, 
                     lower = lower, upper = upper, best = best)

