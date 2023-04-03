#####Problem_01:
#Done

#=========================================================================================================
#####Problem_02:

GMMoneDim <- function(x, C)
{
  
  # Starting values
  pis <- runif(C-1, min = 0.2, max = 0.4)
  pis <- c(pis, 1-sum(pis))
  mu <- seq(min(x), max(x), length.out = C)
  sig2 <- rep(var(x) / C, C)
  
  diff <- 100
  tol <- 1e-5
  iter <- 0
  
  current <- c(pis, mu, sig2)
  
  while(diff>tol)
  {
    iter <- iter + 1
    
    #calculating gamma_ick
    gamma_prob <- matrix(0, nrow = length(x), ncol = C)
    for(c in 1:C)
    {
      gamma_prob[ ,c]  <- dnorm(x, mean = mu[c], sd = sqrt(sig2[c]))* pis[c]
    }
    Ep <- gamma_prob/(rowSums(gamma_prob))
    
    #updating parameters
    pis <- colMeans(Ep)
    mu <- colSums(Ep*x) / colSums(Ep)
    for(c in 1:C)
    {
      sig2[c] <- sum(Ep[,c]*(x - mu[c])^2) / sum(Ep[,c])
    }
    
    update <- c(pis, mu, sig2)
    
    diff <- norm(current - update, type = "2")
    current <- update
  }
  
  #calculating zi's
  gamma_prob <- matrix(0, nrow = length(x), ncol = C)
  for(c in 1:C)
  {
    gamma_prob[ ,c]  <- dnorm(x, mean = mu[c], sd = sqrt(sig2[c]))* pis[c]
  }
  Ep <- gamma_prob/(rowSums(gamma_prob))
  
  Zs <- numeric(length = length(x))
  for (i in 1:length(x))
  {
    Zs[i] <- which(Ep[i,] == max(Ep[i,]))
  }
  
  rtn <- list(Zs, mu, sig2, pis)
}

#rm(list = ls())
data(faithful)
head(faithful)

x <- faithful$eruptions

##for C= 2
erup2 <- GMMoneDim(x,2)
erup2

hist(x, breaks = 30, main = "Eruptions", freq = FALSE)

# add colored points
points(x, rep(0, length(x)), pch = 16, col = erup2[[1]])

##for C= 3
erup3 <- GMMoneDim(x,3)
erup3

hist(x, breaks = 30, main = "Eruptions", freq = FALSE)

# add colored points
points(x, rep(0, length(x)), pch = 16, col = erup3[[1]])


#=========================================================================================================
#####Problem_03:

log.like <- function(x, mu, sig2, pis)
{
  n <- length(x)
  C <- length(mu)
  track <- 0
  for(i in 1:n)
  {
    track2 <- 0
    for (c in 1:C)
    {
      track2 = track2 + pis[c] * dnorm(x[i], mu[c],sqrt(sig2[c]))
    }
    track <- track + log(track2)
  }
  return(track)
}

BIC <- function(x, C)
{
  #C = 2
  output <- GMMoneDim(x,C)
  mu <- output[[2]]
  sig2 <- output[[3]]
  pis <- output[[4]]
  lghd <- log.like(x, mu, sig2, pis)
  p <- 1
  K <- C*p^2 + p*C + C - 1
  n <- length(x)
  ret <- lghd - K * log(n)
  return(ret)
}

BIC(x,2)
BIC(x,3)
BIC(x,4)

#So, BIC is highest for C = 2. So, C = 2 is the best choice.


#=========================================================================================================
#####Problem_04:

#rm(list = ls())
data(faithful)
head(faithful)

x <- faithful$waiting

hist(x, breaks = 30, main = "Waiting")

#from histogram it looks like there are two groups. But we'll check it by BIC.

BIC(x,2)
BIC(x,3)
BIC(x,4)

#So here also C = 2 is best. 

wait <- GMMoneDim(x,2)
wait

hist(x, breaks = 30, main = "Eruptions", freq = FALSE)

# add colored points
points(x, rep(0, length(x)), pch = 16, col = wait[[1]])

#=========================================================================================================
#####Problem_05:
#understood

#=========================================================================================================
#####Problem_06:



GMMforD <- function(x, C, pis, mu, sig2)
{
  p <- dim(x)[2]
  
  diff <- 100
  tol <- 1e-5
  iter <- 0
  
  #setting current as a vector
  current <- c(pis, as.vector(mu), as.vector(sig2))
  
  #while loop
  while(diff>tol)
  {
    iter <- iter + 1
    
    library(mvtnorm)
    
    #calculating gamma_ick
    gamma_prob <- matrix(0, nrow = nrow(x), ncol = C)
    for(c in 1:C)
    {
      gamma_prob[ ,c]  <- dmvnorm(x, mean = mu[,c], sigma = sig2[, , c])* pis[c]
    }
    Ep <- gamma_prob/(rowSums(gamma_prob))
    
    #updating parameters
    pis <- colMeans(Ep)
    
    mu <- matrix(0, ncol = C, nrow = p)
    for (c in 1:C)
    {
      mu[ ,c] <- colSums(Ep[,c]*x) / colSums(Ep)[c]
    }

    sig2 <- array(0, dim = c(p,p,C))
    for (c in 1:C)
    {
      sig2.c <- matrix(0, nrow = 2, ncol = 2)
      for (i in 1:nrow(x))
      {
        sig2.c <- sig2.c + Ep[i,c] * (x[i,] - mu[,c]) %*% t((x[i,] - mu[,c]))
      }
      sig2[, , c] <- sig2.c / colSums(Ep)[c] + + diag(0.01, p)
    }
    
    
    #setting current as a vector
    update <- c(pis, as.vector(mu), as.vector(sig2))
    
    #calculating the gap between current and update
    diff <- norm(current - update, type = "2")
    current <- update
    
    #print(paste0("We are at iteration : ", iter))
    
  }
  
  #calculating zi's
  gamma_prob <- matrix(0, nrow = nrow(x), ncol = C)
  for(c in 1:C)
  {
    gamma_prob[ ,c]  <- dmvnorm(x, mean = mu[,c], sigma = sig2[, , c])* pis[c]
  }
  Ep <- gamma_prob/(rowSums(gamma_prob))
  
  Zs <- numeric(length = nrow(x))
  for (i in 1:nrow(x))
  {
    Zs[i] <- which(Ep[i,] == max(Ep[i,]))
  }
  
  rtn <- list(Zs, mu, sig2, pis, iter)
  return(rtn)
}

#=========================================================================================================
#####Problem_07:

data(faithful)
head(faithful)

x <- as.matrix(faithful)
C = 2

#starting from different starting values
#trial 1:
p <- dim(x)[2]

pis <- runif(C-1, min = 0.2, max = 0.4)
pis <- c(pis, 1-sum(pis))

mu <- matrix(0, ncol = C, nrow = p)
for (i in 1:p)
{
  mu[i,] <- seq(min(x[,i]), max(x[,i]), length.out = C)
}

sig2 <- array(0, dim = c(p,p,C))
for (i in 1:p)
{
  for (j in 1:p)
  {
    sig2[i,j, ] <- rep(cov(x[,i],x[,j]) / C, C)
  }
}

output1 <- GMMforD(x, C, pis, mu, sig2)

color <- output1[[1]] + 1
plot(x, col = color, pch = 20)

#trial 2:
pis_1 <- c(0.4, 0.6)

mu_1 <- matrix(c(1,45,6,90), ncol = 2)

sig2_1 <- array(c(0.6,7,7,95,
                0.7,7.5,7.5,90), dim = c(2,2,2))

output2 <- GMMforD(x, C, pis = pis_1, mu = mu_1, sig2 = sig2_1)

color <- output1[[1]] + 2
plot(x, col = color, pch = 20)

#we need more trials but lyaad lagche.

#=========================================================================================================
#####Problem_08:

library(mvtnorm)

log.like <- function(x, mu, sig2, pis)
{
  n <- nrow(x)
  C <- ncol(mu)
  track <- 0
  for(i in 1:n)
  {
    track2 <- 0
    for (c in 1:C)
    {
      track2 = track2 + pis[c] * dmvnorm(x[i, ], mu[ ,c], sig2[, , c])
    }
    track <- track + log(track2)
  }
  return(track)
}

BIC <- function(x, C)
{
  #C = 2
  output <- GMMforD(x,C)
  mu <- output[[2]]
  sig2 <- output[[3]]
  pis <- output[[4]]
  lghd <- log.like(x, mu, sig2, pis)
  p <- ncol(x)
  K <- C*p^2 + p*C + C - 1
  n <- nrow(x)
  ret <- lghd - K * log(n)
  return(ret)
}

data(faithful)
head(faithful)

x <- as.matrix(faithful)

BIC(x,2)
BIC(x,3)
BIC(x,4)

#so best value of C is 2. 

