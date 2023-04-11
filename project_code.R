library(readr)
library(sparklyr)
library(GGally)
library(nnet)
library(MASS)
library(olsrr)
library(dplyr)
require(mlbench)
require(mvtnorm)
require(coda)
library(tidyr)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
wine <- read_delim("data/winequality-red.csv", 
                   ";", escape_double = FALSE, trim_ws = TRUE)

ggplot(wine) +
  geom_bar(aes(quality ),alpha = 0.4, col='blue',fill= 'blue')+
  xlab("Quality score") + ggtitle("Red wine")+ ylab("Count")+ scale_x_continuous(limits = c(1, 9), breaks=1:9)

set.seed(100)

n <- nrow(wine)
index <- 1:n
train.ind <- sample( index, round(n*3/4) )
test.ind <- setdiff( index, c(train.ind) )

test.wine <- wine [test.ind, ]
train.wine <- wine [train.ind, ]

normalizer <- function(vector) {
  (vector - min(vector))/ ( max(vector) - min(vector) )
}

train.wine[,-12] <- apply( train.wine[,-12], 2, normalizer)
test.wine[,-12] <- apply( test.wine[,-12], 2, normalizer)

# Basic LM model coeffs 
reg_0 <- lm( quality ~ .-1, data=train.wine) 
summary(reg_0)
# 0.638

predictors <- c( "volatile acidity", "chlorides", "total sulfur dioxide", "sulphates", "alcohol" )

var.names <- c('intercept',"fix.acidity","vol.acidity","citric.acid","res.sugar","chlorides",           
               "fr.sulf.diox","tot.sulf.diox","density","pH","sulphates",           
               "alcohol")

#colnames(train.wine) <- c("fix.acidity","vol.acidity","citric.acid","res.sugar","chlorides",           
#                          "fr.sulf.diox","tot.sulf.diox","density","pH","sulphates",           
#                          "alcohol", 'quality')
#ggcorr(train.wine, label=TRUE, label_size = 3, legend.position = "none")

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#   linear model with ALL variables RED
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

y <- as.numeric(train.wine$quality)

# DISIGN MATRIX
X <- train.wine[,-12]
X <- matrix(unlist(X), ncol = dim(X)[2], byrow = FALSE) # MAKE A MATRIX OUT OF DF 
X <- cbind(rep(1,nrow(X)),X) # plus interseption
summary(X)

n <- length(y)
k <- ncol(X)

num.pars <- k + 1
num.chains <- 3

# LOG PRIOR
log.prior <- function(par.vector) {
  beta <- as.matrix(par.vector[1:k])
  sigma2 <- par.vector[num.pars]
  if (sigma2 <= 0)
    return(log(0))
  dens <- (
    sum(dnorm(abs(beta[2:k]),0,sqrt(sigma2),log=TRUE)) +
      dgamma(1/sigma2,4,90,log=TRUE))
  return(dens)
}

# LOG LIKELIHOOD
log.likelihood <- function(y,X,par.vector) {
  beta <- as.matrix(par.vector[1:k])
  sigma2 <- par.vector[num.pars]
  if (sigma2 <= 0)
    return(log(0))
  dens <- sum(dnorm(y,t(X%*%beta),sqrt(sigma2),log=TRUE))
  return(dens)
}

dd <- 300000

sample <- array(0,dim=c(dd,num.pars,num.chains))

for (i in (1:num.chains)) { 
  sample[1,num.pars,i] <- 1/rgamma(1,4,90)       # sigma2
  sample[1,1,i] <- mean(y)                       # beta0
  sample[1,(2:k),i] <- rep(0,k-1)                # beta rest  0
}

acceptance <- rep(1,num.chains)

for (step in 2:dd) {
  for (chain in 1:num.chains) {
    proposed <- rnorm(num.pars,sample[step-1,,chain],0.03)
    r <- log.likelihood(y,X,proposed)+
      log.prior(proposed) -
      log.likelihood(y,X,sample[step-1,,chain]) -
      log.prior(sample[step-1,,chain])
    u <- runif(1)
    if (log(u) < r) {
      sample[step,,chain] <- proposed
      acceptance[chain] <- acceptance[chain]+1
    } else {
      sample[step,,chain] <- sample[step-1,,chain]
    }
  }
}
final.sample <- rbind(sample[(dd-20000):dd,-13,1],sample[(dd-20000),-13,2])
print(acceptance/dd*100)                              # Acceptance rates in percent


plot(sample[,1,1], lty = 1, type='l', xlab="step", ylab="Beta0",col = "red") 
lines(sample[,1,2], col = "blue")
lines(sample[,1,3], col = "green")

plot(sample[,2,1], lty = 1, type='l', xlab="step", ylab="Beta1",col = "red") 
lines(sample[,2,2], col = "blue")
lines(sample[,2,3], col = "green")

plot(sample[,12,1], lty = 1, type='l', xlab="step", ylab="Beta10",col = "red") 
lines(sample[,12,2], col = "blue")
lines(sample[,12,3], col = "green")


plot(density(final.sample[,1]),main = var.names[1], ylab="Density", xlab="value")  
plot(density(final.sample[,3]),main = var.names[3], ylab="Density", xlab="value")   
plot(density(final.sample[,5]),main = var.names[5], ylab="Density", xlab="value")
plot(density(final.sample[,7]),main = var.names[7], ylab="Density", xlab="value")
plot(density(final.sample[,11]),main = var.names[11], ylab="Density", xlab="value")
plot(density(final.sample[,12]),main = var.names[12], ylab="Density", xlab="value")  

boxplot(final.sample, las = 2, names = var.names) + abline(h=0, lty=2)

coefs_0

coefs_0_w <- apply(final.sample,2,mean)
coefs_0_w

# prediction

x_test <- as.matrix(test.wine[,-12])
x_test <- cbind(rep(1,nrow(test.wine)), x_test)
y_test <- as.numeric(unlist(test.wine[,12]))
y_test

predict_set <- final.sample %*% t(x_test)

quant <- apply(predict_set,2,quantile, probs=c(0.025,0.5,0.975) )
quant[2,] <- apply(predict_set,2, mean)

estim_cat <- round(quant)

# model accuracy (mean)
sum(estim_cat[2,] == y_test )/ length(y_test)

# HIT of the interval
hit <- 0
for (i in 1:length(y_test)) {
  if (y_test[i] %in% c(estim_cat[1,i],estim_cat[3,i]) ) {
    hit<- hit+1
  }
  else 
    hit<-hit
}
hit/ length(y_test)


plot_fr <- as.data.frame( cbind(t(estim_cat), y_test  )   )
colnames(plot_fr) <- c("lower",'mean',"upper","true")

ggplot(plot_fr) +
  geom_bar(aes(true),alpha = 0.4, col='blue',fill= 'blue')+
  xlab("Quality score") + ggtitle("True score/Prediction. Mean Beta. Red wine")+ ylab("Count")+
  geom_bar(aes(mean),alpha = 0.2, col='red',fill= 'red')


plot_fr_num <-  as.data.frame( cbind(t(quant), y_test  )   )
plot_fr_num <- apply(plot_fr_num, 2, as.numeric)
colnames(plot_fr_num) <- c("lower",'mean',"upper","true")
plot_fr_num <-  as.data.frame(plot_fr_num)

ggplot(plot_fr_num) +
  geom_density(aes(lower), col='green')+
  geom_density(aes(upper), col='darkgreen')+
  geom_density(aes(mean), col='red')+
  geom_density(aes(true), col='blue')+
  xlab("Quality score") + ggtitle("True score/Prediction. Red wine")+ ylab("Density")+
  xlim(c(2,10))

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  NON linear model with significant coeffs
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

summary(reg_0)
train.wine$quality <- as.factor(train.wine$quality)

# WHITE
plot(train.wine$quality,train.wine$alcohol, ylab="Alcohol", xlab="Quality category") # YES
plot(train.wine$quality,train.wine$pH, ylab="pH", xlab="Quality category") # NO
plot(train.wine$quality,train.wine$`fixed acidity`, ylab="fixed acidity", xlab="Quality category")
plot(train.wine$quality,train.wine$`volatile acidity`, ylab="Volatile acidity", xlab="Quality category")

c("Alcohol","density","pH","Volatile acidity")

plot(train.wine$quality,train.wine$`residual sugar`*train.wine$alcohol, ylab="", xlab="Quality category")

train.wine_log <- train.wine
train.wine_log$quality <- as.factor(train.wine_log$quality)

train.wine_log$alcohol <- log(train.wine$alcohol)
train.wine_log$alcohol[is.infinite(train.wine_log$alcohol)] <- 0


train.wine_log$`residual sugar` <- log(train.wine$`residual sugar`)
train.wine_log$`residual sugar`[is.infinite(train.wine_log$`residual sugar`)] <- 0

train.wine_log$density <- log(train.wine$density)
train.wine_log$density[is.infinite(train.wine_log$density)] <- 0

train.wine_log$`volatile acidity` <- log(train.wine$`volatile acidity`)
train.wine_log$`volatile acidity`[is.infinite(train.wine_log$`volatile acidity`)] <- 0

train.wine_log$`fixed acidity` <- log(train.wine$`fixed acidity`)
train.wine_log$`fixed acidity`[is.infinite(train.wine_log$`fixed acidity`)] <- 0

train.wine_log$pH <- log(train.wine$pH)
train.wine_log$pH[is.infinite(train.wine_log$pH)] <- 0

# FINAL
train.wine_log$quality <- as.numeric(train.wine$quality)
reg_1_w <- lm( quality ~ alcohol*`residual sugar` +`volatile acidity`- 1,
               data=train.wine_log) 
summary(reg_1_w)

# DISIGN MATRIX FOR WHITE

X_n <- train.wine %>% transmute(log_alcohol=log(alcohol), log_sugar=log(`residual sugar`),log_va=log(`volatile acidity`),
                                log_coc= (log(alcohol)*log(`residual sugar`)))

X_n$log_alcohol[is.infinite(X_n$log_alcohol)] <- 0
X_n$log_sugar[is.infinite(X_n$log_sugar)] <- 0
X_n$log_va[is.infinite(X_n$log_va)] <- 0
X_n$log_coc[is.infinite(X_n$log_coc)] <- 0


# RED
train.wine$quality <- as.factor(train.wine$quality)
plot(train.wine$quality,train.wine$alcohol, ylab="Alcohol", xlab="Quality category") 
plot(train.wine$quality,train.wine$`volatile acidity`, ylab="Volatile acidity", xlab="Quality category")
plot(train.wine$quality,train.wine$sulphates, ylab="Sulphates", xlab="Quality category")  

plot(train.wine$quality,train.wine$chlorides/train.wine$`total sulfur dioxide`, ylab="Chlorides/Total sulfur dioxide", xlab="Quality category")

train.wine$quality <- as.numeric(train.wine$quality)

train.wine_log <- train.wine

train.wine_log$alcohol <- log(train.wine$alcohol)
train.wine_log$alcohol[is.infinite(train.wine_log$alcohol)] <- 0

train.wine_log$sulphates <- log(train.wine$sulphates)
train.wine_log$sulphates[is.infinite(train.wine_log$sulphates)] <- 0

train.wine_log$`volatile acidity` <- log(train.wine$`volatile acidity`)
train.wine_log$`volatile acidity`[is.infinite(train.wine_log$`volatile acidity`)] <- 0

train.wine_log$`total sulfur dioxide` <- log(train.wine$`total sulfur dioxide`)
train.wine_log$`total sulfur dioxide`[is.infinite(train.wine_log$`total sulfur dioxide`)] <- 0

train.wine_log$chlorides <- log(train.wine$chlorides)
train.wine_log$chlorides[is.infinite(train.wine_log$chlorides)] <- 0

train.wine_log$`fixed acidity` <- log(train.wine$`fixed acidity`)
train.wine_log$`fixed acidity`[is.infinite(train.wine_log$`fixed acidity`)] <- 0

reg_1_r <- lm( quality ~ alcohol + `volatile acidity` + sulphates + chlorides - 1,
               data=train.wine_log) 
summary(reg_1_r)


# DISIGN MATRIX FOR RED

X_n <- train.wine %>% transmute(log_alcohol=log(alcohol), log_sulphates=log(sulphates),
                                log_va=log(`volatile acidity`),log_chlor=log(chlorides))
X_n$log_chlor[is.infinite(X_n$log_chlor)] <- 0
X_n$log_alcohol[is.infinite(X_n$log_alcohol)] <- 0
X_n$log_sulphates[is.infinite(X_n$log_sulphates)] <- 0
X_n$log_va[is.infinite(X_n$log_va)] <- 0


var.names_n <- colnames(X_n)
X_n <- matrix(unlist(X_n), ncol = dim(X_n)[2], byrow = FALSE) # MAKE A MATRIX OUT OF DF 
summary(X_n)

y_n <- as.numeric(train.wine$quality)


n <- length(y_n)
k <- ncol(X_n)

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  NON linear model with significant coeffs
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# BAYESSIAN RIDGE

num.pars <- k + k + 1
num.chains <- 3

log.prior_r <- function(par.vector) {
  beta <- as.matrix(par.vector[1:k])
  alpha <- par.vector[(k+1):(2*k)]
  sigma2 <- par.vector[num.pars]
  if ( any(alpha <= 0 ))
    return(log(0))
  if ( sigma2 <= 0 )
    return(log(0))
  dens <- (
    sum(dnorm(abs(beta[1:k]),0,sqrt(sigma2/alpha)),log=TRUE) +
      dgamma(1/sigma2,4,60,log=TRUE) +
      sum(dunif(alpha,0,50,log=TRUE))
  )
  return(dens)
}

hod_r <- array(0,dim=c(2000,num.pars,num.chains))
for (i in (1:num.chains)) { 
  hod_r[2000,num.pars,i] <- 1/rgamma(1,4,60)    # sigma2
  hod_r[2000,1:k,i] <- rep(0,k)                 # beta 0
  hod_r[2000,(k+1):(2*k),i] <- rep(1,k)         # alpha     1
}  

hod_r[2000,,1]

stnum <- 0
converged <- FALSE
while (!converged) {
  hod_r[1,,] <- hod_r[2000,,]
  for (step in 2:2000) {
    for (chain in 1:num.chains) {
      proposed <- rnorm(num.pars,hod_r[step-1,,chain],0.03)
      r <- log.likelihood(y_n,X_n,proposed)+
        log.prior_r(proposed) -
        log.likelihood(y_n,X_n,hod_r[step-1,,chain]) -
        log.prior_r(hod_r[step-1,,chain])
      u <- runif(1)
      if (log(u) < r) {
        hod_r[step,,chain] <- proposed
      } else {
        hod_r[step,,chain] <- hod_r[step-1,,chain]
      }
    }
    stnum <- stnum + 1 
  }
  chainlist <- mcmc.list(Chain1=mcmc(hod_r[,,1]), 
                         Chain2=mcmc(hod_r[,,2]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.3)
}

converged
stnum

# 9 million

# Once convergence is achieved,
# we begin collecting the final sample from the posterior
sample_r <- array(0,dim=c(20000,num.pars,num.chains))
sample_r[1,,] <- hod_r[2000,,]
acceptance <- rep(1,num.chains)

for (step in 2:20000) {
  for (chain in 1:num.chains) {
    proposed <- rnorm(num.pars,sample_r[step-1,,chain],0.03)
    r <- log.likelihood(y_n,X_n,proposed)+
      log.prior_r(proposed) -
      log.likelihood(y_n,X_n,sample_r[step-1,,chain]) -
      log.prior_r(sample_r[step-1,,chain])
    u <- runif(1)
    if (log(u) < r) {
      sample_r[step,,chain] <- proposed
      acceptance[chain] <- acceptance[chain]+1
    } else {
      sample_r[step,,chain] <- sample_r[step-1,,chain]
    }
  }
}

final.sample_rr <- rbind(sample_r[,1:k,1],sample_r[,1:k,2])
print(acceptance/200) # Acceptance rates in percent

coefs_r <- apply(final.sample_rr,2,mean)
coefs_r


plot(sample_r[,1,1], lty = 1, type='l', xlab="step", ylab="Beta0",col = "red") 
lines(sample_r[,1,2], col = "blue")
lines(sample_r[,1,3], col = "green")

plot(sample_r[,2,1], lty = 1, type='l', xlab="step", ylab="Beta1",col = "red") 
lines(sample_r[,2,2], col = "blue")


plot(density(final.sample_rr[,1]),main = "log(Alcohol)", ylab="Density", xlab="value")  
plot(density(final.sample_rr[,2]),main = "Log(Sugar)", ylab="Density", xlab="value")   
plot(density(final.sample_rr[,3]),main = "log(Volatile acidity)", ylab="Density", xlab="value")
plot(density(final.sample_rr[,4]),main = "log(Cocktail)", ylab="Density", xlab="value")


boxplot(final.sample_rr, las = 2, names = c("Log(Alcohol)","log(Sugar)","log(Volatile acidity)","log(Cocktail)"), 
        main="Significant Parameters (log transformation)") + abline(h=0, lty=2)

# prediction


#!!!!!!!!!!!!!!!!!!!!   WHITE
x_test <- test.wine %>% transmute(log_alcohol=log(alcohol), log_sugar=log(`residual sugar`),log_va=log(`volatile acidity`),
                                  log_coc= (log(alcohol)*log(`residual sugar`)))

x_test$log_alcohol[is.infinite(x_test$log_alcohol)] <- 0
x_test$log_sugar[is.infinite(x_test$log_sugar)] <- 0
x_test$log_va[is.infinite(x_test$log_va)] <- 0
x_test$log_coc[is.infinite(x_test$log_coc)] <- 0


#!!!!!!!!!!!!!!!!!!!!   RED
x_test  <- test.wine %>% transmute(log_alcohol=log(alcohol),log_pH=log(pH), log_va=log(`volatile acidity`),
                                   log_fix=log(`fixed acidity`), log_coc= (log(alcohol)*log(`residual sugar`)))

x_test $log_alcohol[is.infinite(x_test $log_alcohol)] <- 0
x_test $log_va[is.infinite(x_test $log_va)] <- 0
x_test $log_pH[is.infinite(x_test $log_pH)] <- 0
x_test $log_fix[is.infinite(x_test $log_fix)] <- 0
x_test $log_coc[is.infinite(x_test $log_coc)] <- 0
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!



x_test <- as.matrix(x_test)

y_test <- as.numeric(unlist(test.wine[,12]))
y_test

predict_set_red_n <- final.sample_rr %*% t(x_test)

quant <- apply(predict_set_red_n,2,quantile, probs=c(0.025,0.5,0.975) )
quant[2,] <- apply(predict_set_red_n,2, mean)

estim_cat <- round(quant)

# model accuracy (mean)
sum(estim_cat[2,] == y_test )/ length(y_test)

# HIT of the interval
hit <- 0
for (i in 1:length(y_test)) {
  if (y_test[i] %in% c(estim_cat[1,i],estim_cat[3,i]) ) {
    hit<- hit+1
  }
  else 
    hit<-hit
}
hit/ length(y_test)


plot_fr <- as.data.frame( cbind(t(estim_cat), y_test  )   )
colnames(plot_fr) <- c("lower",'mean',"upper","true")

ggplot(plot_fr) +
  geom_bar(aes(true),alpha = 0.4, col='blue',fill= 'blue')+
  xlab("Quality score") + ggtitle("True score/Prediction. Mean Beta")+ ylab("Count")+
  geom_bar(aes(mean),alpha = 0.2, col='red',fill= 'red')


plot_fr_num <-  as.data.frame( cbind(t(quant), y_test  )   )
plot_fr_num <- apply(plot_fr_num, 2, as.numeric)
colnames(plot_fr_num) <- c("lower",'mean',"upper","true")
plot_fr_num <-  as.data.frame(plot_fr_num)

ggplot(plot_fr_num) +
  geom_density(aes(lower), col='green')+
  geom_density(aes(upper), col='darkgreen')+
  geom_density(aes(mean), col='red')+
  geom_density(aes(true), col='blue')+
  xlab("Quality score") + ggtitle("True score/Prediction")+ ylab("Density")+
  xlim(c(2,10))

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#  MULTINOMIAL REGRESSION  RED WHINE ONLY
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

train.mult_log <- train.wine

reg_m <- lm(quality ~ alcohol+sulphates+`volatile acidity`-1, data=train.mult )
summary(reg_m)

train.mult_log$alcohol <- log(train.mult$alcohol)
train.mult_log$alcohol[is.infinite(train.mult_log$alcohol)] <- 0

train.mult_log$sulphates <- log(train.mult$sulphates)
train.mult_log$sulphates[is.infinite(train.mult_log$sulphates)] <- 0

train.mult_log$`volatile acidity` <- log(train.mult$`volatile acidity`)
train.mult_log$`volatile acidity`[is.infinite(train.mult_log$`volatile acidity`)] <- 0



X_m <- train.mult_log[,c('alcohol',"quality")]          # log transformation on significant predictors
X_m <- matrix(unlist(X_m), ncol = dim(X_m)[2], byrow = FALSE)                           # MAKE A MATRIX OUT OF DF 
X_m <- cbind(rep(1,nrow(X_m)),X_m)                                                      # plus interseption


ncase <- nrow(X_m)
k <- ncol(X_m)-1

category <- c(4,5,6,7,8)                           # minus first one category
category0 <- c(3)

cnum <- length(category)

num.pars <- k*cnum + 1
num.chains <- 2

freq <- train.white %>% group_by(quality) %>% count() %>% transmute( num = n/ncase )

# LOG PRIOR
log.prior_m <- function(par.vector) {
  beta <- as.matrix(par.vector[seq(1,k*cnum,1)])
  sigma2 <- par.vector[num.pars]
  if (sigma2 <= 0)
    return(log(0))
  dens <- (
    sum(dnorm( abs(beta),0,sqrt(sigma2),log=TRUE)) +
      dgamma(1/sigma2,4,60,log=TRUE))
  return(dens)
}

# LOG LIKELIHOOD
log.likelihood_m <- function(M,par) {
  dens <- 0
  beta <- as.matrix(par[ 1:(k*cnum) ])
  link <- c()
  for (i in 1:cnum) {
    D <- M[,-(k+1)]
    link[i] <- exp( t( D %*% beta[((i-1)*k+1):((i-1)*k+k)]))
  }
  for (i in 1:cnum) {
    y <- M[,(k+1)]
    y[ y != category[i]] <- 0
    y[ y == category[i]] <- 1
    p <- link[i]/(sum(link))
    dens <- dens + sum(dbinom(y,1,p,log=TRUE))
  }
  return(dens)
}

hod_m <- array(0,dim=c(2000,num.pars,num.chains))
for (i in (1:num.chains)) { 
  hod_m[2000,num.pars,i] <- 1/rgamma(1,4,60)        # sigma2
  #  hod_m[2000,beta0_ind,i] <- fr$num             # beta0
  hod_m[2000,1:(k*cnum),i] <- rep(0,k*cnum)      # beta rest  0
}

hod_m[2000,,1]

stnum <- 0
converged <- FALSE
while (!converged) {
  hod_m[1,,] <- hod_m[2000,,]
  for (step in 2:2000) {
    for (chain in 1:num.chains) {
      proposed <- rnorm(num.pars,hod_m[step-1,,chain],0.03)
      r <- log.likelihood_m(X_m,proposed)+
        log.prior_m(proposed) -
        log.likelihood_m(X_m,hod_m[step-1,,chain]) -
        log.prior_m(hod_m[step-1,,chain])
      u <- runif(1)
      if (log(u) < r) {
        hod_m[step,,chain] <- proposed
      } else {
        hod_m[step,,chain] <- hod_m[step-1,,chain]
      }
    }
    stnum <- stnum + 1 
  }
  chainlist <- mcmc.list(Chain1=mcmc(hod_m[,,1]), 
                         Chain2=mcmc(hod_m[,,2]))
  converged <- all((gelman.diag(chainlist)$psrf[,2])<1.2)
}

converged
stnum

sample_m <- array(0,dim=c(20000,num.pars,num.chains))
sample_m[1,,] <- hod_m[2000,,]
acceptance <- rep(1,num.chains)

for (step in 2:20000) {
  for (chain in 1:num.chains) {
    proposed <- rnorm(num.pars,sample_m[step-1,,chain],0.03)
    r <- log.likelihood_m(X_m,proposed)+
      log.prior_m(proposed) -
      log.likelihood_m(X_m,sample_m[step-1,,chain]) -
      log.prior_m(sample_m[step-1,,chain])
    u <- runif(1)
    if (log(u) < r) {
      sample_m[step,,chain] <- proposed
      acceptance[chain] <- acceptance[chain]+1
    } else {
      sample_m[step,,chain] <- sample_m[step-1,,chain]
    }
  }
}
final.sample_m <- rbind(sample_m[,,1],sample_m[,,2])
print(acceptance/200) # Acceptance rates in percent

beta <- apply(final.sample_m[,1:(k*cnum)], 2, mean)
beta <- matrix(beta, nrow=cnum, ncol=k, byrow = TRUE)

predictors <- c('alcohol','sulphates',"volatile acidity")


prob_function <- function(case) {
  xx <- as.matrix(log(case[c('alcohol')]))
  xx <- cbind(rep(1,nrow(xx)),xx)  
  y_link <- exp(beta %*% t(xx))/(sum(exp(beta %*% t(xx))))
  y_link <- c((1-sum(y_link)),y_link)
  ind_choice <- match( max(y_link), y_link)
  choice <- ifelse( ind_choice == 1, category0, category[ind_choice-1] )
  return(choice)
}

case <- test.wine

y_case <- unlist(test.wine[,'quality'])

pred.test<- c()
for (i in 1:nrow(case)) { 
  pred.test[i] <- prob_function(case[i,])
}

result <- as.data.frame(cbind(pred.test,as.vector(y_case)))
colnames(result) <- c('predict','true')


ggplot(result) +
  geom_bar(aes(true),alpha = 0.4, col='blue',fill= 'blue')+
  xlab("Quality score") + ggtitle("True score/Prediction. Multinomial model")+ ylab("Count")+
  geom_bar(aes(predict),alpha = 0.2, col='red',fill= 'red')+ scale_x_continuous(limits = c(3, 8), breaks=1:7)

sum( pred.test==as.vector(y_case) )/nrow(case)

