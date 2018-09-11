rm(list=ls())

library(mpath)
library(pscl)
library(zic)
data(docvisits)

source("AMAZonn_2.R")
source("ALasso.R")

#barplot(with(docvisits, table(docvisits)), ylab = "Frequency",
#        xlab = "Doctor office visits")

dt <- docvisits[, -(2:3)]
tmp <- model.matrix(~age30 * health + age35 * health +
                         age40 * health + age45 * health + age50 * health +
                         age55 * health + age60 * health, data = dt)[, -(1:9)]
dat <- cbind(dt, tmp)

# Full ZINB Model With All Covariates

m1 <- zeroinfl(docvisits ~ . | ., data = dat, weights = NULL, dist = "negbin", model = T, y = T, x = F)
summary(m1)
cat("loglik of zero-inflated model", logLik(m1))
cat("BIC of zero-inflated model", AIC(m1, k = log(dim(dat)[1])))
cat("AIC of zero-inflated model", AIC(m1))


# LASSO Estimates

fit.lasso <- zipath(docvisits ~ . | ., data = dat, family = "negbin",
                       nlambda = 100, lambda.zero.min.ratio = 0.001, maxit.em = 300,
                       maxit.theta = 25, theta.fixed = FALSE, trace = FALSE,
                       penalty = "enet", rescale = FALSE)

minBic <- which.min(BIC(fit.lasso))
coef(fit.lasso, minBic)
cat("theta estimate", fit.lasso$theta[minBic])

#Compute standard errors of coecients and theta (the last one for theta).

se(fit.lasso, minBic, log = FALSE)

# Compute AIC, BIC, log-likelihood values of the selected model.

AIC(fit.lasso)[minBic]
BIC(fit.lasso)[minBic]
logLik(fit.lasso)[minBic]


# ADAPTIVE LASSO Estimates

fit.Alasso <- ALasso(docvisits ~ . | ., data = dat, family = "negbin",
                       nlambda = 100, lambda.zero.min.ratio = 0.001, maxit.em = 300,
                       maxit.theta = 25, theta.fixed = FALSE, trace = FALSE,
                       penalty = "enet", rescale = FALSE)

rm(list="param") 

minBic <- which.min(BIC(fit.Alasso))
coef(fit.Alasso, minBic)
cat("theta estimate", fit.Alasso$theta[minBic])

#Compute standard errors of coecients and theta (the last one for theta).

se(fit.Alasso, minBic, log = FALSE)

# Compute AIC, BIC, log-likelihood values of the selected model.

AIC(fit.Alasso)[minBic]
BIC(fit.Alasso)[minBic]
logLik(fit.Alasso)[minBic]


# AMAZonn Estimates

fit.zonn <- AMAZonn(docvisits ~ . | ., data = dat, family = "negbin",
                       nlambda = 100, lambda.zero.min.ratio = 0.001, maxit.em = 300,
                       maxit.theta = 25, theta.fixed = FALSE, trace = FALSE,
                       penalty = "enet", rescale = FALSE)
                       
rm(list="param") 

minBic <- which.min(BIC(fit.zonn))
coef(fit.zonn, minBic)
cat("theta estimate", fit.zonn$theta[minBic])

#Compute standard errors of coecients and theta (the last one for theta).

se(fit.zonn, minBic, log = FALSE)

# Compute AIC, BIC, log-likelihood values of the selected model.

AIC(fit.zonn)[minBic]
BIC(fit.zonn)[minBic]
logLik(fit.zonn)[minBic]