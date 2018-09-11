rm(list=ls())

library(mpath)
library(pscl)
library(zic)
data(docvisits)

source("AMAZonn_2.R")
source("ALasso.R")

#barplot(with(docvisits, table(docvisits)), ylab = "Frequency",
#        xlab = "Doctor office visits")

dat = read.table("racd3_asc.txt")
new.dat = cbind(dat[,c(1,2,4,5,10,11,12)], dat[,13]+dat[,14]) 
colnames(new.dat) = c("Sex","Age","Income","Pr.Insurance","HScore","ChCond1","ChCond2","Dvisits")

# Full ZIP Model With All Covariates

# m1 <- zeroinfl(Dvisits ~ . | ., data = new.dat, weights = NULL, dist = "poisson", model = T, y = T, x = F)
# summary(m1)
# cat("loglik of zero-inflated model", logLik(m1))
# cat("BIC of zero-inflated model", BIC(m1, k = log(dim(dat)[1])))
# cat("AIC of zero-inflated model", AIC(m1))


# LASSO Estimates

val = tab = NULL
t1 = proc.time()

fit.lasso = zipath(Dvisits ~ . | ., data = new.dat, family = "poisson",
                       nlambda = 100, lambda.zero.min.ratio = 0.001, maxit.em = 300,
                       maxit.theta = 25, theta.fixed = FALSE, trace = FALSE,
                       penalty = "enet", rescale = FALSE)

pros.time = (proc.time() - t1)[1]
minBic = which.min(BIC(fit.lasso))
coef(fit.lasso, minBic)
cat("theta estimate", fit.lasso$theta[minBic])

#Compute standard errors of coefficients and theta (the last one for theta).

#se(fit.lasso, minBic, log = FALSE)

# Compute AIC, BIC, log-likelihood values of the selected model.

val = c(val, AIC(fit.lasso)[minBic])
val = c(val, BIC(fit.lasso)[minBic])
val = c(val, logLik(fit.lasso)[minBic])
val = c(val, pros.time)

tab = rbind(tab, val)

# Adaptive LASSO Estimates

val = NULL
t1 = proc.time()

fit.Alasso = ALasso(Dvisits ~ . | ., data = new.dat, family = "poisson",
                       nlambda = 100, lambda.zero.min.ratio = 0.001, maxit.em = 300,
                       maxit.theta = 25, theta.fixed = FALSE, trace = FALSE,
                       penalty = "enet", rescale = FALSE)

pros.time = (proc.time() - t1)[1]
rm(list="param") 

minBic = which.min(BIC(fit.Alasso))
coef(fit.Alasso, minBic)
cat("theta estimate", fit.Alasso$theta[minBic])

#Compute standard errors of coefficients and theta (the last one for theta).

#se(fit.Alasso, minBic, log = FALSE)

# Compute AIC, BIC, log-likelihood values of the selected model.

val = c(val, AIC(fit.Alasso)[minBic])
val = c(val, BIC(fit.Alasso)[minBic])
val = c(val, logLik(fit.Alasso)[minBic])
val = c(val, pros.time)

tab = rbind(tab, val)

# MCP Estimates

val = NULL
t1 = proc.time()

tmp = zipath(Dvisits ~ . | ., data = new.dat, family = "poisson",
gamma.count = 2.7, gamma.zero = 2.7, lambda.zero.min.ratio = 0.1,
maxit = 1, maxit.em = 1, maxit.theta = 2, theta.fixed = FALSE,
penalty = "mnet")

fit.mcp = zipath(Dvisits ~ . | ., data = new.dat, family = "poisson",
gamma.count = 2.7, gamma.zero = 2.7, lambda.count = tmp$lambda.count[1:30],
lambda.zero = tmp$lambda.zero[1:30], maxit.em = 300,
maxit.theta = 25, theta.fixed = FALSE, penalty = "mnet")

pros.time = (proc.time() - t1)[1]
minBic = which.min(BIC(fit.mcp))
coef(fit.mcp, minBic)
cat("theta estimate", fit.mcp$theta[minBic])


# Compute standard errors of coefficients and theta (the last one for theta).
#se(fit.mcp, minBic, log = FALSE)

#Compute AIC, BIC, log-likelihood values of the selected model.
val = c(val, AIC(fit.mcp)[minBic])
val = c(val, BIC(fit.mcp)[minBic])
val = c(val, logLik(fit.mcp)[minBic])
val = c(val, pros.time)

tab = rbind(tab, val)

#Compute SCAD estimates.

val = NULL
t1 = proc.time()

tmp = zipath(Dvisits ~ . | ., data = new.dat, family = "poisson",
gamma.count = 2.5, gamma.zero = 2.5, lambda.zero.min.ratio = 0.01,
maxit = 1, maxit.em = 1, maxit.theta = 2, theta.fixed = FALSE,
penalty = "snet")

fit.scad = zipath(Dvisits ~ . | ., data = new.dat, family = "poisson",
gamma.count = 2.5, gamma.zero = 2.5, lambda.count = tmp$lambda.count[1:30],
lambda.zero = tmp$lambda.zero[1:30], maxit.em = 300,
maxit.theta = 25, theta.fixed = FALSE, penalty = "snet")

#Estimated coefficient parameters with smallest BIC value.

pros.time = (proc.time() - t1)[1]
minBic = which.min(BIC(fit.scad))
coef(fit.scad, minBic)
cat("theta estimate", fit.scad$theta[minBic])

#Compute standard errors of coefficients and theta (the last one for theta).
#se(fit.scad, minBic, log = FALSE)

#Compute AIC, BIC, log-likelihood values of the selected model.
val = c(val, AIC(fit.scad)[minBic])
val = c(val, BIC(fit.scad)[minBic])
val = c(val, logLik(fit.scad)[minBic])
val = c(val, pros.time)

tab = rbind(tab, val)

# AMAZonn Estimates

val = NULL
t1 = proc.time()

fit.zonn = AMAZonn(Dvisits ~ . | ., data = new.dat, family = "poisson",
                       nlambda = 100, lambda.zero.min.ratio = 0.001, maxit.em = 300,
                       maxit.theta = 25, theta.fixed = FALSE, trace = FALSE,
                       penalty = "enet", rescale = FALSE)
                       
pros.time = (proc.time() - t1)[1]
rm(list="param") 

minBic = which.min(BIC(fit.zonn))
coef(fit.zonn, minBic)
cat("theta estimate", fit.zonn$theta[minBic])

#Compute standard errors of coefficients and theta (the last one for theta).

#se(fit.zonn, minBic, log = FALSE)

# Compute AIC, BIC, log-likelihood values of the selected model.

val = c(val, AIC(fit.zonn)[minBic])
val = c(val, BIC(fit.zonn)[minBic])
val = c(val, logLik(fit.zonn)[minBic])
val = c(val, pros.time)

tab = rbind(tab, val)

colnames(tab) = c("AIC","BIC","Loglikelihood","Preocess.Time")
rownames(tab) = c("Lasso","ALasso","MCP","SCAD","AMAZonn")

write.csv(tab,"racd_Results_ZIP.csv")