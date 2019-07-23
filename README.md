# AMAZonn (A Multicollinearity-adjusted Adaptive LASSO for Zero-inflated Count Regression)

# Introduction
Algorithms for fitting standard error adjusted adaptive LASSO for both zero-inflated Poisson (ZIP) and zero-inflated negative binomial (ZINB) regression models. AMAZonn relies on an efficient coordinate descent algorithm embedded within an EM algorithm by imposing standard error adjusted adaptive $L_1$ penalties in both count and zero submodels of the corresponding zero-inflated mixture model. 

# Installation
`AMAZonn` can be installed using the following command (execute from within a fresh R session):
```r
install.packages("devtools")
devtools::install_github("himelmallick/AMAZonn")
library(AMAZonn)
```

# Basic Usage
```r
AMAZonn(formula, data, family, ...)
```
- **formula** symbolic description of the model, similar to that in `zeroinfl` function in the R package `pscl` with | to separate the count and zero submodels
- **data**:	an optional data frame (or object coercible by as.data.frame to a data frame) containing the variables in the model
- **family**: character specification of the count model family with options `poisson`, `negbin`, and `geometric` (a log link is always used)
- **...** other arguments that can be passed from the `zipath` function in the R package `mpath` ([Wang et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/26059498))


The function `AMAZonn` returns a list of following components:
- **coefficients**:	a list with elements `count` and `zero` containing the coefficients from the respective submodels
- **...** other components similar to`zipath` function in the R package `mpath` ([Wang et al., 2015](https://www.ncbi.nlm.nih.gov/pubmed/26059498))

# Examples
```r
library(zic)
data(docvisits)

dt <- docvisits[, -(2:3)]
tmp <- model.matrix(~age30 * health + age35 * health +
                         age40 * health + age45 * health + age50 * health +
                         age55 * health + age60 * health, data = dt)[, -(1:9)]
dat <- cbind(dt, tmp)
AMAZonn Estimates
fit.zonn <- AMAZonn(docvisits ~ . | ., data = dat, family = "negbin")

rm(list="param") 

minBic <- which.min(BIC(fit.zonn))
coef(fit.zonn, minBic)
cat("theta estimate", fit.zonn$theta[minBic])
Compute standard errors of coefficients and theta (the last one for theta).
se(fit.zonn, minBic, log = FALSE)
Compute AIC, BIC, log-likelihood values of the selected model.
AIC(fit.zonn)[minBic]
BIC(fit.zonn)[minBic]
logLik(fit.zonn)[minBic]
```

## References

Wang, Z., Ma, S., and Wang, C.Y. (2015). [Variable Selection for Zero-inflated and Overdispersed Data with Application to Health Care Demand in Germany](https://www.ncbi.nlm.nih.gov/pubmed/26059498). Biometrical Journal 57(5):867-884.

## Citation

Banerjee, P., Garai, B., Mallick, H., Chowdhury, S., Chatterjee, S. (2018). [A Note on the Adaptive LASSO for Zero-inflated Poisson Regression](https://www.hindawi.com/journals/jps/2018/2834183/). Journal of Probability and Statistics, 2834183.


