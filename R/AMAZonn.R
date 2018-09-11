#' A Multicollinearity-adjusted Adaptive LASSO for Zero-inflated Count Regression
#'
#' @description A Multicollinearity-adjusted Adaptive LASSO for Zero-inflated Count Regression fits zero-inflated regression models for count data via penalized maximum likelihood.
#'
#' @param formula symbolic description of the model, see details.
#' @param data,subset,na.action arguments controlling formula processing via \eqn{model.frame}.
#' @param weights optional numeric vector of weights.
#' @param offset optional numeric vector with an a priori known component to be included in the linear predictor of the count model. See below for more information on offsets.
#' @param standardize Logical flag for x variable standardization, prior to fitting the model sequence. The coefficients are always returned on the original scale. Default is \eqn{standardize=TRUE}.
#' @param family character specification of count model family (a log link is always used).
#' @param link character specification of link function in the binary zero-inflation model (a binomial family is always used).
#' @param penalty penalty considered as one of \eqn{enet}, \eqn{mnet}, \eqn{snet}.
#' @param control a list of control arguments specified via \eqn{zeroinfl.control}.
#' @param start starting values for the parameters in the linear predictor.
#' @param model,y,x logicals. If TRUE the corresponding components of the fit (model frame, response, model matrix) are returned.
#' @param nlambda number of \eqn{lambda} value, default value is 100. The sequence may be truncated before \eqn{nlambda} is reached if a close to saturated model for the zero component is fitted.
#' @param lambda.count A user supplied \eqn{lambda.count} sequence. Typical usage is to have the program compute its own \eqn{lambda.count} and \eqn{lambda.zero} sequence based on \eqn{nlambda} and \eqn{lambda.min.ratio}.
#' @param lambda.zero A user supplied \eqn{lambda.zero} sequence.
#' @param penalty.factor.count,penalty.factor.zero These are numeric vectors with the same length as predictor variables. that multiply \eqn{lambda.count},\eqn{lambda.zero}, respectively, to allow differential shrinkage of coefficients. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is same shrinkage for all variables.
#' @param lambda.count.min.ratio,lambda.zero.min.ratio Smallest value for \eqn{lambda.count} and \eqn{lambda.zero}, respectively, as a fraction of \eqn{lambda.max}, the (data derived) entry value (i.e. the smallest value for which all coefficients are zero except the intercepts). Note, there is a closed formula for \eqn{lambda.max} for \eqn{penalty="enet"}. If \eqn{rescale=TRUE}, \eqn{lambda.max} is the same for \eqn{penalty="mnet"} or \eqn{"snet"}. Otherwise, some modifications are required. In the current implementation, for small gamma value, the square root of the computed \eqn{lambda.zero[1]} is used when \eqn{penalty="mnet"} or \eqn{"snet"}.
#' @param alpha.count The elastic net mixing parameter for the count part of model.
#' @param alpha.zero The elastic net mixing parameter for the zero part of model.
#' @param gamma.count The tuning parameter of the snet or mnet penalty for the count part of model.
#' @param gamma.zero The tuning parameter of the snet or mnet penalty for the zero part of model.
#' @param rescale logical value, if TRUE, adaptive rescaling
#' @param init.theta The initial value of theta for \eqn{family="negbin"}.
#' @param theta.fixed Logical value only used for \eqn{family="negbin"}. If TRUE, theta is not updated.
#' @param EM Using EM algorithm. Not implemented otherwise
#' @param convtype convergency type, default is for count component only for speedy computation
#' @param maxit.em Maximum number of EM algorithm
#' @param maxit Maximum number of coordinate descent algorithm
#' @param maxit.theta Maximum number of iterations for estimating theta scaling parameter if \eqn{family=" negbin"}. Default value maxit.theta may be increased, yet may slow the algorithm
#' @param eps.bino a lower bound of probabilities to be claimed as zero, for computing weights and related values when \eqn{family="binomial"}.
#' @param reltol Convergence criteria, default value \eqn{1e-5} may be reduced to make more accurate yet slow
#' @param shortlist logical value, if TRUE, limited results return
#' @param trace If TRUE, progress of algorithm is reported
#'
#' @details The algorithm fits penalized zero-inflated count data regression models using the coordinate descent algorithm within the EM algorithm. The returned fitted model object is of class \eqn{"AMAZonn"} and is similar to fitted \eqn{"glm"} and \eqn{"zeroinfl"} objects. For elements such as \eqn{"coefficients"} a list is returned with elements for the \eqn{zero} and \eqn{count} component, respectively.
#'
#' @return The function AMAZonn returns a list of following components
#' @return \item{coefficients}{a list with elements \eqn{"count"} and \eqn{"zero"} containing the coefficients from the respective models,}
#' @return \item{residuals}{a vector of raw residuals (observed - fitted),}
#' @return \item{fitted.values}{a vector of fitted means,}
#' @return \item{weights}{the case weights used,}
#' @return \item{terms}{a list with elements \eqn{"count"}, \eqn{"zero"} and \eqn{"full"} containing the terms objects for the respective models,}
#' @return \item{theta}{estimate of the additional \eqn{\theta} parameter of the negative binomial model (if a negative binomial regression is used),}
#' @return \item{loglik}{log-likelihood of the fitted model,}
#' @return \item{family}{character string describing the count distribution used,}
#' @return \item{link}{character string describing the link of the zero-inflation model,}
#' @return \item{linkinv}{the inverse link function corresponding to link,}
#' @return \item{converged}{logical value, TRUE indicating successful convergence of \eqn{AMAZonn}, FALSE indicating otherwise}
#' @return \item{call}{the original function call}
#' @return \item{formula}{the original formula}
#' @return \item{levels}{levels of the categorical regressors}
#' @return \item{contrasts}{a list with elements \eqn{"count"} and \eqn{"zero"} containing the contrasts corresponding to levels from the respective models,}
#' @return \item{model}{the full model frame \eqn{(if model = TRUE)},}
#' @return \item{y}{the response count vector \eqn{(if y = TRUE)},}
#' @return \item{x}{a list with elements \eqn{"count"} and \eqn{"zero"} containing the model matrices fromthe respective models \eqn{(if x = TRUE)},}
#'
#' @export
#'
#' @examples
#'
#' data("docvisits", package = "zic")
#' ## with simple inflation (no regressors for zero component)
#' Zonn_zip <- AMAZonn(docvisits ~ . | ., data = docvisits, family = "poisson", nlambda = 10)
#' summary(Zonn_zip)
#' Zonn_zinb <- AMAZonn(docvisits ~ . | ., data = docvisits, family = "negbin", nlambda=10)
#' summary(Zonn_zinb)


AMAZonn <- function(formula, data, weights, subset = NULL, na.action = na.omit, offset = NULL, standardize = TRUE, family = c("poisson", "negbin","geometric"), link = c("logit", "probit", "cloglog", "cauchit", "log"), control = zeroinfl.control(...), penalty = c("enet", "mnet", "snet"), start = NULL, model = TRUE,y = TRUE, x = FALSE, nlambda = 100, lambda.count = NULL, lambda.zero = NULL, penalty.factor.count = NULL, penalty.factor.zero = NULL, lambda.count.min.ratio = .0001, lambda.zero.min.ratio = .1, alpha.count = 1, alpha.zero = alpha.count, gamma.count = 3, gamma.zero = gamma.count, rescale=FALSE, init.theta, theta.fixed=FALSE, EM = TRUE, maxit.em = 200, convtype = c("count", "both"), maxit = 1000, maxit.theta = 1, reltol = 1e-5, eps.bino = 1e-5, shortlist = FALSE, trace = FALSE){

	### Initializing Parameters
	param = NULL

	if(missing(weights)){
		param$weights = as.vector(matrix(1,dim(data)[1]))
	}else{ param$weights = weights}

	if(missing(subset)){
		param$subset = as.vector(matrix(1:dim(data)[1],ncol=1))
	}else{ param$subset = subset}

	if(missing(standardize)){
		param$standardize = TRUE
	}else{ param$standardize = standardize}

	if(missing(family)){
		param$family = "negbin"
	}else{ param$family = family}

	if(missing(link)){
		param$link = "log"
	}else{ param$link = link}

	if(missing(penalty)){
		param$penalty = "enet"
	}else{ param$penalty = penalty}

	if(missing(model)){
		param$model = TRUE
	}else{ param$model = model}

	if(missing(y)){
		param$y = TRUE
	}else{ param$y = y}

	if(missing(x)){
		param$x = FALSE
	}else{ param$x = x}

	if(missing(nlambda)||nlambda < 0){
		param$nlambda = 100
	}else{ param$nlambda = nlambda}

	if(missing(lambda.count.min.ratio)||lambda.count.min.ratio<0){
    	param$lambda.count.min.ratio = 0.0001
  	}else{ param$lambda.count.min.ratio = lambda.count.min.ratio}

	if(missing(lambda.zero.min.ratio)||lambda.zero.min.ratio<0){
    	param$lambda.zero.min.ratio = 0.1
	}else{ param$lambda.zero.min.ratio = lambda.zero.min.ratio}

	if(missing(alpha.count)||alpha.count<0){
    	param$alpha.count = 1
  	}else{ param$alpha.count = alpha.count}

	if(missing(alpha.zero)||alpha.zero<0){
    	param$alpha.zero = alpha.count
  	}else{ param$alpha.zero = alpha.zero}

	if(missing(gamma.count)||gamma.count<0){
    	param$gamma.count = 3
  	}else{ param$gamma.count = gamma.count}

	if(missing(gamma.zero)||gamma.zero<0){
    	param$gamma.zero = gamma.count
  	}else{ param$gamma.zero = gamma.zero}

	if(missing(rescale)){
	   param$rescale = FALSE
	 }else{ param$rescale = rescale}

	if(missing(theta.fixed)){
    	if(family=="negbin"){
    		param$theta.fixed = FALSE
    	}
  	}else{ param$theta.fixed = theta.fixed}

	if(missing(EM)){
    	param$EM =  TRUE
  	}else{ param$EM = EM}

	if(missing(maxit.em)||maxit.em<0){
    	param$maxit.em = 200
  	}else{ param$maxit.em = maxit.em}

	if(missing(maxit)||maxit<0){
    	param$maxit = 1000
  	}else{ param$maxit = maxit}

	if(missing(maxit.theta)||maxit.theta<0){
    	param$maxit.theta = 25
  	}else{ param$maxit.theta = maxit.theta}

	if(missing(reltol)||reltol<0){
    	param$reltol = 1e-5
  	}else{ param$reltol = reltol}

	if(missing(eps.bino)||eps.bino<0){
    	param$eps.bino = 1e-5
  	}else{ param$eps.bino = eps.bino}

	if(missing(shortlist)){
    	param$shortlist = FALSE
  	}else{ param$shortlist = shortlist}

	if(missing(trace)){
    	param$trace = FALSE
  	}else{ param$trace = trace}

	#data = data.frame(data)
	assign("param", param, envir=.GlobalEnv)

	### Fitting of zeroinflated model

	mod.zonn = zeroinfl(formula = formula, data = data, subset = param$subset, dist = param$family, na.action = na.action, weights = param$weights, model = param$model, y = param$y, x = param$x)
	#Sys.getenv("mod.zonn")
	### Checking error and finding weights

	#if(class(mod.zonn)!="try-error" && !is.null(mod.zonn$coefficients)){
		coeff = mod.zonn$coefficients
		count = 1:length(coeff$count)
		coeff = unlist(coeff)
		sd_coeff = sqrt(diag(vcov(mod.zonn)))
		wts = sd_coeff / abs(coeff)
	#}else if(class(mod.zonn)=="try-error"){
	#	mod.zonn = glmregNB(formula, data, weights, nlambda = nlambda, lambda = NULL, lambda.min.ratio = lambda.min.ratio, alpha = 0, gamma = 3, rescale = rescale, standardize = standerdize, penalty.factor = rep(1, nvars), thresh = 0.001, maxit.theta = maxit.theta, maxit = maxit, eps = eps.bino, trace = trace, start = NULL, etastart = NULL, mustart = NULL, theta.est = TRUE, theta0 = NULL, init.theta= ifelse(theta.est, theta0[1],NULL),link = link, penalty = penalty, method = "glmreg_fit", model = model,
 #   x.keep = x, y.keep = y, contrasts = NULL, convex = FALSE, ...)

	#	coeff = unlist(mod.zonn$coefficients)
	#	sd_coeff = sqrt(diag(vcov(mod.zonn)))
	#	wts = sd_coeff / abs(coeff)
	#}

	zero = setdiff(1:length(wts),count)
	wts_count = wts[count][-1]
	wts_zero = wts[zero][-1]

	### Fitting zipath with specific weights

	fit.zonn = zipath(formula, data, na.action = na.action, standardize = param$standardize, family = param$family, penalty = param$penalty, start = start, model = param$model, y = param$y, x = param$x, nlambda = param$nlambda, lambda.count = lambda.count, lambda.zero = lambda.zero, penalty.factor.count = wts_count, penalty.factor.zero = wts_zero, lambda.count.min.ratio = param$lambda.count.min.ratio, lambda.zero.min.ratio = param$lambda.zero.min.ratio, alpha.count = param$alpha.count, alpha.zero = param$alpha.zero, gamma.count = param$gamma.count, gamma.zero = param$gamma.zero, rescale = param$rescale, init.theta = init.theta, theta.fixed = param$theta.fixed, EM = param$EM, maxit.em = param$maxit.em, convtype = convtype, maxit = param$maxit, maxit.theta = param$maxit.theta, reltol = param$reltol, eps.bino = param$eps.bino, shortlist = param$shortlist, trace = param$trace)

}
