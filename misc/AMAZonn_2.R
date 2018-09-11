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
	
	data = data.frame(data)
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