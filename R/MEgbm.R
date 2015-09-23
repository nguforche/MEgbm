#' MEgbm  
#' 
#' Trains a Mixed Effect Boosted regression and  
#'  classification models for longitudinal continuous, binary and count data 
#'  
#' @name MEgbm
#'
#' @param X  data.frame with predictors 
#' @param Y  binary response vector 
#' @param groups character name of the column containing the group identifier
#' @param rand.vars random effect variables 
#' @param gbm.dist gbm loss function   
#' @param para named list of gbm training parameters 
#' @param lme.family glmer control 
#' @param tol convergence tolerance 
#' @param max.iter maximum number of iterations  
#' @param include.RE (logical) to include random effect Zb as predictor in gbm?  
#' @param verbose verbose for lme4 
#' @param likelihoodCheck (logical) to use log likelihood of glmer to check for convergence? 
#' @param type of predictions of gbm to pass to lme4 as population estimates (these will be used as offset) 
#' @param \dots Further arguments passed to or from other methods.
#' @return An object of class MEgbm; a list with items 
#' \item{gbmfit}{fitted gbm model}
#' \item{glmer.fit}{fitted mixed effect logistic regression model}
#' \item{logLik}{log likelihood of mixed effect logistic regression} 
#' \item{random.effects}{random effect parameter estimates}
#' \item{boost.form}{gbm formula for fitted model}
#' \item{glmer.form}{lmer4 formula} 
#' \item{glmer.CI}{estimates of mixed effect logistic regression with 
#'     approximate confidence intervals on the logit scale. More accurate values 
#'     can be obtained by bootstrap}
#' \item{fitted.probs}{fitted probabilites for final model}
#' \item{fitted.class}{fitted class labels for final model}
#' \item{train.perf}{various performance measures for final model on training set}
#' \item{threshold}{classification cut-off}
#'
#' @author  Che Ngufor <Ngufor.Che@@mayo.edu>
#' @import lme4 caret partykit 
NULL 
#'
#'
#' @rdname MEgbm  
#' @export
MEgbm  <- function(X, ...) UseMethod("MEgbm")
#'
#' @rdname MEgbm 
#' @export
#' @examples
#' \dontrun{
#' set.seed(12345)
#' mod <- MEgbm(form, rand.form, data)) 
#'
#' }
#'
MEgbm <- function(X, Y, groups = NULL, rand.vars="1", para = NULL,
                 gbm.dist = "bernoulli", lme.family = binomial,  
                 tol= 1e-5, max.iter =100, include.RE =TRUE, verbose = FALSE, 
                 likelihoodCheck = TRUE, type = "prob", ...){
                
     if(is.null(groups)) stop("please provide grouping variable")
     Y <- as.vector(Y) 
     dat <- cbind.data.frame(response = Y, X)   
     resp.vars <- "response"      
     dat[, resp.vars] <- as.numeric(factor(dat[, resp.vars]))-1  
     X[, groups] <- NULL 
     Y <- factor(Y); levels(Y) <- c("No", "Yes")	 
 
if(!include.RE){     
   if(para$opt.para)
	gbmfit <-  train(X, Y,  method = "gbm", trControl = 
	           trainControl(method = para$method, number =  para$number), 
	           verbose = FALSE,  tuneLength = para$tuneLength, distribution = gbm.dist) 
   else 
	gbmfit <-  train(X, Y,  method = "gbm", trControl = trainControl(method = "none"),  
	           verbose = FALSE, tuneGrid = data.frame(n.trees = para$n.trees, 
               interaction.depth=para$interaction.depth, shrinkage=para$shrinkage,
                n.minobsinnode = para$n.minobsinnode), distribution = gbm.dist)
res <- list(gbmfit = gbmfit, include.RE = include.RE, groups = groups)                
class(res) <- "MEgbm"                
return(res)
}     

        
    old.lik <- -Inf    
	if(likelihoodCheck == FALSE){
	n.re <- sum(rand.vars != 1)+1  
    b.old <-rep(0, n.re*nlevels(factor(dat[, groups])))  ## initial random effect parameters 
    }

### initial random effect component: fit a LME model with no fixed effect, ie 
### a priori known mean of 0 

	dat[, "tree.fit"] <- rep(0, nrow(dat))
    	
	form.glmer <- as.formula(paste0("response ~ ", paste0(c("-1 + offset(tree.fit) + ", 
	"(", paste0(c(rand.vars), collapse = "+"), "|", groups, ")"), collapse = "")))             

	glmer.fit <- glmer(form.glmer, data= dat,family = lme.family, 
	             control = glmerControl(optimizer = "bobyqa"), 
		          nAGQ=  0, verbose = as.numeric(verbose))
		          
	pp = predict(glmer.fit, newdata = dat, type = "response")  
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
	w = pp*(1-pp)    		              
    
    		              
### get the random effect component 
	Zt <-  getME(glmer.fit, name = "Zt")
	b  <-  getME(glmer.fit, name = "b")	
	Zb <-  as.numeric(cprod(x = Zt, y = b))       
	X[, "Zb"] <- Zb 
    
for(ii in 1:max.iter){    	  	

#### fit boosted regression trees 
   if(para$opt.para)
	gbmfit <-  train(X, Y,  method = "gbm", trControl = 
	           trainControl(method = para$method, number =  para$number), 
	           verbose = verbose,  tuneLength = para$tuneLength,
	           weights = w, distribution = gbm.dist) 
   else 
	gbmfit <-  train(X, Y,  method = "gbm", trControl = trainControl(method = "none"),  
	           verbose = verbose, tuneGrid = data.frame(n.trees = para$n.trees, 
               interaction.depth=para$interaction.depth, shrinkage=para$shrinkage,
                n.minobsinnode = para$n.minobsinnode),weights = w, distribution = gbm.dist) 
   if(type == "prob")             
   pp <-  predict(gbmfit, newdata = X, type = "prob")[, 2]
   else 
   pp <-  as.numeric(factor(predict(gbmfit, newdata = X, type = "raw")))-1
   
   dat[, "tree.fit"] <- pp
    
###  Fit mixed effect logistic regression model with nodes as fixed effect predictors
	glmer.fit <- glmer(form.glmer, data= dat,family = lme.family, 
	              control = glmerControl(optimizer = "bobyqa"), 
		          nAGQ=  0, verbose = as.numeric(verbose))


###	 glmer.fit <- glmer(form.glmer, data= dat,family = lme.family, control = glmer.Control, 
###	                 nAGQ=  nAGQ, verbose = as.numeric(verbose))

### compute weights for MOB 	
    pp = predict(glmer.fit, newdata = dat, type = "response")  
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))     
    w = pp*(1-pp)

### test for convergence             
   if(likelihoodCheck){
	new.lik <- as.numeric(logLik(glmer.fit))
	r <- as.numeric(sqrt(t(new.lik - old.lik)%*%(new.lik-old.lik)))
	old.lik <- new.lik	
	} else {
	r <- as.numeric(sqrt(t(b - b.old)%*%(b-b.old)))
	b.old <- b
    } 
	#if(verbose) 
	print(paste("Error: ", r))    
	if( r < tol) break 
	
### compute the adjusted response 
### first get the random effect component 
	Zt <-  getME(glmer.fit, name = "Zt")
	b  <-  getME(glmer.fit, name = "b")	
	Zb <-  as.numeric(cprod(x = Zt, y = b))
	X[,"Zb"] <- Zb 
	
	} ## for loop 
	
	if(r > tol) warning("EM algorithm did not converge")
	
	pp <- ifelse(abs(pp -1) <= 1e-16, pp-1e-16, ifelse(pp <= 1e-16, pp+1e-16, pp))      
        perf <- Performance.measures(pp, Y)
	threshold <- perf$threshold	
	cls <-ifelse(pp >= threshold, 1, 0)	

### get confidence intervals for mixed effect logistic regresion: rough estimates using the SEs
	se <- sqrt(diag(as.matrix(vcov(glmer.fit)))) 
  	tab <- cbind(Est = fixef(glmer.fit), LL = fixef(glmer.fit) - 1.96 * se, 
  	UL = fixef(glmer.fit) + 1.96 *se)

res <- list(gbmfit = gbmfit, glmer.fit = glmer.fit,  groups = groups, 
         rand.vars=rand.vars, logLik=as.numeric(logLik(glmer.fit)), 
         random.effects =ranef(glmer.fit), rand.form = form.glmer, 
         glmer.CI =tab, fitted.probs = pp, 
         fitted.class = cls, train.perf = perf, threshold = threshold, 
         include.RE=include.RE)
class(res) <- "MEgbm"         
return(res)         
}







