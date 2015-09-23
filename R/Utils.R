#' Utility functions.
#' 
#' Some useful miscellaneous functions: plot plot.MECTree plot the 
#'   tree part of MECTree object; heat.tree plot a sort of heat map decision trees 
#'
#' @name Utils
#'
#' @param dat data.frame 
#' @param id group id 
#' @param threshold threshold 
#' @param order.vars variables to order longitudinal data 
#' @param cut.off,propo  cut-off and proportions 
#' @param no.last number of tail observations that goes into test set 
#' @param resp.vars,rhs.vars response and vector of predictor names 
#' @param pred,obs predicted probabilities and true class 
#' @param y,f,nl.n numeric or character vector depending on usage 
#' @param allow.new.levels logical 
#' @param lag lag in train test split 
#' @param formula formula 
#' @param \dots Further arguments passed to or from other methods.
#' 
#' @author  Che Ngufor Ngufor.Che@@mayo.edu
#' @import rpart.plot partykit Matrix
#' 
#' @importFrom PresenceAbsence optimal.thresholds
#' @importFrom PresenceAbsence presence.absence.accuracy
#' @importFrom plyr dlply
#'
NULL 
#'
#' @rdname  Utils
#' @export
#'
###############################################################################
#################  Train-Test for subjects in sample 
insampleTrainTestSplit <- function(dat, id, order.vars = "id", cut.off = 4, propo = 0.25, no.last = 3){ 
tab <- table(dat[, id])
ix <- names(tab[tab >= cut.off])
ixv <- sample(ix, floor(length(ix)* propo))
dd <- dat[dat$id%in%ixv, ]
dd$id <- factor(dd$id)
ixy <- split(1:nrow(dd), f = dd$id)
dd <- dat <- do.call(rbind, lapply(ixy, function(x) dd[x, ]))            

Q <- lapply(ixy, function(x) {
tmp <- dd[x, ]
tmp <- tmp[order(tmp[, order.vars]), ]
train <- head(tmp, nrow(tmp)-no.last) 
test <- tail(tmp, no.last)
list(train = train, test = test)
}
)

train <- do.call(rbind, lapply(Q, function(x) x$train))
dd <- dat[!dat$id%in%ixv, ]
train <- rbind(dd, train)
rownames(train) <- NULL
test <- do.call(rbind, lapply(Q, function(x) x$test))
rownames(test) <- NULL
list(train = train, test = test)
}     

#'
#' @rdname  Utils
#' @export
#'
###############################################################################
#################  Train-Test for subjects in sample 
insampleTrainTestSplit.mtl <- function(dat, id, resp.vars, rhs.vars, cut.off = 4, 
                              propo = 0.25, no.last = 3){ 
                              
r.form <- as.formula(paste0(" ~ ", paste0(c(" -1", rhs.vars), collapse = "+")))
Z <- data.frame(as.matrix(sparse.model.matrix(r.form, dat)))   
dt <- cbind(response = dat[,resp.vars], id =  dat[,id], Z, time = dat$time)
colnames(dt)[1] <- resp.vars
dt[, resp.vars] <- as.numeric(factor(dt[, resp.vars]))-1

rhs.vars <- colnames(Z)

#ixy <- as.character(unique(dt$id))
#dt <- do.call(rbind, lapply(ixy, function(x) cbind(subset(dt, id == x), 
#             time = paste0("t", 1:nrow(subset(dt, id == x)))) ))            
                           
tab <- table(dt[, id])
ix <- names(tab[tab >= cut.off])
ixv <- sample(ix, floor(length(ix)* propo))

### subjects with at least  cut.off number of observations 
dd <- subset(dt, id%in%ixv)
ixy <- as.character(unique(dd$id))

Q <- lapply(ixy, function(x) {
tmp <- subset(dd, id==x)
tmp <- tmp[order(tmp$year), ]
train <- head(tmp, nrow(tmp)-no.last) 
test <- tail(tmp, no.last)
list(train = train, test = test)
}
)

train <- do.call(rbind, lapply(Q, function(x) x$train))
train <- rbind(subset(dt, !id%in%ixv), train)
resp.vars1 <- unique(train$time)
Qt <-  lapply(resp.vars1,  function(x) subset(train, time == x))
names(Qt) <- resp.vars1

Q.train <- lapply(resp.vars1, function(x) {
	colnames(Qt[[x]])[which(colnames(Qt[[x]])==resp.vars)] <- x
	Qt[[x]]})
names(Q.train) <- resp.vars1

test <- do.call(rbind, lapply(Q, function(x) x$test))
resp <- unique(test$time)

Qt <-  lapply(resp,  function(x) subset(test, time == x))
names(Qt) <- resp
Q.test <- lapply(resp, function(x) {
	colnames(Qt[[x]])[which(colnames(Qt[[x]])==resp.vars)] <- x
	Qt[[x]]})
names(Q.test) <- resp
list(Q.train = Q.train, Q.test = Q.test, resp.vars = resp.vars1, test.resp = resp, 
                rhs.vars = rhs.vars)
}     
#'
#' @rdname Utils
# returns t(x).y 
cprod <- function(x,y=NULL){
 if(is.complex(x) | is.complex(y)){
    if(is.null(y)){
     return(crossprod(Conj(x),x))
    } else {
     return(crossprod(Conj(x),y))
    }
 } else {
    return(crossprod(x,y))
 }
}

###  parse and manipulating mixed-model formulas
##' deparse(.) returning \bold{one} string
##' @param x,collapse character vector and how to combine them 
##' @note Protects against the possibility that results from deparse() will be
##'       split after 'width.cutoff' (by default 60, maximally 500)
safeDeparse <- function(x, collapse=" ") paste(deparse(x, 500L), collapse=collapse)

#'
#' @rdname Utils
# Random Effects formula only
reOnly <- function(f) {
    reformulate(paste0("(", vapply(findbars(f), safeDeparse, ""), ")"))
}
#'
#' @rdname Utils
##' @param x a random effect (i.e., data frame with rows equal to levels, columns equal to terms
##' @param n vector of new levels
levelfun <- function(x,nl.n,allow.new.levels=FALSE) {
    ## 1. find and deal with new levels
    if (!all(nl.n %in% rownames(x))) {
        if (!allow.new.levels) stop("new levels detected in newdata")
        ## create an all-zero data frame corresponding to the new set of levels ...
        newx <- as.data.frame(matrix(0, nrow=length(nl.n), ncol=ncol(x),
                                     dimnames=list(nl.n, names(x))))
        ## then paste in the matching RE values from the original fit/set of levels
        newx[rownames(x),] <- x
        x <- newx
    }
    ## 2. find and deal with missing old levels
    ## ... these should have been dropped when making the Z matrices
    ##     etc. in mkReTrms, so we'd better drop them here to match ...
    if (!all(r.inn <- rownames(x) %in% nl.n)) {
        x <- x[r.inn,,drop=FALSE]
    }
    return(x)
}

#'
#' @rdname  Utils
#' @export
#'
###############################################################################
#################  Train-Test for subjects in sample: minimum number of visit time = 3 
trajTrainTestSplit <- function(dat, id, rhs.vars, resp.vars, order.vars = "time", lag=1){ 

test.last <- function(x){
if(nrow(x) < 3) return(NULL) 
x <- x[order(x[, order.vars], decreasing = FALSE), ]
ix1 <- 1:(nrow(x)-1)
ix2 <- 2:nrow(x)
dd <- cbind.data.frame(x[ix1, c(id, rhs.vars, resp.vars)], response = x[ix2, resp.vars]) 
dd[,resp.vars]  <- factor(dd[, resp.vars])
train <- head(dd, nrow(dd)-1)
test <- tail(dd, 1) 
return(list(train = train, test = test))
}

X <- dlply(dat, .variables = id, .fun = test.last)
X[sapply(X, is.null)] <- NULL 
 
train <- do.call(rbind.data.frame, lapply(X, function(x) x$train))
test <- do.call(rbind.data.frame, lapply(X, function(x) x$test))

list(train = train, test = test, resp.vars = "response", rhs.vars = c(rhs.vars, resp.vars))
}     
#'
#' @rdname  Utils
#' @export
#'
######################################################################################
#################  Train-Test split longitdinal data: train data lags test data by lag 
#################  there must be at least lag + 2 visit times for each observation 
LongiLagSplit <- function(dat, id, rhs.vars, resp.vars, order.vars = "time", lag=1){ 

lag.split <- function(x){
if(nrow(x) < (lag+2)) return(NULL) 
x <- x[order(x[, order.vars], decreasing = FALSE), , drop=FALSE]
ix1 <- 1:(nrow(x)-lag)
ix2 <- (lag+1):nrow(x)
dd <- cbind.data.frame(x[ix1, c(id, rhs.vars, resp.vars)], response = x[ix2, resp.vars]) 
dd[,resp.vars]  <- factor(dd[, resp.vars])
train <- head(dd, nrow(dd)-1)
test <- tail(dd, 1) 
return(list(train = train, test = test))
}

X <- dlply(dat, .variables = id, .fun = lag.split)
X[sapply(X, is.null)] <- NULL 
 
train <- do.call(rbind.data.frame, lapply(X, function(x) x$train))
test <- do.call(rbind.data.frame, lapply(X, function(x) x$test))

list(train = train, test = test, resp.vars = "response", rhs.vars = c(rhs.vars, resp.vars))
}     




#' @rdname  Utils 
#' @export
## performance measures 
Performance.measures <- function(pred, obs, threshold=NULL){
obs <- as.numeric(factor(obs))-1 
## get best cut-off 
if(is.null(threshold))
threshold <- opt.thresh(pred, obs)
### get these performance measures
nme = c("PCC", "PCC.sd", "AUC", "AUC.sd", "sensitivity", "sensitivity.sd", 
"specificity", "specificity.sd")
xx = cbind.data.frame(plotID = 1:length(pred), Observed = obs, Predicted = pred)
accuracy <- presence.absence.accuracy(xx, threshold = threshold, st.dev = TRUE)[, nme]
accuracy$G.mean <- sqrt(as.numeric(accuracy$sensitivity)*as.numeric(accuracy$specificity))
accuracy$BER <- 1 - 0.5*(as.numeric(accuracy$sensitivity) + as.numeric(accuracy$specificity)) 
accuracy$threshold = threshold 
return(accuracy)
}
#'
#' @rdname  Utils 
sigmoid <- function(x) 1.0 / (1 + exp(-x))
#' @rdname Utils
# get left hand side of formula
lhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[response] 
}
#' @rdname  Utils
# get right hand side of formula 
rhs.form <- function(formula) {
    tt <- terms(formula)
    vars <- as.character(attr(tt, "variables"))[-1] ## [1] is the list call
    response <- attr(tt, "response") # index of response var
    vars[-response] 
}
#' @rdname  Utils
#' @export
opt.thresh <- function(pred, obs){
thresh = 0.5 
if(length(unique(obs)) > 1){
obs <- as.numeric(as.factor(obs))-1 
SIMDATA = cbind.data.frame(plotID = 1:length(obs), Observed = obs, Predicted = pred)
thresh <- optimal.thresholds(SIMDATA, threshold = 101, which.model = 1, opt.methods = 9)
thresh <- ifelse(length(thresh["Predicted"]) >= 1,as.numeric(thresh["Predicted"]), 0.5)
}
return(thresh)
}





