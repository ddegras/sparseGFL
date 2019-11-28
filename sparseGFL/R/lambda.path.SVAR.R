lambda.path.SVAR <- function(x, p, w=NULL, struct=FALSE,
   intercept=TRUE, alpha=1, nlambda1=20, lambda1.min.ratio=0.01, 
   lambda1=NULL, nlambda2=20, lambda2.min.ratio=0.01,
   control=list(), parallel=FALSE, warm.start=FALSE)
{

#####################
# Data preprocessing
#####################

## Try to load packages glmnet and Matrix, stop if fail
test.glmnet <- require(glmnet)
stopifnot(test.glmnet)
test.Matrix <- require(Matrix)
stopifnot(test.Matrix)

## Data dimensions
stopifnot(p < ncol(x))
N <- nrow(x)
n <- ncol(x)
lag0 <- struct

## Initialize TV weights if needed
if (is.null(w)) { 
   w2 <- rep.int(1,n-p-1L)
} else { 
   stopifnot(all(w > 0)) 
   stopifnot(length(w) == ncol(x)-p-1L || length(w) == 1L)
   w2 <- if (length(w) == 1L) 
      { rep.int(w^2,n-p-1L) } else { w^2 }
}
   
## Set up control parameters 
con <- list(s=.001, decrease=FALSE, maxit=500L, tol=1e-5,
   maxit.proj=100L, tol.proj=1e-6)
   test <- is.element(names(con),names(control))
for (i in seq_along(con)) {
   name.i <- names(con)[i]
   assign(name.i, ifelse(test[i], control[[name.i]], con[[i]]))
}

## Set up predictor (x) and response (y)
y <- x[,(p+1):n]
idx <- mapply(seq.int, 
   from=(p+lag0):(n-1+lag0), to=1:(n-p))
x <- x[,idx]
dim(x) <- c(N*(p+lag0),n-p)
xcopy <- x
if (intercept)
	xcopy <- rbind(xcopy,rep(1,n-p))

# Center if model includes intercept
if (intercept) {
   xbar <- rowMeans(x)
   x <- x - xbar
   ybar <- rowMeans(y)
   y <- y - ybar
}
# Vectorize regression model y = Ax + e
# --> vec(y) = kron(x',I) vec(A) + vec(e)
# with I identity matrix 
I <- Diagonal(N)
x <- as(kronecker(t(x),I), "dgCMatrix")
dim(y) <- NULL
# Remove diagonal coefficients from first regression 
# matrix A(0) if lag 0 used in regression model (SVAR)
if (lag0) {
   diag.idx <- cumsum(c(1,rep.int(N+1,N-1)))
   x <- x[,-diag.idx] }
 

##############
# Elastic net
##############

# In glmnet the sum of squared errors in the objective 
# function is scaled by the length of the response vector 
# Therefore lambda values must be rescaled
if (is.null(lambda1)) {
   lasso <- glmnet(x, y, alpha=alpha, nlambda=nlambda1, 
      lambda.min.ratio=lambda1.min.ratio,
      standardize=FALSE, intercept=FALSE, thres=1e-8) 
      lambda1 <- lasso$lambda * N
} else {
   lasso <- glmnet(x, y, alpha=alpha, lambda=lambda1/N,
      standardize=FALSE, intercept=FALSE, thres=1e-8)
   nlambda1 <- length(lambda1)
}


###########################
# Reshape results and data
###########################

## Regression coefficients
if (!lag0 && !intercept) {
    A <- lasso$beta
} else {
   A <- Matrix(0,(p+lag0)*N^2+(intercept*N),nlambda1) 
   idx <- seq.int(1,(p+lag0)*N^2)
   if (lag0)
      idx <- idx[-diag.idx]
   A[idx,] <- lasso$beta
   if (intercept) {
   last <- seq.int(nrow(A)-N+1,nrow(A))
      A[last,] <- ybar - 
      kronecker(Matrix(t(xbar)),I) %*% A[-last,]
   }
}

## Prepare output
out <- vector("list",nlambda1)
for (i in 1:nlambda1) {
   Ai <- matrix(A[,i],N, N*(p+lag0)+intercept)
   Ai[abs(Ai)<1e-8] <- 0
   out[[i]] <- list(lambda1 = lambda1[i], lambda2 = NULL, 
      alpha=alpha, A = Ai)
}

## Check if lambda2 values are required
if (nlambda2 == 0) 
   return(out)

## Reshape data
# x <- matrix(x@x, ncol=n-p, byrow=TRUE)
# x <- x[seq.int(1,by=N,len=N*(p+lag0)),]
x <- xcopy
dim(y) <- c(N,n-p)


########################
# Projected subgradient 
########################

# Find lambda2 path for each lambda1 value
# --> find maximum lambda2 for which lasso
# solution solves GFL. 
null <- matrix(0,1,0)

## Case: parallel computation
if (parallel) {
   lambda2 <- foreach(i=1:nlambda1, 
   .packages="sparseGFL") %dopar% {  
      lambda1 <- out[[i]]$lambda1
      A <- out[[i]]$A
      dim(A) <- c(N, N*(p+lag0)+intercept)
      lam2 <- 0
      wprev <- wnext <- 0
      Aprev <- Anext <- null
      z <- grad_svar(A, x, y, lambda1*alpha, lam2,
         lambda1*(1-alpha), wprev, Aprev, wnext, Anext, 
         lag0, intercept)
      fixed <- which(A != 0)
      if (intercept) 
         fixed <- union(fixed,last)
      # z <- tcpfun(A %*% x - y, x) 
      # if (alpha < 1)
         # z <- z + (lambda1*(1-alpha)) * as.vector(A)  
      if (length(fixed)) {
         # z[fixed,] <- z[fixed,] + (lambda1*alpha) * sign(A[fixed])
         vfixed <- apply(z[fixed,-(n-p),drop=FALSE],1,cumsum)
     	 c2 <- rowSums(vfixed^2) }

      if (length(fixed) == length(A)) {
      # Case: all lasso regression coefficients are nonzero
         lambda2.max <- sqrt(max(c2/w2))
      } else if (length(fixed) == 0) {
      # Case: all lasso coefficients are zero
         c2 <- numeric(n-p-1L)
         result <-lambda2_max(null, z, lambda1*alpha, c2, w2,
            s, decrease, maxit, tol, maxit.proj, tol.proj)
         lambda2.max <- result$objective
      } else {
   	  # Case: some zero, some nonzero coefficients
         result <-lambda2_max(null, z[-fixed,], lambda1*alpha,
         c2, w2, s, decrease, maxit, tol, maxit.proj, tol.proj)
      lambda2.max <- result$objective
      }
         
      exp(seq(from=log(lambda2.max), len=nlambda2,
         to=log(lambda2.max*lambda2.min.ratio)))
   }
      
   for (i in 1:nlambda1)
      out[[i]]$lambda2 <- lambda2[[i]]
   return(out)
}

## Case: serial computation
u <- null
for (i in 1:nlambda1) { 
   lambda1 <- out[[i]]$lambda1
   A <- out[[i]]$A
   dim(A) <- c(N, N*(p+lag0)+intercept)
   lam2 <- 0
   wprev <- wnext <- 0
   Aprev <- Anext <- null
   z <- grad_svar(A, x, y, lambda1*alpha, lam2,
      lambda1*(1-alpha), wprev, Aprev, wnext, Anext, 
      lag0, intercept)
   fixed <- which(A != 0)
   if (intercept) 
      fixed <- union(fixed,last)
   # z <- tcpfun(A %*% x - y, x)  
   # if (alpha < 1)
      # z <- z + (lambda1*(1-alpha)) * as.vector(A)    
   if (length(fixed)) {
      # z[fixed,] <- z[fixed,] + (lambda1*alpha) * sign(A[fixed])
      vfixed <- apply(z[fixed,-(n-p),drop=FALSE],1,cumsum)
      c2 <- rowSums(vfixed^2) }

   if (length(fixed) == length(A)) {
   # Case: all lasso regression coefficients are nonzero
      lambda2.max <- sqrt(max(c2/w2))
      u <- null
   } else if (length(fixed) == 0) {
   # Case: all lasso coefficients are zero
      c2 <- numeric(n-p-1L)
      result <-lambda2_max(u, z, lambda1*alpha, c2, w2,  
         s, decrease, maxit, tol, maxit.proj, tol.proj)
      lambda2.max <- result$objective   
      u <- if (warm.start) result$u else null
   } else {
   # Case: some zero, some nonzero lasso coefficients	
      if (length(u)) u <- u[-fixed,] 
      result <- lambda2_max(u, z[-fixed,], lambda1*alpha, c2,
      w2, s, decrease, maxit, tol, maxit.proj, tol.proj)
      lambda2.max <- result$objective
      if (warm.start) {
         u <- z ; u[-fixed,] <- result$u 
      } else { u <- null }
   }   
         
   out[[i]]$lambda2 <- exp(seq(from=log(lambda2.max),
      to=log(lambda2.max*lambda2.min.ratio), len=nlambda2))
   }
return(out)
}