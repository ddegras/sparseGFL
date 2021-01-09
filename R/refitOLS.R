######################################
# UTILITY FUNCTIONS FOR DATA ANALYSIS
######################################


## Function for calculating residual sum of squares
## from SGFL fit in model y(t) = A(t) x(t) + e(t)

# Inputs
# A: matrix of vectorized A(t1), ..., A(tK) 
# where t1,...,tK are the segment starts
# x: original matrix of predictor variables 
# (additional row of 1's associated with intercepts 
# will be added if necessary) 
# y: matrix of responses
# cp: change points t1,...,tK

rssfun <- function(A,x,y,cp)
{
	n <- ncol(y); d <- nrow(y); p <- nrow(A)/d
	if (nrow(x) < p) x <- rbind(x,rep(1,n))
	start <- c(1,cp); end <- c(cp-1,n)
	K <- length(start)
	rss <- numeric(K)
	for (k in 1:K) {
		idx <- start[k]:end[k]
		Ak <- matrix(A[,k],d,p)
		rss[k] <- sum((y[,idx] - Ak %*% x[,idx])^2)
	}
	sum(rss)
}





## Function for refitting linear regression model
## after variable selection by sparse group fused 
## lasso/elastic net 
 
# Inputs
# object	- result of a call to SGFL.AXY
# x 		- predictor matrix   
# y			- response matrix
# minlen	- minimum number of time points per segment
# cond		- target condition number for regularized OLS


refit.AXY <- function(object, x, y, minlen=2, cond=1e4)
	# , opts=list())
{
	object$objective <- NULL
	n <- ncol(x)
	d <- nrow(y) # number of response variables
	I <- Diagonal(d)
	A <- as(object[["A"]], "TsparseMatrix")
	m <- nrow(A)
	p <- m/d # number of predictor variables
	intercept <- object$intercept
	if (intercept && nrow(x) < p)
		x <- rbind(x,rep(1,n))
	startc <- c(1,object[["changePoints"]])
	endc <- c(startc[-1]-1L,n)
	nchain <- length(startc)
	seglen <- endc - startc + 1
	long <- which(seglen >= minlen)
	if (length(long) == 0) return(object)
	short <- setdiff(1:nchain,long)
	
	## First pass: refit SVAR model on long segments 
	## (i.e. of sufficient length)
	for (j in long) {
		chain <- startc[j]:endc[j]
		if (intercept && length(chain) == 1) 
			{ A[1:(m-d),chain] <- 0
			A[(m-d+1):m,chain] <- y[,chain] 
			next }
		idx <- (A@i)[(A@j) == j-1] + 1
		if (intercept)
			idx <- union(idx,(m-d+1):m)
		if (length(idx) == 0) next 
		lhs <- kronecker(t(x[,chain]),I)
		lhs <- lhs[,idx]
		rhs <- y[,chain]
		dim(rhs) <- c(length(rhs),1)
		if (ncol(lhs) <= nrow(rhs)) {
			A[idx,j] <- lsfit(lhs,rhs,intercept=FALSE)$coef
		} else {
			rhs <- crossprod(lhs,rhs)
			lhs <- crossprod(lhs)
			dmax <- norm(x[,chain],"2")
			diag(lhs) <- diag(lhs) + dmax^2/(cond-1)
			A[idx,j] <- solve(lhs,rhs)
		}
	}
	if (length(short) == 0) 
		{ object[["A"]] <- A
		return(object) }
	
	## Second pass: fuse short segments with an adjacent
	## long segment
	i <- short[1] 
	start2 <- startc
	# end2 <- endc
	while (i <= nchain) {

		## Find all short segments contiguous to current
		## short segment
		inext <- i
		while (is.element(inext+1,short)) inext <- inext+1
		chain <- startc[i]:endc[inext]	

		## If there exist long segments on both sides of
		## contiguous short segments, determine which long
		## segment the short segments should be fused with
		if (i > 1 && inext < nchain) {
			Aprev <- A[,i-1]
			dim(Aprev) <- c(d,p)
			sse1 <- sum((y[,chain] - Aprev %*% x[,chain])^2) 
			Anext <- A[,inext+1]
			dim(Anext) <- c(d,p)
			sse2 <- sum((y[,chain] - Anext %*% x[,chain])^2) 
		}
		
		## Perform fusion
		if (i == 1 || (inext < nchain && sse1 >= sse2)) {
			start2 <- setdiff(start2,startc[(i+1):(inext+1)])
			# end2 <- setdiff(end2,endc[i:inext])			
		} else {
			start2 <- setdiff(start2,startc[i:inext])
			# end2 <- setdiff(end2,endc[(i-1):(inext-1)])
		}					
		
		## Move to next short segment		
		i <- suppressWarnings(min(short[short>inext]))
	}
			
	object[["A"]] <- A[,long,drop=FALSE]
	object[["changePoints"]] <- start2[-1]
	return(object)	
}

 
 
 




## Quick visualization of time series
plot.mts <- function(x, k=2, xlab="", ylab="")
{
	n <- ncol(x)
	d <- nrow(x)
	xbar <- rowMeans(x)
	x <- x - xbar
	sig <- sqrt(rowMeans(x^2))
	x <- (x/sig) + k*(d:1)
	matplot(t(x), xlab=xlab, ylab=ylab, 
		type="l", yaxt="n", lty=1)
}











# SLIGHTLY SLOWER VERSION
 
# refit.SVAR <- function(object, x, y, minlen=2, cond=10,
	# opts=list(maxitr=5000, ncv=30))
# {
	# if (is.null(object)) return(NULL)
	# RSpectra.flag <- require(RSpectra)
	# n <- ncol(x)
	# N <- nrow(y)
	# I <- Diagonal(N)
	# A <- as(object[["A"]], "TsparseMatrix")
	# m <- nrow(A)
	# intercept <- object$intercept
	# startc <- c(1,object[["changePoints"]])
	# endc <- c(object[["changePoints"]]-1,n)
	# nchain <- length(startc)
	# seglen <- endc - startc + 1
	# long <- which(seglen >= minlen)
	# if (length(long) == 0) return(NULL)
	# short <- (1:nchain)[-long]
	
	# ## First pass: refit SVAR model on long segments 
	# ## (i.e. of sufficient length)
	# for (j in long) {
		# chain <- startc[j]:endc[j]
		# if (intercept && length(chain) == 1) 
			# { A[1:(m-N),chain] <- 0
				# A[(m-N+1):m,chain] <- y[,chain] 
			# next }
		# idx <- (A@i)[(A@j) == j-1] + 1
		# if (length(idx) == 0) next 
		# lhs <- kronecker(t(x[,chain]),I)
		# rhs <- crossprod(lhs[,idx],
			# as.vector(y[,chain]))
		# lhs <- crossprod(lhs[,idx])
		# if (length(idx) <= 100 || !RSpectra.flag) {
			# d <- eigen(lhs,TRUE,TRUE)$values
			# c1 <- (max(d)-cond*min(d)) / (cond-1)
		# if (c1 > 0)
			# diag(lhs) <- diag(lhs) + c1
		# } else {
			# dmax <- eigs(lhs,k=1,opts=opts)$values
			# c1 <- dmax / (cond - 1)
			# diag(lhs) <- diag(lhs) + c1
			# dmin <- eigs(lhs,k=1,which="SM",
				# opts=opts)$values - c1
			# if (dmin >= 1e-3) {
				# c2 <- (dmax - cond * dmin) / (cond - 1)
				# if (c2)
				# diag(lhs) <- diag(lhs) + (c2-c1)
			# } 
		# }
		# A[idx,j] <- as.vector(solve(lhs,rhs))
	# }
	# if (length(short) == 0) 
		# { object[["A"]] <- A
		# return(object) }
	
	# ## Second pass: fuse short segments with an adjacent
	# ## long segment
	# i <- short[1] 
	# start2 <- startc
	# # end2 <- endc
	# while (i < nchain) {

		# ## Find all short segments contiguous to current
		# ## short segment
		# i <- inext <- i
		# while (is.element(inext+1,short)) inext <- inext+1
		# chain <- startc[i]:endc[inext]	

		# ## If there exist long segments on both sides of
		# ## contiguous short segments, determine which long
		# ## segment the short segments should be fused with
		# if (i > 1 && inext < nchain) {
			# Aprev <- A[,i-1]
			# dim(Aprev) <- c(N,m/N)
			# sse1 <- sum((y[,chain] - Aprev %*% x[,chain])^2) 
			# Anext <- A[,inext+1]
			# dim(Anext) <- c(N,m/N)
			# sse2 <- sum((y[,chain] - Anext %*% x[,chain])^2) 
		# }
		
		# ## Perform fusion
		# if (i == 1 || sse1 >= sse2) {
			# start2 <- setdiff(start2,startc[(i+1):(inext+1)])
			# # end2 <- setdiff(end2,endc[i:inext])			
		# } else {
			# start2 <- setdiff(start2,startc[i:inext])
			# # end2 <- setdiff(end2,endc[(i-1):(inext-1)])
		# }					
		
		# ## Move to next short segment		
		# i <- suppressWarnings(min(short[short>inext]))
	# }
			
	# ## Reduce regression coefficient matrix and change points
	# object[["A"]] <- A[,long]
	# object[["changePoints"]] <- start2[-1]
	# return(object)	
# }
