
##############################
# Code for the simulations of 
# "Sparse Group Fused Lasso 
# for Model Segmentation"
# Author: David Degras
# Date: December 2019
##############################


# For this code to run properly, the source 
# package 'sparseGFL' and the C+ file 'allmethods.cpp'
# must be in the working directory 


cat("Start time","\n")
print(Sys.time())

## Install packages as needed
install.packages(c("Rcpp","RcppArmadillo","RSpectra"))
install.packages("sparseGFL", type="source", repos=NULL)

## Clear workspace
rm(list=ls())

## Load packages
library(Rcpp)
library(RcppArmadillo)
library(RSpectra)
library(sparseGFL)

## Source C++ code for other optimization methods
sourceCpp("allmethods.cpp")




########################
# Simulation parameters
########################

## Data dimensions 
d <- 100 # number of response variables
p <- 500 # number of predictor variables
n <- 200 # number of time points

## Signal parameters
s <- 0.9 # sparsity level
K <- 10  # number of segments

## Noise parameters
sigma <- 0 # standard deviation 
rhox <- 0  # correlation 

## Regularization parameters
lam1 <- 20   # lasso penalty parameter
lam2 <- 2000 # total variation penalty parameter

## Number of replications 
nrep <- 100 

## Target for relative distance 
## to minimum of objective function 
tol <- 1e-6

## Set seed for random number generation
set.seed(2019)

## Misc
null <- matrix(,0,1)
verbose <- TRUE # FALSE

## Output table
out <- matrix(,nrep,17)
method <- c("SPG","PD","HYB-C","HYB-R","ADMM","LADMM")
colnames(out) <- c("d","p","n","sigma","rho",
	paste0(method,".time"), paste0(method,".best"))
	
out[,"d"] <- d
out[,"p"] <- p
out[,"n"] <- n
out[,"sigma"] <- sigma
out[,"rho"] <- rhox



############
# MAIN LOOP
############

for (rr in 1:nrep)
{
	
	cat("Replication",rr,"\n")
	
		
	################
	# Generate data
	################
	
	# Generate predictors
	x <- rnorm(d*p*n, sd=sqrt(1-rhox))
	if (rhox > 0) 
		x <- x + rnorm(1, sd=rhox) 
	dim(x) <- c(d,p,n)
	
	# Generate regression coefficients
	beta <- matrix(,p,n)
	start <- seq.int(1, by=floor(n/K), len=K)
	end <- c(start[-1]-1,n)
	for (k in 1:K) {
		val <- rnorm(p)
		val[sample(p,p*s)] <- 0
		beta[,start[k]:end[k]] <- val
	}
	
	# Generate responses
	y <- sapply(1:n, function(i) x[,,i] %*% beta[,i])
	if (sigma > 0) 
		y <- y + rnorm(d*n,sd=sigma)
	
	
	##############
	# Run methods
	##############
	
	L <- apply(x, 3, function(mat) svds(mat,1,0,0)$d[1])^2
	Lmax <- max(L)
	objective0 <- 0.5 * sum(y^2)
	result <- vector("list",length(method))
	names(result) <- method
	
	
	## Smooth proximal gradient
	print("Smooth proximal gradient")
	maxit.spg <- 1e4 #@@@@@@@@
	t1 <- proc.time()[3]
	beta <- null
	objective <- objective0
	o <- seq(floor(log10(objective))-1,-2,-1)
	eps <- c(10^(o[1]+1), rep(c(5,1),length(o)) * rep(10^o, each=2))
	l <- list()
	# Optimization with restarts
	for (i in seq_along(eps)) {
		tmp <- spg_cpp(x, y, lam1, lam2, beta,
			eps[i], Lmax, maxit.spg, verbose)
		l[[i]] <- tmp$trace
		if (tmp$objective <= objective) {
			objective <- tmp$objective
			beta <- tmp$beta } # else { break }
	}	
	t2 <- proc.time()[3]
	result[[1]] <- list(trace=unlist(l), time=t2-t1)
	
	
	## Primal-dual (Condat 2013, Algorithm 3.2)
	print("Primal-dual")
	maxit.pd <- 3000 #@@@@@@@@@
	tau <-10^seq(-6,6)
	objective <- objective0
	l <- list()
	t1 <- proc.time()[3]
	# Select tuning parameter tau
	for (i in seq_along(tau)) {
		tmp <- pd_condat(x, y, lam1, lam2, null, null, Lmax,
			tau[i], sig=0, rho=1.9, tol=-1, maxit=1e2, verbose)
		l[[i]] <- tmp$trace
		if (tmp$objective < objective) {
			beta <- tmp$beta
			ss <- tmp$s
			objective <- tmp$objective
			tau.best <- tmp$tau 
			sig.best <- tmp$sig } else { break }
	}
	# Optimization
	tmp <- pd_condat(x, y, lam1, lam2, beta, ss, Lmax, tau.best, 
		sig.best, rho=1.9, tol=-1, maxit=maxit.pd, verbose)
	t2 <- proc.time()[3]
	result[[2]] <- list(trace=c(unlist(l), tmp$trace), time=t2-t1)
	
	
	## Hybrid
	print("Hybrid")
	con <- list(maxit.bd=1e3, tol.bd=1e-6, tol.fista=1e-6, 
		maxit.pg=500, tol.pg=1e-6, tol.cp=1e-3, 
		maxit.fista=100)
	# Hybrid method with cyclical sweep
	tmp <- SGFL(x, y, lam1, lam2, intercept=FALSE, 
		beta=NULL, L=L, control=con, sweep="cyclical", 
		verbose=verbose)
	# Hybrid method with random sweep (SRSWOR)
	result[[3]] <- list(trace=tmp$trace, time=tmp$time)
	tmp <- SGFL(x, y, lam1, lam2, intercept=FALSE, 
		beta=NULL, L=L, control=con, sweep="srswor", 
		verbose=verbose)
	result[[4]] <- list(trace=tmp$trace, time=tmp$time)
	
	
	## ADMM
	print("ADMM")
	rho <- 10^seq(4,-4,-1)
	maxit.admm <- 2000 #@@@@@@@@
	objective <- objective0
	l <- list()
	# Select augmented Lagrangian parameter rho
	t1 <- proc.time()[3]
	for (i in seq_along(rho)) {
		tmp <- admm_cpp(x, y, null, null,null, lam1, lam2, rho[i], 
			Lmax, maxit_admm=1e2, tol_admm=-1, maxit_fista=8, 
			tol_fista=1e-8, verbose=verbose)
		l[[i]] <- tmp$trace
		if (tmp$objective > objective) break
		objective <- tmp$objective
		rho.best <- rho[i]
		beta <- tmp$beta
		z <- tmp$z
		u <- tmp$u 	
	}
	# Optimization
	tmp <- admm_cpp(x, y, beta, z, u, lam1, lam2, 
		rho.best, Lmax, maxit_admm=maxit.admm, tol_admm=-1, 
		maxit_fista=10, tol_fista=1e-8, verbose)
	t2 <- proc.time()[3]
	result[[5]] <- list(trace = c(unlist(l),tmp$trace), time=t2-t1)
	
	
	## LADMM
	print("LADMM")
	rho <- 10^seq(4,-4,-1)
	maxit.ladmm <- 3000 #@@@@@@@
	objective <- objective0
	l <- list()
	# Select augmented Lagrangian parameter rho
	t1 <- proc.time()[3]
	for (i in seq_along(rho)) {
		tmp <- ladmm_cpp(x, y, null, null, null, lam1, lam2, 
			rho[i], Lmax, maxit=1e2, tol=-1, verbose=verbose)
		l[[i]] <- tmp$trace
		if (tmp$objective > objective) break
		objective <- tmp$objective
		beta <- tmp$beta
		z <- tmp$z
		u <- tmp$u
		rho.best <- tmp$rho 
	}
	# Optimization
	tmp <- ladmm_cpp(x, y, beta, z, u, lam1, lam2, rho.best,
		Lmax, maxit=maxit.ladmm, tol=-1, verbose=verbose)
	t2 <- proc.time()[3]
	result[[6]] <- list(trace = c(unlist(l),tmp$trace), time=t2-t1)
	
	
	
	####################
	# Summarize results
	####################
	
	
	# Best performance
	objective.best <- min(unlist(lapply(result,"[[","trace"))) 
	
	for (i in seq_along(method))
	{
		# Relative distance to minimum 
		rel <-result[[i]]$trace/objective.best - 1
		# Timing 
		time <- result[[i]]$time
		niter <- length(rel)
		if (length(time) == 1)
			time <- seq(0, time, len=niter)
		idx <- min(min(which(rel <= tol)), niter)
		out[rr,paste0(method[i],".time")] <- time[idx]
		out[rr,paste0(method[i],".best")] <- min(rel)
	}

}




cat("Stop time","\n")
print(Sys.time())

