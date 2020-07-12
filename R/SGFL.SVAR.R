
SGFL.SVAR <- function(x, p, lambda1, lambda2, alpha=1, struct=FALSE,
   intercept=TRUE, A=NULL, w=NULL, L=NULL, sweep=c("cyclical",
   "srswr","srswor"), pattern=c(10,5), control=list(), 
   parallel=FALSE, verbose=FALSE)
{

   timing <- proc.time()[3]
   
   ## Rename and/or redefine some parameters
   lam1 <- lambda1 * alpha
   lam2 <- lambda2
   lam3 <- lambda1 * (1-alpha)
   lag0 <- struct
   
   ## Data dimensions 
   if (is.vector(x))  # univariate time series
      x <- t(x) 
   n <- ncol(x) # time series length (original)
   N <- nrow(x) # time series dimension
      
   ## Model lags
   stopifnot(p >= 0)
   p <- as.integer(p)
   lags <- if (lag0) 0:p else 1:p 
   nlags <- p + lag0
   
   ## Check dimensions of initial solution if provided
   dimA <- c(N, N*nlags+intercept, n-p)
   if (is.null(A)) {      
      A <- array(0, dimA)
   } else {
      # Case: univariate time series with only one predictor
      # per time point (A may be specified as a vector)
      test1 <- (dimA[1] == 1 && dimA[2] == 1 && 
         is.vector(A) && length(A) == dimA[3])
      # Case: univariate time series 
      # (A may be specified as a matrix)   
      test2 <- (dimA[1] == 1 && is.matrix(A) && 
         nrow(A) == dimA[2] && ncol(A) == dimA[3])
      # General case
      test3 <- (is.array(A) && length(dim(A)) == 3 &&
         all(dim(A) == dimA))
      if (!(test1 || test2 || test3))
         stop(paste("'A' must be an array of dimensions",
         "N x (N*l+i) x (n-p) with N = number of observed",
         "variables, l = number of lags in model,",
         "i = 1 if model contains intercept, 0 otherwise,",
         "n = number of time points, and p = maximum lag"))
   }
   A.singleton.dim <- any(dimA == 1L)

   # ??? Maybe write faster functions for N=1 ???

   ## Numerical controls
   con <- list(maxit.bd = 1000, tol.bd = 1e-8, 
      maxit.fista = 1000, tol.fista = 1e-8, 
      maxit.fp = 1e4, tol.fp = 1e-8, tol.cp = 1e-3, 
      maxit.pg = 500, tol.pg = 1e-8, tol.equal = 1e-8)
   control.active <- intersect(names(control),names(con))
   con[control.active] <- control[control.active]
   for (name in names(con)) 
      assign(name,con[[name]])
   
   ## Parallel computing 
   PARALLEL <- (parallel && require(foreach) &&
      getDoParRegistered())  
      
   ## Sweep pattern
   sweep <- match.arg(sweep)
   
   ## Set up regression matrices
   if (p == 0 && !intercept)  { 
      # Case: structural model without intercept
      # prediction vector x(t) = response vector y(t) 
      # (with constraint diag(A(t)) = 0)
      y <- x
   } else { # SVAR
      # Response matrix
      y <- x[,(p+1):n,drop=FALSE]
      # Predictor matrix
      idx <- rep((p+1):n,each=nlags) - lags
      x <- x[,idx,drop=FALSE]
      dim(x) <- c(N*nlags,n-p)
      if (intercept)
         x <- rbind(x,rep(1,n-p))
      # sapply((lagmax +1):n, function(i) as.vector(x[,i-lags])) 
   }
   n <- n - p # length of new time series
      
   ## Empty matrix (for C++)   
   null <- matrix(0,0,1)
         
   ## Initialize weights if needed
   if (length(w) == 1) w <- rep_len(w,n-1)
   if (is.null(w)) w <- rep_len(1,n-1)      
         
   ## Lipschitz constants
   if (is.null(L)) L <- colSums(x^2) + lam3
   
   ## Activity status
   active <- rep(TRUE,n)   
   progressBlock <- TRUE 
   progressFusion <- TRUE 
   # Counter for convergence of fusion chains
   chain.counter <- 0
   # Critical value for chain convergence
   chain.counter.c <- min(3,pattern[2]) 
   change <- NULL   
   
   ## Subgradient 
   g <- NULL
   
   ## Cycle counter
   count <- 1

   ## Objective value
   objective <- numeric()
   objective[1] <- 
      objective_global_svar(A,x,y,lam1,lam2,lam3,w,intercept)
   timing[1] <- proc.time()[3]
   ## Display initial objective value
   if (verbose)
      cat("Initial objective ",objective[1],"\n",sep="")
   
   
      
   #############
   # OUTER LOOP
   #############

   repeat
   {

      ################
      # DESCENT CYCLE
      ################
       
      # convergence <- (chain.counter >= chain.counter.c)

      for (j in 1:pattern[1])
      {
            
         if (pattern[1] < 1 || count > maxit.bd) 
            progressBlock <- FALSE
         if (!progressBlock) break
         
         if (verbose) cat("Descent cycle ",count," ",sep="")
                  
         ## Increase counter
         count <- count + 1

         ## Set objective to value from previous cycle
         objective[count] <- objective[count-1]
         
         ## Sweep pattern
         sweep.index <- switch(sweep, cyclical = 1:n, 
         # interlaced = c(seq.int(1L,n,2L),seq.int(2L,n,2L)),
         srswr = sample(n,replace=TRUE), srswor = sample(n))

         for (i in sweep.index) 
         {

            ## Skip if inactive block
            if (!active[i]) next
   
            ## Set predecessor and successor solution
            if (i == 1) {
               Aprev <- null
               wprev <- 0
            } else {
               Aprev <- A[,,i-1]
               if (A.singleton.dim)
                  dim(Aprev) <- dimA[1:2]
               wprev <- w[i-1] }         
            if (i == n) {
               Anext <- null
               wnext <- 0
            } else {
               Anext <- A[,,i+1]
               if (A.singleton.dim)
                  dim(Anext) <- dimA[1:2]
               wnext <- w[i] }                     
   
            ## Initialize activity flag to FALSE
            active[i] <- FALSE

            ## Trivial solution for local objective?   
            A.tmp <- trivial_local_svar(x[,i,drop=FALSE], 
               y[,i,drop=FALSE], lam1, lam2, lam3, wprev, Aprev, 
               wnext, Anext, lag0, intercept, tol.equal)
            if (length(A.tmp) > 0) {
            	   progress <- !all_equal(A[,,i],A.tmp,tol.equal)
            } else {           
            ## If not, perform FISTA
               result <- FISTA_1_svar(A[,,i], x[,i,drop=FALSE],
                  y[,i,drop=FALSE], lam1, lam2, lam3, wprev, Aprev, 
                  wnext, Anext, L[i], lag0, intercept, 
                  maxit.fista,tol.fista, maxit.fp, tol.fp, tol.equal)
               A.tmp <- result[["A"]]
               progress <- result[["progress"]]
            }
            
            ## Update solution and activity status
            if (progress) {   
               A[,,i] <- A.tmp
               active[max(i-1,1):min(i+1,n)] <- TRUE 
            }
 
         } # END for (i) loop   
         
         ## Monitor progress and display objective         
         objective[count] <- 
            objective_global_svar(A,x,y,lam1,lam2,lam3,w,intercept)
   		 timing[count] <- proc.time()[3]

         progressBlock <- 
            (objective[count] <= (1-tol.bd) * objective[count-1])                   
         if (verbose) 
            cat("Objective",objective[count],"\n")
         
         ## Display change points at end of block descent cycle
         if (verbose && (j == pattern[1] || !progressBlock)) {
            fused = compare_slices(A,tol.equal)
            change <- which(c(FALSE,!fused)) 
            cat("Change Points:",change+p,"\n\n") }
                  
         ## Break if insufficient progress  
         # if (!progressBlock) break 
         if (progressBlock && pattern[2] > 0) { #@@@@@@@@@@
         	progressFusion <- TRUE } # else { break } 

      } # END for (j) descent cycle   



      ###############
      # FUSION CYCLE 
      ###############
            

      for (j in 1:pattern[2])
      {
      
         if (pattern[2] < 1 || count > con$maxit.bd)  
         	progressFusion <- FALSE
         if (!progressFusion) break 
         
         ## Change points (first cycle)
         if (j == 1) {
            fused <- compare_slices(A,tol.cp)
            fusedwNext <- c(fused,FALSE)
            fusedwPrev <- c(FALSE,fused)
            change <- which(c(FALSE,!fused)) }
         change.old <- change       
             
         ## Increase counter
         count <- count + 1

         if (verbose) cat("Fusion cycle ",count-1," ", sep="")

         ## Sweep pattern
         sweep.index <- switch(sweep, cyclical = 1:(n-1), 
         # interlaced = c(seq.int(1L,n-1,2L),seq.int(2L,n-1,2L)),
            srswr = sample(n-1,replace=TRUE), srswor = sample(n-1))
                           
         for (i in sweep.index) 
         {
            # Note: according to above fusion rules, all fusions 
            # involving the last block (n) have already been tried 
            # before reaching block n, so there's no use inspecting it.
            # Hence the loop from 1 to (n-1) and not 1 to n. 

            ## Identify chain to optimize
            # Case 1 (interior of chain)      
            if (fusedwPrev[i] & fusedwNext[i]) next    
            # Case 2 (start of chain of length > 1)
            if ((!fusedwPrev[i]) & fusedwNext[i]) {
               # first index in current chain
               startc <- i 
               # last index in current chain
               endc <- i + which.min(fusedwNext[(i+1):n]) } 
            # Case 3 (end of chain or singleton)
            if (fusedwPrev[i] || (!fusedwNext[i])) {
               # first index in extended chain
               startc <- i + 1L - which.min(fusedwPrev[i:1])
               # last index in extended chain 
               endc <- i + which.min(fusedwNext[(i+1):n]) } 
            chain <- startc:endc
            flen <- length(chain) 

            ## Skip singletons
            if (startc == endc) next

            ## Skip chains with inactive neighbors (after 1st cycle)
            if (j > 1) {
               test1 <- (startc>1 && active[startc-1])
               test2 <- (endc<n && active[endc+1])
               test3 <- (startc==1 && endc==n)
               if (!(test1 || test2 || test3)) next }

            ## Lipschitz constant
            # test <- (RSpectra.flag && flen >= 3 && nrow(x) >= 3)   
            # L.chain <- ifelse(test, svds(x[,chain],1,0,0)$d^2, 
               # svd(x,0,0)$d[1]^2) 
			L.chain <- (norm(x[,chain],"2"))^2 + lam3 * flen
            
            ## Predecessor and successor
            if (startc == 1L) {
               Aprev <- null
               wprev <- 0  
            } else {
               Aprev <- A[,,startc-1L]
               if (A.singleton.dim)
                  dim(Aprev) <- dimA[1:2]
               wprev <- w[startc-1L]
            }
            if (endc == n) {
               Anext <- null
               wnext <- 0
            } else {
               Anext <- A[,,endc+1L] 
               if (A.singleton.dim)
                  dim(Anext) <- dimA[1:2]
               wnext <- w[endc]
            }      
   
            ## Objective value for chain (before trying fusion)
            objective.before <- 
               objective_blocks_svar(A[,,chain,drop=FALSE],
               x[,chain,drop=FALSE], y[,chain,drop=FALSE], 
               lam1, lam2, lam3, w[startc:(endc-1L)], wprev, 
               Aprev, wnext, Anext, intercept)
                  
            ## Trivial solution for local objective?
            A.tmp <- trivial_local_svar(x=x[,chain,drop=FALSE],
               y=y[,chain,drop=FALSE], lam1=lam1, lam2=lam2,
               lam3=lam3, wprev=wprev, Aprev=Aprev, wnext=wnext, 
               Anext=Anext, lag0=lag0, intercept=intercept, 
               tol=tol.equal)           
            if (length(A.tmp) > 0) {
               objective.after <- objective_local_svar(A.tmp,
                  x[,chain,drop=FALSE], y[,chain,drop=FALSE], lam1, 
                  lam2, lam3, wprev, Aprev, wnext, Anext, intercept) 
            } else {
            ## If not, FISTA
               result <- FISTA_1_svar(A[,,i], x[,chain,drop=FALSE], 
                  y[,chain,drop=FALSE], lam1, lam2, lam3, wprev, Aprev, 
                  wnext, Anext, L.chain, lag0, intercept, maxit.fista,
                  tol.fista, maxit.fp, tol.fp, tol.equal)
               A.tmp <- result[["A"]]
               objective.after <- result[["objective"]] 
        	    }
      	                   
            ## If reduction in objective, perform fusion and update
            ## fusion & descent status 
            if (objective.after < objective.before) {
               active[max(startc-1,1):min(endc+1,n)] <- TRUE
               fusedwPrev[(startc+1):endc] <- TRUE
               fusedwNext[startc:(endc-1)] <- TRUE
               A[,,chain] <- A.tmp 
            }
         
         } # End for (i) loop
               
         ## Update objective and display progress
         objective[count] <- 
            objective_global_svar(A,x,y,lam1,lam2,lam3,w,intercept)
   		 timing[count] <- proc.time()[3]
         progressFusion <- (objective[count] <= 
            (1 - tol.bd) * objective[count-1]) 
         if (progressFusion && pattern[1] > 0) 
         	progressBlock <- TRUE           
         if (verbose) 
            cat("Objective",objective[count],"\n")

         ## Change points
         change <- which(c(FALSE,!fusedwPrev[-1])) 
         if (length(change) == 0) 
         	progressFusion <- FALSE
     
         ## Are the new change points identical to the previous ones?
         if (identical(change.old,change)) {
	         chain.counter <- chain.counter+1
	     } else { 
            chain.counter <- 0 
            change.old <- change } 

   		 if (chain.counter >= chain.counter.c)  
             progressBlock <- progressFusion <- FALSE 

   
         ## Display exact change points at end of fusion cycle iff required
         last <- (j == pattern[2] || !progressFusion)
         if (verbose && last) {
            fused <- compare_slices(A,tol.equal)
            change <- which(c(FALSE,!fused)) 
            cat("Change Points:",change+p,"\n\n") }
          


      } # END FUSION LOOP 

      


      ##############################################
      # IF NO PROGRESS IN DESCENT AND FUSION CYCLES
      # APPLY FISTA TO ALL CHAINS SIMULTANEOUSLY
      ##############################################

         
      ## Check if algorithm has converged on a set of fusion chains
      # convergence <- (chain.counter >= chain.counter.c)      

      ## Check if there is a single chain
      if (is.null(change)) {
            fused = compare_slices(A,tol.equal)
            change <- which(c(FALSE,!fused)) }
      single.chain <- (length(change) == 0) 

      ## Check on progress of block descent and fusion stages
      progress <- (progressBlock || progressFusion)
         
      if (!progress && !single.chain && count <= maxit.bd) {
   
         if (verbose) 
            cat("Descent (fixed chains) cycle ", count," ", sep="")

         ## Increase counter and temporarily set objective to previous value 
         count <- count + 1
         objective[count] <- objective[count-1]
         
         startc <- c(1L,change) 
         endc <- c(change-1L,n) 
         result <- FISTA_M_svar(A[,,startc,drop=FALSE], 
            x, y, lam1, lam2, lam3, w, lag0, intercept, 
            startc-1L, endc-1L, eta=1.25, maxit.fista, 
            tol.fista, tol.equal)   
            
         ## Update objective and solution if needed   
         if (result$objective < objective[count]) {
            objective[count] <- result$objective 
            startc <- as.integer(result$start) + 1L
            endc <- as.integer(result$end) + 1L 
            change <- startc[-1]  
            for (g in seq_along(startc)) 
               A[,,startc[g]:endc[g]] <- result$A[,,g]
            # progress <- progressBlock <- progressFusion <- TRUE
         } 
         chain.counter <- 0                     
   		 timing[count] <- proc.time()[3]
      
         ## Display progress
         if (verbose) {
            cat(" Objective ",objective[count],"\n",sep="")
            if (!identical(change.old,change)) {
               cat("Change Points:",change+p,"\n\n") } else cat("\n") }
               
      }
      

      #####################################
      # CHECK GLOBAL CONVERGENCE 
      # APPLY SUBGRADIENT METHOD IF NEEDED
      #####################################
      
      if ((!progress) || count > maxit.bd) {
         
         ## Identify fusion chains
         if (verbose) cat("Checking optimality of solution\n")
         fused <- compare_slices(A,tol.equal)
         change <- which(c(FALSE,!fused))
         startc <- c(1,change)
         endc <- c(change-1L,n)
         fchain <- mapply(seq.int,from=startc,to=endc,SIMPLIFY=FALSE)
            
         # Note: computing the subgradient of minimum norm w.r.t. all
         # coordinate blocks reduces to computing the subgradient of 
         # minimum norm for each fusion chain and aggregating the results 
         
         ## Compute subgradient of minimum norm for each chain 
         ## by projected gradient method
         g <- if (PARALLEL) { 
            foreach(chain = fchain, .packages="sparseGFL") %dopar% 
               min.subg.SVAR(chain, A=A, x=x, y=y, lam1=lam1, lam2=lam2,
               lam3=lam3, w=w, lag0=lag0, intercept=intercept, 
               maxit=maxit.pg, tol=tol.pg)
         } else { 
            lapply(fchain, min.subg.SVAR, A=A, x=x, y=y, lam1=lam1,
               lam2=lam2, lam3=lam3, w=w, lag0=lag0, 
               intercept=intercept, maxit=maxit.pg, tol=tol.pg)
         }
         
         ## Combine subgradients
         g <- unlist(g)
         nrm <- mean(abs(g)) # subgradient norm
         dim(g) <- dimA 
         ## If not converged and maximum number of iterations not reached yet, 
         ## apply 1 step of subgradient method (steepest descent)
         if (nrm <= tol.equal || count > maxit.bd) break 
         # if (nrm <= max(tol.equal,1e-5) || count > maxit.bd) break 

         if (verbose) 
            cat("Subgradient method (variable chains) ",
            	count," ",sep="")      
         count <- count + 1
         objective[count] <- objective[count-1]
            
         ## Line search for best step size
         LineSearch <- function(a) objective_global_svar(A-a*g,
            x,y,lam1,lam2,lam3,w,intercept)

         ## Grid search
         grid <- 10^seq.int(-9,1)
         ngrid <- length(grid)
         for (j in 1:ngrid) {
            objective.tmp <- LineSearch(grid[j]) 
            if (objective.tmp < objective[count]) {
               objective[count] <- objective.tmp
            } else break }

         ## Line search
         if (j > 1) {
            opt <- optimize(LineSearch, tol=1e-8,
               interval=c(grid[max(1,j-2)],grid[j]))
            if (opt$objective < objective[count]) {
               A <- A - opt$minimum * g
               objective[count] <- opt$objective
            } # else { A <- A - grid[j-1L] * g }      
         }


         ## Update progress & activity status
         progress <- (objective[count] < 
            (1-tol.bd) * objective[count-1])
         if (progress) {
         	active <- rep(TRUE,n) 
         	progressBlock <- TRUE
         	progressFusion <- TRUE }          
   		 timing[count] <- proc.time()[3]
         
         ## Reset counter for chain convergence
         chain.counter <- 0
         
         ## Enforce sparsity
         A.tmp <- A
         A.tmp[abs(A.tmp) <= tol.equal] <- 0
         objective.tmp <- objective_global_svar(A.tmp,x,y,
            lam1,lam2,lam3,w,intercept)
         if (objective.tmp <= (1+tol.bd)*objective[count]) 
            { A <- A.tmp
            objective[count] <- objective.tmp }         
               
         ## Display progress
         if (verbose) 
            cat("Objective",objective[count],"\n\n")
                              
      } # END CHECK
   
      if ((!progress) || count > maxit.bd) break

   } # END OUTER LOOP
   
   
   ## Change points 
   fused <- compare_slices(A,tol.cp)
   change <- which(c(FALSE,!fused))
   startc <- c(1,change)
   endc <- c(change-1L,n)
   nchain <- length(startc)

   ## Try to simplify solution 
   # Enforce block equality on fusion chains
   A.tmp <- A   
   dim(A.tmp) <- c(dimA[1]*dimA[2],dimA[3])
   for (k in 1:nchain) {
      if (startc[k] < endc[k]) 
          A.tmp[,startc[k]:endc[k]] <- 
             rowMeans(A.tmp[,startc[k]:endc[k],drop=FALSE])
   }

   # Set approximately null coefficients to zero
   A.tmp[abs(A.tmp) <= tol.equal] <- 0  

   ## Accept simplified solution if does not increase 
   ## objective too much 
   dim(A.tmp) <- dimA
   objective.tmp <- objective_global_svar(
      A.tmp,x,y,lam1,lam2,lam3,w,intercept)
   if (objective.tmp <= (1+tol.bd)*objective[count]) {
      A <- A.tmp
      objective[count] <- objective.tmp
   } else {
      fused <- compare_slices(A,tol.equal)
      change <- which(c(FALSE,!fused)) }
   
   names(timing) <- NULL
   
   return(list(A=A, objective=objective[count],
   		changePoints=change+p,trace=objective, g=g, 
   		iters=count-1, lambda1=lambda1, lambda2=lambda2, 
   		alpha=alpha, time=timing-timing[1]))
}   
























###############################
# MINIMIZE SUBGRADIENT NORM
# TO CHECK OPTIMALITY OF SOLUTION
# MIN(U,V) || Z - U - VD' ||_F
###############################

# chain: fusion chain
# A:    3D-array of regression coefficients  
#       Slices of A indexed by 'chain' should be numerically equal
# x:    predictor matrix
# y:    response matrix
# lam1: regularization parameter for lasso penalty
# lam2: global regularization parameter for total variation penalty
# w:    individual weights for total variation penalty
# lag0: logical flag (does model have instantaneous interactions?)
# intercept: logical flag (does regression model have intercept?)
# maxit: maximum number of iterations in projected gradient algorithm
# tol:   tolerance for convergence in project gradient 

min.subg.SVAR <- function(chain, A, x, y, lam1, lam2, lam3, 
   w, lag0, intercept, maxit, tol) 
{

   N <- dim(A)[1] # number of response variables 
   m <- dim(A)[2] # number of predictor variables
   n <- dim(A)[3] # number of time points
   flen <- length(chain) # chain length
   
   if (chain[1] == 1L) {
      Aprev <- matrix(0,1,0)
      wprev <- 0
   } else {
      Aprev <- A[,,chain[1]-1L]
      if (!is.matrix(Aprev))
         dim(Aprev) <- c(N,m) 
      wprev <- w[chain[1]-1L] }

   if (chain[flen] == n) {
      Anext <- matrix(0,1,0)
      wnext <- 0
   } else {
      Anext <- A[,,chain[flen]+1L]
      if (!is.matrix(Anext))
         dim(Anext) <- c(N,m) 
      wnext <- w[chain[flen]] }

   ## Calculate fixed part of subgradient    
   Astart <- A[,,chain[1]]
   if (!is.matrix(Astart))
         dim(Astart) <- c(N,m) 
   z <- grad_svar(Astart, x[,chain,drop=FALSE], 
      y[,chain,drop=FALSE],lam1,lam2,lam3,wprev,Aprev,
      wnext,Anext,lag0,intercept)

   ## Find indices of fixed entries in subgradient 
   fixg <- (Astart != 0)
   if (intercept)
      fixg[,m] <- TRUE
   if (lag0)
      diag(fixg) <- TRUE
   if (flen > 1)
      fixg <- rep.int(fixg,flen)
   fixg <- which(fixg)

   ## Optimize variable part of subgradient       
   if (flen == 1) { # case: single block
      # use soft-thresholding
      g <- numeric(m*N)
      abv <- which(z > lam1)
      blo <- which(z < -lam1)
      g[abv] <- z[abv] - lam1
      g[blo] <- z[blo] + lam1
      g[fixg] <- z[fixg]
   } else { # case: chain of length > 1
      # use projected gradient
      # Indexing starts at 0 in C++
      g <- pg_cpp(z,fixg-1L,lam1,lam2,w,maxit,tol)
   }
         
   return(g)
}   



