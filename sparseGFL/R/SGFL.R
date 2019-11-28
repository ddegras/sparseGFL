
SGFL <- function(x, y, lambda1, lambda2, alpha=1, 
   intercept=TRUE, beta=NULL, w=NULL, L=NULL, 
   sweep=c("cyclical","srswr","srswor"), pattern=c(10,5),
   control=list(), parallel=FALSE, verbose=FALSE)
{

   timing <- numeric()
   timing[1] <- proc.time()[3]
   
   ## Rename and/or redefine some parameters
   lam1 <- lambda1 * alpha
   lam2 <- lambda2
   lam3 <- lam1 * (1-alpha)
   
   ## Check data format, convert to matrices/arrays
   ## if needed [to do: call other function if sparse]
   if (is.list(x)) {
      if (verbose) 
         warning("Converting list 'x' to 3D-array")
      n <- length(x) 
      N <- NROW(x[[1]])
      p <- NCOL(x[[1]]) 
      x <- unlist(x)
      dim(x) <- c(N,p,n)
   } else {
      N <- dim(x)[1]
      p <- dim(x)[2]
      n <- dim(x)[3]   	
   }
   # Expand predictor array to accommodate intercept
   # if required
   if (intercept) {
      newx <- array(,dim=c(N,p+N,n))
      newx[,1:p,] <- x
      newx[,(p+1):(p+N),] <- diag(N)
      x <- newx
      p <- p+N
   }
   if (is.list(y)) {
      if (verbose)
         warning("Converting list 'y' to matrix")
      y <- unlist(y)
      dim(y) <- c(N,n)
   }
     
   ## Numerical controls
   con <- list(maxit.bd = 1000, tol.bd = 1e-7, 
      maxit.fista = 1000, tol.fista = 1e-8, 
      maxit.fp = 1e4, tol.fp = 1e-7, tol.cp = 1e-6,    
      maxit.pg = 100, tol.pg = 1e-5, tol.equal = 1e-8)    
   control.active <- intersect(names(control),names(con))
   con[control.active] <- control[control.active]
   for (name in names(con)) 
      assign(name,con[[name]])
   
   ## Parallel computing 
   PARALLEL <- (parallel && require(foreach) &&
      getDoParRegistered())  
   
   ## Check if package RSpectra is available for 
   ## quickly finding largest eigenvalue of matrix
   RSpectra.flag <- require(RSpectra) && (p >= 3) 
   
   ## Sweep pattern
   sweep <- match.arg(sweep)
    
   ## Initialize regression coefficients if needed
   if (is.null(beta))
      beta <- matrix(0,p,n)
         
   ## Empty matrix (for C++)   
   null <- matrix(0,0,1)
         
   ## Initialize weights if needed
   if (length(w) == 1) w <- rep_len(w,n-1)
   if (is.null(w)) w <- rep_len(1,n-1)      
         
   ## Lipschitz constants
   if (is.null(L)) {
	svdfun <- if (RSpectra.flag && N >= 3) {
   		function(x) (svds(x,1,0,0)$d[1])^2
   	} else { function(x) (svd(x)$d[1])^2 }	
		L <- apply(x,3,svdfun) + lam3
   }
   Lfusion <- matrix(,0,3)
   colnames(Lfusion) <- c("start","end","L")
   
   ## Activity status
   active <- rep(TRUE,n)   
   progressBlock <- TRUE
   progressFusion <- TRUE
   # Counter for convergence of fusion chains
   chain.counter <- 0 #@@@@@@@@@@
   # Critical value for chain convergence 
   chain.counter.c <- min(3,pattern[2]) # 5 #Inf #@@@@@@@@@@
   change <- NULL   
   
   ## Subgradient 
   g <- NULL
   
   ## Cycle counter
   count <- 1

   ## Objective value
   objective <- numeric()
   objective[1] <- objective_global(beta, x, y, 
      lam1, lam2, lam3, w, intercept)
   
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
       
      # convergence <- (chain.counter >= chain.counter.c) #@@@@@@

      for (j in 1:pattern[1])
      {
         if (pattern[1] < 1 || count > maxit.bd) {
         # if (pattern[1] < 1 || count > maxit.bd || convergence) {#@@@@ 
            progressBlock <- FALSE }
         if (!progressBlock) break 
         
         ## Increase counter
         count <- count + 1
         
         if (verbose) 
            cat("Descent cycle ",count-1," ",sep="")
                  
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
               bprev <- null
               wprev <- 0
            } else {
               bprev <- beta[,i-1]
               wprev <- w[i-1] }         
            if (i == n) {
               bnext <- null
               wnext <- 0
            } else {
               bnext <- beta[,i+1]
               wnext <- w[i] }                     
   
            ## Initialize activity flag to FALSE
            active[i] <- FALSE

            ## Simple solution for local objective?   
            beta.tmp <- trivial_local(beta[,i],
               x[,,i,drop=FALSE], y[,i,drop=FALSE],
               lam1, lam2, lam3, wprev, bprev, 
               wnext, bnext, intercept, tol.equal)
               
            if (length(beta.tmp) > 0) {
            	progress <- !all_equal(beta[,i,drop=FALSE],
            	   		beta.tmp,tol.equal)
            } else {           
            ## If not, perform FISTA
               result <- FISTA_1(beta[,i,drop=FALSE], 
                  x[,,i,drop=FALSE], y[,i,drop=FALSE], lam1, 
                  lam2, lam3, wprev, bprev, wnext, bnext,
                  L[i], intercept, maxit.fista, tol.fista,
                  maxit.fp, tol.fp, tol.equal)                
               beta.tmp <- result[["beta"]]
               progress <- result[["progress"]]
            }
            
            ## Update solution and activity status
            if (progress) {   
               beta[,i] <- beta.tmp
               active[max(i-1,1):min(i+1,n)] <- TRUE 
            }
 
         } # END for (i) loop   
         
         ## Monitor progress and display objective         
         objective[count] <- objective_global(
            beta, x, y, lam1, lam2, lam3, w, intercept)
         timing[count] <- proc.time()[3]
         # level[count] <- 1
         progressBlock <- (objective[count] <= 
           (1-tol.bd) * objective[count-1])                              
         if (verbose) 
            cat("Objective",objective[count],"\n")
         
         ## Display change points at end of block descent cycle
         if (verbose && (j == pattern[1] || !progressBlock)) {
            fused = compare_cols(beta,tol.equal)
            change <- which(c(FALSE,!fused)) 
            cat("Change Points:",change,"\n\n") }
                  
         if (progressBlock && (pattern[2] > 0)) 
         	{ progressFusion <- TRUE } # else break #@@@@@@@

        	

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
            fused <- compare_cols(beta,tol.cp)
            fusedwNext <- c(fused,FALSE)
            fusedwPrev <- c(FALSE,fused)
            change <- which(c(FALSE,!fused)) }
         change.old <- change       
             
         ## Increase counter
         count <- count + 1
        if (verbose) cat("Fusion cycle ",count-1," ", sep="")

         ## Sweep pattern
         sweep.index <- switch(sweep, cyclical = 1:(n-1), 
            srswr = sample(n-1,replace=TRUE), 
            srswor = sample(n-1))
         # interlaced = c(seq.int(1L,n-1,2L),
         # seq.int(2L,n-1,2L)),
                           
         for (i in sweep.index) 
         {
         # Note: according to above fusion rules, all fusions 
         # involving the last block (n) have already been 
         # tried before reaching block n, so there's no use 
         # inspecting it, hence the loop from 1 to (n-1) and
         # not from 1 to n 

            ## Determine which chain to optimize over
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

            ## Skip chains with inactive neighbors 
            ## (after 1st cycle)
            if (j > 1) {
               test1 <- (startc>1 && active[startc-1])
               test2 <- (endc<n && active[endc+1])
               test3 <- (startc==1 && endc==n)
               if (!(test1 || test2 || test3)) next }

            ## Lipschitz constant
            idx <- which(Lfusion[,"start"] == startc & 
               Lfusion[,"end"] == endc)
            if (length(idx)) {
            	   Lchain <- Lfusion[idx,"L"]
            } else {
               xchain <- aperm(x[,,chain,drop=FALSE],
                  perm=c(1,3,2))
               dim(xchain) <- c(N*flen,p)
               test <- (RSpectra.flag && N*flen >= 3)   
               Lchain <- ifelse(test, svds(xchain,1,0,0)$d^2, 
                  svd(xchain,0,0)$d[1]^2) + (lam3 * flen)
               Lfusion <- rbind(Lfusion,
                  t(c(startc,endc,Lchain)))
               rm(xchain)
            }
            
            ## Predecessor and successor
            if (startc == 1L) {
               bprev <- null
               wprev <- 0  
            } else {
               bprev <- beta[,startc-1L]
               wprev <- w[startc-1L]
            }
            if (endc == n) {
               bnext <- null
               wnext <- 0
            } else {
               bnext <- beta[,endc+1L] 
               wnext <- w[endc]
            }      
   
            ## Objective value for chain (before trying fusion)
            objective.before <- objective_blocks(
               beta[,chain,drop=FALSE], x[,,chain,drop=FALSE],
               y[,chain,drop=FALSE], lam1, lam2, lam3,
               w[startc:(endc-1L)], wprev, bprev, wnext, 
               bnext, intercept)
                  
            ## Trivial solution for local objective?
            beta.tmp <- trivial_local(beta=null,
               x=x[,,chain,drop=FALSE], y=y[,chain,drop=FALSE],
               lam1=lam1, lam2=lam2, lam3=lam3, wprev=wprev,
               bprev=bprev, wnext=wnext, bnext=bnext,
               intercept=intercept, tol=tol.equal)           
            if (length(beta.tmp) > 0) {
               objective.after <- objective_local(
                  beta.tmp, x[,,chain,drop=FALSE],
                  y[,chain,drop=FALSE], lam1, lam2, lam3,
                  wprev, bprev, wnext, bnext, intercept) 
            } else {
            ## If not, FISTA
               result <- FISTA_1(
                  beta[,i,drop=FALSE], x[,,chain,drop=FALSE], 
                  y[,chain,drop=FALSE], lam1, lam2, lam3,
                  wprev, bprev, wnext, bnext, Lchain,intercept,
                  maxit.fista, tol.fista, maxit.fp, tol.fp, 
                  tol.equal)
                 
               beta.tmp <- result[["beta"]]
               objective.after <- result[["objective"]] 
        	    }
      	                   
            ## If reduction in objective, perform fusion
            ## and update fusion & descent status 
            if (objective.after < objective.before) {
               active[max(startc-1,1):min(endc+1,n)] <- TRUE
               fusedwPrev[(startc+1):endc] <- TRUE
               fusedwNext[startc:(endc-1)] <- TRUE
               beta[,chain] <- beta.tmp 
            }
         
         } # End for (i) loop
               
         ## Update objective and display progress
         objective[count] <- objective_global(
            beta, x, y, lam1, lam2, lam3, w, intercept)
         timing[count] <- proc.time()[3]
         # level[count] <- 2 #@@@@@@@
         progressFusion <- (objective[count] <= 
            (1 - tol.bd) * objective[count-1]) 
         if (progressFusion && pattern[1] > 0) 
         	progressBlock <- TRUE           
         if (verbose) 
            cat("Objective",objective[count],"\n")

         ## Change points
         change <- which(c(FALSE,!fusedwPrev[-1])) 
         # Interrupt fusion cycle if no change point
         if (length(change) == 0) 
         	progressFusion <- FALSE

         ## Are the new change points identical 
         ## to the previous ones?
         if (identical(change.old,change)) { 
            chain.counter <- chain.counter + 1 
            # # progressFusion <- FALSE   #@@@@@@
         } else { 
            chain.counter <- 0  
            change.old <- change }

   		 if (chain.counter >= chain.counter.c) { #@@@@@@
             # convergence <- TRUE
             progressBlock <- progressFusion <- FALSE 
         } # else { convergence <- FALSE }
   
         ## Display exact change points at end of 
         ## fusion cycle iff required
         # last <- (j == pattern[2] || !progressFusion || convergence)
         last <- (j == pattern[2] || !progressFusion) #@@@@
         if (verbose && last) {
            fused <- compare_cols(beta,tol.equal)
            change <- which(c(FALSE,!fused)) 
            cat("Change Points:",change,"\n\n")            
            }
                            
         ## Break if insufficient progress 
         # if (last) break #@@@@@@@
         

      } # END FUSION LOOP 

      


      ##############################################
      # IF NO PROGRESS IN DESCENT AND FUSION CYCLES
      # APPLY FISTA TO ALL CHAINS SIMULTANEOUSLY
      ##############################################

        
      ## Check if algorithm has converged on a set of fusion chains
      # convergence <- (chain.counter >= chain.counter.c) 

      ## Check if there is a single chain
     if (is.null(change)) {
            fused = compare_cols(beta,tol.equal)
            change <- which(c(FALSE,!fused)) }
      single.chain <- (length(change) == 0) 

      ## Check progress of block descent and fusion stages
      progress <- (progressBlock || progressFusion)
         
      # if ((!progress || convergence) && !single.chain && 
         # count <= maxit.bd) { 
      if (!progress && !single.chain && count <= maxit.bd) {
   
         if (verbose) 
            cat("Descent (fixed chains) cycle ", count," ", 
               sep="")

         ## Increase counter and temporarily set objective 
         ## to previous value 
         objective[count+1] <- objective[count]
         count <- count + 1
         
         change.old <- change
         startc <- c(1L,change) 
         endc <- c(change-1L,n) 
         result <- FISTA_M(beta[,startc,drop=FALSE], 
            x, y, lam1, lam2, lam3, w, intercept,
            startc-1L, endc-1L, eta=1.25, maxit.fista,
            tol.fista, tol.equal)   
            
         ## Update objective and solution if needed
         if (result$objective < objective[count]) {
            startc <- as.integer(result$start) + 1L
            endc <- as.integer(result$end) + 1L 
            change <- startc[-1]  
            for (g in seq_along(startc)) 
               beta[,startc[g]:endc[g]] <- result$beta[,g]
            if (result$objective <= 
               ((1-tol.bd) * objective[count])) {
               # if (pattern[1] > 0)         #@@@@@@@ 	  
                   # progressBlock <- TRUE  #@@@@@@@
               # if (!identical(change.old,change)) #@@@@@@@
	              # progressFusion <- TRUE #@@@@@@@
	           # progress <- (progressBlock || progressFusion)
            }
            objective[count] <- result$objective 
         } #@@@@@@@@@
         chain.counter <- 0                     
         timing[count] <- proc.time()[3]
         # level[count] <- 3#@@@@@@@
         
         ## Display progress
         if (verbose) {
            cat(" Objective ",
               objective[count],"\n",sep="")
            if (!identical(change.old,change)) {
               cat("Change Points:",change,"\n\n") 
            } else cat("\n") 
         }
               
      }
      

      #####################################
      # CHECK GLOBAL CONVERGENCE 
      # APPLY SUBGRADIENT METHOD IF NEEDED
      #####################################
      
      # if ((!progress) || convergence || count > maxit.bd) {
      if ((!progress) || count > maxit.bd) {   #@@@@@@@
         ## Identify fusion chains
         if (verbose) cat("Checking optimality of solution\n")
         fused <- compare_cols(beta,tol.equal)
         change <- which(c(FALSE,!fused))
         startc <- c(1,change)
         endc <- c(change-1L,n)
         fchain <- mapply(seq.int, from=startc, to=endc,
            SIMPLIFY=FALSE)
            
         # Note: computing the subgradient of minimum norm 
         # w.r.t.all coordinate blocks reduces to computing 
         # the subgradient of minimum norm for each fusion 
         # chain and aggregating the results 
         
         ## Compute subgradient of minimum norm for each 
         ## chain by projected gradient method
         g <- if (PARALLEL) { 
            foreach(chain = fchain, .packages="sparseGFL") %dopar% 
               min.subg(chain, beta=beta, x=x, y=y,
                  lam1=lam1, lam2=lam2, lam3=lam3, w=w,
                  intercept=intercept, maxit=maxit.pg, 
                  tol=tol.pg)
         } else { 
            lapply(fchain, min.subg, beta=beta, x=x, 
               y=y, lam1=lam1, lam2=lam2, lam3=lam3, w=w,
               intercept=intercept, maxit=maxit.pg, 
               tol=tol.pg)
         }

         ## Combine subgradients
         g <- unlist(g)
         nrm <- mean(abs(g)) # subgradient norm
         dim(g) <- c(p,n)

         ## If not converged and maximum number of iterations 
         ## not reached yet, apply 1 step of subgradient
         ## method (steepest descent)
         if (nrm <= max(tol.equal,1e-5) || count > maxit.bd) 
            break 

         if (verbose) 
            cat("Subgradient method (variable chains) ",
               count," ",sep="")      
         objective[count+1] <- objective[count]
         count <- count + 1
            
         ## Line search for best step size
         LineSearch <- function(a) objective_global(
            beta-a*g, x, y, lam1, lam2, lam3, w, intercept)

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
            opt <- optimize(LineSearch, 
               c(grid[max(1,j-2)],grid[j]), tol=1e-8)
            if (opt$objective < objective[count]) {
               beta <- beta - opt$minimum * g
               objective[count] <- opt$objective
            } # else { beta <- beta - grid[j-1L] * g }  @@@@    
         }

         ## Update progress & activity status
         progress <- (objective[count] < (1-tol.bd) * 
              objective[count-1])
         if (progress) 
          { active <- rep(TRUE,n)             
            progressBlock <- TRUE
            progressFusion <- TRUE }
         timing[count] <- proc.time()[3]
         # level[count] <- 4 #@@@@@@@
            
         ## Enforce sparsity 
         beta.tmp <- beta
         beta.tmp[abs(beta.tmp) <= tol.equal] <- 0
         objective.tmp <- objective_global(
            beta.tmp, x, y, lam1, lam2, lam3, w, intercept)
         if (objective.tmp <= (1+tol.bd) * 
           objective[count]) {
            beta <- beta.tmp
            objective[count] <- objective.tmp
         }         
                                            
         ## Reset counter for chain convergence
         chain.counter <- 0 #@@@@@@
               
         ## Display progress
         if (verbose) 
            cat("Objective",objective[count],"\n\n")
                              
      } # END CHECK
   
      if (!progress || count > maxit.bd) break

   } # END OUTER LOOP
   
   
   ## Change points 
   fused <- compare_cols(beta,tol.cp)
   change <- which(c(FALSE,!fused))
   startc <- c(1,change)
   endc <- c(change-1L,n)
   nchain <- length(startc)

   ## Try to simplify solution 
   # Enforce strict block equality on fusion chains
   beta.tmp <- beta   
   for (k in 1:nchain) {
      if (startc[k] < endc[k]) 
          beta.tmp[,startc[k]:endc[k]] <- 
             rowMeans(beta.tmp[,startc[k]:endc[k],drop=FALSE])
   }

   # Set approximately null coefficients to zero
   beta.tmp[abs(beta.tmp) <= tol.equal] <- 0  

   ## Accept simplified solution if does not increase 
   ## objective too much 
   ##@@ Comment out for performance analysis @@##
   objective.tmp <- objective_global(
      beta.tmp, x, y, lam1, lam2, lam3, w, intercept)
   if (objective.tmp <= (1+tol.bd)*objective[count]) {
      beta <- beta.tmp
      objective[count] <- objective.tmp
   } else {
      fused <- compare_cols(beta,tol.equal)
      change <- which(c(FALSE,!fused)) }
   
   return(list(beta=beta, objective=objective[count], 
      changePoints=change, trace=objective, g=g, 
      iters=count-1, lambda1=lambda1, lambda2=lambda2,
      alpha=alpha, time=timing-timing[1]))
     
}   










###############################
# MINIMIZE SUBGRADIENT NORM
# TO CHECK OPTIMALITY OF SOLUTION
# MIN(U,V) || Z - U - VD' ||_F
###############################

# chain: fusion chain
# beta: matrix of regression coefficients  
#       columns of beta indexed by 'chain' should be
#       numerically equal
# x:    predictor array
# y:    response matrix
# lam1: regularization parameter for lasso penalty
# lam2: regularization parameter for total variation penalty
# lam1: regularization parameter for elastic net penalty
# w:    individual weights for total variation penalty
# intercept: logical flag (does regression model have intercept?)
# maxit: maximum number of iterations in projected gradient algorithm
# tol:   tolerance for convergence in project gradient 

min.subg <- function(chain, beta, x, y, lam1, lam2, lam3, 
   w, intercept, maxit, tol) 
{

   N <- nrow(y) # number of response variables 
   n <- ncol(y) # number of time points
   p <- nrow(beta) # number of predictor variables
   flen <- length(chain) # chain length
   
   if (chain[1] == 1L) {
      bprev <- matrix(0,1,0)
      wprev <- 0
   } else {
      bprev <- beta[,chain[1]-1L]
      wprev <- w[chain[1]-1L] }

   if (chain[flen] == n) {
      bnext <- matrix(0,1,0)
      wnext <- 0
   } else {
      bnext <- beta[,chain[flen]+1L]
      wnext <- w[chain[flen]] }

   ## Calculate fixed part of subgradient    
   bstart <- beta[,chain[1]]
   z <- grad_cpp(bstart, x[,,chain,drop=FALSE], 
      y[,chain,drop=FALSE], lam1, lam2, lam3,
      wprev, bprev, wnext, bnext, intercept)

   ## Find indices of fixed entries in subgradient 
   fixg <- (bstart != 0)
   if (intercept)
      fixg[(p-N+1):p] <- TRUE
   if (flen > 1)
      fixg <- rep.int(fixg,flen)
   fixg <- which(fixg)

   ## Optimize variable part of subgradient       
   if (flen == 1) { # case: single block
      # use soft-thresholding
      g <- numeric(p)
      over <- which(z > lam1)
      undr <- which(z < -lam1)
      g[over] <- z[over] - lam1
      g[undr] <- z[undr] + lam1
      g[fixg] <- z[fixg]
   } else { # case: chain of length > 1
      # use projected gradient
      # Indexing starts at 0 in C++
      g <- pg_cpp(z,fixg-1L,lam1,lam2,w,maxit,tol)
   }
         
   return(g)
}   



