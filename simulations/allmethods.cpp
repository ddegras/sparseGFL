// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppArmadillo.h>
#include <progress.hpp>
using namespace std;
using namespace Rcpp;
using namespace arma;





///////////////////////////////////////////////////////////////////////////////////




// GLOBAL OBJECTIVE FUNCTION
// [[Rcpp::export]]
double objective_global (const mat & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2) 
{
   
   uword i, n = y.n_cols;
   double sse = 0.0, lasso, tv = 0.0;

   // Sum of squared errors
   for (i=0; i<n; i++) 
      sse += pow(arma::norm(y.col(i)-x.slice(i)*beta.col(i)),2); 
   sse *= 0.5;

   // l1 penalty
   lasso = lam1 * arma::accu(arma::abs(beta)); 

   // TV penalty
   for(i=0; i<(n-1); i++)
      tv += arma::norm(beta.col(i+1)-beta.col(i));
   tv *= lam2; 

   return (sse+lasso+tv);
}

   

// SOFT-THRESHOLDING
// [[Rcpp::export]]
mat soft_cpp(const mat & x, double lam)
{
   mat out = zeros<mat>(size(x));
   uword i, n = x.n_elem;
   double nlam = -lam;
   for (i=0; i<n; i++) {
      if (x[i] > lam) { out[i] = x[i] - lam; }
      else if (x[i] < nlam) { out[i] = x[i] + lam; } 
   }
   return out; 
}



///////////////////////////////////////////////////////////////////////////////////



    
// PARTIAL OBJECTIVE FUNCTION FOR ADMM (FOR LASSO)
// [[Rcpp::export]]
double objective_admm(const mat & beta, const cube & x, const mat & y, 
   const mat & zpu, const double rho, const double lam1)
{

   int i, n = y.n_cols;
   double primal = 0.0, dual = 0.0;

   for (i=0; i<n; i++) {
      primal += accu(pow(y.col(i)-x.slice(i)*beta.col(i),2));
      if (i < n-1)
         dual += accu(pow(zpu.col(i)+beta.col(i)-beta.col(i+1),2));
   }

   primal = 0.5 * primal + lam1 * accu(abs(beta));
   dual *= 0.5 * rho;
   return (primal+dual); 
} 






// FISTA FOR LASSO PROBLEM (ADMM)
// [[Rcpp::export]]
mat lasso_fista(const cube & x, const mat & y, mat beta, const mat & z, 
   const mat & u, const double lam1, const double rho, const double L, 
   const int maxit, const double tol)
{
   
   uword n = beta.n_cols; // number time points
   uword p = beta.n_rows; // number variables
   uword i, j; // counters
   vec grad(p); // gradient
   mat zpu = z + u;
   vec Delta(p), Delta_old(p); // increments of solution
   mat c(p,n); // shift for gradient  
   c.col(0) = zpu.col(0);
   for (j=1; j<(n-1); j++)
      c.col(j) = zpu.col(j)-zpu.col(j-1);
   c.col(n-1) = - zpu.col(n-2);
   mat beta_prev(p,n); // previous solution
   mat beta_adj = beta; // adjusted solution
   double theta = 1.0, theta_prev; // acceleration parameter   
   double objective = objective_admm(beta,x,y,zpu,rho,lam1), objective_prev; 

   // FISTA
   for (i=0; i<maxit; i++) {
      
      // Store previous results
      beta_prev = beta;
      theta_prev = theta;
      objective_prev = objective;  
 
      // Gradient step
      Delta.zeros();
      for (j=0; j<n; j++) {
         Delta_old = Delta;
         if (j < (n-1)) {
            Delta = beta_adj.col(j+1) - beta_adj.col(j);
         } else { Delta.zeros(); }
         grad = trans(x.slice(j)) * (x.slice(j) * beta_adj.col(j)  - y.col(j))  + 
            rho * (c.col(j) + Delta_old - Delta);
         beta.col(j) = beta_adj.col(j) - grad / L;
      }         

      // Proximal step
      beta = soft_cpp(beta,lam1/L);
    
      // Nesterov acceleration
      theta = 0.5 * (1.0 + sqrt(1.0 + 4.0 * theta*theta));
      beta_adj = beta + ((theta_prev-1.0) / theta) * (beta - beta_prev);
            
      // Update objective and check convergence 
      objective = objective_admm(beta,x,y,zpu,rho,lam1);
      if (objective >= (1-tol) * objective_prev) break;
   }

   return beta;   
}



// ADMM: MAIN FUNCTION
// [[Rcpp::export]]
Rcpp::List admm_cpp(const cube & x, const mat & y, mat beta, mat z, 
   mat u, const double lam1, const double lam2, double rho, double Lx, 
   const int maxit_admm, const double tol_admm, const int maxit_fista, 
   const double tol_fista, const bool verbose)
{
   uword i, j; 
   uword n = y.n_cols;
   uword p = x.n_cols;

   // Lipschitz constant
   double nrm;   
   if  (Lx == 0.0) {
      for (j=0; j<n; j++) 
         {nrm = max(svd(x.slice(j)));
         if (nrm > Lx) Lx = nrm;}
      Lx = Lx*Lx;
   } 
   double L = Lx + 4 * rho;

   // Initialization
   if (beta.is_empty()) beta = zeros<mat>(p,n);
   mat beta_best = beta;
   if (z.is_empty()) z = diff(beta,1,1);
   mat z_best = z;
   if (u.is_empty()) u = zeros<mat>(p,n-1);
  mat u_best = u;
   vec objective(maxit_admm+1);
   objective(0) = objective_global(beta,x,y,lam1,lam2);
   double objective_best = objective(0);
   if (verbose) 
      Rprintf("Iteration 0 Objective = %f\n",objective(0));
   double lam2_rho = lam2/rho; 
   vec zj(p); 
   vec Delta(p), Delta_prev; 
// double nrm_prim, nrm_dual;
   Progress prog(0, false); 

   for (i=1; i<=maxit_admm; i++) {

      // Update primal variable beta
      beta = lasso_fista(x,y,beta,z,u,lam1,rho,L,
         maxit_fista,tol_fista);
    
      // Update primal variable z and dual variable u
   // if (adaptive) z_prev = z;
      z.zeros();
      for (j=0; j<(n-1); j++) {
         Delta = beta.col(j+1) - beta.col(j);
         zj = Delta - u.col(j);
         nrm = norm(zj);
         if (nrm > lam2_rho) 
            z.col(j) = (1.0 - (lam2_rho/nrm)) * zj; 
         u.col(j) += z.col(j) - Delta;
      }    

      // Adjust rho if needed
      // if (adaptive && i <= 100) { 
         // // Calculate residual norms
         // nrm_prim = norm(z-diff(beta,1,1),"fro");
         // nrm_dual = 0.0;
         // Delta.zeros();
         // for (j=0; j<n; j++) {
            // Delta_prev = Delta;
            // if (j<n-1) {
               // Delta = z.col(j) - z_prev.col(j);
            // } else { Delta.zeros(); }
            // nrm_dual += accu(pow(Delta-Delta_prev,2));
         // }
         // nrm_dual = rho * sqrt(nrm_dual);   
         // if (nrm_prim > 10 * nrm_dual) { rho *= 2; }
         // else if (nrm_dual > 10 * nrm_prim) { rho *= 0.5; }
      // }

      // Update objective
      objective(i) = objective_global(beta,x,y,lam1,lam2); 
      if (objective(i) < objective_best)
         { objective_best = objective(i);
         beta_best = beta; u_best = u; z_best = z; } 
      if (verbose && i % 100 == 0) 
            Rprintf("Iteration %d Objective = %f\n",i,objective(i));

      // Check convergence and user interruption
      if (i % 10 == 0 && Progress::check_abort()) 
         break;
      if (abs(objective(i)-objective(i-1)) <= tol_admm * objective(i-1)) 
         break;
   }

   if (i > maxit_admm) i = maxit_admm;

   // Return primal solution, dual solution, and objective
   return List::create(Named("beta") = beta_best,
                        Named("z") = z_best,
                        Named("u") = u_best,
                        Named("objective") = objective_best,
                        Named("trace") = objective.head(i+1),
                        Named("rho") = rho);

}





///////////////////////////////////////////////////////////////////////////////////




// LADMM
// [[Rcpp::export]]
List ladmm_cpp(const cube & x, const mat & y, mat beta, mat z, mat u, 
   const double lam1, const double lam2, double rho, double Lx, 
   const int maxit, const double tol, const bool verbose)
{
   uword i, j; 
   uword n = y.n_cols;
   uword p = x.n_cols;

   // Lipschitz constant 
   double nrm;
   if (Lx == 0.0) {
      for (j=0; j<n; j++) 
         {nrm = max(svd(x.slice(j)));
      if (nrm > Lx) Lx = nrm;}
      Lx = Lx*Lx;
   }
   double L = Lx + 4 * rho;

   // Initialization
   if (beta.is_empty()) beta = zeros<mat>(p,n);
   if (z.is_empty()) z = diff(beta,1,1);
   if (u.is_empty()) u = zeros<mat>(p,n-1);
   mat beta_best = beta, z_best = z, u_best = u;
   vec objective(maxit+1);
   objective(0) = objective_global(beta,x,y,lam1,lam2);
   if (verbose) 
      Rprintf("Iteration 0 Objective = %f\n",objective(0));
   double objective_best = objective(0);
   vec grad(p), Delta(p), Delta_prev(p);
   vec zu(p), zu_prev(p);
   double alpha;   
   vec zj(p); 
   Progress prog(0, false); 

   for (i=1; i<=maxit; i++) {

      // Update primal variable beta 
      Delta.zeros();
      zu.zeros();
      for (j=0; j<n; j++) {
         Delta_prev = Delta;
         zu_prev = zu;
         if (j < (n-1)) {
            Delta = beta.col(j+1) - beta.col(j);
            zu = z.col(j) + u.col(j);
         } else { Delta.zeros(); zu.zeros(); }
         grad = trans(x.slice(j)) * 
            (x.slice(j) * beta.col(j)  - y.col(j)) + 
            rho * (zu - Delta - zu_prev + Delta_prev);
         beta.col(j) -= grad/L;
      } 
      beta = soft_cpp(beta,lam1/L);        
    
      // Update primal variable z and dual variable u
      z.zeros();
      for (j=0; j<(n-1); j++) {
         Delta = beta.col(j+1) - beta.col(j);
         zj = Delta - u.col(j);
         nrm = norm(zj);
         if (nrm > 0.0) {
            alpha = 1 - lam2 / (rho * nrm);
            if (alpha > 0.0) z.col(j) = alpha * zj; 
         }
         u.col(j) += z.col(j) - Delta;
      }    

      // Update objective
      objective(i) = objective_global(beta,x,y,lam1,lam2); 
      if (objective(i) < objective_best)
         { objective_best = objective(i);
         beta_best = beta; 
         u_best = u; z_best = z;} 
//    prog.increment();
      if (verbose && i % 100 == 0) 
            Rprintf("Iteration %d Objective = %f\n",i,objective(i));
      // Check convergence (objective) and user interruption
      if (i % 100 == 0 && Progress::check_abort()) 
         break;
      if (abs(objective(i)-objective(i-1)) <= tol * objective(i-1)) 
         break;
   }

   if (i > maxit) i = maxit;

   return List::create(Named("beta") = beta_best,
                        Named("z") = z_best,
                        Named("u") = u_best,
                        Named("objective") = objective_best,
                        Named("trace") = objective.head(i+1),
                        Named("rho") = rho);

}
 




///////////////////////////////////////////////////////////////////////////////////




// PRIMAL-DUAL (CONDAT 2013)
// [[Rcpp::export]]
Rcpp::List pd_condat(const cube & x, const mat & y, const double lam1,
   const double lam2, mat beta, mat s, double Lx, double tau, double sig, 
   const double rho, const double tol, const int maxit, const bool verbose)
{
   uword n = y.n_cols, p = x.n_cols;
   uword i,j;
   double lam1tau = lam1*tau;
   double nrm;
   vec objective(maxit+1);
   if (beta.is_empty()) beta = zeros<mat>(p,n);
   if (s.is_empty()) s = zeros<mat>(p,n-1);
   mat beta_old, s_old;
   mat grad(p,n);
   sp_mat D(n,n-1);
   D.diag(0).fill(-1.0);
   D.diag(-1).ones();
   sp_mat tD = trans(D);
   if (sig <= 0.0) {
      if (Lx == 0.0) {
         for (j=0; j<n; j++) {
            nrm = max(svd(x.slice(j)));
            if (nrm > Lx) Lx = nrm; }
         Lx = Lx*Lx; 
         }
      sig = 0.25 * (1/tau - Lx); 
   }

   objective(0) = objective_global(beta,x,y,lam1,lam2);
   double objective_best = objective(0);
   mat beta_best = beta, s_best = s;
   if (verbose) 
      Rprintf("Initial objective %f \n",0,objective_best);
   Progress prog(0,false);

   for (i=1;i<=maxit;i++)
   {
      // Store current beta and s
      beta_old = beta;
      s_old = s;

      // Update s by projecting on (cartesian product of) l2 balls
      s += sig * beta * D; 
      for (j=0; j<(n-1); j++) {
         nrm = arma::norm(s.col(j));
         if (nrm > lam2) s.col(j) *= (lam2/nrm);
      }
      
      // Update b by soft-thresholding
      for (j=0; j<n; j++) 
         grad.col(j) = trans(x.slice(j)) * (x.slice(j) * beta.col(j) - y.col(j));
      beta -=  tau * (grad + (2.0 * s - s_old) * tD);
      beta = soft_cpp(beta,lam1tau);
            
      // Relaxation
      beta = rho * beta + (1-rho) * beta_old;
      s = rho * s + (1-rho) * s_old;

      // Update objective
      objective(i) = objective_global(beta,x,y,lam1,lam2);
      if (verbose && i % 100 == 0) 
         Rprintf("Iteration %d Objective %f\n",i,objective(i)); 
      if (objective(i) < objective_best) {
         objective_best = objective(i);
         beta_best = beta;
         s_best = s; }
      if (abs(objective(i)-objective(i-1)) <= tol * objective(i-1)) 
         break;

      if (i % 100 == 0 && Progress::check_abort()) break;
    
   }
   if (i > maxit) i = maxit;
   return List::create(Named("beta") = beta_best,
                          Named("s") = s_best,
                          Named("objective") = objective_best, 
                          Named("trace") = objective.head(i+1),
                          Named("tau") = tau,
                          Named("sig") = sig,
                          Named("rho") = rho);

}




///////////////////////////////////////////////////////////////////////////////////




Rcpp::List pd3o(const cube & x, const mat & y,
   const double lam1, const double lam2, mat beta, double tau, 
   double sig, const int maxit, const bool verbose)
{
   uword n = y.n_cols, p = x.n_cols;
   uword i,j;
   double lam1tau = lam1*tau;
   double nrm;
   vec objective(maxit+1); 
   if (beta.is_empty()) beta = zeros<mat>(p,n);
   mat grad(p,n), s(p,n-1);
   s.zeros();
   sp_mat D(n,n-1);
   D.diag(0).fill(-1.0);
   D.diag(-1).ones();
   sp_mat tD = trans(D);
   double nrmx, b = 0.0;
   for (j=0; j<n; j++) {
      nrmx = max(svd(x.slice(j)));
      if (nrmx > b) b = nrmx; }
   b = b * b;
   if (tau >= 2.0 / b) tau = 1.99 / b;
   if (sig <= 0.0) sig = 0.25 / tau;
   
   objective(0) = objective_global(beta,x,y,lam1,lam2);
   double objective_best = objective(0);
   mat beta_best = beta, beta_bar = beta, beta_h;
   if (verbose) Rprintf("Initial objective %f \n",objective_best);

   for (i=1;i<=maxit;i++)
   {
 
      // Update s by projecting on (cartesian product of) l2 balls
      s += sig * beta_bar * D; 
      for (j=0; j<(n-1); j++) {
         nrm = arma::norm(s.col(j));
         if (nrm > lam2) s.col(j) *= (lam2/nrm);
      }
      
      // Update beta by soft-thresholding
      for (j=0; j<n; j++) 
         grad.col(j) = trans(x.slice(j)) * (x.slice(j) * beta.col(j) - y.col(j));
      beta_h = beta - tau * grad;
      beta = soft_cpp(beta_h - tau * s * tD,lam1tau);
            
      // Update beta_bar
      for (j=0; j<n; j++) 
         grad.col(j) = trans(x.slice(j)) * (x.slice(j) * beta.col(j) - y.col(j));
      beta_bar = 2 * beta - beta_h - tau * grad;

      // Update objective
      objective(i) = objective_global(beta,x,y,lam1,lam2);
      if (objective(i) < objective_best) {
         objective_best = objective(i);
         beta_best = beta; }

      if (verbose && i % 100 == 0) 
         Rprintf("Iteration %d Objective %f\n",i,objective(i)); 

   }

   return List::create(Named("beta") = beta_best, 
                          Named("objective") = objective_best, 
                          Named("trace") = objective,
                          Named("tau") = tau,
                          Named("sig") = sig);

}





///////////////////////////////////////////////////////////////////////////////////



// SMOOTH PROXIMAL GRADIENT (GRADIENT CALCULATION)
// [[Rcpp::export]]
mat grad_spg(const cube & x, const mat & y, const mat & beta, 
   double mu, double lam2)
{

   uword i, n = beta.n_cols, p = beta.n_rows;
   double nrm;
   vec grad_g(p), grad_fmu(p);
   vec alpha(p), alpha_prev(p);
   mat grad_h(p,n);

   // Gradient calculation
   alpha.zeros();
   for (i=0; i<n; i++) {
      alpha_prev = alpha;
      grad_g = trans(x.slice(i)) * (x.slice(i) * beta.col(i) - y.col(i));
      if (i < n-1) {
         alpha = (beta.col(i+1) - beta.col(i)) * (lam2 / mu);
         nrm = norm(alpha);
         if (nrm > 1.0) alpha /= nrm;
      } else { alpha.zeros(); }
      grad_fmu = lam2 * (alpha_prev - alpha);
      grad_h.col(i) = grad_g + grad_fmu;
   }

   return grad_h;
}



// SMOOTH PROXIMAL GRADIENT (MAIN FUNCTION)
//[[Rcpp::export]]
Rcpp::List spg_cpp(const cube & x, const mat & y, const double lam1, 
   const double lam2, mat beta, const double eps, double Lx, 
   const int maxit, const bool verbose)
{
   // Initialization
   uword i, n = y.n_cols, p = x.n_cols;
   double D = (n - 1.0) / 2.0; 
   double mu = eps / (2.0 * D); 
   double nrmC = 2 * lam2; 
   if (Lx == 0.0) {
      double nrmx;
      for (i=0; i<n; i++) {
         nrmx = max(svd(x.slice(i)));
         if (nrmx > Lx) Lx = nrmx; }
      Lx = Lx*Lx;
   }
   double Lmu = nrmC*nrmC/mu;
   double L = Lx + Lmu; 
   if (beta.is_empty()) beta = zeros<mat>(p,n);
   mat w = beta, beta_old(size(beta)), grad(size(beta));
   double theta_old, theta = 1;
   vec objective(maxit+1);
   objective(0) = objective_global(beta,x,y,lam1,lam2);
   double objective_best = objective(0);
   mat beta_best = beta;
   Progress prog(0,false);
   int count = 0;
   if (verbose) 
      Rprintf("Initial objective %f\n", objective_best);

   // Loop
   for (i=1; i<=maxit; i++) {  
      beta_old = beta;
      theta_old = theta;
      grad = grad_spg(x,y,w,mu,lam2);
      beta = soft_cpp(w - grad/L, lam1/L);            
      theta = 2.0/(i+2.0);
      w = beta + ((1.0-theta_old)/theta_old*theta) * (beta-beta_old);
      objective(i) = objective_global(beta,x,y,lam1,lam2); 
      if (objective(i) < objective_best) {
         count = 0;
         objective_best = objective(i);
         beta_best = beta; } else {count++;}
      if (verbose && i % 100 == 0)
         Rprintf("Iteration %d Objective %f\n",i,objective(i));
//      if (objective(i) > objective(i-1)) break;
//      if (abs(objective(i)-objective(i-1)) <= tol * objective(i-1))
//         break;
      if (count == 100) break;
      if (i % 100 == 0 && Progress::check_abort()) 
         break; 
   }
   if (i > maxit) i = maxit;
   return List::create(Named("beta") = beta_best, 
                          Named("objective") = objective_best, 
                          Named("trace") = objective.head(i+1));
}




