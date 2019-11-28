// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace std;
using namespace Rcpp;
using namespace arma;



//***************************************************************
// Problem: minimize f(u) = max_{1<=t<T} f_t(u)
// where 
// f_t(u) = 1/w_t^2 (|| z_1 + ... + z_t - u_1 - ... - u_t||_2^2 + c_t^2)
// with u = (u_1,...,u_T) 
// subject to 
// 1) box constraints: || u_t ||_∞ <= lambda_1, 1 <= t <= T
// 2) linear constraint:  u_1 + ... + u_T = z_1 + ... + z_T
// Solution: projected subgradient
// (In fact, objective function is sqrt(f(u)) )
//***************************************************************





// Projection step with method of alternating projections
//[[Rcpp::export]]
List lambda2_max(mat u, const mat & z, const double lam1, 
   const vec & c2, const vec & w2, const double s,
   const bool decrease, const uword maxit, 
   const double tol, const uword maxit_proj, 
   const double tol_proj)
{

   double sw = s; // effective step size
   double nlam1 = -lam1;
   uword i,j,t; // counters
   uword a; // active function in objective
   const uword n = z.n_cols;   
   vec g(z.n_rows); // negative subgradient (unscaled)
   vec xbar(z.n_rows);
   const vec zbar = arma::mean(z,1);
   uvec idx; 
   
   // Initialize u
   if (u.is_empty()) 
      { u = repmat(zbar,1,n); }
   else { 
      for (j=0; j<maxit_proj; j++) 
      {
         u = clamp(u,nlam1,lam1);
         xbar = zbar - arma::mean(u,1);
         idx = arma::find(arma::abs(xbar) > tol_proj);
         if (idx.is_empty())
            break;
         u.rows(idx) += repmat(xbar(idx),1,n); 
      }   
   } 

   // Objective function values
   double objective_prev, objective, objective_best;
   double delta;
   vec cumx = zeros<vec>(z.n_rows);
   double ftu, fu = -1.0;
   for (t=0; t<(n-1); t++) {
      cumx += z.col(t) - u.col(t);
      ftu = (accu(pow(cumx,2))+c2(t))/w2(t);
      if (ftu > fu) {
         a = t;   
         g = cumx;
         fu = ftu; }
   }
   objective = sqrt(fu);
   objective_best = objective;
   mat u_best = u; 
   
   // MAIN LOOP
   for (i=0; i<maxit; i++)
   {
      objective_prev = objective;
      
      // Subgradient step 
      if (decrease) 
         sw = s / sqrt(i+1); 
      u.cols(0,a).each_col() += (2.0/w2(a)*sw) * g;

      // Projection step (alternating projections)
      for (j=0; j<maxit_proj; j++)
      {
         // Clamping
         u = clamp(u,nlam1,lam1);
         // Centering
         xbar = zbar - arma::mean(u,1);
         idx = arma::find(arma::abs(xbar) > tol_proj);
         if (idx.is_empty())
            break;
         u.rows(idx) += repmat(xbar(idx),1,n); 
      } 

      // Calculate objective value
      cumx.zeros();
      fu = -1.0;
      for (t=0; t<(n-1); t++) {
         cumx += z.col(t) - u.col(t);
         ftu = (accu(pow(cumx,2))+c2(t))/w2(t);
         if (ftu > fu) {
            a = t;   
            g = cumx;
            fu = ftu; }
      }
      objective = sqrt(fu);
      if (objective < objective_best && j < maxit_proj) {
         objective_best = objective;
         u_best = u; }
               
      // Check stopping criterion
      if (i>0) {
         delta = abs(objective-objective_prev);
         if (delta <= tol || 
            delta <= tol * objective_prev)     
            break; }
   }

   return List::create(Named("objective") = objective_best,
         Named("u") = u_best);
}


