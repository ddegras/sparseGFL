// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace std;
using namespace arma;
using namespace Rcpp;


// FUNCTION PROTOTYPES

double objective_global (const mat &, const cube &, const mat &, 
   const double, const double, const double, const vec &, const bool); 

double objective_local (const vec &, const cube &, const mat &, 
   const double, const double, const double, const double, const vec &, 
   const double, const vec &, const bool);

double objective_blocks (const mat &, const cube &, const mat &, 
   const double, const double, const double, const vec &, const double,
   const vec &, const double, const vec &, const bool); 

double objective_chains (const mat &, const cube &, const mat &, 
   const double, const double, const double, const vec &,
   const bool, const uvec &, const uvec &); 

vec trivial_local (const vec &, const cube &, const mat &, 
   const double, const double, const double, const double, const vec &, 
   const double, const vec &, const bool, const double);

vec trivial_mm (const double, const double, const double, 
   const vec &, const double, const vec &, const double, 
   const vec &, const bool, const uword, const double);

vec fp1_vec(vec, const double, const vec &, const double, const double, 
   const double, const vec &, const bool, const bool, const int,
   const int, const double);

vec fp2_vec(vec, const double, const vec &, const double, const double, 
   const double, const vec &, const double, const vec &, const bool, 
   const bool, const int, const int, const double);

List FISTA_1(vec, const cube &, const mat &, const double, 
   const double, const double, const double, const vec &, const double, 
  const vec &,  const double, const bool, const uword, const double, 
   const uword, const double, const double);

List FISTA_M(mat, const cube &, const mat &, const double, 
   const double, const double, const vec &, const bool, uvec, uvec, 
   double, const uword, const double, const double);
 
mat grad_cpp(const vec &, const cube &, const mat &, const double, 
   const double, const double, const double, const vec &, const double, 
   const vec &, const bool);

mat pg_cpp(const mat &, const uvec &, double, double, vec, uword, double);

LogicalVector compare_cols(const mat &, const double);








// OBJECTIVE FUNCTION (GLOBAL)
// [[Rcpp::export]]
double objective_global (const mat & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const vec & w, 
   const bool intercept) 
{
   
   uword i, N = y.n_rows, p = beta.n_rows, n = y.n_cols;
   double sse = 0.0, lasso = 0.0, tv = 0.0, el = 0.0;

   // Sum of squared errors
   for (i=0; i<n; i++) 
      sse += pow(arma::norm(y.col(i)-x.slice(i)*beta.col(i)),2); 
   sse *= 0.5;

   // l1 penalty
   if (lam1 > 0.0) {
      if (intercept) 
         { lasso = lam1 * arma::accu(arma::abs(beta.rows(0,p-N-1))); }
      else 
         { lasso = lam1 * arma::accu(arma::abs(beta)); }
   }

   // TV penalty
   if (lam2 > 0.0) {
      for(i=0; i<(n-1); i++)
         tv += w(i) * arma::norm(beta.col(i+1)-beta.col(i));
      tv *= lam2; 
   }

   // Elastic net penalty 
   if (lam3 > 0.0) {
      if (intercept) 
         { el = 0.5 * lam3 * arma::accu(arma::square(beta.rows(0,p-N-1))); }
      else 
         { el = 0.5 * lam3 * arma::accu(arma::square(beta)); }
   }

   return (sse+lasso+tv+el);
}

   


//****************************************************************************




// OBJECTIVE FUNCTION (SINGLE BLOCK OR SINGLE FUSION CHAIN)
// [[Rcpp::export]]
double objective_local (const vec & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const double wprev, 
   const vec & bprev, const double wnext, const vec & bnext, const bool intercept) 
{
   
   uword i, N = y.n_rows, p = beta.n_elem, n = y.n_cols;

   // Sum of squared errors
   double sse = 0.0;
   for (i=0; i<n; i++)
      sse += arma::accu(arma::square(y.col(i)-x.slice(i)*beta));
   sse *= 0.5; 

   // l1 penalty
   double lasso = 0.0;
   if (lam1 > 0.0) {
      if (intercept) 
         { lasso = lam1 * n * arma::accu(arma::abs(beta.head(p-N))); }
      else 
         { lasso = lam1 * n * arma::accu(arma::abs(beta)); }
   }

   // TV penalty
   double tv = 0.0; 
   if (lam2 > 0.0) {
      if (!bprev.is_empty())
         tv += lam2 * wprev * arma::norm(beta-bprev);
      if (!bnext.is_empty())
         tv += lam2 * wnext * arma::norm(beta-bnext);
   }

   // Elastic net penalty 
   double el = 0.0;
   if (lam3 > 0.0) {
      if (intercept) 
         { el = 0.5 * lam3 * n * arma::accu(arma::square(beta.head(p-N))); }
      else 
         { el = 0.5 * lam3 * n * arma::accu(arma::square(beta)); }
   }

   return (sse+lasso+tv+el);
}




//****************************************************************************




// OBJECTIVE FUNCTION (MULTIPLE CONSECUTIVE BLOCKS, NOT NECESSARILY FUSED)
// [[Rcpp::export]]
double objective_blocks (const mat & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const vec & w, 
   const double wprev, const vec & bprev, const double wnext, const vec & bnext, 
   const bool intercept) 
{
   
   uword i, N = y.n_rows, p = beta.n_rows, n = y.n_cols;

   // Sum of squared errors
   double sse = 0.0;
   for (i=0; i<n; i++) 
      sse += arma::accu(arma::square(y.col(i)-x.slice(i)*beta.col(i))); 
   sse *= 0.5;

   // l1 penalty
   double lasso = 0.0;
   if (lam1 > 0.0) {
      if (intercept) 
         { lasso = lam1 * arma::accu(arma::abs(beta.rows(0,p-N-1))); }
      else 
         { lasso = lam1 * arma::accu(arma::abs(beta)); }
   }

   // TV penalty
   double tv = 0.0; 
   if (lam2 > 0.0) {
      for(i=0; i<(n-1); i++)
         tv += w(i) * arma::norm(beta.col(i+1)-beta.col(i));
      if (!bprev.is_empty())
         tv += wprev * arma::norm(beta.col(0)-bprev);
      if (!bnext.is_empty())
         tv += wnext * arma::norm(beta.col(n-1)-bnext);
      tv *= lam2; }

   // Elastic net penalty
   double el = 0.0;
   if (lam3 > 0.0) {
      if (intercept) 
         { el = 0.5 * lam3 * arma::accu(arma::square(beta.rows(0,p-N-1))); }
      else  
         { el = 0.5 * lam3 * arma::accu(arma::square(beta)); }
   }

   return (sse+lasso+tv+el);
}

    


//****************************************************************************




// OBJECTIVE FUNCTION (ALL FUSION CHAINS)
double objective_chains (const mat & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const vec & w,
   const bool intercept, const uvec & startc, const uvec & endc) 
{
   uword i,g,G = startc.n_elem;
   uvec len = endc - startc + ones<uvec>(G);

   // Sum of squared errors
   double sse = 0.0;
   for (g=0; g<G; g++) 
      for (i=startc(g); i<=endc(g); i++)
         sse += arma::accu(arma::square(y.col(i)-x.slice(i)*beta.col(g)));  
   sse *= 0.5;

   // l1 penalty
   double lasso = 0.0; 
   uword N = y.n_rows, p = beta.n_rows;
   if (lam1 > 0.0) {
      if (intercept) {
         for (g=0; g<G; g++) 
            lasso += len(g) * 
               arma::accu(arma::abs(beta(span(0,p-N-1),g))); }
      else {
         for (g=0; g<G; g++) 
            lasso += len(g) * arma::accu(arma::abs(beta.col(g))); }
      lasso *= lam1; 
   } 

   // TV penalty
   double tv = 0.0; 
   if (lam2 > 0.0) {
      for(g=0; g<(G-1); g++)
         tv += w(endc(g)) * arma::norm(beta.col(g+1)-beta.col(g));
      tv *= lam2; }

   // Elastic net penalty
   double el = 0.0;
   if (lam3 > 0.0) {
      if (intercept) 
         { for (g=0; g<G; g++) 
            el += len(g) * 
               arma::accu(arma::square(beta(span(0,p-N-1),g))); }
      else 
         { for (g=0; g<G; g++)
            el += len(g) * arma::accu(arma::square(beta.col(g))); }
      el *= 0.5 * lam3; 
   }
   
   return (sse+lasso+tv+el);
}
   





//****************************************************************************




// CHECK FOR TRIVIAL SOLUTIONS IN LOCAL OBJECTIVE 
// [[Rcpp::export]]
vec trivial_local (const vec & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, 
   const double wprev, const vec & bprev, const double wnext,
   const vec & bnext, const bool intercept, const double tol)
{

   vec null_vec;
   null_vec.reset();
   bool is_null_bprev = (bprev.is_empty() || wprev == 0.0);
   bool is_null_bnext = (bnext.is_empty() || wnext == 0.0);
   double d; // distance b/w two candidate solutions
   mat Delta;
   uword i, N = y.n_rows, p = x.n_cols, n = y.n_cols;
   double lam1n = lam1 * n, lam3n = lam3 * n;
   vec grad(p); // subgradient 
   uvec idx; // index of zero values in candidate solution

   // If current coefficient block is provided and different from previous
   // and next blocks, see if it is simple solution to local problem
   if (!(beta.is_empty() || approx_equal(beta,bprev,"both",tol,tol) ||  
       approx_equal(beta,bnext,"both",tol,tol)))
      { grad.zeros();
      for (i=0; i<n; i++) 
         grad += arma::trans(x.slice(i)) * (x.slice(i)*beta-y.col(i));
      if (lam1 > 0.0) {
         if (intercept) {
            grad.head(p-N) += lam1n * arma::sign(beta.head(p-N));
         } else {
            grad += lam1n * arma::sign(beta); }
      }
      if (!is_null_bprev) {
         Delta = beta - bprev;
         grad += (lam2*wprev/norm(Delta)) * Delta; }
      if (!is_null_bnext) {
         Delta = beta - bnext;
         grad += (lam2*wnext/norm(Delta)) * Delta; }
      if (lam3 > 0.0) {
         if (intercept) {
            grad.head(p-N) += lam3n * beta.head(p-N);
         } else {
            grad += lam3n * beta; } 
      } 

      idx = find(beta);
      if (arma::norm(grad(idx),"inf") <= tol) {
         if (intercept) {
            idx = arma::find(beta.head(p-N) == 0.0);
         } else { idx = arma::find(beta == 0.0); }
         if (norm(grad(idx),"inf") <= lam1n)
            return beta; 
      }
   } 


   // Return empty matrix if no candidate for trivial solution
   if (is_null_bprev && is_null_bnext) return null_vec;
      
   bool single_candidate = (is_null_bprev || is_null_bnext || 
      approx_equal(bprev,bnext,"both",tol,tol));
   if (!single_candidate) {
      Delta = bnext - bprev;
      d = norm(Delta); 
   }

   // Check whether previous coefficient block minimizes local
   // objective (if it exists)
   if (!is_null_bprev) {
      grad.zeros();
      for (i=0; i<n; i++) 
         grad += trans(x.slice(i)) * (x.slice(i) * bprev - y.col(i));
      if (lam1 > 0.0) {
         if (intercept) {
            idx = arma::find(bprev.head(p-N) == 0.0);      
            grad.head(p-N) += lam1n * sign(bprev.head(p-N)); 
         } else {
            idx = arma::find(bprev == 0.0);            
            grad += lam1n * sign(bprev); }
      }
      if (!single_candidate) 
         grad -= (lam2*wnext/d) * Delta; 
      if (lam3 > 0.0) {
         if (intercept) {
            grad.head(p-N) += lam3n * bprev.head(p-N); }
         else { grad += lam3n * bprev; }
      }
      if (idx.n_elem > 0)
         grad(idx) -= clamp(grad(idx),-lam1n,lam1n);

      if (single_candidate) {
         if (arma::norm(grad) <= (lam2*(wprev+wnext))) 
            {return bprev;} else {return null_vec;}   
      } else if (arma::norm(grad) <= (lam2*wprev)) 
         {return bprev;}     
   }

   // Check whether next coefficient block minimizes local 
   // objective (if it exists)
   if (!is_null_bnext) {
      grad.zeros();
      for (i=0; i<n; i++) 
         grad += trans(x.slice(i)) * (x.slice(i) * bnext - y.col(i));
      if (lam1 > 0.0) {
         if (intercept) {
            grad.head(p-N) += lam1n * sign(bnext.head(p-N)); 
            idx = arma::find(bnext.head(p-N) == 0.0);      
         } else {
            grad += lam1n * sign(bnext); 
            idx = arma::find(bnext == 0.0); }           
      }
      if (!is_null_bprev) 
         grad += (lam2*wprev/d) * Delta; 
      if (lam3 > 0.0) {
         if (intercept) {
            grad.head(p-N) += lam3n * bnext.head(p-N); }
         else { grad += lam3n * bnext; }
      }
      if (idx.n_elem > 0)
         grad(idx) -= clamp(grad(idx),-lam1n,lam1n);

      if (arma::norm(grad) <= wnext*lam2) return bnext;
   } 

   return null_vec;
}




//****************************************************************************





// CHECK FOR TRIVIAL SOLUTIONS IN MAJORIZATION-MINIMIZATION 
// [[Rcpp::export]]
vec trivial_mm (const double lam1, const double lam2, const double L, 
   const vec & z, const double wprev, const vec & bprev, const double wnext, 
   const vec & bnext, const bool intercept, const uword N, const double tol)
{

   vec null_vec;
   null_vec.reset();
   bool is_null_bprev = (bprev.is_empty() || wprev == 0.0);
   bool is_null_bnext = (bnext.is_empty() || wnext == 0.0);
   bool single_candidate = (is_null_bprev || is_null_bnext || 
      arma::approx_equal(bprev,bnext,"both",tol,tol));
   double d; // distance b/w two candidate solutions
   vec Delta;
   if (!single_candidate) 
      {Delta = bnext - bprev;
      d = arma::norm(Delta);} 
   uword p = z.n_elem;
   vec grad(p); // subgradient 
   uvec idx; // index of null values in candidate solution
   
  
   // Check whether previous coefficient block minimizes quadratic 
   // approximation to local objective (if it exists)
   if (!is_null_bprev) {
      grad = L * (bprev - z);
      if (!single_candidate) 
         grad -= (lam2*wnext/d) * Delta; 
      if (lam1 > 0.0) {
         if (intercept) {
            idx = find(bprev.head(p-N) == 0.0);      
            grad.head(p-N) += lam1 * sign(bprev.head(p-N));
         } else { 
            idx = find(bprev == 0.0);
            grad += lam1 * sign(bprev); }
         if (idx.n_elem > 0)
            grad(idx) -= clamp(grad(idx),-lam1,lam1);
      }
      if (single_candidate) {
         if (arma::norm(grad) <= (lam2*(wprev+wnext))) 
            {return bprev;} else {return null_vec;}   
      } else if (arma::norm(grad) <= lam2*wprev) {return bprev;}
   }

   // Check whether next coefficient block minimizes quadratic
   // approximation to local objective (if it exists)
   if (!is_null_bnext) {
      grad = L * (bnext - z);
      if (!is_null_bprev) 
         grad += (lam2 * wprev / d) * Delta; 
      if (lam1 > 0.0) {
         if (intercept)
            {idx = arma::find(bnext.head(p-N) == 0.0);
            grad.head(p-N) += 
               lam1 * arma::sign(bnext.head(p-N));}
         else
            {idx = find(bnext == 0.0);
            grad += lam1 * sign(bnext);}
         if (idx.n_elem > 0)
            grad(idx) -= arma::clamp(grad(idx),-lam1,lam1); 
      }
      if (arma::norm(grad) <= lam2*wnext) return bnext;
   } 

   return null_vec;
}





//****************************************************************************





// FIXED POINT ITERATION (1 NEIGHBOR)
// [[Rcpp::export]]
vec fp1_vec(vec beta, const double L, const vec & z, const double lam1,
   const double lam2, const double w, const vec & bprev,  
   const bool intercept, const int N, const int maxit, const double tol)
{

   uword i, p = beta.n_elem;
   double tol_equal = tol * sqrt(p);
   double nlam1 = -lam1;
   double dbu;
   const mat Lz = L * z;     
   vec beta_old(p); 

   for (i=1; i<=maxit; i++) {

         beta_old = beta; 
      
      // Distance between current solution and neighbor
      dbu = arma::norm(beta-bprev);
      if (dbu <= tol_equal) return bprev;

      // Gradient step 
      beta = Lz + (lam2*w/dbu) * bprev;

      // Soft-thresholding 
      if (lam1 > 0.0) {
         if (intercept) {
            beta.head(p-N) -= 
               arma::clamp(beta.head(p-N),nlam1,lam1);  
         } else {
            beta -= arma::clamp(beta,nlam1,lam1); }
      }

      // Scaling
      beta /=  L + lam2*w/dbu;

      // Check convergence
      if (arma::approx_equal(beta,beta_old,"both",tol,tol)) 
         break;

   } 

   return beta;

}




// FIXED POINT ITERATION (2 NEIGHBORS)
// [[Rcpp::export]]
vec fp2_vec(vec beta, const double L, const vec & z, const double lam1, 
   const double lam2, const double wprev, const vec & bprev, 
   const double wnext, const vec & bnext, const bool intercept,
   const int N, const int maxit, const double tol)
{

   uword i, p = beta.n_elem;
   double tol_equal = tol * sqrt(p);
   double nlam1 = -lam1;
   vec beta_old(p); 
   double dbu, dbv;
   const vec Lz = L * z;

   for (i = 1; i <= maxit; i++) {

      beta_old = beta;

      // Distance from neighboring blocks
      dbu = arma::norm(beta-bprev);
      if (dbu <= tol_equal) return bprev;
      dbv = arma::norm(beta-bnext);
      if (dbv <= tol_equal) return bnext;    

      // Gradient step 
      beta = Lz + (lam2*wprev/dbu) * bprev 
         + (lam2*wnext/dbv) * bnext;

      // Soft-thresholding 
      if (lam1 > 0.0) {
         if (intercept) 
            { beta.head(p-N) -= 
            	   arma::clamp(beta.head(p-N),nlam1,lam1);}
         else 
            { beta -= arma::clamp(beta,nlam1,lam1); }
      }

      // Scaling
      beta /=  L + (lam2*wprev/dbu) + (lam2*wnext/dbv);

      // Check convergence
      if (arma::approx_equal(beta,beta_old,"both",tol,tol)) break;
   }     

   return beta;

}





//****************************************************************************



//#################
// FISTA FUNCTIONS
//#################


// [[Rcpp::export]]
Rcpp::List FISTA_1(vec beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const double wprev,
   const vec & bprev, const double wnext, const vec & bnext,  const double L,  
   const bool intercept, const uword maxit_fista, const double tol_fista, 
   const uword maxit_fp, const double tol_fp, const double tol_equal)
{
      
   bool first = bprev.is_empty(), last = bnext.is_empty();
   uword nbounds;
   if ((first && last) || lam2 == 0.0) { nbounds = 0; }
   else if (first || last) { nbounds = 1; }
   else { nbounds = 2; }
   uword i, j, N = y.n_rows, p = beta.n_elem, n = y.n_cols;
   double objective = objective_local(beta, x, y, lam1, lam2, lam3,
      wprev, bprev, wnext, bnext, intercept);
   vec beta_best = beta;
   double objective_best = objective, objective_init = objective, 
      objective_old;   
   vec beta_adj = beta; // adjusted solution
   double theta = 1.0, theta_old; // acceleration parameter
   vec beta_old(p), z(p);   
   bool progress;
   if (first && last && L == 0.0) 
     {beta.zeros(); 
      objective_best = 0.5 * accu(pow(y,2));
      progress = (objective_best < objective_init);
      return List::create(Named("beta") = beta, 
                          Named("objective") = objective, 
                          Named("progress") = progress); }

   mat xtx;
   vec xty;   
   bool store_cp = (p <= n);
   if (store_cp) {
      xtx.set_size(p,p);
      xtx.zeros();
      xty.set_size(p);
      xty.zeros();
      for (i=0; i<n; i++) {
         xtx += arma::trans(x.slice(i)) * x.slice(i);
         xty += arma::trans(x.slice(i)) * y.col(i); }
      xtx /= L; xty /= L;
   }
   const double lam1n = lam1 * n, lam3n = lam3 * n;
   const bool L_positive = (L > 0.0), lam3_positive = (lam3 > 0.0);
   vec grad(p);

   // FISTA
   for (i=1; i<=maxit_fista; i++) {
      
      // Store previous results
      beta_old = beta;
      objective_old = objective;
      theta_old = theta;
 
      if (i % 10 == 0)
         Rcpp::checkUserInterrupt();

      // Gradient step
      if (L_positive) {    
         if (store_cp) {
            grad = xtx*beta_adj-xty;
         } else { 
         	grad.zeros();
         	for (j=0; j<n; j++)
         	   grad += arma::trans(x.slice(j)) * 
         	   (x.slice(j)*beta_adj-y.col(j));
         	grad /= L; }
         z = beta_adj - grad;     
      } else {
         z = beta_adj;
      }
      if (lam3_positive) {
         if (intercept) 
            {z.head(p-N) -= 
               (lam3n/L) * beta_adj.head(p-N);} 
         else 
            {z -= (lam3n/L) * beta_adj;}
      }

      // Check for trivial solutions 
      beta = trivial_mm(lam1n, lam2, L, z, wprev, bprev, 
         wnext, bnext, intercept, N, tol_equal);
         
      // If no trivial solution, apply fixed point iteration 
      // (or soft-thresholding if no neighbors)
      if (beta.is_empty()) {
         if (approx_equal(beta_adj, bprev, "both", tol_equal, tol_equal) 
            || approx_equal(beta_adj, bnext, "both", tol_equal, tol_equal)) {
               beta = z; } else { beta = beta_adj; }

         switch (nbounds) {
            case 0:
               beta = z - clamp(z,-lam1n/L,lam1n/L);
               break;
            case 1:
               if (first) {
               beta = fp1_vec(beta, L, z, lam1n, lam2, wnext, bnext, 
                  intercept, N, maxit_fp, tol_fp); } else {
               beta = fp1_vec(beta, L, z, lam1n, lam2, wprev, bprev, 
                  intercept, N, maxit_fp, tol_fp); }
               break;
            case 2:
               beta = fp2_vec(beta, L, z, lam1n, lam2, wprev, bprev, wnext,
               bnext, intercept, N, maxit_fp, tol_fp);
               break;
            }
      }
                  
      // Nesterov method
      theta = (1.0 + sqrt(1.0+4.0*theta*theta)) / 2.0;
      beta_adj = beta + ((theta_old-1.0) / theta) * (beta - beta_old);
      
      // Update objective (every iteration? every 10 iterations?)
      objective = objective_local(beta, x, y, lam1, lam2, lam3,
         wprev, bprev, wnext, bnext, intercept);         
      if (objective < objective_best) {
         beta_best = beta;
         objective_best = objective;}
      
      // Check convergence
      if (abs(objective - objective_old) <= tol_fista * objective_old ||
         objective <= tol_fista) break;
   }
   
   progress = (objective_best < objective_init); 

   return List::create(Named("beta") = beta_best,
                       Named("objective") = objective_best,
                       Named("progress") = progress);

}
   

// FISTA OPTIMIZATION OVER FIXED CHAINS
// [[Rcpp::export]]
List FISTA_M(mat beta, const cube & x, const mat & y, const double lam1, 
   const double lam2, const double lam3, const vec & w, 
   const bool intercept, uvec startc, uvec endc, double eta, 
   const uword maxit, const double tol_fista, const double tol_equal)
{
   uword g,i,j; 
   uword nchain = startc.n_elem; // number of chains
   uvec len = endc - startc + ones<uvec>(nchain); // chains lengths
   uword N = y.n_rows, p = x.n_cols; 

   if (beta.n_cols > nchain) // reduce beta if needed
      beta = beta.cols(startc);
        
   // Initial objective value
   double objective = objective_chains(beta, x, y, lam1, lam2, lam3,
      w, intercept, startc, endc);   

   // Lipschitz constants for gradient step (estimates)
   vec Lvec(y.n_cols), d;  
   double val;
   mat sumxtx(p,p);
   for (g=0; g<nchain; g++) {
      sumxtx.zeros();
      for (j=startc(g); j<=endc(g); j++) 
         sumxtx += arma::trans(x.slice(j)) * x.slice(j);   
      d = arma::eig_sym(sumxtx);
      val = d(p-1) + (lam3 * len(g));
      Lvec.subvec(startc(g),endc(g)).fill(val); }
   double L = Lvec.max();    

   // Gradient quantities
   mat grad(size(beta));   
   vec Delta;
   
   // Misc
   double objective_fixed, objective_tmp, nrm, quad; // for backtracking 
   uvec fused_w_prev = zeros<uvec>(nchain);
   uvec fused_w_next = zeros<uvec>(nchain);
   uvec start_best = startc, end_best = endc;
   mat beta_best(beta), beta_old, beta_tmp;
   mat beta_adj(beta); // adjusted solution
   double theta = 1.0, theta_old; // acceleration parameter   
   double objective_best(objective), objective_old; 
   bool backtrack, lam3_positive = (lam3 > 0.0);        
   uword count;

   // FISTA
   for (i=1; i<=maxit; i++) {
      if (i % 20 == 0)
         Rcpp::checkUserInterrupt();

      // Store previous results
      beta_old = beta;
      theta_old = theta;

      // Gradient
      for (g=0; g<nchain; g++) {
         grad.col(g).zeros();
         for (j=startc(g); j<=endc(g); j++)
            grad.col(g) += arma::trans(x.slice(j)) * 
              (x.slice(j) * beta.col(g) - y.col(j)); 
         if (g > 0) {
            Delta = beta_adj.col(g) - beta_adj.col(g-1); 
            grad.col(g) += (lam2*w(endc(g-1))/arma::norm(Delta)) * Delta; }
         if (g < (nchain-1)) {
            Delta = beta_adj.col(g+1) - beta_adj.col(g); 
            grad.col(g) -= (lam2*w(endc(g))/arma::norm(Delta)) * Delta; }
         if (lam3_positive) {
            if (intercept) {
               grad(span(0,p-N),g) += 
                  (lam3*len(g)) * beta_adj(span(0,p-N),g); }
            else { grad.col(g) += (lam3*len(g)) * beta_adj.col(g); }
         }
      }
      nrm = arma::accu(arma::square(grad));    

      // Backtracking
      objective_fixed = objective_chains(beta_adj, x, y, 0.0, lam2, 
         lam3, w, intercept, startc, endc);
      count = 0;
       do {
         // Gradient step
         beta_tmp = beta_adj -  grad / L;
         // Proximal step
         beta = beta_tmp;
         if (lam1 > 0.0) {
            for (g=0; g<nchain; g++) 
               beta.col(g) -= clamp(beta_tmp.col(g), -lam1*len(g)/L,
                  lam1*len(g)/L); 
            if (intercept) 
               beta.tail_rows(N) = beta_tmp.tail_rows(N); 
         } 

         // Objective at proximal point (omit l1 penalty)
         objective_tmp = objective_chains(beta, x, y, 0.0, lam2, lam3,
            w, intercept, startc, endc);

         // Quadratic approximation (omit l1 penalty)
         quad = objective_fixed + 0.5 * L * arma::accu(arma::square(beta-beta_tmp))
            - (0.5 / L) * nrm;

         // Check if objective no greater than approximation
         backtrack = (objective_tmp > quad); 
         if (backtrack) { L *= eta; }

      count++;
      } while (backtrack && count < 100);
//    if (count == 100) Rprintf("More than 100 backtracking updates\n");
            
      // Nesterov acceleration method
      theta = 0.5 * (1.0 + sqrt(1.0 + 4.0 * theta * theta));
      beta_adj = beta + ((theta_old - 1.0) / theta) * (beta - beta_old);
      
      // Check if new fusions have occurred  
      for (g=1; g<nchain; g++) {
         if (approx_equal(beta_adj.col(g-1),beta_adj.col(g),
            "both",tol_equal,tol_equal)) 
            {fused_w_prev(g) = 1; fused_w_next(g-1) = 1;}
      }
      // If yes, update fusion chains
      if (any(fused_w_prev)) {
         uvec start_tmp = startc, end_tmp = endc;  
         uvec idx = find(fused_w_prev == 0);
         beta = beta.cols(idx);
         beta_adj = beta_adj.cols(idx);
         startc = startc(idx);
         idx = find(fused_w_next == 0);
         endc = endc(idx);
         nchain = startc.n_elem;
         len = endc - startc + ones<uvec>(nchain);      
         fused_w_prev.set_size(nchain); fused_w_prev.zeros();
         fused_w_next.set_size(nchain); fused_w_next.zeros();
         grad.set_size(size(beta));

         // Recalculate Lipschitz constants
         for (g=0; g<nchain; g++) {
            idx = find(start_tmp == startc(g)); 
            if (end_tmp(idx(0)) == endc(g)) continue;        
            sumxtx.zeros();
            for (j=startc(g); j<=endc(g); j++) 
               sumxtx += trans(x.slice(j)) * x.slice(j);   
            d = arma::eig_sym(sumxtx);
            val = d(p-1) + (lam3*len(g));
            Lvec.subvec(startc(g),endc(g)).fill(val);
         }
         L = Lvec.max();

         // Restart Nesterov acceleration sequence
         theta = 1.0;
      }  

      // Update objective and check convergence
      if (i % 10 == 0) {
         objective_old = objective;
         objective = objective_chains(beta, x, y,
            lam1, lam2, lam3, w, intercept, startc, endc);   
         if (objective < objective_best) {
               beta_best = beta;
               start_best = startc;
               end_best = endc;
               objective_best = objective; }
         if (objective_best <= tol_fista || 
            abs(objective-objective_old) <= tol_fista * objective_old)
               break;
      }
      
   } 

   return List::create(Named("beta") = beta_best,
      Named("objective") = objective_best,
      Named("start") = start_best,
      Named("end") = end_best);  
}
   
   

// FIXED PART OF SUBGRADIENT OF OBJECTIVE OVER FUSION CHAIN
// 0.5 * sum_t || y_t - X_t beta ||_2^2 + lam1 * sum_t || beta ||_1 
//  + lam2 * wprev * || beta - bprev ||_2 + lam2 * wnext * || beta - bnext ||_2
// + 0.5 * lam3 * sum_t || beta ||_2^2
// [[Rcpp::export]]
mat grad_cpp(const vec & beta, const cube & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, 
   const double wprev, const vec & bprev, const double wnext, 
   const vec & bnext, const bool intercept)
{
   uword i, N = y.n_rows, n = y.n_cols, p = beta.n_elem;
   vec s = lam1 * sign(beta), Delta;
   if (lam3 > 0.0) s += lam3 * beta;
   if (intercept) s.tail(N).zeros();
   mat grad(p,n);  
   for (i=0; i<n; i++) 
      grad.col(i) = trans(x.slice(i)) * (x.slice(i) * beta - y.col(i)) + s;
   if (!bprev.is_empty()) {
      Delta = beta - bprev; 
      grad.col(0) += (lam2*wprev/arma::norm(Delta)) * Delta; } 
   if (!bnext.is_empty()) {
      Delta = bnext - beta; 
      grad.col(n-1) -= (lam2*wnext/arma::norm(Delta)) * Delta; } 
   
   return grad;
}


   


// FIND SUBGRADIENT OF MINIMUM NORM (MIN_u,v || Z + U + V D' ||_F)
// BY PROJECTED GRADIENT (BOX CONSTRAINTS ON U, L2-NORM CONSTRAINTS ON V)
// [[Rcpp::export]]
mat pg_cpp(const mat & z, const uvec & fixg, double lam1, double lam2, 
   vec w, int maxit, double tol) 
{
uword i,j;
uword n = z.n_cols;
bool nge3 = (n >= 3);
double nrm;
const double nlam1 = -lam1;
sp_mat D(n,n-1); // differencing matrix
D.diag(0).fill(-1.0);
D.diag(-1).ones(); 
sp_mat Dt = trans(D);
w *= lam2;

// Initialization 
mat u = -clamp(z,nlam1,lam1);
u(fixg).zeros();
mat v = -trans(solve(mat(Dt*D),trans((z+u)*D)));
for (j=0; j<(n-1); j++) {   
   nrm = arma::norm(v.col(j));
   if (nrm > w(j)) v.col(j) *= w(j)/nrm; }
mat u_old, v_old;
mat u_adj(u), v_adj(v);
mat v_best(v);
double theta_old, theta = 1;
mat grad_u = u + z;
grad_u.col(0) -= v.col(0);
grad_u.col(n-1) += v.col(n-2);
if (nge3) grad_u.cols(1,n-2) -= diff(v,1,1);  
double objective = arma::norm(grad_u,"fro");
if (objective <= tol) return (grad_u); 
double objective_best = objective;
double objective_old;


for (i=1; i<=maxit; i++) {
   u_old = u;
   v_old = v;
   theta_old = theta;
   if (i % 50 == 0)
      Rcpp::checkUserInterrupt();

   // Gradient step
   grad_u = u_adj + z;
   grad_u.col(0) -= v_adj.col(0);
   grad_u.col(n-1) += v_adj.col(n-2);
   if (nge3) grad_u.cols(1,n-2) -= diff(v_adj,1,1);  
   u = u_adj - 0.2 * grad_u;
   v = v_adj - 0.2 * diff(grad_u,1,1);

   // Projection step
   u(fixg).zeros();
   u = clamp(u,nlam1,lam1);
   for (j=0; j<n-1; j++) {   
      nrm = arma::norm(v.col(j));
      if (nrm > w(j)) v.col(j) *= w(j)/nrm; }

   // Acceleration step
   theta = 0.5 * (1.0 + sqrt(1.0+4.0*theta*theta));
   u_adj = u + ((theta_old-1.0)/theta) * (u - u_old);
   v_adj = v + ((theta_old-1.0)/theta) * (v - v_old);

   // Update objective and check convergence every 10 iterations
   if (i % 10 == 0 || i == maxit) {
      objective_old = objective;
      grad_u = u + z;
      grad_u.col(0) -= v.col(0);
      grad_u.col(n-1) += v.col(n-2);
      if (nge3) grad_u.cols(1,n-2) -= diff(v,1,1);    
      objective = arma::norm(grad_u,"fro");
      if (objective < objective_best) {
         objective_best = objective;
         v_best = v;}
      if (abs(objective-objective_old) 
           <= tol * objective_old) break;  
   }
} // end for

// Assign z+vD' to grad_u 
grad_u = z;
grad_u.col(0) -= v_best.col(0);
grad_u.col(n-1) += v_best.col(n-2);
if (nge3) grad_u.cols(1,n-2) -= diff(v_best,1,1);
u(fixg).zeros();
u = clamp(grad_u,nlam1,lam1);
return (grad_u-u);
}

 


// [[Rcpp::export]]
LogicalVector compare_cols(const mat & x, const double tol)
{
   uword i, n = x.n_cols;
   LogicalVector out(n-1);
   for (i=0; i<n-1; i++)
      out[i] = approx_equal(x.col(i),x.col(i+1),"both",tol,tol);
   return out;	
}


// [[Rcpp::export]]
bool all_equal(const mat & x, const mat & y, const double tol)
{return approx_equal(x,y,"both",tol,tol);}







