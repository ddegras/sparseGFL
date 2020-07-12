// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace std;
using namespace arma;
using namespace Rcpp;


// FUNCTION PROTOTYPES

double objective_global_svar (const cube &, const mat &, const mat &, 
   const double, const double, const vec &, const bool); 

double objective_local_svar (const mat &, const mat &, const mat &, 
   const double, const double, const double, const double, const mat &, 
   const double, const mat &, const bool);

double objective_blocks_svar (const cube &, const mat &, const mat &, 
   const double, const double, const vec &, const double,
   const mat &, const double, const mat &, const bool); 

double objective_chains_svar (const cube &, const mat &, const mat &, 
   const double, const double, const vec &,
   const bool intercept, const uvec &, const uvec &); 

mat trivial_local_svar (const mat &, const mat &, const double, 
   const double, const double, const mat &, const double, 
   const mat &, const bool, const bool, const double);

mat trivial_mm_svar (const double, const double, const double, 
   const mat &, const double, const mat &, const double, 
   const mat &, const bool, const bool, const double);

mat fp1_svar(mat, const double, const mat &, const double, const double, 
   const double, const mat &, const bool, const bool, const int, 
   const double);

mat fp2_svar(mat, const double, const mat &, const double, const double, 
   const double, const mat &, const double, const mat &, const bool, 
   const bool, const int, const double);

List FISTA_1_svar(mat, const mat &, const mat &, const double, 
   const double, const double, double, mat &, double, mat &,  double, 
   const bool, const bool, const uword, const double, const uword, 
   const double, const double);

List FISTA_M_svar(cube, const mat &, const mat &, const double, 
   const double, const vec &, const bool, const bool, uvec, uvec, 
   double, const uword, const double, const double);
 
mat grad_svar(const mat &, const mat &, const mat &, const double, 
   const double, const double, const mat &, const double, const mat &, 
   const bool, const bool);

LogicalVector compare_slices(const cube &, const double);







// OBJECTIVE FUNCTION (GLOBAL)
// [[Rcpp::export]]
double objective_global_svar (const cube & A, const mat & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const vec & w, 
   const bool intercept) 
{
   
   uword i, m = x.n_rows, n = x.n_cols;
   double sse = 0.0, lasso = 0.0, tv = 0.0, el = 0.0;

   // Sum of squared errors
   for (i=0; i<n; i++) 
      sse += pow(norm(y.col(i)-A.slice(i)*x.col(i)),2); 
   sse *= 0.5;

   // l1 penalty
   if (lam1 > 0.0) {
      if (intercept) 
         { lasso = lam1 * accu(abs(A.cols(0,m-2))); }
      else 
         { lasso = lam1 * accu(abs(A)); }
   }

   // TV penalty
   if (lam2 > 0.0) {
      for(i=0; i<(n-1); i++)
         tv += w(i) * norm(A.slice(i+1)-A.slice(i),"fro");
      tv *= lam2; 
   }

   // Elastic net penalty 
   if (lam3 > 0.0) {
      if (intercept) 
         { el = 0.5 * lam3 * accu(pow(A.cols(0,m-2),2)); }
      else 
         { el = 0.5 * lam3 * accu(pow(A,2)); }
   }

   return (sse+lasso+tv+el);
}

   


//****************************************************************************




// OBJECTIVE FUNCTION (SINGLE BLOCK OR SINGLE FUSION CHAIN)
// [[Rcpp::export]]
double objective_local_svar (const mat & A, const mat & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const double wprev, 
   const mat & Aprev, const double wnext, const mat & Anext, const bool intercept) 
{
   
   uword m = x.n_rows, n = x.n_cols;

   // Sum of squared errors
   double sse = 0.5 * accu(pow(A*x-y,2)); 

   // l1 penalty
   double lasso = 0.0;
   if (lam1 > 0.0) {
      if (intercept) 
         { lasso = lam1 * n * accu(abs(A.head_cols(m-1))); }
      else 
         { lasso = lam1 * n * accu(abs(A)); }
   }

   // TV penalty
   double tv = 0.0; 
   if (lam2 > 0.0) {
      if (!Aprev.is_empty())
         tv += lam2 * wprev * norm(A-Aprev,"fro");
      if (!Anext.is_empty())
         tv += lam2 * wnext * norm(A-Anext,"fro");
   }

   // Elastic net penalty 
   double el = 0.0;
   if (lam3 > 0.0) {
      if (intercept) 
         { el = 0.5 * lam3 * n * accu(pow(A.head_cols(m-1),2)); }
      else 
         { el = 0.5 * lam3 * n * accu(pow(A,2)); }
   }

   return (sse+lasso+tv+el);
}




//****************************************************************************




// OBJECTIVE FUNCTION (MULTIPLE CONSECUTIVE BLOCKS, NOT NECESSARILY FUSED)
// [[Rcpp::export]]
double objective_blocks_svar (const cube & A, const mat & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const vec & w, 
   const double wprev, const mat & Aprev, const double wnext, const mat & Anext, 
   const bool intercept) 
{
   
   uword i, m = x.n_rows, n = x.n_cols;

   // Sum of squared errors
   double sse = 0.0;
   for (i=0; i<n; i++) 
      sse += accu(pow(A.slice(i)*x.col(i)-y.col(i),2)); 
   sse *= 0.5;

   // l1 penalty
   double lasso = 0.0;
   if (lam1 > 0.0) {
      if (intercept) 
         { lasso = lam1 * accu(abs(A.cols(0,m-2))); }
      else 
         { lasso = lam1 * accu(abs(A)); }
   }

   // TV penalty
   double tv = 0.0; 
   if (lam2 > 0.0) {
      for(i=0; i<(n-1); i++)
         tv += w(i) * norm(A.slice(i+1)-A.slice(i),"fro");
      if (!Aprev.is_empty())
         tv += wprev * norm(A.slice(0)-Aprev,"fro");
      if (!Anext.is_empty())
         tv += wnext * norm(A.slice(n-1)-Anext,"fro");
      tv *= lam2; }

   // Elastic net penalty
   double el = 0.0;
   if (lam3 > 0.0) {
      if (intercept) 
         { el = 0.5 * lam3 * accu(pow(A.cols(0,m-2),2)); }
      else  
         { el = 0.5 * lam3 * accu(pow(A,2)); }
   }

   return (sse+lasso+tv+el);
}

    


//****************************************************************************




// OBJECTIVE FUNCTION (ALL FUSION CHAINS)
double objective_chains_svar (const cube & A, const mat & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, const vec & w,
   const bool intercept, const uvec & startc, const uvec & endc) 
{
   uword g, G = startc.n_elem;
   uvec len = endc - startc + ones<uvec>(G);

   // Sum of squared errors
   double sse = 0.0;
   for (g=0; g<G; g++) 
      sse += pow(norm(A.slice(g) * x.cols(startc(g),endc(g)) - 
         y.cols(startc(g),endc(g)),"fro"),2);
   sse *= 0.5;

   // l1 penalty
   double lasso = 0.0; 
   uword m = A.n_cols;
   if (lam1 > 0.0) {
      if (intercept) {
         for (g=0; g<G; g++) 
            lasso += len(g) * 
               accu(abs(A(span(),span(0,m-2),span(g)))); }
      else if (!intercept) {
         for (g=0; g<G; g++) 
            lasso += len(g) * accu(abs(A.slice(g))); }
      lasso *= lam1; 
   } 

   // TV penalty
   double tv = 0.0; 
   if (lam2 > 0.0) {
      for(g=0; g<(G-1); g++)
         tv += w(endc(g)) * norm(A.slice(g+1)-A.slice(g),"fro");
      tv *= lam2; }

   // Elastic net penalty
   double el = 0.0;
   if (lam3 > 0.0) {
      if (intercept) 
         { for (g=0; g<G; g++) 
            el += len(g) * 
               accu(pow(A(span(),span(0,m-2),span(g)),2)); }
      else 
         { for (g=0; g<G; g++)
            el += len(g) * accu(pow(A.slice(g),2)); }
      el *= 0.5 * lam3; 
   }
   
   return (sse+lasso+tv+el);
}
   





//****************************************************************************




// CHECK FOR TRIVIAL SOLUTIONS IN LOCAL OBJECTIVE 
// [[Rcpp::export]]
mat trivial_local_svar (const mat & x, const mat & y, const double lam1, 
   const double lam2, const double lam3, const double wprev, 
   const mat & Aprev, const double wnext, const mat & Anext, 
   const bool lag0, const bool intercept, const double tol)
{

   mat null_mat(1,0);
   bool is_null_Aprev = (Aprev.is_empty() || wprev == 0.0);
   bool is_null_Anext = (Anext.is_empty() || wnext == 0.0);

   // Return empty matrix if no candidate for trivial solution
   if (is_null_Aprev && is_null_Anext) return null_mat;
      
   bool single_candidate = (is_null_Aprev || is_null_Anext || 
      approx_equal(Aprev,Anext,"both",tol,tol));
   double d; // distance b/w two candidate solutions
   mat Delta;
   if (!single_candidate) {
      Delta = Anext - Aprev;
      d = norm(Delta,"fro"); 
   }
   uword N = y.n_rows, m = x.n_rows, n = x.n_cols;
   double lam1new = lam1 * n, lam3new = lam3 * n;
   mat grad(N,m); // subgradient 
   uvec idx; // index of zero values in candidate solution

   // Check whether previous coefficient block minimizes local
   // objective (if it exists)
   if (!is_null_Aprev) {
      grad = (Aprev*x-y) * trans(x);
      if (lam3 > 0.0) {
         if (intercept)
            { grad.head_cols(m-1) += 
               lam3new * Aprev.head_cols(m-1); }
         else 
            { grad += lam3new * Aprev; }
      }
      if (!single_candidate) 
         grad -= (lam2*wnext/d) * Delta; 
      if (lam1 > 0.0) {
         idx = find(Aprev == 0.0);
         if (intercept) 
            { idx = idx.elem(find(idx < (m-1)*N));      
            grad.head_cols(m-1) += 
               lam1new * sign(Aprev.head_cols(m-1)); }
         else
            { grad += lam1new * sign(Aprev); }
         if (idx.n_elem > 0)
            grad(idx) -= clamp(grad(idx),-lam1new,lam1new);
      }
      if (lag0) 
         grad.diag().zeros();
      if (single_candidate) {
         if (norm(grad,"fro") <= (lam2*(wprev+wnext))) 
            {return Aprev;} else {return null_mat;}   
      } else if (norm(grad,"fro") <= (lam2*wprev)) 
         {return Aprev;}     
   }

   // Check whether next coefficient block minimizes local 
   // objective (if it exists)
   if (!is_null_Anext) {
      grad = (Anext*x-y) * trans(x);
      if (lam3 > 0.0) {
         if (!intercept)
            { grad += lam3new * Anext; }
         else if (intercept)
            { grad.head_cols(m-1) += 
               lam3new * Anext.head_cols(m-1); }
      }
      if (!is_null_Aprev) 
         grad += (lam2*wprev/d) * Delta; 
      if (lam1 > 0.0) {
         idx = find(Anext == 0.0);
         if (intercept)
            {idx = idx.elem(find(idx < (m-1)*N));      
            grad.head_cols(m-1) += 
               lam1new * sign(Anext.head_cols(m-1));}
         else 
            {grad += lam1new * sign(Anext);}
         if (idx.n_elem > 0)
            grad(idx) -= clamp(grad(idx),-lam1new,lam1new);
      }
      if (lag0) 
         grad.diag().zeros();
      if (norm(grad,"fro") <= wnext*lam2) return Anext;
   } 

   return null_mat;
}




//****************************************************************************






// CHECK FOR TRIVIAL SOLUTIONS IN MAJORIZATION-MINIMIZATION 
// [[Rcpp::export]]
mat trivial_mm_svar (const double lam1, const double lam2, const double L, 
   const mat & Z, const double wprev, const mat & Aprev, const double wnext, 
   const mat & Anext, const bool lag0, const bool intercept, const double tol)
{

   mat null_mat(1,0);
   bool is_null_Aprev = (Aprev.is_empty() || wprev == 0.0);
   bool is_null_Anext = (Anext.is_empty() || wnext == 0.0);
   bool single_candidate = (is_null_Aprev || is_null_Anext || 
      approx_equal(Aprev,Anext,"both",tol,tol));
   double d; // distance b/w two candidate solutions
   mat Delta;
   if (!single_candidate) 
      {Delta = Anext - Aprev;
      d = norm(Delta,"fro");} 
   uword m = Z.n_cols;
   mat grad(size(Z)); // subgradient 
   uvec idx; // index of null values in candidate solution
   
  
   // Check whether previous coefficient block minimizes quadratic 
   // approximation to local objective (if it exists)
   if (!is_null_Aprev) {
      grad = L * (Aprev - Z);
      if (!single_candidate) 
         { grad -= (lam2*wnext/d) * Delta; } 
      if (intercept) 
         { idx = find(Aprev.head_cols(m-1) == 0.0);      
         grad.head_cols(m-1) += lam1 * sign(Aprev.head_cols(m-1)); }
      else 
         { idx = find(Aprev == 0.0);
         grad += lam1 * sign(Aprev); }
      if (idx.n_elem > 0)
         grad(idx) -= clamp(grad(idx),-lam1,lam1);
      
      if (lag0) grad.diag().zeros();
      if (single_candidate) {
         if (norm(grad,"fro") <= (lam2*(wprev+wnext))) 
            { return Aprev; } else { return null_mat; } }   
      else if (norm(grad,"fro") <= lam2*wprev) 
         { return Aprev; }
   }

   // Check whether next coefficient block minimizes quadratic
   // approximation to local objective (if it exists)
   if (!is_null_Anext) {
      grad = L * (Anext - Z);
      if (!is_null_Aprev) 
         grad += (lam2 * wprev / d) * Delta; 
      if (intercept)
         { idx = find(Anext.head_cols(m-1) == 0.0);
         grad.head_cols(m-1) += 
            lam1 * sign(Anext.head_cols(m-1)); }
      else
         { idx = find(Anext == 0.0);
         grad += lam1 * sign(Anext); }
      if (idx.n_elem > 0)
         grad(idx) -= clamp(grad(idx),-lam1,lam1); 

      if (lag0) grad.diag().zeros();
      if (norm(grad,"fro") <= lam2*wnext) return Anext;
   } 

   return null_mat;
}




//****************************************************************************





// FIXED POINT ITERATION (1 NEIGHBOR)
// [[Rcpp::export]]
mat fp1_svar(mat A, const double L, const mat & Z, const double lam1,
   const double lam2, const double w, const mat & Aprev, const bool lag0, 
   const bool intercept, const int maxit, const double tol)
{

   uword i, m = A.n_cols;
   double tol_equal = tol * sqrt(A.n_elem);
   double nlam1 = -lam1;
   double dAU = norm(A-Aprev,"fro");
   const mat LZ = L * Z;
   double objective = 0.5 * L * accu(pow(A-Z,2)) + lam2 * w * dAU; 
   if (intercept)  
      { objective += lam1 * accu(abs(A.head_cols(m-1))); }
   else if (!intercept)
      { objective += lam1 * accu(abs(A)); }
   double objective_best = objective;
     
   mat A_old(size(A)); 

   for (i=1; i<=maxit; i++) {

     if (i % 1000 == 0)
         Rcpp::checkUserInterrupt();

      if (i % 10 == 0)
         A_old = A; 
      
      // Check if neighboring block is trivial solution
      if (dAU <= tol_equal) return Aprev;

      // If not, fixed point iteration 
      A = LZ + (lam2*w/dAU) * Aprev;

      // Soft-thresholding 
      if (lam1 > 0.0)
         { if (intercept)
            { A.head_cols(m-1) -= clamp(A.head_cols(m-1),nlam1,lam1); } 
         else 
            { A -= clamp(A,nlam1,lam1); }
      }

//      if (lam1 > 0.0)
//         A = soft_svar(A,lam1,intercept);

      // Scaling
      A /=  L + (lam2*w/dAU);

      // Apply diagonal constraints if required
      if (lag0) A.diag().zeros();

      // Distance between solution at current block and neighbor
      dAU = norm(A-Aprev,"fro");

      // Check convergence
      if (i % 10 == 0) 
         { if (approx_equal(A,A_old,"both",tol,tol)) 
            break;
         objective = (0.5 * L * accu(pow(A-Z,2))) + (lam2 * w * dAU); 
         if (intercept)  
            { objective += lam1 * accu(abs(A.head_cols(m-1))); }
         else 
            { objective += lam1 * accu(abs(A)); }
         if (objective >= (1-tol) * objective_best)
            break;
         if (objective < objective_best)
            objective_best = objective;
      }

   } // END (for i)    

   return A;

}






//****************************************************************************






// FIXED POINT ITERATION (2 NEIGHBORS)
// [[Rcpp::export]]
mat fp2_svar(mat A, const double L, const mat & Z, const double lam1, 
   const double lam2, const double wprev, const mat & Aprev, 
   const double wnext, const mat & Anext, const bool lag0, 
   const bool intercept, const int maxit, const double tol)
{

   uword i, m = A.n_cols;
   double tol_equal = tol * sqrt(A.n_elem);
   double nlam1 = -lam1;
   mat A_old(size(A)); 
   if (lag0) A.diag().zeros();
   double dAU = norm(A-Aprev,"fro");
   double dAV = norm(A-Anext,"fro");
   const mat LZ = L * Z;
   double objective = 0.5 * L * accu(pow(A-Z,2)) 
      + lam2 * (wprev * dAU + wnext * dAV); 
   if (intercept)  
      { objective += lam1 * accu(abs(A.head_cols(m-1))); }
   else 
      { objective += lam1 * accu(abs(A)); }
   double objective_best = objective;


   for (i = 1; i <= maxit; i++) {

      if (i % 1000 == 0)
         Rcpp::checkUserInterrupt();

      if (i % 10 == 0)
         A_old = A;

      // Check for trivial solutions
      if (dAU <= tol_equal) return Aprev;
      if (dAV <= tol_equal) return Anext;    

      // Fixed point iteration 
      A = LZ + ((lam2*wprev/dAU) * Aprev) 
         + ((lam2*wnext/dAV) * Anext);

      // Soft-thresholding 
         if (intercept) 
            { A.head_cols(m-1) -= 
            clamp(A.head_cols(m-1),nlam1,lam1); }
         else 
            { A -= clamp(A,nlam1,lam1); }
      
      // Scaling
      A /=  L + (lam2*wprev/dAU) + (lam2*wnext/dAV);

      // Apply diagonal constraints if required
      if (lag0) A.diag().zeros();

      // Distance between solution and neighboring blocks
      dAU = norm(A-Aprev,"fro");
      dAV = norm(A-Anext,"fro");
 
      // Check convergence
       if (i % 10 == 0) 
         { if (approx_equal(A,A_old,"both",tol,tol)) 
            break;
         objective = (0.5 * L * accu(pow(A-Z,2))) 
            + lam2 * (wprev * dAU + wnext * dAV); 
         if (intercept)  
            { objective += lam1 * accu(abs(A.head_cols(m-1))); }
         else 
            { objective += lam1 * accu(abs(A)); }
         if (objective >= (1-tol) * objective_best)
            break;
         if (objective < objective_best)
            objective_best = objective;
      }

   }     

   return A;

}







//****************************************************************************



//#################
// FISTA FUNCTIONS
//#################


// [[Rcpp::export]]
Rcpp::List FISTA_1_svar(mat A, const mat & x, const mat & y, const double lam1, 
   const double lam2, const double lam3, double wprev, mat & Aprev, 
   double wnext, mat & Anext,  double L, const bool lag0, 
   const bool intercept, const uword maxit_fista, const double tol_fista, 
   const uword maxit_fp, const double tol_fp, const double tol_equal)
{
      
   bool first = Aprev.is_empty(), last = Anext.is_empty();
   uword nbounds;
   if ((first && last) || lam2 == 0.0) { nbounds = 0; }
   else if (first || last) { nbounds = 1; }
   else { nbounds = 2; }
   uword i;
   double objective = objective_local_svar(A, x, y, lam1, lam2, lam3,
      wprev, Aprev, wnext, Anext, intercept);
   mat A_best = A;
   double objective_best = objective, objective_init = objective, 
      objective_old;   
   mat A_adj = A; // adjusted solution
   double theta = 1.0, theta_old; // acceleration parameter
   mat A_old(size(A)), Z(size(A));   
   bool progress;
   if (nbounds == 0 && L == 0.0) 
     { A.zeros(); 
      objective_best = 0.5 * accu(pow(y,2));
      progress = (objective_best < objective_init);
      return List::create(Named("A") = A, 
                          Named("objective") = objective, 
                          Named("progress") = progress); }

   mat xtx, ytx;   
   uword m = x.n_rows, n = x.n_cols;
   const bool L_positive = (L > 0.0), lam3_positive = (lam3 > 0.0);
   bool store_cp = (m <= n && L_positive);
   if (store_cp) 
      { xtx = x * trans(x) / L;
        ytx = y * trans(x) / L; }
   const double lam1new = lam1 * n, lam3new = lam3 * n;


   // FISTA
   for (i=1; i<=maxit_fista; i++) {
      
      // Store previous results
      A_old = A;
      objective_old = objective;
      theta_old = theta;
 
      if (i % 100 == 0)
         Rcpp::checkUserInterrupt();

      // Gradient step
      if (L_positive) {
         if (store_cp) 
            { Z = A_adj - (A_adj*xtx-ytx); }
         else
            { Z =  A_adj - (A_adj*x-y) * trans(x/L); }
         if (lam3_positive) {
            if (intercept) {
               Z.head_cols(m-1) -= (lam3new/L) * A_adj.head_cols(m-1); } 
            else 
               { Z -= (lam3new/L) * A_adj; }
         }
      }
      else { Z = A_adj; }
      if (lag0) Z.diag().zeros();

      // Check for trivial solutions 
      A = trivial_mm_svar(lam1new, lam2, L, Z, wprev, Aprev, 
         wnext, Anext, lag0, intercept, tol_equal);
         
      // If no trivial solution, apply fixed point iteration 
      // (or soft-thresholding if no neighbors)
      if (A.is_empty()) {
         if (approx_equal(A_adj, Aprev, "both", tol_equal, tol_equal) 
            || approx_equal(A_adj, Anext, "both", tol_equal, tol_equal)) {
               A = Z; } else { A = A_adj; }

         if (lag0) A.diag().zeros();
         switch (nbounds) {
            case 0:
               if (intercept) {
                  A.head_cols(m-1) = Z.head_cols(m-1) - 
                     clamp(Z.head_cols(m-1),-lam1new/L,lam1new/L);
                  A.col(m-1) = Z.col(m-1);
               } else {
                  A = Z - clamp(Z,-lam1new/L,lam1new/L); }
               break;
            case 1:
               if (first) {
               A = fp1_svar(A, L, Z, lam1new, lam2, wnext, Anext, 
                  lag0, intercept, maxit_fp, tol_fp); } else {
               A = fp1_svar(A, L, Z, lam1new, lam2, wprev, Aprev, lag0, 
                  intercept, maxit_fp, tol_fp); }
               break;
            case 2:
               A = fp2_svar(A, L, Z, lam1new, lam2, wprev, Aprev, wnext,
               Anext, lag0, intercept, maxit_fp, tol_fp);
               break;
            }
      }
                  
      // Nesterov method
      theta = (1.0 + sqrt(1.0+4.0*theta*theta)) / 2.0;
      A_adj = A + ((theta_old-1.0) / theta) * (A - A_old);
      
      // Update objective (every iteration? every 10 iterations?)
      objective = objective_local_svar(A, x, y, lam1, lam2, lam3,
         wprev, Aprev, wnext, Anext, intercept);         
//  Rcout << "Iteration " << i << " Objective " << objective << std::endl;  
       if (std::isinf(objective)) 
         Rcpp::stop("Unexpected error: infinite objective");
      if (std::isnan(objective)) 
         Rcpp::stop("Unexpected error: NaN objective");
      if (objective < objective_best) {
         A_best = A;
         objective_best = objective;}
      
      // Check convergence
      if (abs(objective - objective_old) <= tol_fista * objective_old ||
         objective <= tol_fista) break;
   }
   
   progress = (objective_best < objective_init); 

   return List::create(Named("A") = A_best,
                       Named("objective") = objective_best,
                       Named("progress") = progress);

}
   





//****************************************************************************







// FISTA OPTIMIZATION OVER FIXED CHAINS
// [[Rcpp::export]]
List FISTA_M_svar(cube A, const mat & x, const mat & y, const double lam1, 
   const double lam2, const double lam3, const vec & w, const bool lag0, 
   const bool intercept, uvec startc, uvec endc, double eta, 
   const uword maxit, const double tol_fista, const double tol_equal)
{
   uword g,i; 
   uword nchain = startc.n_elem; // number of chains
   uvec len = endc - startc + ones<uvec>(nchain); // chains lengths
   uword m = A.n_cols; // number of predictors per block

   if (A.n_slices > nchain) // reduce A if needed
      A = A.slices(startc);
        
   // Initial objective value
   double objective = objective_chains_svar(A,x,y,lam1,lam2,lam3,
      w,intercept,startc,endc);   

   // Lipschitz constants for gradient step (estimates)
   vec Lvec(x.n_cols);  
   double d, val;
   for (g=0; g<nchain; g++) {
      d = norm(x.cols(startc(g),endc(g)));
      val = (d * d) + (lam3 * len(g));
      Lvec.subvec(startc(g),endc(g)).fill(val); }
   double L = Lvec.max();    
 
   // Gradient quantities
   cube grad(size(A));   
   mat Delta;
   
   // Misc
   double objective_fixed, objective_tmp, nrm, quad; // for backtracking 
   uvec fused_w_prev = zeros<uvec>(nchain);
   uvec fused_w_next = zeros<uvec>(nchain);
   uvec start_best = startc, end_best = endc;
   cube A_best(A), A_old, A_tmp;
   cube A_adj(A); // adjusted solution
   double theta = 1.0, theta_old; // acceleration parameter   
   double objective_best(objective), objective_old; 
   bool backtrack, lam3_positive = (lam3 > 0.0);        

   // FISTA
   for (i=1; i<=maxit; i++) {
      
      if (i % 50 == 0) 
         Rcpp::checkUserInterrupt();

      // Store previous results
      A_old = A;
      theta_old = theta;

      // Gradient
      for (g=0; g<nchain; g++) {
         grad.slice(g) = 
            (A_adj.slice(g) * x.cols(startc(g),endc(g)) 
               - y.cols(startc(g),endc(g))) 
               * trans(x.cols(startc(g),endc(g)));
         if (g > 0) {
            Delta = A_adj.slice(g) - A_adj.slice(g-1); 
            grad.slice(g) += (lam2*w(endc(g-1))/norm(Delta,"fro")) * Delta; }
         if (g < (nchain-1)) {
            Delta = A_adj.slice(g) - A_adj.slice(g+1); 
            grad.slice(g) += (lam2*w(endc(g))/norm(Delta,"fro")) * Delta; }
         if (lam3_positive) {
            if (intercept)
               { grad(span(),span(0,m-2),span(g)) += (lam3*len(g)) *
                    A_adj(span(),span(0,m-2),span(g)); }
            else { grad.slice(g) += (lam3*len(g)) * A_adj.slice(g); }
         }
         if (lag0) 
            grad.slice(g).diag().zeros();
      }
      nrm = accu(pow(grad,2)); // squared Frobenius norm    
 
      // Backtracking
      objective_fixed = objective_chains_svar(A_adj,x,y,0.0,lam2,lam3,
         w,intercept,startc,endc);
      do {
         // Gradient step
         A_tmp = A_adj -  grad / L;
         // Proximal step
         A = A_tmp;
         if (lam1 > 0.0) {
            for (g=0; g<nchain; g++) {
               if (intercept)
                  { A(span(),span(0,m-2),span(g)) -= 
                    clamp(A_tmp(span(),span(0,m-2),span(g)),
                    -lam1*len(g)/L,lam1*len(g)/L); }
               else
               { A.slice(g) -= clamp(A_tmp.slice(g),
               -lam1*len(g)/L,lam1*len(g)/L); }}
         } 

         // Objective at proximal point (omit l1 penalty)
         objective_tmp = objective_chains_svar(A, x, y, 0.0, lam2, lam3,
            w, intercept, startc, endc);

         // Quadratic approximation (omit l1 penalty)
         quad = objective_fixed + 0.5 * L * accu(pow(A-A_tmp,2))
            - (0.5 / L) * nrm;
         // Check if objective no greater than approximation
         backtrack = (objective_tmp > quad); 
         if (backtrack) { L *= eta; }
      } while (backtrack);
                  
      // Nesterov acceleration method
      theta = 0.5 * (1.0 + sqrt(1.0 + 4.0 * theta * theta));
      A_adj = A + ((theta_old - 1.0) / theta) * (A - A_old);
      
      // Check if new fusions have occurred  
      for (g=1; g<nchain; g++) {
         if (approx_equal(A_adj.slice(g-1),A_adj.slice(g),
            "both",tol_equal,tol_equal)) 
            {fused_w_prev(g) = 1; fused_w_next(g-1) = 1;}
      }
      // If yes, update fusion chains
      if (any(fused_w_prev)) {
         uvec start_tmp = startc, end_tmp = endc;  
         uvec idx = find(fused_w_prev == 0);
         A = A.slices(idx);
         A_adj = A_adj.slices(idx);
         startc = startc(idx);
         idx = find(fused_w_next == 0);
         endc = endc(idx);
         nchain = startc.n_elem;
         len = endc - startc + ones<uvec>(nchain);      
         fused_w_prev.set_size(nchain); fused_w_prev.zeros();
         fused_w_next.set_size(nchain); fused_w_next.zeros();
         grad.set_size(size(A));

         // Recalculate Lipschitz constants
         for (g=0; g<nchain; g++) {
            idx = find(start_tmp == startc(g)); 
            if (end_tmp(idx(0)) == endc(g)) 
               { continue ; } 
            else {
               d = norm(x.cols(startc(g),endc(g)));
               val = (d*d) + (lam3*len(g));
               Lvec.subvec(startc(g),endc(g)).fill(val); }}
         L = Lvec.max();

         // Restart Nesterov acceleration sequence
         theta = 1.0;
      }  

      // Update objective and check convergence
      if (i % 10 == 0) {
         objective_old = objective;
         objective = objective_chains_svar(A,x,y,
            lam1,lam2,lam3,w,intercept,startc,endc);   
         if (objective < objective_best) {
               A_best = A;
               start_best = startc;
               end_best = endc;
               objective_best = objective; }
         if (objective_best <= tol_fista || 
            abs(objective-objective_old) <= tol_fista * objective_old)
               break;
      }
      
   } 

   return List::create(Named("A") = A_best,
      Named("objective") = objective_best,
      Named("start") = start_best,
      Named("end") = end_best);  
}
   
   




//****************************************************************************






// FIXED PART OF SUBGRADIENT OF OBJECTIVE OVER FUSION CHAIN
// 0.5 * sum_t || y_t - A x_t ||_2^2 + lam1 * sum_t || A ||_1 
//  + lam2 * wprev * || A - Aprev ||_F + lam2 * wnext * || A - Anext ||_F
// + 0.5 * lam3 * sum_t || A ||_F^2
// [[Rcpp::export]]
mat grad_svar(const mat & A, const mat & x, const mat & y, 
   const double lam1, const double lam2, const double lam3, 
   const double wprev, const mat & Aprev, const double wnext, 
   const mat & Anext, const bool lag0, const bool intercept)
{
   uword i, N = A.n_rows, m = A.n_cols, n = x.n_cols;
   mat s = lam1 * sign(A);
   if (lam3 > 0.0)
      s += lam3 * A;
   if (intercept) 
      s.col(m-1).zeros();
   cube z(N,m,n);  
   for (i=0; i<n; i++) {
      z.slice(i) = (A * x.col(i) - y.col(i)) * trans(x.col(i)) + s;
      if (lag0) z.slice(i).diag().zeros(); }

   mat Delta;
   if (!Aprev.is_empty())  
      { Delta = A - Aprev; 
      z.slice(0) += (lam2*wprev/norm(Delta,"fro")) * Delta; } 
   if (!Anext.is_empty())  
      { Delta = Anext - A; 
      z.slice(n-1) -= (lam2*wnext/norm(Delta,"fro")) * Delta; } 

   z.reshape(m*N,n,1);
   return z.slice(0);
}





//****************************************************************************






// [[Rcpp::export]]
LogicalVector compare_slices(const cube & x, const double tol)
{
	uword i, n = x.n_slices;
	LogicalVector out(n-1);
	for (i=0; i<n-1; i++)
	   out[i] = approx_equal(x.slice(i),x.slice(i+1),"both",tol,tol);
   return out;	
}




//****************************************************************************



// Elastic net: 
// min 0.5 || Y - AX ||_F^2 + lam ( alpha || A ||_1 + (1-alpha)/2 || A ||_F^2 )
// (l1 and l2 penalties do not apply to intercept)



//[[Rcpp::export]]
List AXY_elnet(mat x, mat y, double lambda, double alpha,
   bool intercept, bool diag0, double eps, uword maxit)
{
   vec xbar, ybar;
   if (intercept) {
      xbar = mean(x,1);
      ybar = mean(y,1);
      x.each_col() -= xbar;
      y.each_col() -= ybar; } 

   const double lam1 = lambda * alpha, lam2 = lambda * (1.0-alpha);
   double objective = datum::inf, objective_old;
   double theta = 1.0, theta_old;
   const double L = pow(norm(x,2),2) + lam2;
   mat A = zeros<mat>(y.n_rows, x.n_rows);
   mat A_old(size(A)), A_adj = A, grad(size(A)); 
   uword i;
   bool precalc = (y.n_rows <= y.n_cols && x.n_rows <= x.n_cols);     
   mat xxt, yxt;
   double yty, xty;
   if (precalc) {
      xxt = x * x.t(); 
      yxt = y * x.t(); 
      xty = trace(yxt);
      yty = pow(norm(y,"fro"),2); }

   for (i=0; i<maxit; i++) {

      A_old = A;
      theta_old = theta;
      objective_old = objective;

      // Gradient step
      if (precalc) 
         { grad = A_adj * xxt - yxt + lam2 * A_adj; }
      else  
         { grad = (A_adj * x - y) * x.t() + lam2 * A_adj; }
      A = A_adj - grad / L;

      if (diag0) A.diag().zeros();

      // Proximal step
      A -= clamp(A,-lam1/L,lam1/L);

      // Acceleration step
      theta = (1.0 + sqrt(1.0+4.0*theta*theta)) / 2.0;
      A_adj = A + ((theta_old-1.0) / theta) * (A - A_old);

      // Objective
      if (precalc) 
         { objective = 0.5 * (trace(A * xxt) - 2.0 * xty + yty) 
            + lam1 * accu(abs(A)) + lam2 * 0.5 * pow(norm(A,"fro"),2); }
      else 
         { objective = 0.5 * pow(norm(A*x-y,"fro"),2) + lam1 * accu(abs(A))
         + 0.5 * lam2 * pow(norm(A,"fro"),2); }

      // Check convergence
      if (abs(objective-objective_old) <= eps * objective || 
         approx_equal(A,A_old,"both",eps,eps)) break;
   }

   mat out; 
   if (intercept) { out = join_rows(A, ybar - A * xbar); } else { out = A; }
   return List::create(Named("A") = out, Named("objective") = objective); 

} 

