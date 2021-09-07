# sparseGFL
Hybrid Approach to Sparse Group Fused Lasso

Fast hybrid algorithm for joint model segmentation and sparse estimation in high-dimensional piecewise linear regression models.   
The algorithm finds global solutions to the sparse group fused lasso (SGFL) problem whose objective function is the sum of a squared loss function (to control fit to the data), an elastic net penalty (to promote sparsity in regression coefficients), and a total variation penalty (to promote parsimony in detected change points or segments).   

## Package installation
``` 
install.packages("Rcpp", "RcppArmadillo", "Matrix", "foreach") # package dependencies
install.packages("devtools")
library(devtools)  
install_github("ddegras/sparseGFL") 
```
## Main functions
- `SGFL`: solve SGFL for piecewise regression model y(t) = X(t) b(t) + e(t) 
(matrix predictors X, vector regression coefficients b)
- `SGFL.AXY`: solve SGFL for piecewise regression model y(t) = A(t) x(t) + e(t) 
(vector predictors x, matrix regression coefficients A)
- `SGFL.VAR`: solve SGFL for piecewise vector autoregressive (VAR) model y(t) = A1(t) y(t-1) + ... + Ap(t) y(t-p) + e(t)

### Reference
Degras, D. "Sparse group fused lasso for model segmentation: a hybrid approach." *Advances in Data Analysis and Classification* **15**, 625–671 (2021). 
