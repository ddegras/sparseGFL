\name{SGFL.AXY}
\alias{SGFL.AXY}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse Group Fused Lasso for Piecewise Linear Regression 
with Vector Predictor
}
\description{
Fits piecewise linear regression model \eqn{y(t) = A(t) x(t) + e(t), 1 \le t \le T} by Sparse Group Fused Lasso (SGFL) where \eqn{y(t)} is a response vector, \eqn{x(t)} a predictor vector, and \eqn{A(t)} an unknown (sparse) regression matrix.  
}
\usage{
SGFL.AXY(x, y, lambda1, lambda2, alpha = 1, intercept = TRUE, A = NULL, 
	w = NULL, sweep = c("cyclical", "srswr", "srswor"), pattern = c(10,5),
	control = list(), parallel = FALSE, verbose = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
matrix of predictors (dimension \eqn{p}-by-\eqn{T} with \eqn{p} number of predictor variables and \eqn{T} time series length)
}
  \item{y}{
matrix of responses (dimension \eqn{N}-by-\eqn{T} with \eqn{N} number of response variables)
}
  \item{lambda1}{
regularization parameter for lasso penalty
}
  \item{lambda2}{
regularization parameter for total variation (TV) penalty
}
  \item{alpha}{
mixing weight between lasso and TV penalty (\code{alpha=1} is pure lasso, \code{alpha=0} is pure TV)
}
  \item{intercept}{
logical: does the model contain a time-varying intercept?
}
  \item{A}{
optional initial regression matrices as starting point for optimization (3D array, \eqn{N %*% p %*% T}) 
}
  \item{w}{
optional individual weights for TV penalty components (default = 1) 
}
  \item{sweep}{
sweeping pattern: \code{cyclical}, \code{srswr} (simple random sampling with replacement), or  \code{srswor} (simple random sampling without replacement)
}
  \item{pattern}{
optional vector of length 2: maximum number of block coordinate descent cycles 
and maximum number of fusion cycles per epoch 
}
  \item{control}{
optional list of tuning parameters for optimization
}
  \item{parallel}{
logical: use parallel computing? (if \code{TRUE}, parallel setup must be registered beforehand)
}
  \item{verbose}{
logical: display optimization progress?
}
}
\details{
If \code{intercept=TRUE}, the regression model is \eqn{y(t) = A(t) x(t) + b(t) + e(t)} wih \eqn{b(t)} a piecewise constant intercept vector. 
}
\value{
A list with components 
\item{A}{fitted regression matrices as 3D array, the third dimension being time. The \eqn{t}-th slice of \code{A} contains the matrix of coefficient estimates \eqn{[A_k(t),...,A_p(t),\mu(t)]}} 
\item{objective}{best objective value} 
\item{changePoints}{detected change points or \code{NULL}. 
	The first change point is at least \eqn{p+2}}
\item{trace}{objective values at each iteration}	
\item{g}{subgradient of minimum norm at solution \code{A}. This matrix should have most of its values very close to zero.}
\item{iters}{number of iterations performed}
\item{lambda1}{input \code{lambda1}}
\item{lambda2}{input \code{lambda2}}
\item{alpha}{input \code{alpha}}
\item{time}{timing information for all optimization steps}
}


\references{
Degras (2020). Sparse group fused lasso for model segmentation: a hybrid approach. \emph{Advances in Data Analysis and Classification.} \url{https://doi.org/10.1007/s11634-020-00424-5}
}


\seealso{
\code{\link{SGFL}}, \code{\link{SGFL.SVAR}}
}
\examples{
In progress
}


