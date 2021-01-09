\name{SGFL}
\alias{SGFL}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Sparse Group Fused Lasso For Piecewise Linear Regression with Matrix Predictor}

\description{Sparse group fused lasso for the general piecewise constant linear regression model \eqn{y(t) = X(t) \beta(t) + e(t)} 
}

\usage{SGFL(x, y, lambda1, lambda2, alpha = 1, intercept = TRUE, 
  beta = NULL, w = NULL, L = NULL, sweep = c("cyclical", "srswr", "srswor"), 
  pattern = c(10, 5), control = list(), parallel = FALSE, verbose = FALSE)
}

\arguments{
  \item{x}{A 3D array of predictor variables. The third dimension should be time.}
  \item{y}{A matrix of response variables. Rows correspond to variables, columns to time. The dimensions of \code{x} and \code{y} must satisfy: \eqn{nrow(y)=dim(x)[1]} and \eqn{ncol(y)=dim(x)[3]}}
  \item{lambda1}{
Regularization parameter for the lasso penalty.}
  \item{lambda2}{
Regularization parameter for the total variation penalty.}
  \item{alpha}{
Mixing parameter for the elastic net penalty. See Details. 
}
  \item{intercept}{
Logical: should the regression model include an intercept? Default = TRUE.}
  \item{beta}{
Optional argument: starting point for the optimization algorithm. If provided, this should be a matrix of dimensions \eqn{dim(x)[2]} by \eqn{dim(x)[3]}.
}
  \item{w}{
Optional argument: individual weights for the total variation penalty. If specified, this should be a numerical vector of length \eqn{ncol(y)}. By default, all weights are set to 1. 
}
  \item{L}{
Optional vector of Lipschitz constants for the gradient of the smooth part of the objective function. If specified, should be a numerical vector of length \eqn{ncol(y)}
}
  \item{sweep}{
Optional argument: sweeping pattern for the block coordinate descent and fusion stages of the optimization. If specified, must one of the characters \code{"cyclical"}, \code{"srswr"} (simple sampling with replacement), \code{"srswor"} (simple sampling without replacement). Default = \code{"cyclical"}.}
  \item{pattern}{
Optional argument: vector of length 2 containing the maximum number of consecutive block coordinate descent cycles and the maximum number of consecutive fusion cycles. Default = \eqn{c(10,5)}. 
}
  \item{control}{
Optional list of tuning parameters for the optimization algorithm. 
}
  \item{parallel}{
Logical: should parallel computing be used in the algorithm? Default = \code{FALSE}. If set to \code{TRUE}, the user should register a parallel cluster before executing the function. See Details.
}
  \item{verbose}{
Logical: should the algorithm progress be displayed during the optimization? Default = \code{FALSE}.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with components
\item{beta}{regression vectors in columns}
\item{objective}{best objective value} 
\item{changePoints}{detected change points}
\item{trace}{objective values in all iterations}
\item{g}{subgradient at last iteration}
\item{iters}{number of iterations performed}
\item{lambda1}{input \code{lambda1}}
\item{lambda2}{input \code{lambda2}}
\item{alpha}{input \code{alpha}}
\item{time}{timing information for all optimization steps}
}

\references{
Degras (2020). Sparse group fused lasso for model segmentation: a hybrid approach. \emph{Advances in Data Analysis and Classification.} \url{https://doi.org/10.1007/s11634-020-00424-5}
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SGFL.SVAR}}
}
\examples{
m <- 5   # number of (time-varying) observation variables
n <- 100 # time series length
p <- 10  # number of (time-varying) predictor variables 

## Array of predictors, regression matrix, and response matrix
x <- array(rnorm(m*n*p),c(m,p,n)) # random coefficients
beta <- matrix(0,p,n) # random coefficients with 50\% sparsity  
beta[1:5,1:50] <- runif(5) # and one change point at time 51
beta[6:p,51:n] <- runif(5)
y <- matrix(,m,n)
for (i in 1:n)
  y[,i] <- x[,,i] \%*\% beta[,i]) + rnorm(m,sd=.01)

lambda1 <- .05
lambda2 <- 1.5
alpha <- .9

## Not run
# result <- SGFL(x, y, lambda1, lambda2, alpha, 
#	intercept=FALSE, verbose=TRUE)

}
