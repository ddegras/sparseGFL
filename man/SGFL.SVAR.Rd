\name{SGFL.SVAR}
\alias{SGFL.SVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit piecewise (Structural) VAR model by the sparse group fused lasso}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SGFL.SVAR(x, p, lambda1, lambda2, alpha=1, struct=FALSE,
   intercept=TRUE, A=NULL, w=NULL, L=NULL, sweep=c("cyclical",
   "srswr","srswor"), pattern=c(10,5), control=list(), 
   parallel=FALSE, verbose=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The data matrix with variables in rows and time in columns.
 (rows=variables, columns=time)}
  \item{p}{The order of the VAR model (see details).}
  \item{lambda1}{The regularization parameter for lasso penalty in objective function.}
  \item{lambda2}{The 
regularization parameter for total variation penalty in objective function.}
\item{alpha}{real number in [0,1]: tradeoff between lasso penalty and ridge regression  penalty. 1 = pure lasso, 0 = pure ridge. Default = 1}
  \item{struct}{
logical: does the model have a structural (lag 0) component? Default = FALSE}
  \item{intercept}{
logical: does the model include an intercept? Default = TRUE}
  \item{A}{
Initial value for model parameters (optional). 
If provided, this should be an array of dimensions (N,N*(p+struct)+intercept,T-p), 
where N and T are the respective numbers of rows and columns of \code{x},
 and struct and intercept are recoded as 0 (FALSE) or 1 (TRUE)  }
  \item{w}{
individual weights for the total variation penalty (optional). Can be specified as a single number or as a vector of length T-p-1 
}
  \item{L}{Vector of length T (optional). Lipschitz constants of the gradients of the functions that make up the loss function}
  \item{sweep}{
Character string specifying the sweep method for the block coordinate descent component of the fitting algorithm (optional). This can be one of "cyclical", "srswr", or "srswor". Default = "cyclical" 
}
  \item{pattern}{
integer of vector of length 2 specifying the maximum number of consecutive epochs spent in block coordinate descent, resp. fusion cycles (optional). Defaults to c(10,5) 
}
  \item{control}{
optional list for fine tuning the optimization algorithm (optional)
}
  \item{ncores}{
number of cpu cores used in optimization (optional). If more than 1, parallel computing will be used. Default = 1
}
  \item{verbose}{
logical: should the optimization progress be displayed? 
}
}
\details{
The time-varying vector autoregressive (VAR) model to be estimated is  
\deqn{x(t) = \mu(t) +  A_{1}(t) x(t-1) + ... + A_{p}(t) x(t-p) + e(t)}
where \eqn{x(t)} is the observed data vector at time \eqn{t}, \eqn{\mu(t)} is an unknown intercept vector, \eqn{A_{0}(t),...,A_{p}(t)} are time-varying transition matrices, and \eqn{e(t)} is a random error vector. The time-varying model parameter \eqn{\theta(t) = (\mu(t),A_{0}(t),...,A_{p}(t))} is assumed to be piecewise constant.A structural component (\eqn{A_{0}(t) x(t)}) can be added to the right-hand side of the equation to get a full time-varying structural VAR (SVAR) model. In this case \eqn{A_{0}(t)} is constrained to have diagonal coefficients equal to zero. 

The intercept \eqn{\mu(t)} is included by default in the model. To remove it, set \code{intercept=FALSE}. Also, by default, the time-varying VAR model is used. Set {struct=TRUE} to estimate a SVAR model.  

The objective function to minimize is the sum of a loss term (sum of squared errors), an elastic net penalty (=lasso + ridge regression penalty), and a total variation penalty. Specifically, the loss term is \deqn{(1/2) \sum [ x(t) - \mu(t) -  A_{0}(t) x(t) - ... - A_{p}(t) x(t-p) ]^2.}
The elastic net penalty is 
\deqn{ \lambda_{1} \alpha \sum ||A_{k}(t) ||_{1} + \lambda_{1} (1-\alpha) (1/2) \sum ||A_{k}(t)||_{F}^2 }
where the sums are over all lags \eqn{k} and times \eqn{t}, \eqn{|| ||_1} is the \eqn{l_1} norm and \eqn{|| ||_F} is the Frobenius norm (Euclidean norm for matrices). The total variation penalty is   
\deqn{\lambda_{2} \sum w_{t} [ || \mu(t+1) - \mu(t) ||_{2}^2 + \sum ||A_{k}(t+1)-A_{k}(t)||_{F}^2 ]^{1/2}} where the first sum ranges over \eqn{t}, the second over \eqn{k}, and \eqn{|| ||_2} is the Euclidean (vector) norm. 
}
\value{
A list with components 
\item{A}{The gitted model parameters given as a 3D array, the third dimension being time. The \eqn{t}-th slice of \code{A} contains the matrix of coefficient estimates \eqn{[A_1(t),...,A_p(t),\mu(t)]}} with \eqn{A_0(t)} to the far left for the SVAR model. 

\item{objective}{The attained minimum of the objective function.} 
\item{changePoints}{The estimated change points range or NULL if no change point. The first change point is at least \eqn{p+2}.}
\item{trace}{The history of objective values at each iteration of the algorithm.}	
\item{subgradient}{The subgradient of minimum norm at solution \code{A}. This matrix should have most of its values very close to zero.}
\item{iters}{The number of iterations/epochs performed in the optimization.}
\item{lambda1}{The regularization parameter \code{lambda1}.}
\item{lambda2}{The regularization parameter \code{lambda2}.}
\item{alpha}{The elastic net mixing parameter \code{alpha}.}
\item{time}{Time elapsed.}
}
\references{
%% ~put references to the literature/web site here ~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ optimize }% use one of  RShowDoc("KEYWORDS")
\keyword{ ts }% __ONLY ONE__ keyword per line
