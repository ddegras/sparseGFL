\name{lambda.path.SVAR}
\alias{lambda.path.SVAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Path Selection for Regularization Parameters in Sparse Group Fused Lasso
}
\description{
Select path for regularization parameters in sparse group fused lasso for piecewise SVAR model}
\usage{
lambda.path.SVAR(x, p, w=NULL, struct=FALSE,
   intercept=TRUE, nlambda1=20, lambda1.min.ratio=0.01, 
   lambda1=NULL, nlambda2=20, lambda2.min.ratio=0.01,
   control=list(), parallel=FALSE, warm.start=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{multivariate time series matrix (variables in rows, times in columns)}
  \item{p}{order of the VAR model}
  \item{w}{vector of individual weights for the total variation penalty. Length must be 1 or number of time points - p - 1. Weights default to 1}
  \item{struct}{Logical: does the model include a structural component? Default = FALSE}
\item{intercept}{Logical: does the model include an intercept? Default = TRUE}
\item{nlambda1}{Number of values for the regularization parameter lambda1 
associated with the lasso penalty. Default = 20}
\item{lambda1.min.ratio}{Ratio of largest to smallest value of lambda1. Default = 0.01}
\item{lambda1}{Values of lambda1 (optional)}
\item{nlambda2}{Number of values for the regularization parameter lambda2 associated with the total variation penalty. Default = 20}
\item{lambda2.min.ratio}{Ratio of largest to smallest lambda2. Default = 0.01}
\item{control}{List of optional arguments for the path search for lambda2.}
\item{parallel}{Should calculations be performed in parallel? Default = FALSE. 
To use parallel computing, first register a cluster and set this argument to TRUE}
\item{warm.start}{For sequential computation, should warm starts be used to find lambda2 paths? Default = FALSE. Argument ignored if \code{parallel=TRUE}}
}

\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with \code{nlambda1} components, each containing the fields
\item{lambda1}{A value of lambda1}
\item{lambda2}{Associated path of \code{nlambda2} values for lambda2}
\item{A}{Associated lasso solution}
The components are sorted by decreasing order of lambda1
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

