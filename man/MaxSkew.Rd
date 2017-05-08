\name{MaxSkew}
\alias{MaxSkew}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MaxSkew: skewness-based projection pursuit
}
\description{Finds Orthogonal Data Projections with Maximal Skewness
}
\usage{
MaxSkew(data, iterations, components, plot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Data matrix where rows and columns represent units and variables.
%%     ~~Describe \code{data} here~~
}
  \item{iterations}{It is a positive integer
%%     ~~Describe \code{iterations} here~~
}
  \item{components}{Number of orthogonal projections maximizing skewness. It is a positive integer smaller than the number of variables.
%%     ~~Describe \code{components} here~~
}
  \item{plot}{Dichotomous variable: TRUE/FALSE. If plot is set equal to TRUE (FALSE) the scatterplot appears (does not appear) in the output.
%%     ~~Describe \code{plot} here~~
}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{projectionmatrix}{Matrix of projected data. The i-th row represents the i-th unit, while the j-th column represents the j-th projection.}
\item{pairs(projectionmatrix[,2:i],labels=values,main="Projections") }{It is the multiple scatterplot of the projections maximizing skewness.}
\item{.projectionBIV}{Vector of projected data when the original data are bivariate.The user can obtain a scatterplot of the projection by writing plot(.projectionBIV)} 
}
\references{de Lathauwer L., de Moor B.and Vandewalle J. (2000). Onthebestrank-1andrank-(R_{1},R_{2},...R_{N}) approximation of high-order tensors. SIAM Jour. Matrix Ana. Appl. 21, 1324-1342.

Loperfido, N. (2010). Canonical Transformations of Skew-Normal Variates. Test 19, 146-165.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics and Probability Letters 83, 93-99.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179
%% ~put references to the literature/web site here ~
}
\author{Cinzia Franceschini and Nicola Loperfido
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
