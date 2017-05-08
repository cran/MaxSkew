\name{.MaxSkewThree}
\alias{.MaxSkewThree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MaxSkew: skewness-based projection pursuit
}
\description{Finds Orthogonal Data Projections with Maximal Skewness for Trivariate Random Vectors
}
\usage{
.MaxSkewThree(data, iterations)
}

\arguments{
  \item{data}{Data matrix where rows and columns represent units and variables.
}
  \item{iterations}{Number of required iterations.
}
}
\details{It is an internal function called by MaxSkew
}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
\references{de Lathauwer L., de Moor B.and Vandewalle J. (2000). Onthebestrank-1andrank-(R_{1},R_{2},...R_{N}) approximation of high-order tensors. SIAM Jour. Matrix Ana. Appl. 21, 1324-1342.

Loperfido, N. (2010). Canonical Transformations of Skew-Normal Variates. Test 19, 146-165.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics and Probability Letters 83, 93-99.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179
%% ~put references to the literature/web site here ~
}
\author{Cinzia Franceschini and Nicola Loperfido
}
%%\note{
%%  ~~further notes~~
%%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%%}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
