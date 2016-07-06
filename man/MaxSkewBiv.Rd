\name{MaxSkewBiv}
\alias{MaxSkewBiv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{MaxSkewBiv: skewness-based projection pursuit for bivariate data
}
\description{Finds Orthogonal Data Projections with Maximal Skewness for Bivariate Random Vectors
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
MaxSkewBiv(x,y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{it is a numerical variable}
  \item{y}{it is a numerical variable
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{.projectionBIV }{Vector of projected data when the original data are bivariate. The user can obtain it by writing ".projectionBIV", an he can obtain a scatterplot of the projeciton by writing plot(.projectionBIV).}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{de LathauwerL., de Moor B.and Vandewalle J. (2000). Onthebestrank-1andrank-(R_{1},R_{2},...R_{N}) approximation of high-order tensors. SIAM Jour. Matrix Ana. Appl. 21, 1324-1342.

Loperfido, N. (2010). Canonical Transformations of Skew-Normal Variates. Test 19, 146-165.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics and Probability Letters 83, 93-99.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179.
%% ~put references to the literature/web site here ~
}
\author{Cinzia Franceschini and Nicola Loperfido

}



