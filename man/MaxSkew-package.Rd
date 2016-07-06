\name{MaxSkew-package}
\alias{MaxSkew-package}
%\alias{MaxSkew}
\docType{package}
\title{MaxSkew: skewness-based projection pursuit
}
\description{
Finds Orthogonal Data Projections with Maximal Skewness
}
\details{
\tabular{ll}{
Package: \tab MaxSkew\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2016-06-17\cr
License: \tab GPL-2\cr
}

}
\author{Cinzia Franceschini and Nicola Loperfido


Maintainer: Cinzia Franceschini cinziafranceschini@msn.com
}
\references{de Lathauwer L.,de Moor B.and Vandewalle, J. (2000). On the best rank-1 and rank-(R_{1},R_{2},...R_{N}) approximation of high-order tensors. SIAM Jour. Matrix Ana. Appl. 21, 1324-1342.

Loperfido, N. (2010). Canonical Transformations of Skew-Normal Variates. Test 19, 146-165.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics and Probability Letters 83, 93-99.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179.
}

\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## Run MaxSkew on the iris data
data(iris)
iris<-data.matrix(iris) #returns the matrix obtained by converting the data frame to numeric mode
MaxSkew(iris[,1:3],5,2) # the matrix whose columns are the two projections with maximal skewness
MaxSkew(iris[,1:2],5,1) #obtains the projection with maximal skewness of the first two variables
MaxSkewBiv(iris[,1],iris[,2]) #obtains the same of MaxSkew(iris[,1:2],5,1)
}
