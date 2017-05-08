\name{MaxSkew-package}
\alias{MaxSkew-package}
%\alias{MaxSkew}
\docType{package}
\title{MaxSkew: skewness-based projection pursuit
%\packageTitle{MaxSkew}
}
\description{Finds Orthogonal Data Projections with Maximal Skewness
%\packageDescription{"MaxSkew"}
}
\details{Package: MaxSkew

Type: Package

Title: Orthogonal Data Projections with Maximal Skewness

Version: 1.1

Date: 2017-05-02

Author: Cinzia Franceschini, Nicola Loperfido

Maintainer: Cinzia Franceschini <cinziafranceschini@msn.com>

Description: It finds Orthogonal Data Projections with Maximal Skewness.
The first data projection in the output is the most skewed among all linear data projections.
The second data projection in the output is the most skewed among all data projections orthogonal to the first one, and so on.

License: GPL-2

%The DESCRIPTION file:
%\packageDESCRIPTION{MaxSkew}
%\packageIndices{MaxSkew}
%%~~ An overview of how to use the package, including the most important functions ~~
}
\author{Cinzia Franceschini and Nicola Loperfido
%\packageAuthor{MaxSkew}

%\Maintainer{Cinzia Franceschini cinziafranceschini@msn.com}
%\packageMaintainer{MaxSkew}
}
\references{de Lathauwer L., de Moor B.and Vandewalle J. (2000). Onthebestrank-1andrank-(R_{1},R_{2},...R_{N}) approximation of high-order tensors. SIAM Jour. Matrix Ana. Appl. 21, 1324-1342.

Loperfido, N. (2010). Canonical Transformations of Skew-Normal Variates. Test 19, 146-165.

Loperfido, N. (2013). Skewness and the Linear Discriminant Function. Statistics and Probability Letters 83, 93-99.

Malkovich, J.F. and Afifi, A.A. (1973). On Tests for Multivariate Normality. J. Amer. Statist. Ass. 68, 176-179
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
\keyword{ package }
%\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
%}
\examples{
## Example 1. Run MaxSkew on the iris data 
data(iris)
iris<-data.matrix(iris) #returns the matrix obtained by converting the data frame to numeric mode
MaxSkew(iris[,1:3],5,2,FALSE) # matrix whose columns are the two projections with maximal skewness
MaxSkew(iris[,1:2],5,1,FALSE) #projection with maximal skewness of the first two variables
#MaxSkewBiv(iris[,1],iris[,2]) #obtains the same of MaxSkew(iris[,1:2],5,1)

## Example 2. Run MaxSkew on the OLYMPIC_DECATHLON_2016 data 
data(OLYMPIC_DECATHLON_2016)
OLYMPIC_DECATHLON_2016_matrix<-data.matrix(OLYMPIC_DECATHLON_2016) #returns a data matrix
MaxSkew(OLYMPIC_DECATHLON_2016_matrix[,4:13],10,2,TRUE) #it returns also the scatterplot
MaxSkew(OLYMPIC_DECATHLON_2016_matrix[,4:13],10,2,FALSE)#as in example 1

OLYMPIC_DECATHLON_2016_projections<-MaxSkew(OLYMPIC_DECATHLON_2016_matrix[,4:13],10,2,FALSE)
plot(OLYMPIC_DECATHLON_2016_projections) #scatterplot of the first two projections
##install.packages("calibrate")
##library(calibrate)
##textxy(OLYMPIC_DECATHLON_2016_projections[,1],OLYMPIC_DECATHLON_2016_projections[,2],
##OLYMPIC_DECATHLON_2016$ATHLETE,offset=0.5)
}
