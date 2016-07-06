MaxSkew <-
function(data,iterations,components){
  
    #Data projection with maximal skewness. The skewness of a distribution is 
    #the third standardized moment of the distribution itself.
    
    #skew= vector containing the projections' skewnesses 
    #linear= vector determining the linear projection
    #projection= vector of projected data
    #projectionmatrix= matrix of projected data
    #data= data matrix where rows and columns represent units and variables
    #iterations= number of required iterations
    #components= number of orthogonal projections maximizing skewness
    
    #PRELIMINARIES
    
    n<-nrow(data)#number of rows of the data matrix
    d<-ncol(data)#number of columns of the data matrix
    .projectionBIV=NULL
    rm(.projectionBIV)
    
  if(d==2){
           MaxSkewBiv(data[,1],data[,2])
           return(.projectionBIV)
           plot(.projectionBIV)
  } 
   
   #INITIALIZATION
    
    values<-as.list(seq(1:(d-1)))
    
    projectionmatrix<-matrix(c(0),nrow=n,ncol=1)#initialization of the matrix of the projections
    linearmatrix<-matrix(c(0),nrow=1,ncol=1)#initialization of the matrix linearmatrix
    Skewmatrix<-matrix(c(0),nrow=1,ncol=1)#initialization of the matrix skewmatrix
    .Skew<-c() #initialization of the matrix of .Skew
   
    
    for (i in 1:(d-1)){
        h<-d-i+1
        mx<-colMeans(data) # vector mean
        MaxSkewThree(data,iterations)
        projectionmatrix<-cbind(projectionmatrix,.projection)#where projection is taken from MaxSkewThree
        linearmatrix<-rbind(linearmatrix,.linear)
        
       Skewmatrix<-rbind(Skewmatrix,.Skew)
        
        mp<-mean(.projection) #mean of projection
        ssquarep<-c(var(.projection)*(n-1)/n) #variance of projection
        spx<-cov(.projection,data)*(n-1)/n #covariance (projection and data)
        pen<-spx/ssquarep #slope
        intercept<-mx-pen*mp #intercept
        teo<-matrix(1,n,1)%*%intercept+.projection%*%pen #matrix of predicted values
        res<-data-teo #matrix of  the residuals
        covres<-cov(res)*(n-1)/n #covariance of the residuals
        eigen(covres)#spectral decomposition  of covres
        o <- order(eigen(covres)$values, decreasing=FALSE)#we reorder 
        eigen(covres)$values[o]#we reorder
        V<-eigen(covres)$vectors[,o]# spectral decomposition of covres,eigenvector ordered in ascending order
        proiezione<-res%*%V[,2:h]#data projection orthogonal to skewed components
        data<-proiezione
    }
  if(d>2){
          
            for (i in 3:(components+1)){
              dev.new()
        pairs(projectionmatrix[,2:i],labels=values,main="Projections")#scatterplot of the projectionmatrix
    }
    
    projectionmatrix<-projectionmatrix[,2:(components+1)]
    linearmatrix<-as.matrix(linearmatrix[2:sum(seq((components+1):1)),])
    Skewmatrix<-as.matrix(Skewmatrix[2:(components+1)])
    return(projectionmatrix)
    .projection<-.projection
    .linear<-.linear
  }
    
}

