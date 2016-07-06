MaxSkewThree <-function(data,iterations){
    #Data projection with maximal skewness. The skewness of a distribution is 
    #the third standardized moment of the distribution itself.
    #MaxSkewThree calls the function "MaxSkewBiv".
    #skew= vector containing the projections' skewnesses 
    #linear= vector determining the linear projection.
    #projection= vector of projected data.
    #data= data matrix where rows and columns represent units and variables.
    #iterations= number of required iterations
 
  
    #PRELIMINARIES
    n<-nrow(data)#number of rows of the data matrix
    d<-ncol(data)#number of columns of the data matrix
    .Skew<-c()
    .skewdesideredBIV<-c(0)
    
  .projection=NULL
  rm(.projection)
  .linear=NULL
  rm(.linear)
    
    provaconcentrationmatrix<-solve(cov(data)*(n-1)/n)
    provaSS<-provaconcentrationmatrix
    provaaut<-eigen(provaSS)
    sort(provaaut$values,decreasing=FALSE)
    provaVV<-matrix(c(0),nrow=d)
    provaDD<-diag(sort(provaaut$values,decreasing=FALSE))
    
    for(i in ncol(provaaut$vectors):1){
        provaVV<-cbind(provaVV,provaaut$vectors[,i])
    }
    
    provaVV<-provaVV[,2:ncol(provaVV)]
    provaDD.sqrt<-solve(sqrt(provaDD)) 
    provaAA<-provaVV%*%provaDD.sqrt%*%t(provaVV)
    provaQ<-solve(provaAA)
    
    x.mean <- colMeans(data)#vector mean
    m<-sweep(data, 2, x.mean)#centered data
    Z<-m%*%provaQ
    
    A<-matrix(c(0),nrow=d*d,ncol=d)#initialization of the matrix which sums the tensor products
    
    for(i in 1:n){
        A<-A+kronecker(kronecker(Z[i,],t(Z[i,])),Z[i,])#update of the matrix summing the tensor products
    }
    
    #INITIALIZATION
    s<-svd(A/n)#singular value decomposition of the third standardized cumulant	
    c<-matrix(-s$v[,1])#first right singular vector of the third multivariate cumulant
    v<-provaQ%*%c #starting value of "linear"
    .linear<<-v
  
    .projection<<-data%*%.linear #starting value of "projection"
    
    s3<-sqrt(var(.projection)*(n-1)/n)^3 #we compute the skewness of projection...
    mx<-mean(.projection)
    sk<-sum((.projection-mx)^3)/s3
    M<-sk/n# we have computed the skewness of projection
    .Skew<<-M
    
    if (M<0){
        .Skew<<--M #change the sign of the skewness
        v<--v
        .projection<<--.projection
    }

    
    
    #ITERATIONS
    for (i in 1:iterations){
        for (j in 1:d){
            y<-data[,j]#j-th column of the data matrix
            v[j]<-c(0) #removes the j-th variable from the linear combination
            
           MaxSkewBiv(data%*%v,y)
            s<-.skewdesideredBIV #skewness of the new linear combination
            v[j]<-.linearBIV[2]/.linearBIV[1]#updated vector of coefficients
            if (s>M) {#if the new skewness is greater than the current maximum...
                M<-s #...make the new skewness the current maximum
                .Skew<<-s #set the new skewness equal to Skew...
                .linear<<-v#...set the vector of the coefficients equal to linear...
                .projection<<-data%*%v #project the data onto the direction of v
                
            }
            
        }
    }
    .linearBIV<-.linearBIV
  
}
