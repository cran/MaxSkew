MaxSkewBiv <-function(x,y){
    #Skewness-based projection pursuit for bivariate random vectors.
    #linear=coefficients of the linear projection.
    #projection=vector of projected data.
    #skewdesidered= third standardized moment of the projected data.
    
    data<-cbind(x,y)# We made the data matrix called "dati" by  means of the two variables x and y
   
    n<-length(x)#number of the elements of the vector x
    reale<-matrix(c(0),3)#initialization of the matrix "reale" 
    k<-matrix(c(0),3)#initialization of the k matrix
    
    .Skew<-c()
    .skewdesideredBIV<-c(0)
    .linearBIV<-matrix(c(0),nrow=2)  
    .projectionBIV<-matrix(c(0),nrow=n)
    
    concentrationmatrix<-solve(cov(data)*(n-1)/n)#the inverse of the covariance matrix of dati
    
    SS<-concentrationmatrix
    aut<-eigen(SS)#computes the eigen vectors  and the eigenvalues 
    sort(aut$values,decreasing=FALSE)#sorts the eigenvalues not in an increasing order
    
    VV<-matrix(0,nrow=ncol(data))
    DD<-diag(sort(aut$values,decreasing=FALSE))
    for(i in ncol(aut$vectors):1){
        VV<-cbind(VV,aut$vectors[,i])
    }
    
    VV<-VV[,2:ncol(VV)]
    
    DD.sqrt<-solve(sqrt(DD))
    AA<-VV%*%DD.sqrt%*%t(VV)
    
    Q<-solve(AA)#it's equal to sqrtm()
    
    x.mean <- colMeans(data)#mean vector
    m<-sweep(data, 2, x.mean)#centered data
    Y<-m%*%Q #standardized data
    
    z1<-Y[,1]# first standardized variable
    z2<-Y[,2]#second standardized variable
    m111<-mean(z1*z1*z1) #Expectation of z1 cubed
    m112<-mean(z1*z1*z2) #Expectation of the product of z1 squared and z2 
    m122<-mean(z1*z2*z2) #Expectation of the product z1 and z2 squared
    m222<-mean(z2*z2*z2) #Expectation of z2 cubed
    
    z<-c(-m122,m222-2*m112,2*m122-m111,m112)#cubic function
    v<-polyroot(z)#roots of a cubic function
    v<-matrix(c(v[3],v[2],v[1]))#roots of a cubic function in a reversed order
    
    for(i in 1:3){
        reale[i]<-Re(v[i])# we take the real part of the roots
        h<-rbind(reale[i],1)
        linearp<-(Y%*%h)
        s3<-sqrt(var(linearp)*(n-1)/n)^3 #we calculate the skewness...
        mx<-mean(linearp)
        sk<-sum((linearp-mx)^3)/s3
        M<-sk/n #we have calculated the skewness
        .Skew<<-M
        k[i,1]<-M #skewness of the linear projection
        reale
    }
    
    M<-max(abs(k))#maximum absolute skewness of the projections
    
    for (i in 1:length(reale)){
    
    if (abs(k[i])==M)#if the i-th absolute skewness attains the maximum...
        {   
        .skewdesideredBIV<<-k[i,1]#it is the desired skewness
            .linearBIV<<-Q%*%rbind(reale[i],1)#compute the corresponding lnear function
        .projectionBIV<<-data%*%(Q%*%rbind(reale[i],1))
        
       }
    }
    
  
    if (.skewdesideredBIV<0){ #if skewness is negative...
        .linearBIV<<--.linearBIV #change the sign of the linear function
        .skewdesideredBIV<--.skewdesideredBIV #change the sign of the skewness
        .projectionBIV<<--.projectionBIV #change the sign of the linear projection
        
        }
    

}

