
BicMix <- function(Y_TMP_param,nrow_param, ncol_param, a_param,b_param, nf_param, itr_param, LAM_out, EX_out, Z_out, O_out, EXX_out, nf_out){
    Y_TMP_param <- as.numeric(as.character(Y_TMP_param))   
    LAM_out <- rep(0,nrow_param*nf_param)
    EX_out <- rep(0,nf_param*ncol_param)
    EXX_out <- rep(0,nf_param*nf_param)
    Z_out <- rep(0,nf_param)
    O_out <- rep(0,nf_param)
    nf_out <- rep(0,1)
    
    result <- .C ("BicMix",
                  as.double(Y_TMP_param),as.integer(nrow_param), as.integer(ncol_param), as.double(a_param),as.double(b_param), as.integer(nf_param), as.integer(itr_param), LAM=as.double(LAM_out), EX=as.double(EX_out), Z=as.double(Z_out), O=as.double(O_out),EXX=as.double(EXX_out), nf=as.integer(nf_out))

    nf <- result[['nf']][1]
    
    LAM=result[['LAM']]
    LAM <- LAM[1:(nrow_param*nf)]
    LAM <- matrix(LAM,nrow=nrow_param,ncol=nf)

    EX=result[['EX']]
    EX <- EX[1:(ncol_param*nf)]
    EX <- matrix(EX,nrow=nf,ncol=ncol_param)
    
    EXX=result[['EXX']]
    EXX <- EXX[1:(nf*nf)]
    EXX <- matrix(EXX,nrow=nf,ncol=nf)
    
    Z <- 1- result[["Z"]][1:nf]
    O <- 1- result[["O"]][1:nf]
    
    #return(result)
    return(list(lam=LAM,ex=EX,z=Z,o=O,exx=EXX,nf=nf))
    ##return(list(LAM=result$LAM_out,EX=result$EX_out))
}

#' An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a sparse factor matrix.

#' @author Chuan Gao <chuan.gao.cornell@@gmail.com>

#' @param y matrix to be decmoposed, no missing values are allowed
#' @param nf the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50
#' @param a paramater one for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe
#' @param b paramater two for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe
#' @param itr The maximum number of iterations the algorithm is allowed to run, default to 500

#' @return lam: the sparse loading matrix
#' @return ex: the factor matrix
#' @return z: a vector indicating whether the corresponding loading is sparse (value of 1)
#' @return o: a vector indicating whether the corresponding factor is sparse (value of 1)
#' @return nf: the number of factors learned by the model
#' @return exx: the expected value of the covarance matrix, E(XX^T)

#' @examples
#' library(BicMix)
#' ## simulate data, the parameter std specifies the standard error of non-zero entries in the
#' ## loading and factor matrices, where a normal distribution of mean zero
#' ## is assumed for these values.
#' data = gen_BicMix_data(std=2)
#' ## Visulize the loading matrix
#' image(t(data$lam),x=1:ncol(data$lam),y=1:nrow(data$lam),xlab="Loadings",ylab="Samples")
#' ## Visulize the factor matrix
#' image(t(data$ex),x=1:ncol(data$ex),y=1:nrow(data$ex),xlab="Samples",ylab="Factors")

#' ## run algorithm on the simulated data
#' result = BicMixR(data$y,nf=50,a=0.5,b=0.5,itr=1000)


#' ## calculate a correlation matrix of the estimated loading matrix 
#' ## and the true loading matrix. Ideally, there should be one and 
#' ## only one big correlation value for a given row and column of the 
#' ## correlation matrix if the recovered sparse loadings and the true sparse loadings
#' cor.est.real = cor(result$lam[,result$z==1],data$lams)
#' ## visulize the correlation matrix
#' image(cor.est.real,x=1:nrow(cor.est.real),y=1:ncol(cor.est.real),
#' xlab="Recovered loadings",ylab="True loadings")
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004791}

BicMixR <- function(y=y,nf=50,a=0.5,b=0.5,itr=500){

    sn = nrow(y)
    dy = ncol(y)
    
    LAM_out <- c()
    EX_out <- c()
    EXX_out <- c()
    Z_out <- c()
    O_out <- c()
    nf_out <- c()

    result <- BicMix(y,sn,dy,a,b,nf,itr,LAM.out,EX.out,Z_out,O_out,EXX_out, nf.out)
    return(result)
}

#' Simulate matrix with dimension of 500 x 20. Number of loadings and factors is set to 15, where 10 loadings and 10 factors are sparse. The sparse loadings and factors cotain mostly zeros, and random blocks of nonzero values generated from N(0,std). The dense loadings and factors are also generated from N(0,std). The error matrix is generated from N(0,1).

#' @param std standard deviation for the normal distribution of the non-zero entries of the sparse components

#' @return a list containing the following
#' @return lams: the sparse loadings
#' @return lamd: the dense loadings
#' @return lam: the loading matrix combining both the sparse and dense loading

#' @return exs: the sparse factors matrix
#' @return exd: the dense factors matrix
#' @return ex: the factors matrix combining both the sparse and dense factors
#' @return y: the y matrix calculated as y = lam * ex + err
 
gen_BicMix_data <- function(std=2){
    nf.s <- 10
    nf.d <- 5

    ng <- 500
    ns <- 200

    nf <- nf.s+nf.d


    lams <- matrix(0,nrow=ng,ncol=nf.s)
    lamd <- matrix(rnorm(ng*nf.d,0,std),nrow=ng,ncol=nf.d)

    block <- ng/nf.s
    for(i in 1:nf.s){
        ne <- sample(20:40,1)
        lams[(i-1)*block + sample(1:block,ne,replace=F),i] = rnorm(ne,0,std)
    }
    

    lam <- as.matrix(cbind(lams,lamd))

    lam <- lam[,sample(1:nf,nf,replace=F)]
    
    
    exs <- matrix(0,ncol=ns,nrow=nf.s)
    exd <- matrix(rnorm(ns*nf.d,0,std),nrow=nf.d,ncol=ns)
    block <- ns/nf.s
    for(i in 1:nf.s){
        ne <- sample(10:20,1)
        exs[i,(i-1)*block + sample(1:block,ne,replace=F)] = rnorm(ne,0,std)
    }
    
    ex <- as.matrix(rbind(exs,exd))
    ex <- ex[sample(1:nf,nf,replace=F),]
    
    err <- matrix(rnorm(ng*ns),nrow=ng,ncol=ns)

    y <- lam %*% ex + err
    
    return(list(y=y,lams=lams,lamd=lamd,lam=lam,exs=exs,exd=exd,ex=ex))
}

#y = gen_BicMix_data(std=1)$y
#BicMixR(y=y,nf=50,a=0.5,b=0.5,itr=500)


