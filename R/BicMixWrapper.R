
BicMix <- function(Y_TMP_param,nrow_param, ncol_param, a_param,b_param, nf_param, itr_param, LAM_out, EX_out, Z_out, O_out, EXX_out, PSI_out, nf_out, out_itr, out_dir,rsd, x_method, tol){
    
    Y_TMP_param <- as.numeric(as.character(Y_TMP_param))   
    LAM_out <- rep(0,nrow_param*nf_param)
    EX_out <- rep(0,nf_param*ncol_param)
    EXX_out <- rep(0,nf_param*nf_param)
    PSI_out <- rep(0,nrow_param)
    Z_out <- rep(0,nf_param)
    O_out <- rep(0,nf_param)
    nf_out <- rep(0,1)
    
    result <- .C ("BicMix",
                  as.double(Y_TMP_param),as.integer(nrow_param), as.integer(ncol_param), as.double(a_param),as.double(b_param), as.integer(nf_param), as.integer(itr_param), LAM=as.double(LAM_out), EX=as.double(EX_out), Z=as.double(Z_out), O=as.double(O_out),EXX=as.double(EXX_out),PSI=as.double(PSI_out), nf=as.integer(nf_out), as.integer(out_itr), as.character(out_dir),as.integer(rsd),as.character(x_method), as.double(tol),PACKAGE="BicMix")

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
    
    PSI <- result[["PSI"]][1:nrow_param]
    
    #return(result)
    return(list(lam=LAM,ex=EX,z=Z,o=O,exx=EXX,psi=PSI,nf=nf,rsd=rsd))
    ##return(list(LAM=result$LAM_out,EX=result$EX_out))
}

#' An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a sparse factor matrix.

#' @author Chuan Gao <chuan.gao.cornell@@gmail.com>

#' @param y matrix to be decomposed, no missing values are allowed
#' @param nf the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 100
#' @param a paramater one for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe
#' @param b paramater two for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe
#' @param itr The maximum number of iterations the algorithm is allowed to run, default to 5000
#' @param out_itr Iteration number out_itr, the algorithm will write temporary results into the specified directory (see below) every out_itr number of iterations. default to 500
#' @param out_dir Directory where the algorithm will write temporary results into at the specified iteration number(see above)
#' @param rsd random seed for initializing the parameter values, default to be randomly drawn
#' @param x_method whether induce sparsity on the X matrix, take values either "sparse" or "dense". default to "sparse"
#' @param tol tolerance threshold for convergence, default to 1e-5
#' @param qnorm whether to qq-normalize the gene expression matrix, default to TRUE

#' @return lam: the sparse loading matrix
#' @return ex: the factor matrix
#' @return z: a vector indicating whether the corresponding loading is sparse (value of 1)
#' @return o: a vector indicating whether the corresponding factor is sparse (value of 1)
#' @return nf: the number of factors learned by the model
#' @return exx: the expected value of the covariance matrix, E(XX^T)

#' @examples
#' library(BicMix)
#' ## simulate data, the parameter std specifies the standard error of non-zero entries in the
#' ## loading and factor matrices, where a normal distribution of mean zero
#' ## is assumed for these values.
#' data = gen_BicMix_data(std=2)
#' ## Visualize the loading matrix
#' image(t(data$lam),x=1:ncol(data$lam),y=1:nrow(data$lam),xlab="Loadings",ylab="Samples")
#' ## Visualize the factor matrix
#' image(t(data$ex),x=1:ncol(data$ex),y=1:nrow(data$ex),xlab="Samples",ylab="Factors")

#' ## run algorithm on the simulated data
#' system("mkdir results")
#' result = BicMixR(data$y,nf=100,a=0.5,b=0.5,itr=5000,out_dir="results",tol=1e-5,x_method="sparse",rsd=123)


#' ## calculate a correlation matrix of the estimated loading matrix 
#' ## and the true loading matrix. Ideally, there should be one and 
#' ## only one big correlation value for a given row and column of the 
#' ## correlation matrix if the recovered sparse loadings and the true sparse loadings
#' cor.est.real = cor(result$lam[,result$z==1],data$lams)
#' ## visualize the correlation matrix
#' image(cor.est.real,x=1:nrow(cor.est.real),y=1:ncol(cor.est.real),
#' xlab="Recovered loadings",ylab="True loadings")
#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004791}

BicMixR <- function(y=y,nf=100,a=0.5,b=0.5,itr=5001,rsd=NULL,out_itr=500,out_dir=NULL, x_method=NULL, tol=NULL, qnorm = TRUE){
    
	#if(! "preprocessCore" %in% rownames(installed.packages())){
	#	source("https://bioconductor.org/biocLite.R")
   	#	biocLite("preprocessCore")
	#}
	#library("preprocessCore")

    if(qnorm){
	#y <- normalize.quantiles(y,copy=TRUE)
        y <- t(apply(y,1,function(i){return(qqnorm(i,plot=F)$x)}))
	#y <- t(apply(y,1,function(i){return(scale(i))}))
    }
    
    if(is.null(rsd)){
        rsd = sample(1:1000000,1)
        cat ("You didn't specify a random seed, one is randomly generated.\n")
    }
    
    #set.seed(rsd)
    
    if(is.null(x_method)){
        x_method = "sparse"
    }
    
    if(x_method != "sparse" && x_method != "dense"){
        cat("x_method can only be sparse or dense\n")
        stop()
    }
    
    if(is.null(tol)){
        tol = 1e-5
    }
    out_dir2 = out_dir
    out_dir2 = gsub("/","%",out_dir2)
    if(missing(y)){
        stop("Please check documentation for the correct usage of this methods, you are missing an input matrix!")
    }
    if(is.null(out_dir)){
        out_dir2 = "NULL"
    }else{
        if(!file.exists(out_dir)){
            stop("Your specified output directory is not found!")
        }
    }
    out_dir2 <- gsub("/","",out_dir2)

    sn = nrow(y)
    dy = ncol(y)
    
    LAM_out <- c()
    EX_out <- c()
    EXX_out <- c()
    PSI_out <- c()
    Z_out <- c()
    O_out <- c()
    nf_out <- c()

    result <- BicMix(y,sn,dy,a,b,nf,itr,LAM_out,EX_out,Z_out,O_out,EXX_out, PSI_out, nf_out, out_itr, out_dir2, rsd, x_method, tol)
    return(result)
}

#' Simulate matrix with dimension of 1000 x 200. Number of loadings and factors is set to 30, where 20 loadings and 20 factors are sparse. The sparse loadings and factors cotain mostly zeros, and random blocks of nonzero values generated from N(0,std). The dense loadings and factors are also generated from N(0,std). The error matrix is generated from N(0,1).

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
    nf.s <- 20
    nf.d <- 10

    ng <- 1000
    ns <- 200

    nf <- nf.s+nf.d

    lams <- matrix(0,nrow=ng,ncol=nf.s)
    lamd <- matrix(rnorm(ng*nf.d,0,std),nrow=ng,ncol=nf.d)

    #block <- ng/nf.s
    block <- 50
    for(i in 1:nf.s){
        #ne <- sample(20:40,1)
	start <- sample(1:(ng-50),1)
        lams[start + sample(1:block,block,replace=T),i] = rnorm(block,0,std)
    }
    

    lam <- as.matrix(cbind(lams,lamd))

    lam <- lam[,sample(1:nf,nf,replace=F)]
    
    
    exs <- matrix(0,ncol=ns,nrow=nf.s)
    exd <- matrix(rnorm(ns*nf.d,0,std),nrow=nf.d,ncol=ns)
    #block <- ns/nf.s
	block = 30
    for(i in 1:nf.s){
        #ne <- sample(20:30,1)
	start <- sample(1:(ns-30),1)
        exs[i,start + sample(1:block,block,replace=T)] = rnorm(block,0,std)
    }
    
    ex <- as.matrix(rbind(exs,exd))
    ex <- ex[sample(1:nf,nf,replace=F),]
    
    err <- matrix(rnorm(ng*ns),nrow=ng,ncol=ns)

    y <- lam %*% ex + err
    
    return(list(y=y,lams=lams,lamd=lamd,lam=lam,exs=exs,exd=exd,ex=ex))
}


