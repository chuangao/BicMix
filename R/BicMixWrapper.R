BicMix <- function(Y_TMP_param,nrow_param, ncol_param, a_param,b_param, c_param, d_param, e_param, f_param, nf_param, itr_param, LAM_out, EX_out, Z_out, O_out, EXX_out, PSI_out, nf_out, out_itr, out_dir,rsd, lam_method, x_method, tol, nf_min){

    Y_TMP_param <- as.numeric(as.character(Y_TMP_param))   
    LAM_out <- rep(0,nrow_param*nf_param)
    EX_out <- rep(0,nf_param*ncol_param)
    EXX_out <- rep(0,nf_param*nf_param)
    PSI_out <- rep(0,nrow_param)
    Z_out <- rep(0,nf_param)
    O_out <- rep(0,nf_param)
    nf_out <- rep(0,1)
    
    result <- .C ("BicMix",
                  as.double(Y_TMP_param),as.integer(nrow_param), as.integer(ncol_param), as.double(a_param),as.double(b_param), as.double(c_param),as.double(d_param),as.double(e_param),as.double(f_param), as.integer(nf_param), as.integer(itr_param), LAM=as.double(LAM_out), EX=as.double(EX_out), Z=as.double(Z_out), O=as.double(O_out),EXX=as.double(EXX_out),PSI=as.double(PSI_out), nf=as.integer(nf_out), as.integer(out_itr), as.character(out_dir),as.integer(rsd),as.character(lam_method), as.character(x_method), as.double(tol),as.integer(nf_min), PACKAGE="BicMix")

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
    
    cat("finished running c\n")
    #return(result)
    
    return(list(lam=LAM,ex=EX,z=Z,o=O,exx=EXX,psi=PSI,nf=nf,rsd=rsd))
    
    
    ##return(list(LAM=result$LAM_out,EX=result$EX_out))
}

#' An algorithm for decomposing a high dimensional matrix into the product of a sparse loading matrix, and a sparse factor matrix.

#' @author Chuan Gao <chuan.gao.cornell@@gmail.com>

#' @param y matrix to be decomposed, no missing values are allowed, each row is a feature and each column is a sample
#' @param nf the number of factors to start with, default to 100
#' @param a first shape parameter for the tpb distribution at the local level of the hierarchy, default to 0.5 to recapitulate horseshoe
#' @param b second shape parameter for the tpb distribution at the local level of the hierarchy, default to 0.5 to recapitulate horseshoe
#' @param c first shape parameter for the tpb distribution at the component specific level of the hierarchy, default to 0.5 to recapitulate horseshoe
#' @param d second shape parameter for the tpb distribution at the component specific level of the hierarchy, default to 0.5 to recapitulate horseshoe
#' @param e first shape parameter for the tpb distribution at the gloabal level of the hierarchy, default to 0.5 to recapitulate horseshoe
#' @param f second shape parameter for the tpb distribution at the global level of the hierarchy, default to 0.5 to recapitulate horseshoe
#' @param itr the maximum number of iterations the algorithm is allowed to run, default to 5000
#' @param out_itr iteration number at which the algorithm writes temporary output, default to 200
#' @param out_dir directory where the algorithm will write temporary output
#' @param rsd random seed used for initializing the parameter values, default to a randomly drawn number
#' @param lam_method the method used to update the loading matrix, take values either "matrix" or "element".  if "matrix", then all component are updated simultaneously (slower but more stable, don't need as many iterations to converge); if "element", each component is updated sequentially (faster but could beless stable, and need more iterations to converge), default to "matrix"
#' @param x_method whether induce sparsity on the X matrix, take values either "sparse" or "dense". default to "dense"
#' @param tol tolerance threshold for convergence, default to 1e-5
#' @param qnorm whether to quantile-normalize the gene expression matrix, default to TRUE

#' @return lam: the sparse loading matrix
#' @return ex: the factor matrix
#' @return z: a vector indicating whether the loadings are sparse (1 indicate sparse)
#' @return o: a vector indicating whether the factors are sparse (1 indicate sparse)
#' @return nf: the number of factors learned by the model
#' @return exx: the expected value of the covariance matrix, E(XX^T)

#' @examples
#' library(BicMix)
#' ## The following is an example on how to use BicMix to obtain biclusters (sparsity is induced on both the loading and factor matrices)
#' ## simulate data
#' data = gen_BicMix_data(std=2, type.factor="mixture",rsd=123)
#' ## Visualize the loading matrix
#' image(t(data$lam),x=1:ncol(data$lam),y=1:nrow(data$lam),xlab="Loadings",ylab="Samples")
#' ## Visualize the factor matrix
#' image(t(data$ex),x=1:ncol(data$ex),y=1:nrow(data$ex),xlab="Samples",ylab="Factors")

#' ## run BicMix on the simulated data
#' dir.create("results")
#' result = BicMixR(data$y,nf=100,out_dir="results",tol=1e-10,x_method="sparse",rsd=123)

#' ## calculate a correlation matrix of the estimated loading matrix 
#' ## and the true loading matrix. Ideally, there should be one and 
#' ## only one big correlation value for a given row and column of the 
#' ## correlation matrix
#' cor.est.real = cor(result$lam[,result$z==1],data$lams)
#' ## visualize the correlation matrix
#' image(cor.est.real,x=1:nrow(cor.est.real),y=1:ncol(cor.est.real),
#' xlab="Recovered loadings",ylab="True loadings")
#' ## calculate similarity score of the recovered loading matrix and the true loading matrix
#' cal_score_sparse(result$lam[,result$z>0.9],data$lams)

#' ## The following is an example on how to use BicMix to obtain clustering of the genes, while keeping the X matrix sparse (sparsity is induced on the loading matrix, not on the factor matrix)
#' ## simulate data
#' data = gen_BicMix_data(std=2, type.factor="dense", rsd = 123)
#' ## perform analysis
#' result = BicMixR(data$y,nf=100,out_dir="results",tol=1e-10,x_method="dense",rsd=123)
#' ## calculate similarity score of the recovered loading matrix and the true loading matrix
#' cal_score_sparse(result$lam[,result$z>0.9],data$lams)

#' @references \url{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004791}

#' @export BicMixR
BicMixR <- function(y=NULL,nf=100,a=0.5,b=0.5,c=0.5,d=0.5, e=0.5,f=0.5, itr=5001,rsd=NULL,out_itr=200,out_dir=NULL, lam_method="matrix", x_method="dense", tol=1e-10, qnorm = TRUE, nf_min = 1){
    
	#if(! "preprocessCore" %in% rownames(installed.packages())){
	#	source("https://bioconductor.org/biocLite.R")
   	#	biocLite("preprocessCore")
	#}
	#library("preprocessCore")
   if(is.null(y)){
        stop("Please check documentation for the correct usage of this methods, you are missing an input matrix!")
    }
	
    if(qnorm){
	#y <- normalize.quantiles(y,copy=TRUE)
        y <- t(apply(y,1,function(i){return(qqnorm(i,plot=F)$x)}))
	#y <- t(apply(y,1,function(i){return(scale(i))}))
    }
    
    if(is.null(rsd)){
        rsd = sample(1:1000000,1)
        warning("You didn't specify a random seed, one is randomly generated.\n")
    }
    
    if(x_method != "sparse" && x_method != "dense"){
        warning("x_method can only be sparse or dense\n")
        stop()
    }

    if(lam_method != "matrix" && lam_method != "element"){
        ##cat("lam_method can only be matrix or element\n")
        stop("lam_method can only be matrix or element\n")
    }

    if(lam_method == "element" && x_method == "sparse"){
        warning("lam_method can't be element when x_method is sparse\n")
        stop()
    }

    out_dir2 = out_dir
    out_dir2 = gsub("/","%",out_dir2)
 
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

    result <- BicMix(y,sn,dy,a,b,c, d, e, f, nf,itr,LAM_out,EX_out,Z_out,O_out,EXX_out, PSI_out, nf_out, out_itr, out_dir2, rsd, lam_method, x_method, tol, nf_min)

    cat("results returned from c\n")
    return(result)
}

gen_sparse_matrix <- function(std=2, rsd = 123, nf = 10, ng = 500){
  set.seed(rsd)
  
  lam <- matrix(0,nrow=ng,ncol=nf)
  
  ########## simulate lam
  for(i in 1:nf){
    block.size <- sample(min(20,ng/2):min(50,ng),1)
    index <- sample(1:ng,block.size,replace=F)
    lam[index,i] = rnorm(length(index),0,std)
  }
  #index <- unlist(apply(lam,2,function(x){return(which(x!=0))}))
  lam
}


gen_dense_matrix <- function(std=2, rsd = 123, nf = 5, ng = 500){
  set.seed(rsd)
  lam <- matrix(0,nrow=ng,ncol=nf)
  lam <- matrix(rnorm(ng*nf,0,std),nrow=ng,ncol=nf)
}

#' Simulate a matrix y = lam * ex + err

#' @param std standard deviation for the normal distribution of the non-zero entries of the sparse loading
#' @param rsd random seed
#' @param std.err standard deviation of the error term
#' @param nf total number of the loadings
#' @param nfs number of the sparse loadings
#' @param ng number of genes
#' @param ns number of samples
#' @param type.loading type of loading matrix wanted, "mixture" for a matrix mixed with sparse and dense components, "sparse" for a matrix with only sparse components
#' @param type.factor type of factor matrix wanted, "mixture" for a matrix mixed with sparse and dense components, "dense" for a matrix with only dense components

#' @return a list containing the following
#' @return lams: the sparse loadings
#' @return lamd: the dense loadings
#' @return lam: the loading matrix combining both the sparse and dense loading

#' @return ex: the factors matrix combining both the sparse and dense factors
#' @return err: matrix of the error term
#' @return y: the y matrix calculated as y = lam * ex + err
#' @export

gen_BicMix_data <- function(std=2, rsd = NULL, std.err=1, nf = 15, nfs = 10, ng = 500, ns=200, type.loading="mixture",type.factor="dense"){
  
  if(is.null(rsd)){
    rsd = sample(1:1000000,1)
    warning("You didn't specify a random seed, one is randomly generated.\n")
  }
  
  set.seed(rsd)
  
  lam <- NULL
  lams <- NULL
  lamd <- NULL
  
  ex <- NULL
  
  if(type.loading == "sparse"){
    lam <- gen_sparse_matrix(std=std, rsd = rsd, nf = nf, ng = ng)
    lams <- lam
  }else if(type.loading=="mixture"){
    lams <- gen_sparse_matrix(std=std, rsd = rsd, nf = nfs, ng = ng)
    lamd <- gen_dense_matrix(std=std, rsd = rsd, nf = nf - nfs, ng = ng)
    lam <- cbind(lamd,lams)
  }else{stop("type.loading can only be sparse, or mixture")}
  
  if(type.factor == "mixture"){
    xs <- gen_sparse_matrix(std=std, rsd = rsd, nf = nfs, ng = ns)
    xd <- gen_dense_matrix(std=std, rsd = rsd, nf = nf - nfs, ng = ns)
    ex <- rbind(t(xs),t(xd))
    ex <- ex[sample(1:nrow(ex),nrow(ex),replace=F),]
  }else if(type.factor == "dense"){
    ex <- t(gen_dense_matrix(std=std, rsd = rsd, nf = nf, ng = ns))
  }else{stop("type.factor can only be dense, or mixture")}
  ###################### simulate error
  err <- matrix(rnorm(ng*ns,0,std.err),nrow=ng,ncol=ns)
  
  y <- as.matrix(lam) %*% as.matrix(ex) + err
  
  return(list(y=y,lams=lams,lamd=lamd,lam=lam,ex=ex, err=err))
}


#' calculate the similarity score of two sparse matrices, explained in the SFAmix paper

#' @param m1 sparse matrix 1
#' @param m2 sparse matrix 2
#' @param precis when TRUE, only use the most correlated components in the two
#' matrices to calculate the score, this is to account for the scenario where
#' the two matrices have different dimensions, but all components in the
#' smaller matrix have high correlation with a subset of big matrix
#' @return a similarity score
#' @export cal_score_sparse
cal_score_sparse <- function(m1, m2, precis = FALSE){
  
  sigma <- abs(cor(m1,m2))
  nr <- nrow(sigma)
  nc <- ncol(sigma)
  
  r1 <- apply(sigma,1,function(x){
    xm <- max(x)
    x <- x[x!=max(x)]
    return(xm - mean(x[x>=mean(x)]))
  })
  r1 <- as.numeric(r1)
  
  r2 <- 0
  if(nrow(sigma) > 1){
    r2 <- apply(sigma,2,function(x){
      xm <- max(x)
      x <- x[x!=max(x)]
      return(xm - mean(x[x>=mean(x)]))  
    })
  }
  r2 <- as.numeric(r2)
  
  nf <- min(nc,nr)
  o1 <- order(r1,decreasing=T)    
  o2 <- order(r2,decreasing=T)
  
  r <- (mean(r1)+mean(r2))/2
  if(precis){
    r <- (mean(r1[o1[1:nf]]) + mean(r2[o2[1:nf]]))/2
  }
  return(r)
}

#' calculate the similarity score of two dense matrices, explained in the SFAmix paper
#' @param m1 dense matrix 1
#' @param m2 dense matrix 2
#' @param precis when TRUE, only use the most correlated components in the two
#' matrices to calculate the score, this is to account for the scenario where
#' the two matrices have different dimensions, but all components in the
#' smaller matrix have high correlation with a subset of big matrix
#' @return a similarity score
#' @export cal_score_dense
cal_score_dense <- function(m1,m2,precis=FALSE){
  
  m1 <- apply(m1,2,function(x){scale(x)})
  m2 <- apply(m2,2,function(x){scale(x)})
  
  if(precis){
    cr <- cor(m1,m2)
    cr <- abs(cr)
    nr <- nrow(cr)
    nc <- ncol(cr)
    nf <- min(nr,nc)
    maxr <- apply(cr,1,max)
    maxc <- apply(cr,2,max)
    o.r <- order(maxr,decreasing=T)
    o.c <- order(maxc,decreasing=T)
    if(nr > nc){
      m1 <- m1[,o.r[1:nc]]
    }else{
      m2 <- m2[,o.c[1:nr]]
    }
  }
  sigma <- m1 %*% t(m1) - m2 %*% t(m2)
  r <- sum(diag(sigma * sigma))/nrow(sigma)/nrow(sigma)
  return(r)
}
