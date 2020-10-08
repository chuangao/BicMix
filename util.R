
run_sim <- function(param.config, itr=2001, inputDir=NULL, outputDir = NULL, script.path = NULL, bfrm.path=NULL, matlab.path = NULL, nfs=20, nf=30, ng= 1000, ns=500, min_fac=20,step=2, mc.cores=5){
    results<- mclapply(((1:nrow(param.config))),function(i){
        #cat(i,"\n")
        #write.table(i,file.path(output.dir,"itr.txt"))
        
        nfs <- nfs
        nf <- nf

        nfd <- nf - nfs

        ng <- ng
        ns <- ns
        
        n.effects <- param.config[i,"n.effects"]
        std.err <- param.config[i,"std.err"]
        i.seed <- param.config[i,"seed"]
        method <- as.character(param.config[i,"method"])
        std.effect <- param.config[i,"std.effect"]
        b <- param.config[i,"b"]
        dense <- as.logical(as.character(param.config[i,"dense"]))

        param <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed, method = method, std.effect = std.effect, b=b, dense=dense,stringsAsFactors =F)
        
        print(i)
        print(param)

        data=gen_SFA_data(std=2, rsd = i.seed, std.err=std.err, n.effects = nf, nfs=nfs, ng=ng,ns=ns,dense=dense)
                  #function(std=2, rsd = 123, std.err=1, n.effects = 30, nfs = 20, ng = 1000, ns=200, dense=T)
        data.file.name <- paste(paste0(names(param)[names(param) != "method"],param[names(param) != "method"]),collapse="_")
        data.file.name=paste0(data.file.name,".txt")
        write.table(data$y,file.path(inputDir,data.file.name),col.names=F,row.names=F)

        #param2 <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed, std.effect = std.effect, dense=dense)
        
        #file.name <-  paste(paste0(names(param2),param2[1,]),collapse="_")    
     
        #c=a;g=a;
        #d=round(b/2);h=round(b/2);

        #b = 10000000
        lams <- NULL
        lamd <- NULL
        z <- NULL
        lam <- NULL
        if(method == "SFAmix"){
            outputDir.sfamix=file.path(outputDir,"sfamix")
            dir.create(outputDir.sfamix)
            
            results <- BicMixR(y=data$y, nf = 50, out_dir=outputDir.sfamix, rsd = 123, tol=1e-10, itr = 2001)
            lam <- results$lam
            z <- results$z
            lams <- lam[,z>0.8,drop=F]
            lamd <- lam[,z<=0.8,drop=F]
        }else if(method == "SPCA"){
            if(dense){
                sp<- SPC(t(data$y),sumabsv=4,K=nf,niter=100)
                lamd <- sp$v[,1:nfd,drop=F]
                lams <- sp$v[,(nfd+1):nf,drop=F]
            }else{
                sp <- SPC(t(data$y),sumabsv=4,K=nfs,niter=100)
                lams <- sp$v[,1:(nfs),drop=F]
            }
            
        }else if(method == "KSVD"){
            outputDir.ksvd=file.path(outputDir,"KSVD")
            dir.create(outputDir.ksvd)
            res <- c()
            if(dense){
                res <- ksvd(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.ksvd,scriptFile=file.path(script.path, 'KSVD/mine_sim.m'),matlabWhere=matlab.path)
                ksvd(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.ksvd,scriptFile=file.path(script.path, 'KSVD/mine_sim.m'),matlabWhere=matlab.path)
                #inputDir=NULL,inputFileName=NULL,n=NA,p=NA,k=NA, outputDir=NULL,scriptFile=NULL,matlabWhere
            }else{
                res <- ksvd(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nfs,outputDir=outputDir.ksvd,scriptFile=file.path(script.path, 'KSVD/mine_sim.m'),matlabWhere=matlab.path)
            }
            
            lam <- res$lam
            index.lam <- order(apply(lam,2,var))

            lams <- lam[,index.lam[1:(nfs)],drop=F]
            if(dense){  
                lamd <- lam[,index.lam[(nfs+1):nf],drop=F]
            }
            
        }else if(method == "IFA"){ 
            outputDir.ifa=file.path(outputDir,'IFA')
            dir.create(outputDir.ifa)
            
            res <- c()
            if(dense){
                res <- ifa(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.ifa,scriptFile=file.path(script.path, 'IFA/IFA_chuan.m'),matlabWhere=matlab.path)
            }else{
                res <- ifa(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nfs,outputDir=outputDir.ifa,scriptFile=file.path(script.path,'IFA/IFA_chuan.m'),matlabWhere=matlab.path)
            }
            
            lam <- res$lam
            lam <- lam[,ncol(lam):1,drop=F]
            count <- apply(lam,2,function(x){return(sum(x!=0))})
            lam <- lam[,count>0,drop=F]
            ncol <- ncol(lam)
            if(ncol <= nfs){
                lams <- lam[,1:ncol,drop=F]
            }else{
                lams <- lam[,1:nfs,drop=F]
            }
            if(dense){  
                if(ncol > nfs){
                    lamd <- lam[,(nfs+1):ncol,drop=FALSE]
                }else{
                    lamd <- NULL
                }
            }
        }else if(method == "BFRM"){  
            outputDir.bfrm=file.path(outputDir, 'BFRM')
            dir.create(outputDir.bfrm)
            
            #system(paste0("mkdir -p ",outputDir))
            if(dense){
                res <- bfrm(inputDir=inputDir, inputFileName=data.file.name,n=ns,p=ng,k=nf, outputDir=outputDir.bfrm,scriptFile=file.path(script.path, 'BFRM/parameters.txt'),bfrmWhere=bfrm.path)
            }else{
                res <- bfrm(inputDir=inputDir, inputFileName=data.file.name,n=ns,p=ng,k=nfs, outputDir=outputDir.bfrm,scriptFile=file.path(script.path, 'BFRM/parameters.txt'),bfrmWhere=brfm.path)
            }
            
            lam <- res$lam
            index.lam <- order(apply(lam,2,var))
    
            lams <- lam[,index.lam[1:(nfs)],drop=F]
            if(dense){  
                lamd <- lam[,index.lam[(nfs+1):nf],drop=F]
            }

        }
       
        nfs.o <- ifelse(is.null(lams),NA,ncol(lams))
        nfd.o <- ifelse(is.null(lamd),NA,ncol(lamd))

        score.sparse <- NA
        score.sparse.precis <- NA

        if(length(lams) == 0){
            score.sparse <- NA
        }else{
            corr.sparse <- cor(lams,data$lams)
            score.sparse <- cal_score_sparse(corr.sparse)
            score.sparse.precis <- cal_score_sparse(corr.sparse,precis=TRUE)
        }
        score.sparse.all <- data.frame(param,score.sparse = score.sparse,score.sparse.precis = score.sparse.precis, nfs=nfs.o,nfd=nfd.o)

        score.dense <- NA
        score.dense.precis <- NA

        if(dense){
            if(length(lamd) == 0){
                score.dense <- NA
            }else{
                score.dense <- cal_score_dense(lamd,data$lamd)
                score.dense.precis <- cal_score_dense(lamd,data$lamd, precis=TRUE)

            }
        }
        score.dense.all <- data.frame(param,score.dense = score.dense,score.dense.precis = score.dense.precis, nfs=nfs.o,nfd=nfd.o)
       
        if(dense){
            return(list(score.sparse=score.sparse.all, score.dense=score.dense.all))
        }else{
            return(list(score.sparse=score.sparse.all, score.dense=data.frame(param,score.dense = NA,score.dense.precis = NA, nfs=nfs.o,nfd=nfd.o)))
        }
        
        
    },mc.cores = mc.cores)
    return(results)
}


cal_score_sparse <- function(sigma, precis = FALSE){

    sigma <- abs(sigma)
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

gen_bulk_data <- function(param.config, output.dir = "data"){
    mclapply(1:nrow(param.config),function(i){
        #cat(i,"\n")
        #write.table(i,file.path(output.dir,"itr.txt"))
        
        n.effects <- param.config[i,"n.effects"]
        std.err <- param.config[i,"std.err"]
        i.seed <- param.config[i,"seed"]
        method <- param.config[i,"method"]
        std.effect <- param.config[i,"std.effect"]
        b <- param.config[i,"b"]
        dense <- as.logical(as.character(param.config[i,"dense"]))

        param <- data.frame(n.effects=n.effects, std.err=std.err,seed=i.seed, std.effect = std.effect,dense=dense)
        
        data=gen_SFA_data(std=2, rsd = i.seed, std.err=std.err, n.effects = 30, nfs = 20, ng=1000,ns=500)
        file.name <- paste(paste0(names(param),param[1,]),collapse="_")
        
        write.table(data$y,file.path(output.dir,paste0(file.name,".txt")),col.names=F,row.names=F)
 },mc.cores = 10)
}

extract_res <- function(results){
    
    ##### extract betas
    score.sparse <- do.call(rbind,lapply(results,function(x){
        return(x$score.sparse)     
    }))

    ##### extract predicted values
    score.dense <- do.call(rbind,lapply(results,function(x){
        return(x$score.dense)     
    }))

    # lams.recover <- do.call(rbind,lapply(results,function(x){
    #     return(x$lams.recover)     
    # }))

    # lams.true <- do.call(rbind,lapply(results,function(x){
    #     return(x$lams.true)     
    # }))

    # return(list(score.sparse = score.sparse, score.dense = score.dense, lams.recover = lams.recover, lams.true = lams.true))
    return(list(score.sparse = score.sparse, score.dense = score.dense))
}


cal_roc <- function(df){
    obj <- roc(df$prob.true,df$prob)
    return(data.frame(fp=1-obj$specificities,tp=obj$sensitivities))
}

cal_quant <- function(df){
    #df <- df[order(df$fp),]
    #print(df[1,])
    prob=seq(0,1,0.01)
    fp.map.to.prob <- cut(df$fp,prob,include.lowest=TRUE)
    index.non.dup <- !duplicated(fp.map.to.prob)
    tp.map.to.prob <- df$tp[index.non.dup]
    fp.tmp = fp.map.to.prob[index.non.dup]

    d = data.frame(
        fp = prob[fp.tmp],
        tp = tp.map.to.prob
    )
    return(d)
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
#' @export


gen_SFA_data <- function(std=2, rsd = 123, std.err=1, n.effects = 30, nfs = 20, ng = 1000, ns=200, dense=T){
  set.seed(rsd)
  
  nf <- n.effects
  nfd <- nf - nfs
  
  lam <- c()
  lams <- c()
  lamd <- c()
  
  lams <- matrix(0,nrow=ng,ncol=nfs)
  
  ########## simulate lam
  for(i in 1:nfs){
    block <- sample(20:60,1)
    start <- sample(1:(ng-block),1)
    lams[start + sample(1:block,block,replace=T),i] = rnorm(block,0,std)
  }
  
  if(dense){
    lamd <- matrix(rnorm(ng*nfd,0,std),nrow=ng,ncol=nfd)
    lam <- cbind(lams,lamd)
    #lam <- lam[,sample(1:nf,nf,replace=F)]
  }else{
    lam <- lams
    #lam <- lam[,sample(1:nfs,nfs,replace=F)]
  }
  
  ############### simulate ex
  if(!dense){
    nf <- nfs
  }
  
  ex <- matrix(rnorm(ns*nf,0,std),nrow=nf,ncol=ns)
  
  ###################### simulate error
  err <- matrix(rnorm(ng*ns,0,std.err),nrow=ng,ncol=ns)
  
  y <- lam %*% ex + err
  
  if(dense){
    return(list(y=y,lams=lams,lamd=lamd,lam=lam,ex=ex, err=err))
  }else{
    return(list(y=y,lams=lams,lam=lam,ex=ex, err=err))
  }
}


mine_heatmap <- function(x,title=NA){
    m <- melt(x)
    p <- ggplot(m,aes(x=Var1,y=Var2)) + geom_tile(aes(fill=value)) + coord_flip() +theme_minimal() + theme(legend.position="none",axis.title = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),plot.title = element_text(face = "bold", hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) + scale_fill_gradient2(low="blue",mid="white",high="red") 
    #scale_fill_gradient2(low="green",mid="white",high="red") 
    #panel.border = element_rect(colour = "black", fill=NA, size=1)
}

sfa_scheme <- function(){
    #data=gen_SFA_data2(std=2, rsd = 123, std.err=1, n.effects = 5, nfs=3, ng=10,ns=5,dense=TRUE)
    data=gen_SFA_data3(std=0.5, rsd = 123, std.err=1, n.effects = 15, nfs=10, ng=50,ns=20,dense=TRUE)

    plam <- mine_heatmap(data$lam,"Loading")
    pex <- mine_heatmap(data$ex,"Factor")
    py <- mine_heatmap(data$y,"Y")
    perr <- mine_heatmap(data$err,"Error")
    return(list(plam = plam, pex = pex, py = py, perr = perr))
}
