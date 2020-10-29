
run_sim <- function(data.config, method.config, itr=2001, inputDir=NULL, outputDir = NULL, script.path = NULL, bfrm.path=NULL, matlab.path = NULL, nfs=10, nf=15, ng= 500, ns=300, mc.cores=5){
    results<- do.call(rbind,lapply((1:nrow(data.config)),function(i){
        results<- do.call(rbind,lapply((1:nrow(method.config)),function(j){
            
            nfd <- nf - nfs
            
            std.err <- data.config[i,"std.err"]
            i.seed <- data.config[i,"seed"]
            std.effect <- data.config[i,"std.effect"]
            dense <- as.logical(as.character(data.config[i,"dense"]))
            
            method <- as.character(method.config[i,"method"])
            local <- as.character(method.config[i,"local"])
            component <- as.character(method.config[i,"component"])
            global <- as.character(method.config[i,"global"])
            
            param <- data.frame(nf=nf, std.err=std.err,seed=i.seed, method = method, std.effect = std.effect, 
                                local=local,component=component,global=global, dense=dense,stringsAsFactors =F)
            
            print(i)
            print(param)
            
            type.loading="mixture"
            if(!dense){type.loading="sparse"}
            
            data=gen_BicMix_data(std=2, rsd = i.seed, std.err=std.err, nf = nf, 
                                 nfs=nfs, ng=ng,ns=ns,type.loading=type.loading,type.factor="dense")
            
            data.file.name <- paste(paste0(names(data.config),data.config[i,]),collapse="_")
            data.file.name=paste0(data.file.name,".txt")
            write.table(data$y,file.path(inputDir,data.file.name),col.names=F,row.names=F)
            
            shapes <- data.frame(Horseshoe=c(0.5,0.5), Strawderman=c(1,0.5), Uniform=c(1,1))
          
            
            lams <- NULL
            lamd <- NULL
            z <- NULL
            lam <- NULL
            
            if(method == "SFAmix"){
                outputDir.sfamix=file.path(outputDir,"sfamix")
                dir.create(outputDir.sfamix)
                
                a <- shapes[,method.config[i,"local"]][1]; b <- shapes[,method.config[i,"local"]][2]
                c <- shapes[,method.config[i,"component"]][1]; d <- shapes[,method.config[i,"component"]][2]
                e <- shapes[,method.config[i,"global"]][1]; f <- shapes[,method.config[i,"global"]][2]
                
                results <- BicMixR(y=data$y, nf = 50, a=a,b=b,c=c,d=d,e=e,f=f, out_dir=outputDir.sfamix, 
                                   rsd = 123, tol=1e-10, itr = 2001, lam_method="matrix",x_method="dense")
                lam <- results$lam
                z <- results$z
                lams <- lam[,z>0.8,drop=F]
                lamd <- lam[,z<=0.8,drop=F]
                
            }else if(method == "SBIF"){ 
                outputDir.sbif=file.path(outputDir,'SBIF')
                dir.create(outputDir.sbif)
                res <- sbif(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.sbif,scriptFile=file.path(script.path, 'SBIF/SBIF_chuan.m'),matlabWhere=matlab.path)
                
                lam <- res$lam
                lam <- lam[,ncol(lam):1,drop=F]
                count <- apply(lam,2,function(x){return(sum(x!=0))})
                lam <- lam[,count>0,drop=F]
                lams <- lam[,1:min(ncol(lam),nfs),drop=F]
                
                if(dense){  
                    lamd <- ifelse(ncol(lam) > nfs, lam[,(nfs+1):ncol(lam),drop=FALSE],NULL)
                }
            }else{
                if(method == "SPCA"){
                    lam <- SPC(t(data$y),sumabsv=4,K=nf,niter=100)
                }else if(method == "KSVD"){
                    outputDir.ksvd=file.path(outputDir,"KSVD")
                    dir.create(outputDir.ksvd)
                    res <- ksvd(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.ksvd,scriptFile=file.path(script.path, 'KSVD/mine_sim.m'),matlabWhere=matlab.path)
                    lam <- res$lam
                }else if(method == "BFRM"){  
                    outputDir.bfrm=file.path(outputDir, 'BFRM')
                    dir.create(outputDir.bfrm)
                    res <- bfrm(inputDir=inputDir, inputFileName=data.file.name,n=ns,p=ng,k=nf, outputDir=outputDir.bfrm,scriptFile=file.path(script.path, 'BFRM/parameters.txt'),bfrmWhere=bfrm.path)
                    lam <- res$lam
                }
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
            
            if(length(lams) != 0){
                score.sparse <- cal_score_sparse(lams,data$lams)
                score.sparse.precis <- cal_score_sparse(lams,data$lams, precis=TRUE)
            }
            
            score.dense <- NA
            score.dense.precis <- NA
            
            
            if(length(lamd) != 0){
                score.dense <- cal_score_dense(lamd,data$lamd)
                score.dense.precis <- cal_score_dense(lamd,data$lamd, precis=TRUE)
            }
            
            score.all <- data.frame(param,score.sparse = score.sparse,score.sparse.precis = score.sparse.precis, score.dense = score.dense, score.dense.precis = score.dense.precis, nfs=nfs.o,nfd=nfd.o)
            
            return(score.all)
        }))
        return(results)
    }))
    return(results)
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
