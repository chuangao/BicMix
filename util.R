
run_sim <- function(data.config, method.config, itr=2001, inputDir=NULL, outputDir = NULL, script.path = NULL, bfrm.path=NULL, matlab.path = NULL, nfs=10, nf=15, ng= 500, ns=200, mc.cores=5){
    results<- do.call(rbind,mclapply((1:nrow(data.config)),function(i){
        results<- do.call(rbind,mclapply(1:length(method.config),function(j){
            
            nfd <- nf - nfs
            
            std.err <- data.config[i,"std.err"]
            i.seed <- data.config[i,"seed"]
            std.effect <- data.config[i,"std.effect"]
            dense <- as.logical(as.character(data.config[i,"dense"]))
            
            method <- as.character(method.config[j])
        
            
            param <- data.frame(nf=nf, std.err=std.err,seed=i.seed, method = method, std.effect = std.effect,
                                dense=dense,stringsAsFactors =F)
            
            print(i)
            print(param)
            
            type.loading="mixture"
            if(!dense){type.loading="sparse"}
            
            data=gen_BicMix_data(std=std.effect, rsd = i.seed, std.err=std.err, nf = 15, 
                                 nfs=10, ng=ng,ns=ns,type.loading=type.loading,type.factor="dense")
            
            data.file.name <- paste(paste0(names(data.config),data.config[i,]),collapse="_")
            data.file.name=paste0(data.file.name,".txt")
            write.table(data$y,file.path(inputDir,data.file.name),col.names=F,row.names=F)
          
            
            lams <- NULL
            lamd <- NULL
            z <- NULL
            lam <- NULL
            
            a <- 0.5; b <- 0.5; c <- 0.5; d <- 0.5; e <- 0.5; f <- 0.5
            if(method == "SFAMix"){
                outputDir.sfamix=file.path(outputDir,"sfamix",data.file.name)
                dir.create(outputDir.sfamix, recursive = T)
                
                results <- BicMixR(y=data$y, nf = 50, a=a,b=b,c=c,d=d,e=e,f=f, out_dir=outputDir.sfamix, out_itr = 500, 
                                   rsd = 123, tol=1e-10, itr = itr, lam_method="element",x_method="dense",nf_min=1)
                lam <- results$lam
                z <- results$z
                lams <- lam[,z>0.8,drop=F]
                lamd <- lam[,z<=0.8,drop=F]
       
            }else if(method == "SBIF"){ 
                outputDir.sbif=file.path(outputDir,'SBIF')
                dir.create(outputDir.sbif, recursive = T)
                res <- sbif(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.sbif,scriptFile=file.path(script.path, 'SBIF/SBIF_chuan.m'),matlabWhere=matlab.path)
                
                lam <- res$lam
                #lam <- lam[,ncol(lam):1,drop=F]
                count <- apply(lam,2,function(x){return(sum(x!=0))})
                lam <- lam[,count>0,drop=F]

            }else if(method == "SPCA"){
                outputDir.spca=file.path(outputDir,'SPCA',data.file.name)
                dir.create(outputDir.spca, recursive = T)
                lam <- SPC(t(data$y),sumabsv=4,K=nf,niter=100)
                lam <- lam$v
                write.table(lam,file.path(outputDir.spca,"lam.txt"))
            }else if(method == "KSVD"){
                outputDir.ksvd=file.path(outputDir,"KSVD")
                dir.create(outputDir.ksvd, recursive = T)
                res <- ksvd(inputDir=inputDir,inputFileName=data.file.name,n=ns,p=ng,k=nf,outputDir=outputDir.ksvd,scriptFile=file.path(script.path, 'KSVD/mine_sim.m'),matlabWhere=matlab.path)
                lam <- res$lam
            }else if(method == "BFRM"){  
                outputDir.bfrm=file.path(outputDir, 'BFRM')
                dir.create(outputDir.bfrm, recursive = T)
                res <- bfrm(inputDir=inputDir, inputFileName=data.file.name,n=ns,p=ng,k=nf, outputDir=outputDir.bfrm,scriptFile=file.path(script.path, 'BFRM/parameters.txt'),bfrmWhere=bfrm.path)
                lam <- res$lam
            }
            
            if(method != "SFAMix"){
                index.lam <- order(apply(lam,2,var))
                lam <- lam[,index.lam]
                
                if(dense){  
                    lams <- lam[,1:min(ncol(lam),nfs),drop=F]
                    if(ncol(lam) > nfs){
                      lamd <- lam[,(nfs+1):ncol(lam),drop=FALSE]
                    }
                }else{
                    lams <- lam[,1:min(ncol(lam),nf),drop=F]
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
            
            
            if(length(lamd) != 0 && length(data$lamd) != 0){
                score.dense <- cal_score_dense(lamd,data$lamd)
                score.dense.precis <- cal_score_dense(lamd,data$lamd, precis=TRUE)
            }
            
            score.all <- data.frame(param,score.sparse = score.sparse,score.sparse.precis = score.sparse.precis, score.dense = score.dense, score.dense.precis = score.dense.precis, nfs=nfs.o,nfd=nfd.o)
            
            return(score.all)
        },mc.cores=10))
        return(results)
    }))
    return(results)
}


format_font <- function(p, size=14){
    p <- ggpar(p,
               font.main=size,
               font.x=size,
               font.y=size,
               font.legend=size,
               font.tickslab = size-4)
    p <- p + grids()
}

format_font_ggplot <- function(p, size=16){
    p <- p + theme(axis.text=element_text(size=size), text = element_text(size=size),axis.text.x = element_text(angle=45, hjust=1))
    p <- p + theme(strip.text = element_text(size = size), legend.text=element_text(size=size))
    p
}


gen_schema <- function(){
    data <- gen_BicMix_data(rsd=123, std=1)
    index <- unlist(apply(data$lams,2,function(x){return(which(x!=0))}))
    yS <- data$lams[index,] %*% data$ex[1:ncol(data$lams),]
    corYS <- abs(cor(t(yS)))
    corYSM <- melt(corYS)
    
    yD <- data$lamd[index,] %*% data$ex[1:ncol(data$lamd),]
    corYD <- abs(cor(t(yD)))
    corYDM <- melt(corYD)
    
    corY <- abs(cor(t(data$y[index,])))
    corYM <- melt(corY)
    
    corM <- rbind(cbind(corYSM,type="Cor(sparse)"),cbind(corYDM,type="Cor(dense)"),cbind(corYM,type="Cor(Y)"))
    names(corM) <- c("Gene1","Gene2","Corr","Type")
    corM$Type <- factor(corM$Type,levels=c("Cor(Y)","Cor(sparse)","Cor(dense)"))
    p <- ggplot(corM,aes(x=Gene1,y=Gene2,fill=Corr)) + geom_tile() + facet_grid(.~Type) + theme_minimal()
    p <- format_font_ggplot(p, size=14) + xlab("Genes") + ylab("Genes") + scale_fill_gradient2(low = "white",high = "red")
    p
}

plot_comparison <- function(res, precis=F){
    
    if(precis){
        names(res) <- c("K","Std.err","Seed","Method","Std.effect","Dense","Score.sparse","Score Sparse","Score.dense","Score Dense","Ks","Kd")
    }else{
        names(res) <- c("K","Std.err","Seed","Method","Std.effect","Dense","Score Sparse","Score.sparse.precis","Score Dense","Score.dense.precis","Ks","Kd")
    }
  
    #res <- res[as.character(res$alg) == "element",]
    #names(res) <- c("K","Std.err","Seed","Method","alg","Std.effect","Local","Component","Global","Dense","Score Sparse","Score.sparse.precis","Score Dense","Score.dense.precis","Ks","Kd")
    res <- melt(res,id.vars = c("Std.err","Seed","Method","Dense"), measure.vars = c("Score Sparse","Score Dense"),variable.name = "Score.type", value.name = "Score")
    res <- res[!(!res$Dense & res$Score.type=="Score Dense"),] ## sparse simulations have no dense scores
    #res <- res[res$Dense,]
    
    res$Data <- ifelse(res$Dense,"Sim2","Sim1")
    res$Type <- paste(res$Score.type, res$Data)
    res$Std.err <- paste0("Std_err=",res$Std.err)
    
    res$Type <- factor(res$Type, levels=c("Score Sparse Sim1","Score Sparse Sim2","Score Dense Sim2"))
    
    #res$Score = as.numeric(res$Score)
    p <- ggboxplot(res,x="Method",y="Score",color="Method",palette="jco",add="jitter",shape="Method")
    p <- facet(p,facet.by = c("Type","Std.err"), panel.labs.font = list(size=14),scales = "free_y",)
    p <- format_font(p, size=18) + rotate_x_text(angle=45) 
    p <- ggpar(p, legend.title = NULL) + panel_border() + background_grid()
    
    p
}

# 
# extract_res <- function(results){
#     
#     ##### extract betas
#     score.sparse <- do.call(rbind,lapply(results,function(x){
#         return(x$score.sparse)     
#     }))
# 
#     ##### extract predicted values
#     score.dense <- do.call(rbind,lapply(results,function(x){
#         return(x$score.dense)     
#     }))
# 
#     return(list(score.sparse = score.sparse, score.dense = score.dense))
# }
# 
# cal_roc <- function(df){
#     obj <- roc(df$prob.true,df$prob)
#     return(data.frame(fp=1-obj$specificities,tp=obj$sensitivities))
# }
# 
# cal_quant <- function(df){
#     #df <- df[order(df$fp),]
#     #print(df[1,])
#     prob=seq(0,1,0.01)
#     fp.map.to.prob <- cut(df$fp,prob,include.lowest=TRUE)
#     index.non.dup <- !duplicated(fp.map.to.prob)
#     tp.map.to.prob <- df$tp[index.non.dup]
#     fp.tmp = fp.map.to.prob[index.non.dup]
# 
#     d = data.frame(
#         fp = prob[fp.tmp],
#         tp = tp.map.to.prob
#     )
#     return(d)
# }
# 
# 
# mine_heatmap <- function(x,title=NA){
#     m <- melt(x)
#     p <- ggplot(m,aes(x=Var1,y=Var2)) + geom_tile(aes(fill=value)) + coord_flip() +theme_minimal() + theme(legend.position="none",axis.title = element_blank(),axis.ticks=element_blank(),axis.text=element_blank(),plot.title = element_text(face = "bold", hjust = 0.5),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + ggtitle(title) + scale_fill_gradient2(low="blue",mid="white",high="red") 
#     #scale_fill_gradient2(low="green",mid="white",high="red") 
#     #panel.border = element_rect(colour = "black", fill=NA, size=1)
# }
# 
# sfa_scheme <- function(){
#     #data=gen_SFA_data2(std=2, rsd = 123, std.err=1, n.effects = 5, nfs=3, ng=10,ns=5,dense=TRUE)
#     data=gen_SFA_data3(std=0.5, rsd = 123, std.err=1, n.effects = 15, nfs=10, ng=50,ns=20,dense=TRUE)
# 
#     plam <- mine_heatmap(data$lam,"Loading")
#     pex <- mine_heatmap(data$ex,"Factor")
#     py <- mine_heatmap(data$y,"Y")
#     perr <- mine_heatmap(data$err,"Error")
#     return(list(plam = plam, pex = pex, py = py, perr = perr))
# }

# }else if(method == "SFAMix2"){
#     outputDir.sfamix2=file.path(outputDir,"sfamix2",data.file.name)
#     dir.create(outputDir.sfamix2, recursive = T)
#     
#     prc <- prcomp(data$y, center = TRUE, scale = FALSE)
#     y <- lm(data$y ~ prc$x[,1:5])$residual
#     
#     results <- BicMixR(y=y, nf = 50, a=a,b=b,c=c,d=d,e=e,f=f, out_dir=outputDir.sfamix2, 
#                        rsd = 123, tol=1e-10, itr = itr, lam_method="element",x_method="dense",nf_min=1)
#     lam <- results$lam
#     z <- results$z
#     lams <- lam[,z>0.8,drop=F]
#     lamd <- lam[,z<=0.8,drop=F]
#  
