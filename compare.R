
rm(list=ls())

options(scipen=10)

setwd("/Users/cg253/BicMix")

#source("./TPBayes/fastBVSR-2019v1/R-mac/bvsr.R")
source("./util.R")

set.seed(123)
library(ggplot2)
#library(TPBayes)
library(parallel)
library(glmnet)
library(gridExtra)
library(SSLASSO)
library(dplyr)
library(ggpubr)
library(reshape2)
library("ggsci")
library(BicMix)

#library(BicMix)
#library(BicMix2)
library(PMA)
numCores <- 10

#library(knitr)
#require(markdown)

#setwd("./code/R")

## file directories
results.path <- "./AOAS"

## testing directory creations
table.path <- file.path(results.path, "table")
#stashDirCreate(table.path, otherrx=FALSE, grouprx=TRUE)
plot.path <- file.path(results.path, "plot")
#stashDirCreate(plot.path, otherrx=FALSE, grouprx=TRUE)
data.path <- file.path(results.path, "data")

dir.create(table.path)
dir.create(plot.path)
dir.create(data.path)

#stashDirCreate(results.path, otherrx=FALSE, grouprx=TRUE)

n.effects.list <- c(30)
#std.err.list <- c(1, 2, 3, 4, 5)
std.err.list <- c(1,3,5)

dense.list <- c("TRUE","FALSE")

method.list <- c("SFAmix","IFA","BFRM","SPCA","KSVD")

seed.list <- 1:10
std.effect.list <- 2
#b.list <- c(10) * 100000
b.list <- 0.5
#b.list <- c(500000)
sn <- 500
dy <- 200

ng=sn
ns=dy

nfs = 10
nf = 15

a=0.5
#b=1000000
#c=0.1; d=1000; g=0.1; h=100000


param.config <- expand.grid(n.effects.list, std.err.list, method.list, seed.list, std.effect.list, b.list,dense.list)
names(param.config) <- c("n.effects", "std.err", "method", "seed", "std.effect","b","dense")

source("./KSVD/ksvd.r")
source("./IFA/ifa.r")
source("./BFRM/bfrm.r")
source("./util.R")

#method <- "BicMix2"

script.path = "/Users/cg253/BicMix/"
matlab.path = "/Applications/Matlab_R2019b.app/bin/matlab"
bfrm.path="/Users/cg253/BicMix/BFRM/bfrm"
outputDir = '/Users/cg253//BicMix/AOAS/table'
inputDir= '/Users/cg253/BicMix/AOAS/data'
dir.create(inputDir)


itr <- 2001
#res <- run_sim(param.config[param.config$method == method & param.config$dense == TRUE,][1,], nfs=nfs, nf=nf, ng= ng, ns=ns, itr = itr)
itr = itr
inputDir=inputDir
outputDir
script.path = script.path
bfrm.path=bfrm.path
matlab.path = matlab.path
nfs=nfs
ng=sn
ns=dy
mc.cores = 10

res <- run_sim(param.config[1,], itr = itr, inputDir=inputDir, outputDir = outputDir, script.path = script.path, bfrm.path=bfrm.path, matlab.path = matlab.path, nfs=nfs, nf=nf, ng=sn, ns=dy, mc.cores = 10)
#results <- res[[1]]
#count.prob <- apply(results$z,2,function(x){return(sum(x > 0.5)/sn)})

res2 <- extract_res(res)
#res2.bak <- res2
library("ggsci")

score.sparse.all <- res2$score.sparse
score.dense.all <- res2$score.dense


score.sparse.all$std.err <- factor(score.sparse.all$std.err,levels=unique(score.sparse.all$std.err))
score.dense.all$std.err <- factor(score.dense.all$std.err,levels=unique(score.dense.all$std.err))

score.sparse.all$method <- factor(score.sparse.all$method,levels=unique(score.sparse.all$method))
score.dense.all$method <- factor(score.dense.all$method,levels=unique(score.dense.all$method))

names(score.sparse.all) <- c("n.effects", "Std.err", "seed",  "Method", "Std.effect","b", "dense", "Score")
names(score.dense.all) <- c("n.effects", "Std.err", "seed",  "Method", "Std.effect","b", "dense", "Score")

score.dense.all <- score.dense.all[score.dense.all$dense == TRUE,]

score.sparse.all$Data <- ifelse(score.sparse.all$dense,"Sparse-sim1","Sparse-sim2")
score.dense.all$Data <- ifelse(score.dense.all$dense,"Dense-sim1")

score.all <- rbind(score.sparse.all,score.dense.all)
score.all$Data <- factor(score.all$Data,levels=unique(score.all$Data))


file.csv <- file.path(table.path,paste0("comparison_stdErr1_7_stdEff2_nf",nf,"_nfs",nfs,"_p",500,"_n",200,"_a0.1_bstart",1000000,"_nfGrt_10_bfGrt1000.csv"))
#write.csv(score.all,file.csv,row.names=F)

score.all <- read.csv(file.csv,as.is=T) 

score.all$Std.err <- factor(score.all$Std.err,levels=unique(score.all$Std.err))
score.all$Std.err <- paste0("Sd.err=",score.all$Std.err)
score.all$Method[score.all$Method == "BicMix2"] <- "BARF"
score.all$Method <- factor(score.all$Method,levels=c("BARF","BicMix","BFRM","IFA","KSVD","SPCA"))
score.all$Data <- factor(score.all$Data,levels=c("Sparse-sim1","Sparse-sim2","Dense-sim1"))

#p <- ggplot(score.all, aes(x=Std.err,y=Score, color=Method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + facet_grid(Data ~ . , scale="free_y")

p <- ggplot(score.all, aes(x=Method,y=Score, color=Method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + facet_grid(Data ~ Std.err , scale="free_y")+guides(colour=guide_legend(nrow=1,byrow=TRUE))



file.pdf <- file.path(plot.path,paste0("comparison_stdErr1_7_stdEff2_nf",nf,"_nfs",nfs,"_p",500,"_n",200,"_a0.1_bstart",1000000,"_nfGrt_10_bfGrt1000.pdf"))
pdf(file.pdf,width=8,height=6)
print(p)
dev.off()

########### plot up the score for data with dense components
# p.sparse1 <- ggplot(score.sparse.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + ylim(0,1)
# #ggthemes::geom_tufteboxplot()
# p.sparse2 <- ggplot(score.sparse.all[score.sparse.all$dense == FALSE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none",legend.title=element_blank()) + scale_color_jco() + ylim(0,1)

# lg <- as_ggplot(get_legend(p.sparse1)) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
   
# p.sparse1 <- p.sparse1 + theme(legend.position="none")

# p.dense <- ggplot(score.dense.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.dense, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none") + scale_color_jco() + ylim(0,1)


# #lay <- matrix(c(rep(c(1,2),6),3,3),nrow=2)
# #lay <- matrix(rep(1:4,c(1,4,4,4)),ncol=1)
# lay <- matrix(c(c(1,2,2,2,2),c(1,3,3,3,3),c(1,4,4,4,4)),ncol=3)


# file.png <- file.path(plot.path,paste0("comparison_self_stdErr1_7_stdEff2_p500_n200_a0.1_bstart1000000_nfGrt_sqrt(dy)_bfGrt1000_mixture.png"))
# png(file.png,width=2400,height=1600, res = 300)
# grid.arrange(lg, p.sparse1,p.sparse2,p.dense,layout_matrix = lay) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
# dev.off()

############# plot up results with different b values
p.sparse1 <- ggplot(score.sparse.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="top",legend.title=element_blank()) + scale_color_jco() + ylim(0,1) + facet_grid(. ~ b)
#ggthemes::geom_tufteboxplot()
p.sparse2 <- ggplot(score.sparse.all[score.sparse.all$dense == FALSE,], aes(x=std.err,y=score.sparse, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none",legend.title=element_blank()) + scale_color_jco() + ylim(0,1) + facet_grid(. ~ b)

lg <- as_ggplot(get_legend(p.sparse1)) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
   
p.sparse1 <- p.sparse1 + theme(legend.position="none")

p.dense <- ggplot(score.dense.all[score.sparse.all$dense == TRUE,], aes(x=std.err,y=score.dense, color=method)) + geom_boxplot(width=0.5,outlier.size=0.5) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme(legend.pos="none") + scale_color_jco() + facet_grid(. ~ b)


#lay <- matrix(c(rep(c(1,2),6),3,3),nrow=2)
#lay <- matrix(rep(1:4,c(1,4,4,4)),ncol=1)
#lay <- matrix(c(c(1,2,2,2,2),c(1,3,3,3,3),c(1,4,4,4,4)),ncol=3)
lay <- matrix(c(
    rep(1:3,6),
    rep(4,3)
),nrow=3)


file.png <- file.path(plot.path,paste0("comparison_self_stdErr1_7_stdEff2_p500_n200_a0.1_bstart1-3-5-7-10E5_nfGrt_1_bfGrt1000.png"))
png(file.png,width=2400,height=1600, res = 300)
grid.arrange(p.sparse1,p.sparse2,p.dense,ncol=1) + theme(plot.margin = unit(c(0,0,0,0), "cm"))
dev.off()

