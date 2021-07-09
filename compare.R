
rm(list=ls())

options(scipen=10)

#install_github("chuangao/BicMix",force=T)

## do this to install BicMix from command line
#git clone https://github.com/chuangao/BicMix 
#R CMD INSTALL BicMix

library(BicMix)
library(ggplot2)
library(parallel)
library(glmnet)
library(gridExtra)
library(SSLASSO)
library(dplyr)
library(ggpubr)
library(reshape2)
library(ggsci)
library(BicMix)
library(PMA)
library(cowplot)
library(devtools)

setwd("~/BicMix")

source("./util.R")
source("./KSVD/ksvd.r")
source("./SBIF/sbif.r")
source("./BFRM/bfrm.r")
source("./util.R")

script.path = "~/BicMix/"
matlab.path = "/Applications/Polyspace/R2020b.app/bin/matlab"
bfrm.path="~/BicMix/BFRM/bfrm"

################################################################### file directories
results.path <- "./AOAS"; dir.create(results.path)

table.path <- file.path(results.path, "table"); dir.create(table.path)
plot.path <- file.path(results.path, "plot"); dir.create(plot.path)
data.path <- file.path(results.path, "data"); dir.create(data.path)

outputDir = table.path
inputDir= data.path

##########################################################################
std.err.list <- c(0.5, 1, 2)
dense.list <- c("TRUE","FALSE") ## whether loading has dense components
seed.list <- 1:10
std.effect.list <- 1

data.config <- expand.grid(std.err.list, seed.list, std.effect.list, dense.list)
names(data.config) <- c("std.err", "seed", "std.effect","dense")

method.list <- c("SFAMix","SBIF","BFRM","SPCA","KSVD")
method.config <- method.list

ng <- 500
ns <- 200

nfs = 10
nf = 15

itr <- 2001

i=1

res <- run_sim(data.config,method.config, itr = itr, inputDir=inputDir, outputDir = outputDir, script.path = script.path, bfrm.path=bfrm.path, matlab.path = matlab.path, nfs=nfs, nf=nf, ng=ng, ns=ns, mc.cores = 10)
res.bak <- res
write.csv(res,file.path(table.path,"res.csv"),row.names=F)


################################### plot comparing to other method, element

res <- read.csv(file.path(table.path,"res.csv"))
#res <- res[(!is.na(res$alg) & as.character(res$alg) == "element") | is.na(res$alg),]
#res$method[res$method=="SFAmix"] <- "SFAMix"
#res$method[!is.na(res$alg)] <- paste0(res$method[!is.na(res$alg)], "_", res$alg[!is.na(res$alg)])
res$method[res$method=="SFAMix"] <- "BicMix"
res$method <- factor(res$method, levels=c("BicMix", "SBIF","BFRM","SPCA","KSVD"))

p.comp <- plot_comparison(res, precis=F)

file.pdf <- file.path(plot.path,"score_compare_methods_element.pdf")
pdf(file.pdf,width=8,height=8)
print(p.comp)
dev.off()


################################################################# 

nf.init = 20
nfs.init = 15

file.name <- paste0("score_compare_methods_element_nf_unknown_nfs",nfs.init,"_nf",nf.init)

data.config2 <- data.config
res <- run_sim(data.config2,method.config, itr = itr, inputDir=inputDir, outputDir = outputDir, script.path = script.path, bfrm.path=bfrm.path, matlab.path = matlab.path, nfs=nfs.init, nf=nf.init, ng=ng, ns=ns, mc.cores = 10)
res.bak <- res
write.csv(res,file.path(table.path,paste0(file.name,".csv")),row.names=F)


################################### plot comparing to other method, element

res <- read.csv(file.path(table.path,paste0(file.name,".csv")))
#res <- res[(!is.na(res$alg) & as.character(res$alg) == "element") | is.na(res$alg),]
#res$method[res$method=="SFAmix"] <- "SFAMix"
#res$method[!is.na(res$alg)] <- paste0(res$method[!is.na(res$alg)], "_", res$alg[!is.na(res$alg)])
res$method[res$method=="SFAMix"] <- "BicMix"
res$method <- factor(res$method, levels=c("BicMix", "SBIF","BFRM","SPCA","KSVD"))

p.comp <- plot_comparison(res)

file.pdf <- file.path(plot.path,paste0(file.name,".pdf"))
pdf(file.pdf,width=8,height=8)
print(p.comp)
dev.off()


################################################################################ 

data.config2 <- data.config[data.config$std.err==2 & data.config$dense == "FALSE",]
res <- run_sim(data.config2,method.config[1], itr = itr, inputDir=inputDir, outputDir = outputDir, script.path = script.path, bfrm.path=bfrm.path, matlab.path = matlab.path, nfs=20, nf=30, ng=ng, ns=ns, mc.cores = 10)
#res <- run_sim(data.config2,method.config, itr = itr, inputDir=inputDir, outputDir = outputDir, script.path = script.path, bfrm.path=bfrm.path, matlab.path = matlab.path, nfs=nfs.init, nf=nf.init, ng=ng, ns=ns, mc.cores = 10)

res.bak <- res
write.csv(res,file.path(table.path,"res_nf_unknown20.csv"),row.names=F)


################################### plot comparing to other method, element

res <- read.csv(file.path(table.path,"res_nf_unknown20.csv"))
#res <- res[(!is.na(res$alg) & as.character(res$alg) == "element") | is.na(res$alg),]
#res$method[res$method=="SFAmix"] <- "SFAMix"
#res$method[!is.na(res$alg)] <- paste0(res$method[!is.na(res$alg)], "_", res$alg[!is.na(res$alg)])
res$method[res$method=="SFAMix"] <- "BicMix"
res$method <- factor(res$method, levels=c("BicMix", "SBIF","BFRM","SPCA","KSVD"))

p.comp <- plot_comparison(res, precis=T)

file.pdf <- file.path(plot.path,"score_compare_methodselement_nf_unknown_other_methods.pdf")
pdf(file.pdf,width=8,height=8)
print(p.comp)
dev.off()









##########################################################
# 
# ab.config <- c("Horseshoe","Strawderman","Uniform")
# #alg.config <- c("matrix","element")
# 
# sfa.config <- cbind("SFAmix",expand.grid(ab.config, ab.config, ab.config))
# names(sfa.config) <- c("method","local","component","global")
# 
# method.config <- rbind(sfa.config,
#                        cbind(method=method.list,local="NA",component="NA",global="NA"))
# 
# local <- as.character(method.config[j,"local"])
# component <- as.character(method.config[j,"component"])
# global <- as.character(method.config[j,"global"])
# 
# shapes <- data.frame(Horseshoe=c(0.5,0.5), Strawderman=c(1,0.5), Uniform=c(1,1))
# 
# a <- shapes[,method.config[j,"local"]][1]; b <- shapes[,method.config[j,"local"]][2]
# c <- shapes[,method.config[j,"component"]][1]; d <- shapes[,method.config[j,"component"]][2]
# e <- shapes[,method.config[j,"global"]][1]; f <- shapes[,method.config[j,"global"]][2]