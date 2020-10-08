bfrm <- function(inputDir=NULL, inputFileName=NULL,n=NA,p=NA,k=NA, outputDir=NULL,scriptFile=NULL,bfrmWhere){

    script <- readLines(scriptFile)

    inputFile2 <- file.path(inputDir,paste0(inputFileName))
    outputDir2 <- file.path(outputDir,inputFileName)

    system(paste0("mkdir -p ",outputDir2))

    itemToReplace <- c("inputFileToReplace", "outputDirToReplace","nsToReplace","ngToReplace","nfToReplace")

    #inputFile <- paste0(inputFileName,".txt")
    itemToReplaceWith <- c(inputFile2,outputDir2,n,p,k)

    for(i in 1:length(itemToReplace)){
        index <- grepl(itemToReplace[i],script)
        script[index] <- gsub(itemToReplace[i],itemToReplaceWith[i],script[index])
    }
    

    tmpScript <- file.path(outputDir2,paste0("param.txt"))
    writeLines(script,tmpScript)

    #cmd <- paste0('/Applications/Matlab_R2019b.app/bin/matlab -nodesktop /r ','"inputDir=','"',inputDir,'"',",inputFile=",'"',inputFile,'"',"outputDir=", '"',outputDir,'"','"',";./BicMix/KSVD/mine_sim.m")
    cmd <- paste0("cd ",outputDir2,";",bfrmWhere, " ", tmpScript) 
    #/Applications/Matlab_R2019b.app/bin/matlab -nodesktop -r mine_sim '/Users/cg253/data_bicmix' 'n.effects30_std.err1_seed1_std.effect2.txt' '/Users/cg253/results_bicmix/KSVD'
    system(cmd)
    lam <- as.matrix(read.table(file.path(outputDir2,"mA.txt")))
    return(list(lam=lam[,-1]))
 }
 