ifa <- function(inputDir=NULL,inputFileName=NULL,n=NA,p=NA,k=NA, outputDir=NULL,scriptFile=NULL,matlabWhere){

    system(paste0("mkdir -p ",outputDir))

    script <- readLines(scriptFile)

    itemToReplace <- c("inputDirToReplace", "inputFileToReplace", "outputDirToReplace","nToReplace","pToReplace","kinitToReplace")

    #inputFile <- paste0(inputFileName,".txt")
    itemToReplaceWith <- c(inputDir,inputFileName,outputDir,n,p,k)

    for(i in 1:length(itemToReplace)){
        index <- grepl(itemToReplace[i],script)
        script[index] <- gsub(itemToReplace[i],itemToReplaceWith[i],script[index])
    }
    

    tmpScript <- file.path(outputDir,paste0(inputFileName,".m"))
    writeLines(script,tmpScript)

    #cmd <- paste0('/Applications/Matlab_R2019b.app/bin/matlab -nodesktop /r ','"inputDir=','"',inputDir,'"',",inputFile=",'"',inputFile,'"',"outputDir=", '"',outputDir,'"','"',";./BicMix/KSVD/mine_sim.m")
    cmd <- paste0(matlabWhere,' -nodesktop -nodisplay < ', tmpScript) 
    #/Applications/Matlab_R2019b.app/bin/matlab -nodesktop -r mine_sim '/Users/cg253/data_bicmix' 'n.effects30_std.err1_seed1_std.effect2.txt' '/Users/cg253/results_bicmix/KSVD'
    system(cmd)
    lam <- as.matrix(read.table(file.path(outputDir,paste0(inputFileName,"_lam_ifa.txt"))))
    return(list(lam=lam))
 }
 
 