% KSVD running file
% in this file a synthetic test of the K-SVD algorithm is performed. First,
% a random dictionary with normalized columns is being generated, and then
% a set of data signals, each as a linear combination of 3 dictionary
% element is created, with noise level of 20SNR. this set is given as input
% to the K-SVD algorithm.

% a different mode for activating the K-SVD algorithm is until a fixed
% error is reached in the Sparse coding stage, instead until a fixed number of coefficients is found
% (it was used by us for the
% denoising experiments). in order to switch between those two modes just
% change the param.errorFlag (0 - for fixed number of coefficients, 1 -
% until a certain error is reached).

%function res = mine_sim(inputDir, inputFile, outputDir)
%    disp(0.5)
%    disp(inputDir);

inputDir='inputDirToReplace'
inputFile='inputFileToReplace'
outputDir='outputDirToReplace'


%inputDir='/Users/cg253/data_bicmix'
%inputFile='n.effects30_std.err1_seed1_methodKSVD_std.effect2_b1000000_denseFALSE.txt'
%outputDir='/Users/cg253/results_bicmix/KSVD'


dataFile=strcat(inputDir,'/',inputFile);

%param.K=20;


%factor=[ 10 20 50 ];
%atom=[ 1 3 10 ];

outputFile = inputFile;

%file=strcat();

param.L = 1;   % number of elements in each linear combination.

param.K=FactorNumberToReplace;


param.numIteration = 1000; % number of iteration to execute the K-SVD algorithm.

param.errorFlag = 0; % decompose signals until a certain error is reached. do not use fix number of coefficients.
                     %param.errorGoal = sigma;
param.preserveDCAtom = 0;

%param.errorGoal=1
%%%%%%% creating the data to train on %%%%%%%%
%N = 1500; % number of signals to generate
%n = 20;   % dimension of each data
%SNRdB = 20; % level of noise to be added
%[param.TrueDictionary, D, x] = gererateSyntheticDictionaryAndData(N, param.L, n, param.K, SNRdB);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%% initial dictionary: Dictionary elements %%%%%%%%
param.InitializationMethod =  'DataElements';

param.displayProgress = 0;
disp('Starting to  train the dictionary');


fid = fopen(dataFile);
A = fscanf(fid, '%g', [nrowToReplace ncolToReplace]); %nrowToReplace is the number of samples. ncolToReplace is the number of genes
fclose(fid);

addpath('./KSVD/')

[Dictionary,output]  = KSVD(A,param);

size(Dictionary)

dlmwrite(strcat(outputDir,'/',outputFile,'.dic'), Dictionary, 'delimiter', '\t','precision', 6);
dlmwrite(strcat(outputDir,'/',outputFile,'.coef'), full(output.CoefMatrix), 'delimiter', '\t','precision', 6);

exit(0);

%dlmwrite(strcat('~/SFA/KSVD/result/Dic_file_',file,'_fac_',num2str(param.K),'_atom_',num2str(atom(k)),'.txt'), Dictionary, 'delimiter', '\t','precision', 6);
%dlmwrite(strcat('~/SFA/KSVD/result/Coef_file_',file,'_fac_',num2str(param.K),'_atom_',num2str(atom(k)),'.txt'), full(output.CoefMatrix), 'delimiter', '\t','precision', 6);

