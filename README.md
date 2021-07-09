## BicMix (now integrated both SFAmix and BicMix)

BicMix is a sparse matrix decomposition tool. Given a matrix Y with dimension of P by N, BicMix decompose it into the product of two sparse matrices LAM and X

This is the C++ implementation of BicMixC wrapped in R.

## Use devtools to install in R

`library(devtools)` <br/>
`install_github("chuangao/BicMix")` <br/>

I found that install_github does not always work. Please install from source (see below) if this is the case.  

## Install from source

`git clone https://github.com/chuangao/BicMix` <br/>
`R CMD INSTALL BicMix` <br/>

## Usage

BicMixR(y=your_input_matrix, nf=100, a=0.5, b=0.5, c=0.5, d=0.5, e=0.5, f=0.5, itr=5000, out_itr=200,  out_dir=director_to_save_results, rsd=NULL, lam_method="matrix", x_method="dense", tol=1e-10, qnorm = TRUE, nf_min = 5) <br/>

**The input file y should have no headers, no missing values, just pure numbers** <br/>
**Also no corrections of confounding beforehand, BicMix will handle that in the dense components** <br/>
**For a gene expression matrix, each gene is a row and each sample is a column** <br/> 

### Arguments
**y** matrix to be decomposed, no missing values are allowed, no headers, space or tab delimited <br/>
**nf** the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50 <br/>
**a** parameter one for the three parameter beta distribution at local level, default to 0.5 to recapitulate horseshoe <br/>
**b** parameter two for the three parameter beta distribution at local level, default to 0.5 to recapitulate horseshoe <br/>
**c** parameter one for the three parameter beta distribution at component level, default to 0.5 to recapitulate horseshoe <br/>
**d** parameter one for the three parameter beta distribution at component level, default to 0.5 to recapitulate horseshoe <br/>
**e** parameter one for the three parameter beta distribution at global level, default to 0.5 to recapitulate horseshoe <br/>
**f** parameter one for the three parameter beta distribution at global level, default to 0.5 to recapitulate horseshoe <b
**itr** the maximum number of iterations the algorithm is allowed to run, default to 5000 <br/>
**out_itr** number of iterations at which temporary output will be written into the specified directory (see below) <br/>
**out_dir** directory where the algorithm will write temporary results (see above) <br/>
**rsd** random seed for initializing the parameter values, default to be randomly drawn <br/>
**lam_method** the method used to update the loading matrix, take values either "matrix" or "element".  if "matrix", then all component are updated simultaneously (slower but more stable, don't need as many iterations to converge); if "element", each component is updated sequentially (faster but less stable, and need more iterations to converge), default to "matrix" <br/>
**x_method** whether induce sparsity on the X matrix, take values either "sparse" or "dense". default to "sparse" <br/>
**tol** tolerance threshold for convergence, default to 1e-10 <br/>
**qnorm** whether to qq-normalize the gene expression matrix, default to TRUE <br/>
**nf_min** the minimum number of factors that needed to be kept (when the signals in the data are small, the default shrinkage parameters in BicMix can be too aggressive that zero factors are left. nf_min make sure at lease some factors are kept, default to 5) <br/>

### Output
**lam** the sparse loading matrix <br/>
**ex** the factor matrix <br/>
**z** a vector indicating whether the corresponding loading is sparse (value of 1) <br/> 
**o** a vector indicating whether the corresponding factor is sparse (value of 1) <br/> 
**nf** the number of factors learned by the model <br/>
**exx** the expected value of the covariance matrix, E(XX^T) <br/>
**itr** the number of iterations for the algorithm to converge <br/>

### Example
library(BicMix)<br/>
#### simulate data where the loading is a mixture of sparse and dense components, and factor is dense <br/>
data = gen_BicMix_data(std=2) <br/>
#### Visulize the loading matrix  <br/>
image(t(data$lam),x=1:ncol(data$lam),y=1:nrow(data$lam),xlab="Loadings",ylab="Samples") <br/>
#### Visulize the factor matrix <br/>
image(t(data$ex),x=1:ncol(data$ex),y=1:nrow(data$ex),xlab="Samples",ylab="Factors") <br/>
#### run algorithm on the simulated data <br/>
result = BicMixR(data$y,nf=100,a=0.5,b=0.5,itr=5000,out_dir="results",tol=1e-10,x_method="sparse",rsd=123) <br/>
#### calculate a correlation matrix of the estimated loading matrix and the true loading matrix. Ideally, there should be one and only one big correlation value for a given row and column of the correlation matrix <br/>
cor.est.real = cor(result$lam[,result$z>0.9],data$lams) <br/>

#### visulize the correlation matrix <br/>
image(cor.est.real,x=1:nrow(cor.est.real),y=1:ncol(cor.est.real),xlab="Recovered loadings",ylab="True loadings") <br/>

#### calculate similarity score of the recovered sparse loading components and the true sparse loading components
cal_score_sparse(result$lam[,result$z>0.9],data$lams)
#### calculate similarity score of the recovered dense loading components and the true dense loading components
cal_score_dense(result$lam[,result$z<=0.9],data$lamd)

#### simulate data where the loading is a mixture of sparse and dense components, and factor is dense <br/>
data = gen_BicMix_data(std=2, type.factor="dense", rsd = 123)
#### perform analysis <br/>
result = BicMixR(data$y,nf=100,out_dir="results",tol=1e-10,x_method="dense",rsd=123)
#### calculate similarity score of the recovered sparse loading components and the true sparse loading components <br/>
cal_score_sparse(result$lam[,result$z>0.9],data$lams)
#### calculate similarity score of the recovered dense loading components and the true dense loading components <br/>
cal_score_dense(result$lam[,result$z<=0.9],data$lamd)

### Documentation
Refer to BicMix.pdf for more usage details <br/>

### Reference
Context Specific and Differential Gene Co-expression Networks via Bayesian Biclustering <br/>
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004791



