## BicMix

###### BicMix is a sparse matrix decomposition tool. Given a matrix Y with dimension of P by N, BicMix decompose it into the product of two sparse matrices LAM and X

###### This is the C++ implementation of BicMixC wrapped in R. It is compiled single threaded. If you want to run multiple threaded, please check the BicMixC package 


## Install from source

`git clone https://github.com/chuangao/BicMix` <br/>
`R CMD INSTALL BicMix` <br/>

## Use devtools to install in R

I found that install_github does not always work. Please install from source (see above) is this is the case.  

`library(devtools)` <br/>
`install_github("chuangao/BicMix")` <br/>

## Usage

BicMixR(data$y,nf=100,a=0.5,b=0.5,itr=5000,out_dir="results",tol=1e-5,x_method="dense",rsd=123) <br/>

**Please no headers in the input matrix, no missing values, just pure numbers, ideally quantile normalized** <br/>
**Also no corrections of confounding beforehand, BicMix will handle that in the dense components** <br/>
**For a gene expression matrix, it is prefered that each gene is a row and each sample is a column** <br/> 

### Arguments
**y** matrix to be decmoposed, no missing values are allowed <br/>
**nf** the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50 <br/>
**a** paramater one for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe <br/>
**b** paramater two for the three parameter beta distribution, default to 0.5 to recapitulate horseshoe <br/>
**itr** The maximum number of iterations the algorithm is allowed to run, default to 5000 <br/>
**out_itr** the algorithm will write temporary results into the specified directory (see below) every out_itr number of iterations <br/>
**out_dir** (Optional) Directory where the algorithm will write temporary results into at the specified iteration number(see above) <br/>
**rsd** random seed for initializing the parameter values, default to be randomly drawn <br/>
**x_method** whether induce sparsity on the X matrix, take values either "sparse" or "dense". default to "sparse" <br/>
**tol** tolerance threshold for convergence, default to 1e-5 <br/>


### Output
**lam** the sparse loading matrix <br/>
**ex** the factor matrix <br/>
**z** a vector indicating whether the corresponding loading is sparse (value of 1) <br/> 
**o** a vector indicating whether the corresponding factor is sparse (value of 1) <br/> 
**nf** the number of factors learned by the model <br/>
**exx** the expected value of the covarance matrix, E(XX^T) <br/>
**itr** the number of iterations for the algorithm to converge <br/>

### Example
library(BicMix)<br/>
\## simulate data, the parameter std specifies the standard error of non-zero entries in the ## loading and factor matrices, where a normal distribution of mean zero is assumed for these values. <br/>
data = gen_BicMix_data(std=2) <br/>
\## Visulize the loading matrix  <br/>
image(t(data$lam),x=1:ncol(data$lam),y=1:nrow(data$lam),xlab="Loadings",ylab="Samples") <br/>
\## Visulize the factor matrix <br/>
image(t(data$ex),x=1:ncol(data$ex),y=1:nrow(data$ex),xlab="Samples",ylab="Factors") <br/>
\## run algorithm on the simulated data <br/>
result = BicMixR(data$y,nf=100,a=0.5,b=0.5,itr=5000,out_dir="results",tol=1e-5,x_method="dense",rsd=123) <br/>
\## calculate a correlation matrix of the estimated loading matrix and the true loading matrix. Ideally, there should be one and only one big correlation value for a given row and column of the correlation matrix <br/>
cor.est.real = cor(result$lam[,result$z==1],data$lams) <br/>
\## visulize the correlation matrix <br/>
image(cor.est.real,x=1:nrow(cor.est.real),y=1:ncol(cor.est.real),xlab="Recovered loadings",ylab="True loadings") <br/>

### Documentation
Please refer to BicMix.pdf for more usage details <br/>

### Reference
Context Specific and Differential Gene Co-expression Networks via Bayesian Biclustering <br/>
http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004791



