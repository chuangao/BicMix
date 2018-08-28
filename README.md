## BicMix

###### BicMix is a sparse matrix decomposition tool. Given a matrix Y with dimension of P by N, BicMix decompose it into the product of two sparse matrices LAM and X

###### This is the C++ implementation of BicMixC wrapped in R. It is compiled single threaded. If you want to run multiple threaded, please check the BicMixC package 

## Use devtools to install in R,

`library(devtools)` <br/>
`install_github("chuangao/BicMix")` <br/>

### If install_github command fails, then clone library into one of your local directory, then install from source

`git clone https://github.com/chuangao/BicMix` <br/>
`R CMD INSTALL BicMix` <br/>

## Usage

BicMixR(y = y, nf = 100, itr = 5000) <br/>

**Please no headers in the input matrix, no missing values, just pure numbers, ideally quantile normalized** <br/>
**Also no corrections of confounding beforehand, BicMix will handle that in the dense components** <br/>
**For a gene expression matrix, it is prefered that each gene is a row and each sample is a column** <br/> 

### Arguments
**y** matrix to be decmoposed, no missing values are allowed <br/>
**nf** the number of factors for the algorithm to start with, will be shrank to a smaller number reflecting the number of factors needed to explain the variance, default to 50 <br/>
**itr** The maximum number of iterations the algorithm is allowed to run, default to 500 <br/>
**out_itr** (Optional) the algorithm will write temporary results into the specified directory (see below) every out_itr number of iterations <br/>
**out_dir** (Optional) Directory where the algorithm will write temporary results into at the specified iteration number(see above) <br/>

### Value
**lam** the sparse loading matrix <br/>
**ex** the factor matrix <br/>
**z** a vector indicating whether the corresponding loading is sparse (value of 1) nf: the number of factors learned by the model <br/>
**exx** the expected value of the covarance matrix, E(XX^T) <br/>
**itr** the number of iterations for the algorithm to converge <br/>

### Examples
library(BicMix) <br/>
\# simulate data <br/>
data = gen_BicMix_data(std=2) <br/>
\# run algorithm on the simulated  <br/>
result = BicMixR(data$y,nf=50,itr=5000) <br/>
\# calculate a correlation matrix of the estimated loading matrix <br/>
\# and the true loading matrix. Ideally, there should be one and only one big correlation value for a given row and column of the correlation matrix <br/>
cor.est.real = cor(result$lam[,result$z==1],data$lams) <br/>
\# visulize the correlation matrix <br/>
image(cor.est.real) <br/>
    
Please refer to BicMix.pdf for more usage details
This work is detailed in https://arxiv.org/abs/1310.4792



