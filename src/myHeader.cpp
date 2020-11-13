#include <stdlib.h>
#include <math.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <limits>
#include <stdio.h>

#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <iostream>

#include <ctime>
#include <time.h>
#include <chrono>


//#include "mkl.h"


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <float.h>
#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_statistics.h>

#include <unistd.h>

#include <dirent.h>
#include <errno.h>

#include <numeric>      // std::iota
#include <algorithm>    // std::sort
#include <functional>   // std::bind

//using namespace std;
using namespace std;
using namespace std::placeholders;

using namespace Eigen;


using Eigen::MatrixXd;



const double PI_mine  =3.141592653589793238462;

void convert_results(MatrixXd &LAM, MatrixXd &EX,MatrixXd &EXX,MatrixXd &Z,MatrixXd &O,VectorXd &PSI,int nf, double *LAM_out, double *EX_out, double *EXX_out, double *Z_out, double *O_out, double *PSI_out, int *nf_out, int s_n, int d_y){
    
    for(int i = 0; i < nf; i++){
        for(int j = 0; j < s_n; j++){
            LAM_out[i*s_n + j] = LAM(j,i);
        }
    }
    
    for(int i = 0; i < d_y; i++){
        for(int j = 0; j < nf; j++){
            EX_out[i*nf + j] = EX(j,i);
        }
    }
    
    for(int i = 0; i < nf; i++){
        for(int j = 0; j < nf; j++){
            EXX_out[i*nf + j] = EXX(j,i);
        }
    }
    
    for(int i = 0; i < nf; i++){
        Z_out[i] = Z(1,i);
    }
    
    for(int i = 0; i < nf; i++){
        O_out[i] = O(1,i);
    }
    
    for(int i =0; i < s_n; i++){
        PSI_out[i] = PSI(i);
    }
    
    nf_out[0] = nf;
}


template <class T>

//Eigen::initParallel();

void write_file(T M, string sinput){
    
    stringstream ss;
    ss.str("");
    ss.clear();
    
    ss << sinput;
    
    ofstream f_file (ss.str().c_str());
    if (f_file.is_open()){
        f_file << M << endl;
    }
    f_file.close();
    
}

bool compare(int a, int b, double* data)
{
    return data[a]<data[b];
}

void convert_row_to_double(MatrixXd& m, double * d, int i_row, int d_y){
    for(int j=0;j<d_y;j++){
        d[j] = m(i_row,j);
    }
}

/*
void quant_norm(MatrixXd& D, MatrixXd& DN, int s_n, int d_y){
    for(int i=0; i < s_n; i++){
        double* din = new double[d_y];
        //double din[10];
        convert_row_to_double(D,din, i,d_y);
        
        //int index[10];
        int* index = new int[d_y];
        
        //for(int j =0; j< d_y; j++){index[j] = j;}
        std::iota(index, index + d_y, 0); // fill index with {0,1,2,...} This only needs to happen once
        
        //std::sort(std::begin(index), std::end(index), std::bind(compare,  _1, _2, din));
        std::sort(index,index+d_y,std::bind(compare,  _1, _2, din));
        
        double* din_sorted = new double[d_y];
        for(int j = 0; j < d_y; j++){
            din_sorted[j] = din[index[j]];
        }
        
        //double qn = gsl_stats_quantile_from_sorted_data (din_sorted, size_t stride, size_t n, double f)
        
        double* p = new double[d_y];
        for(int j = 0; j < d_y; j++){
            p[j] = double(j)/(d_y+1);
            DN(i,index[j]) = gsl_stats_quantile_from_sorted_data(din_sorted,1, d_y,p[j]);
        }
        delete [] din;
        delete [] index;
        delete [] din_sorted;
        delete [] p;
    }
}
*/

// count the number of samples and the number of genes
void cal_y_dimensions(string file_y, string sep, int &s_n, int &d_y){
    string line;
    string field;
    ifstream f;
    
    f.open(file_y.c_str());
    if (! f.is_open()){
        printf("Gene expression or covariate file open failed\n");exit(0);
    }
    
    getline(f,line);
    s_n++;
    istringstream iss(line);
    if(sep.compare("space")==0){
        while(getline(iss,field,' ')){d_y++;}
    }else if(sep.compare("tab")==0){
        while(getline(iss,field,'\t')){d_y++;}
    }else{
        cout << "Please specify a valid separator." << endl << endl;
    }
    while(getline(f,line)){s_n++;}
    f.close();
}

// read in the gene expression matrix
void read_y(string file_y, string sep, MatrixXd &Y_TMP){
    
    string line;
    string field;
    ifstream f;
    
    //f.clear();
    //f.seekg (0, ios_base::beg);
    
    f.open(file_y.c_str());
    if (! f.is_open()){
        printf("Gene expression or covariate file open failed\n");exit(0);
    }
    
    int i=0,j=0;
    while(getline(f,line)){
        istringstream iss(line);
        j=0;
        if(sep.compare("space")==0){
            while(getline(iss,field,' ')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }else{
            while(getline(iss,field,'\t')){
                Y_TMP(i,j)=atof(field.c_str());
                j++;
            }
            i++;
        }
    }
    f.close();
}

void cal_mu(MatrixXd &Y, string dir_out, int s_n, int d_y){
    VectorXd mu = VectorXd::Constant(s_n,0);
    for(int i=0;i<d_y;i++){
        mu += Y.col(i);
    }
    mu *= double(1)/d_y;
    
    for(int i=0;i<d_y;i++){
        Y.col(i) -= mu;
    }
    
    string sin = dir_out + "/mu";
    write_file <VectorXd> (mu,sin);
    
}

void init_lam(MatrixXd &LAM, int s_n, int nf, long seed){
    
    gsl_rng *r;  // random number generator
    r=gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (r, seed);
    
    for (int i=0; i<s_n; i++) {
        for(int j=0;j<nf;j++){
            LAM(i,j)=gsl_ran_gaussian(r,1);
        }
    }
}

void init_ex(MatrixXd &EX, int d_y, int nf, long seed){
    
    gsl_rng *r;  // random number generator
    r=gsl_rng_alloc(gsl_rng_mt19937);
    gsl_rng_set (r, seed);
    
    for (int i=0; i<nf; i++) {
        for(int j=0;j<d_y;j++){
            EX(i,j)=gsl_ran_gaussian(r,1);
        }
    }
}

void inv_psi_vec(VectorXd& psi,VectorXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i)=double(1)/psi(i);
    }
}
void cpy_col_matrix(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int col1,int col2){
    for(int i=0;i<n;i++){
        M1(i,col1)=M2(index(i),col2);
    }
}
void cpy_col_matrix_bak(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int col1,int col2){
    for(int i=0;i<n;i++){
        M1(index(i),col1)=M2(i,col2);
    }
}
void cpy_row_matrix(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int row){
    for(int i=0;i<n;i++){
        M1(0,i)=M2(row,index(i));
    }
}
void cpy_row_matrix_bak(MatrixXd& M1,MatrixXd& M2,VectorXd& index, int n,int row){
    for(int i=0;i<n;i++){
        M1(row,index(i))=M2(0,i);
    }
}

double log_norm(double x,double e, double v){
    //const double PI_mine  =3.141592653589793238462;
    double a = std::numeric_limits<double>::infinity();
    //cout << "a " << a << endl;
    if(v==0){
        return a;
    }else{
        return -0.5*(x-e)*(x-e)/v-0.5*log(v)-0.5*log(2*PI_mine);
    }
}

double E_log_norm(double x,double vx, double v){
    //const double PI  =3.141592653589793238462;
    double a = std::numeric_limits<double>::infinity();
    //cout << "a " << a << endl;
    if(v==0){
        return a;
    }else{
        return -0.5*(x*x+vx)/v-0.5*log(v)-0.5*log(2*PI_mine);
    }
}

double log_gamma(double x, double a, double beta){
    //return a*log(beta)-log(gsl_sf_gamma(a))+(a-1)*log(x)-beta*x;
    if(a==1&&x==0){
        return(beta);
    }else{
        return a*log(beta)-lgamma(a)+(a-1)*log(x)-beta*x;
    }
    //return(log(gsl_ran_gamma_pdf(x,a,double(1)/beta)));
    //boost::lambda
}

double E_log_gamma(double x, double lnx, double a, double beta){
    //return a*log(beta)-log(gsl_sf_gamma(a))+(a-1)*log(x)-beta*x;
    return a*log(beta)-lgamma(a)+(a-1)*lnx-beta*x;
    //boost::lambda
}

void cumsum(VectorXd& S, VectorXd& D, int n){
    for(int i=0;i<n;i++){
        if(i==0){
            D(i)=S(i);
        }
        D(i)=D(i-1)+S(i);
    }
}

double fx(double x,double c){
    return -1*log(x)+0.5/x+c;
}

double dfx(double x){
    return -1*(double(1)/x+0.5/x/x);
}

double NR(double c){
    double x=1e-10;
    for(int i=0;i<500;i++){
        x=x-fx(x,c)/dfx(x);
    }
    return x;
}

int count_total_non_zero(MatrixXd& M,int r, int c){
    int count=0;
    for(int i=0;i<r;i++){
        for(int j=0;j<c;j++){
            if(M(i,j)!=0){
                count++;
            }
        }
    }
    return count;
}


void cal_tau(VectorXd &KAPPA, VectorXd &LAMX, double OMEGA, int nf, double c, double d){
    //Eigen::initParallel();
    //#pragma omp parallel for
    for(int i=0;i<nf;i++){
        KAPPA(i)=double(c+d)/(LAMX(i)+OMEGA);
    }
    
}

/*
 void cal_lamx(VectorXd& LAMX, MatrixXd& O, VectorXd& KAPPA, MatrixXd& RHO,MatrixXd& EX, int nf, int d_y, double b, double c){
 Eigen::initParallel();
 //#pragma omp parallel for
 for(int i=0;i<nf;i++){
 double x_sum=0;
 double sum_c=d_y*b*O(0,i)+c-1-0.5*d_y*O(1,i);
 double at = 2*(KAPPA(i)+O(0,i)*(RHO.row(i).sum()));
 
 x_sum=EX.row(i).dot(EX.row(i));
 
 double bt = O(1,i)*x_sum;
 LAMX(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
 
 if(LAMX(i)<1e-50){
 LAMX(i)=1e-50;
 }
 
 }
 }
 */

void red_dim(MatrixXd& EX,MatrixXd& EXX, MatrixXd& LAM, MatrixXd& THETA, MatrixXd& DELTA, VectorXd& PHI, VectorXd& TAU, MatrixXd& Z, MatrixXd& logZ, VectorXd& count_lam, VectorXd& index, VectorXd& KAPPA, VectorXd& LAMX, MatrixXd& RHO, MatrixXd& SIGMA, VectorXd& count_x, MatrixXd& O, MatrixXd& logO, MatrixXd& LPL, int nf, int nt, int s_n, int d_y, int nmix, double zi, bool include_x){
    
    MatrixXd EX2=MatrixXd::Constant(nt,d_y,0);
    //MatrixXd TEX2=MatrixXd::Constant(d_y,nt,0);
    //MatrixXd VX2=MatrixXd::Constant(nt,nt,0);
    MatrixXd EXX2=MatrixXd::Constant(nt,nt,0);
    
    MatrixXd LAM2=MatrixXd::Constant(s_n,nt,0);
    
    MatrixXd THETA2=MatrixXd::Constant(s_n,nf,0);
    MatrixXd DELTA2=MatrixXd::Constant(s_n,nf,0);
    VectorXd PHI2 = VectorXd::Constant(nf,0);
    VectorXd TAU2 = VectorXd::Constant(nf,0);
    MatrixXd Z2 = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ2 = MatrixXd::Constant(nmix,nf,log(zi));
    VectorXd count_lam2 = VectorXd::Constant(nf,0);
    VectorXd index2 = VectorXd::Constant(nf,0);
    //ID2.diagonal()=id_v2;
    
    // for EX related
    VectorXd KAPPA2,LAMX2;
    MatrixXd RHO2,SIGMA2,O2,logO2;
    
    if(include_x){
        KAPPA2 = VectorXd::Constant(nf,0);
        LAMX2 = VectorXd::Constant(nf,0);
        RHO2=MatrixXd::Constant(nf,d_y,0);
        SIGMA2=MatrixXd::Constant(nf,d_y,0);
        O2 = MatrixXd::Constant(nmix,nf,zi);
        logO2 = MatrixXd::Constant(nmix,nf,log(zi));
    }
    
    VectorXd count_x2 = VectorXd::Constant(nf,0);
    
    //MatrixXd LOGVO2 = MatrixXd::Constant(nmix,1,log(0.5));
    
    MatrixXd LPL2 = MatrixXd::Constant(nf,nf,0);
    
    for(int i=0;i<nf;i++){
        EX2.row(i)=EX.row(index(i));
        //TEX2=EX2.transpose();
        for(int j=0;j<nf;j++){
            EXX2(i,j)=EXX(index(i),index(j));
            LPL2(i,j)=LPL(index(i),index(j));
            
        }
        //LP2.row(i)=LP.row(index(i));
        LAM2.col(i)=LAM.col(index(i));
        
        THETA2.col(i)=THETA.col(index(i));
        DELTA2.col(i)=DELTA.col(index(i));
        PHI2(i)=PHI(index(i));
        TAU2(i)=TAU(index(i));
        Z2.col(i)=Z.col(index(i));
        logZ2.col(i)=logZ.col(index(i));
        
        count_lam2(i)=count_lam(index(i));
        index2(i)=index(i);
        count_x2(i)=count_x(index(i));
        // EX related
        if(include_x){
            KAPPA2(i)=KAPPA(index(i));
            LAMX2(i)=LAMX(index(i));
            RHO2.row(i)=RHO.row(index(i));
            SIGMA2.row(i)=SIGMA.row(index(i));
            
            O2.col(i)=O.col(index(i));
            logO2.col(i)=logO.col(index(i));
        }
    }
    
    // Assign the new parameters back
    EX=EX2;
    //TEX=TEX2;
    //VX=VX2;
    EXX=EXX2;
    //LP=LP2;
    //ID=ID2;
    LAM=LAM2;
    //ELL=ELL2;
    //LAM_BAK=LAM_BAK2;
    THETA=THETA2;
    DELTA=DELTA2;
    PHI=PHI2;
    TAU=TAU2;
    Z=Z2;
    logZ=logZ2;
    //LAM_TOP=LAM_TOP2;
    //LAM_BOT=LAM_BOT2;
    count_lam=count_lam2;
    index=index2;
    count_x=count_x2;
    
    // EX related
    if(include_x){
        KAPPA=KAPPA2;
        LAMX=LAMX2;
        RHO=RHO2;
        SIGMA=SIGMA2;
        
        O=O2;
        logO=logO2;
    }
    LPL=LPL2;
    //partV=partV2;
    
}



/*
void cal_lamx(VectorXd& LAMX, MatrixXd& O, VectorXd& KAPPA, MatrixXd& RHO,MatrixXd& EX, int nf, int d_y, double b, double c){
    Eigen::initParallel();
    //#pragma omp parallel for
    for(int i=0;i<nf;i++){
        double x_sum=0;
        double sum_c=d_y*b*O(0,i)+c-1-0.5*d_y*O(1,i);
        double at = 2*(KAPPA(i)+O(0,i)*(RHO.row(i).sum()));
        
        //x_sum=EXX(i,i);
        
        x_sum=EX.row(i).dot(EX.row(i));
        
        double bt = O(1,i)*x_sum;
        LAMX(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
        
        if(LAMX(i)<1e-50){
            LAMX(i)=1e-50;
        }  
        
    }
}
*/

void cal_phi(VectorXd& PHI, MatrixXd& Z, VectorXd& TAU, MatrixXd& DELTA,MatrixXd& LAM, int nf, int s_n, double b, double c, bool cln){
    //Eigen::initParallel();
    
    //#pragma omp parallel for
    for(int i=0;i<nf;i++){
        double lam_sum=0;
        double sum_c=s_n*b*Z(0,i)+c-1-0.5*s_n*Z(1,i);
        double at = 0;
        
        if(cln){
            at = 2*(TAU(i)+Z(0,i)*(DELTA.col(i).sum()));
            lam_sum=LAM.col(i).dot(LAM.col(i));
        }else{
            at = 2*(TAU(i)+Z(0,i)*(DELTA.row(i).sum()));
            lam_sum=LAM.row(i).dot(LAM.row(i));
        }
        
        double bt = Z(1,i)*lam_sum;
        PHI(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
        
        if(PHI(i)<1e-300){
            PHI(i)=1e-300;
        }
        
    }
    
}

/*
void cal_rho(MatrixXd& RHO, MatrixXd& SIGMA, VectorXd& LAMX, double a, double b, int d_y, int nf){
    Eigen::initParallel();
    //#pragma omp parallel for
    for(int i=0;i<d_y;i++){
        //#pragma omp parallel for
        for(int j=0;j<nf;j++){
            RHO(j,i)=double((a+b))/(SIGMA(j,i)+LAMX(j));
        }
    }
}
*/

void cal_delta(MatrixXd& DELTA,MatrixXd& THETA,VectorXd& PHI,double a, double b,int nf, int s_n, bool cln){
    
    //Eigen::initParallel();
    //#pragma omp parallel for collapse(2)
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            if(cln){
                DELTA(i,j)=double((a+b))/(THETA(i,j)+PHI(j));
            }else{
                DELTA(j,i)=double((a+b))/(THETA(j,i)+PHI(j));
            }
        }
    }
    
}

/*
void cal_sigma(MatrixXd& SIGMA,MatrixXd& EX, MatrixXd& RHO, double a, int d_y, int nf){
    double a23=(2*a-3);
    Eigen::initParallel();
    //#pragma omp parallel for
    for(int i=0;i<d_y;i++){
        //#pragma omp parallel for
        for(int j=0;j<nf;j++){
            SIGMA(j,i)=double(a23+sqrt(a23*a23+8*EX(j,i)*EX(j,i)*RHO(j,i)))/4/RHO(j,i);
        }
    }
}
*/

void cal_theta(MatrixXd& THETA,MatrixXd& LAM,MatrixXd& DELTA,double a, int s_n, int nf, bool cln){
    
    double a23=(2*a-3);
    //Eigen::initParallel();
    //#pragma omp parallel for collapse(2)
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            if(cln){
                THETA(i,j)=double(a23+sqrt(a23*a23+8*LAM(i,j)*LAM(i,j)*DELTA(i,j)))/4/DELTA(i,j);
            }else{
                THETA(j,i)=double(a23+sqrt(a23*a23+8*LAM(j,i)*LAM(j,i)*DELTA(j,i)))/4/DELTA(j,i);
            }
        }
    }
    
}

int count_nonzero_lam_ex(VectorXd &count_x, VectorXd &count_lam, VectorXd &index, MatrixXd &LAM, MatrixXd &EX, VectorXd &PHI, VectorXd &LAMX, int nf, int s_n, int d_y){
    
    count_x.setZero();
    for(int i=0;i<nf;i++){
        for(int j=0;j<d_y;j++){
            if(EX(i,j)!=0){
                count_x(i) +=  1;
            }
        }
    }
    
    // count the number of non-zero values in each column of the lambda matrix
    count_lam.setZero();
    for(int i=0;i<nf;i++){
        for(int j=0;j<s_n;j++){
            if(LAM(j,i)!=0){
                count_lam(i) +=  1;
            }
        }
    }
    
    // Count the number of loadings or factors that have either at least one non-zero values (phi is non-zero too).
    int count_all = 0;
    for(int i=0;i<nf;i++){
        if(count_lam(i)>1&&PHI(i)!=0&&count_x(i)>1&&LAMX(i)!=0){
            index(count_all)=i;
            count_all ++;
        }
    }
    return count_all;
}


int count_nonzero_lam_cov(VectorXd& count_lam, VectorXd& index, VectorXd& index_cov_pos, MatrixXd& LAM, VectorXd& PHI, int nf, int s_n){
    
    // count the number of non-zero values in each column of the lambda matrix
    count_lam.setZero();
    for(int i=0;i<nf;i++){
        for(int j=0;j<s_n;j++){
            if(LAM(j,i)!=0){
                count_lam(i) +=  1;
            }
        }
    }
    
    // Count the number of loadings or factors that have either at least one non-zero values (phi is non-zero too).
    int count_all = 0;
    for(int i=0;i<nf;i++){
        if(count_lam(i)>1&&PHI(i)!=0){
            index(count_all)=i;
            count_all ++;
        }
    }
    
    VectorXd index_tmp = VectorXd::Constant(count_all,0);
    for(int i=0;i<count_all;i++){
        index_tmp(i) = index_cov_pos(index(i));
    }
    
    index_cov_pos = index_tmp;
    return count_all;
}
/*
void cal_ex_single_thread(MatrixXd& EX, MatrixXd& LAM,MatrixXd& Y,VectorXd& PSI_INV, MatrixXd& EXX,MatrixXd& SIGMA,VectorXd& LAMX,MatrixXd& O,MatrixXd& LPL,int d_y, int s_n, int nf){
    
    
    MatrixXd LAM_T=LAM.transpose();
    MatrixXd partR = MatrixXd::Constant(nf,d_y,0);
    for(int i=0;i<s_n;i++){
        partR = partR + LAM_T.col(i)*PSI_INV(i)*Y.row(i);
    }
    
    EXX.setZero();
    
    //Eigen::initParallel();

    for(int j=0;j<d_y;j++){
        int count_indexALL=0;
        VectorXd indexALL = VectorXd::Constant(nf,0);
        for(int i=0;i<nf;i++){
            if(SIGMA(i,j)!=0&&LAMX(i)!=0){
                indexALL(count_indexALL)=i;
                count_indexALL++;
            }
        }
        
        if(count_indexALL==0){
            EX.col(j).setZero();
            continue;
        }
        
        //partV.setZero();
        MatrixXd partV = MatrixXd::Constant(nf,nf,0);
        //#pragma omp parallel for
        for(int i1=0;i1<count_indexALL;i1++){
            for(int i2=0;i2<count_indexALL;i2++){
                partV(indexALL(i1),indexALL(i2)) = LPL(indexALL(i1),indexALL(i2));
            }
            partV(indexALL(i1),indexALL(i1)) += O(0,indexALL(i1))/SIGMA(indexALL(i1),j)+O(1,indexALL(i1))/LAMX(indexALL(i1));
        }
        
        MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
 	    //#pragma omp parallel for collapse(2)
        for(int i1=0;i1<count_indexALL;i1++){
            for(int i2=0;i2<count_indexALL;i2++){
                partVI(i1,i2)=partV(indexALL(i1),indexALL(i2));
            }
        }
        MatrixXd partRI=MatrixXd::Constant(count_indexALL,1,0);
        MatrixXd EXI=MatrixXd::Constant(count_indexALL,1,0);
        MatrixXd EXXI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
        cpy_col_matrix(partRI,partR,indexALL,count_indexALL,0,j);
        
        MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);
        MatrixXd vx=MatrixXd::Constant(count_indexALL,count_indexALL,0);
        vx = partVI.lu().solve(IDNF);
        EXI=vx*partRI;
        for(int i=0;i<count_indexALL;i++){
            if(SIGMA(indexALL(i),j)==0){
                EXI(i)=0;
            }
        }
        
        cpy_col_matrix_bak(EX,EXI,indexALL,count_indexALL,j,0);
        EXXI=EXI*EXI.transpose();
 	    //#pragma omp parallel for collapse(2)
        for(int i1=0;i1<count_indexALL;i1++){
            for(int i2=0;i2<count_indexALL;i2++){
                EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2)+vx(i1,i2);
            }
        }
    }
    
    for(int i=0;i<nf;i++){
        for(int j=0;j<d_y;j++){
            if(SIGMA(i,j)==0){
                EX(i,j)=0;
            }
        }
    }
}
*/

void cal_ex(MatrixXd& EX, MatrixXd& LAM,MatrixXd& Y,VectorXd& PSI_INV, MatrixXd& EXX,MatrixXd& SIGMA,VectorXd& LAMX,MatrixXd& O,MatrixXd& LPL,int d_y, int s_n, int nf){
    
    MatrixXd LAM_T=LAM.transpose();
    MatrixXd partR = MatrixXd::Constant(nf,d_y,0);
    for(int i=0;i<s_n;i++){
        partR = partR + LAM_T.col(i)*PSI_INV(i)*Y.row(i);
    }
    
    EXX.setZero();
    
    //Eigen::initParallel();

    #pragma omp parallel
    {
        MatrixXd EXX_tmp = MatrixXd::Constant(nf,nf,0); 
        #pragma omp for nowait
        for(int j=0;j<d_y;j++){
            int count_indexALL=0;
            VectorXd indexALL = VectorXd::Constant(nf,0);
            for(int i=0;i<nf;i++){
                if(SIGMA(i,j)!=0&&LAMX(i)!=0){
                    indexALL(count_indexALL)=i;
                    count_indexALL++;
                }
            }
            
            if(count_indexALL==0){
                EX.col(j).setZero();
                continue;
            }
            
            //partV.setZero();
            MatrixXd partV = MatrixXd::Constant(nf,nf,0);
            //#pragma omp parallel for
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    partV(indexALL(i1),indexALL(i2)) = LPL(indexALL(i1),indexALL(i2));
                }
                partV(indexALL(i1),indexALL(i1)) += O(0,indexALL(i1))/SIGMA(indexALL(i1),j)+O(1,indexALL(i1))/LAMX(indexALL(i1));
            }
            
            MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            //#pragma omp parallel for collapse(2)
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    partVI(i1,i2)=partV(indexALL(i1),indexALL(i2));
                }
            }
            MatrixXd partRI=MatrixXd::Constant(count_indexALL,1,0);
            MatrixXd EXI=MatrixXd::Constant(count_indexALL,1,0);
            MatrixXd EXXI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            cpy_col_matrix(partRI,partR,indexALL,count_indexALL,0,j);
            
            MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);
            MatrixXd vx=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            vx = partVI.lu().solve(IDNF);
            EXI=vx*partRI;
            for(int i=0;i<count_indexALL;i++){
                if(SIGMA(indexALL(i),j)==0){
                    EXI(i)=0;
                }
            }
            
            cpy_col_matrix_bak(EX,EXI,indexALL,count_indexALL,j,0);
            EXXI=EXI*EXI.transpose();
            //#pragma omp parallel for collapse(2)
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    EXX_tmp(indexALL(i1),indexALL(i2)) += EXXI(i1,i2)+vx(i1,i2);
                }
            }
        }
        #pragma omp critical
        {
           EXX += EXX_tmp;
        }
    }

    for(int i=0;i<nf;i++){
        for(int j=0;j<d_y;j++){
            if(SIGMA(i,j)==0){
                EX(i,j)=0;
            }
        }
    }
}

void cal_vl(MatrixXd& vl, VectorXd& indexALL,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV, int count_indexALL, int j){
    
    //Eigen::initParallel();
    
    partV.setZero();
    //#pragma omp parallel for
    for(int i1=0;i1<count_indexALL;i1++){
        for(int i2=0;i2<count_indexALL;i2++){
            partV(indexALL(i1),indexALL(i2)) = PSI_INV(j)*EXX(indexALL(i1),indexALL(i2));
        }
        partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
    }
    
    MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
    //#pragma omp parallel for collapse(2)
    for(int i1=0;i1<count_indexALL;i1++){
        for(int i2=0;i2<count_indexALL;i2++){
            partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
        }
    }
    
    MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);
    vl = partVI.ldlt().solve(IDNF);
    //vl = partVI.lu().solve(IDNF);
}
/*
void cal_lam_single_thread(MatrixXd& LAM, MatrixXd& Y,MatrixXd& EX,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI, int s_n, int d_y, int nf){
    
    MatrixXd PSIY=MatrixXd::Constant(s_n,d_y,0);
    for(int i=0;i<s_n;i++){
        PSIY.row(i) = PSI_INV(i)*Y.row(i);
    }
    
    MatrixXd partL = PSIY*EX.transpose();
    
    LPL.setZero();
    
    //Eigen::initParallel();
    //#pragma omp parallel for
    for(int j=0;j<s_n;j++){
        int count_indexALL=0;
        VectorXd indexALL = VectorXd::Constant(nf,0);
        for(int i=0;i<nf;i++){
            if(THETA(j,i)!=0 && PHI(i)!=0){
                indexALL(count_indexALL)=i;
                count_indexALL++;
            }
        }
        
        if(count_indexALL==0){
            LAM.row(j).setZero();
            continue;
        }
        MatrixXd vl = MatrixXd::Constant(count_indexALL, count_indexALL,0);
        
        MatrixXd partV = MatrixXd::Constant(nf,nf,0);
        cal_vl(vl,indexALL,PSI_INV,EXX,Z,THETA,PHI,partV,count_indexALL,j);
        
        MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
        MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
        MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
        cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
        
        LAMI=partLI*vl;
        for(int i=0;i<count_indexALL;i++){
            if(THETA(j,indexALL(i))==0){
                LAMI(i)=0;
            }
        }
        
        cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
        LLI=LAMI.transpose()*LAMI;
        
        //#pragma omp parallel for collapse(2)
        for(int i1=0;i1<count_indexALL;i1++){
            for(int i2=0;i2<count_indexALL;i2++){
                LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j)*(LLI(i1,i2)+vl(i1,i2));
                //vLXL(j,j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
            }
        }
        
    }
    
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            if(THETA(i,j)==0){
                LAM(i,j)=0;
            }
        }
    }
    
    
}
*/

void cal_lam(MatrixXd& LAM, MatrixXd& Y,MatrixXd& EX,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI, int s_n, int d_y, int nf){
  
    MatrixXd PSIY=MatrixXd::Constant(s_n,d_y,0);
    for(int i=0;i<s_n;i++){
        PSIY.row(i) = PSI_INV(i)*Y.row(i);
    }
    
    MatrixXd partL = PSIY*EX.transpose();
    
    LPL.setZero();

    
    //cout << "nf " << nf << endl;
    //cout << LAM.block(0,0,5,5) << endl;
    //cout << LPL.block(0,0,5,5) << endl;

    #pragma omp parallel
    {
        MatrixXd LPL_tmp = MatrixXd::Constant(nf, nf,0);

        #pragma omp for nowait
        for(int j=0;j<s_n;j++){
               
                int count_indexALL=0;
                VectorXd indexALL = VectorXd::Constant(nf,0);
                for(int i=0;i<nf;i++){
                    if(THETA(j,i)!=0 && PHI(i)!=0){
                        indexALL(count_indexALL)=i;
                        count_indexALL++;
                    }
                }
                
                if(count_indexALL==0){
                    LAM.row(j).setZero();
                    //continue;
                }else{
                    MatrixXd vl = MatrixXd::Constant(count_indexALL, count_indexALL,0);
                    
                    MatrixXd partV = MatrixXd::Constant(nf,nf,0);
                    cal_vl(vl,indexALL,PSI_INV,EXX,Z,THETA,PHI,partV,count_indexALL,j);
                    
                    MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
                    MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
                    MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
                    cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
                    
                    LAMI=partLI*vl;
                    for(int i=0;i<count_indexALL;i++){
                        if(THETA(j,indexALL(i))==0){
                            LAMI(i)=0;
                        }
                    }
                    
                    cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
                    LLI=LAMI.transpose()*LAMI;
                    
                    for(int i1=0;i1<count_indexALL;i1++){
                        for(int i2=0;i2<count_indexALL;i2++){
                            LPL_tmp(indexALL(i1),indexALL(i2)) = LPL_tmp(indexALL(i1),indexALL(i2)) + PSI_INV(j)*(LLI(i1,i2)+vl(i1,i2));
                            //vLXL(j,j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
                        }
                    }
                }
                
            
        }
        #pragma omp critical
        {
           LPL += LPL_tmp;
        }
    }
    
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            if(THETA(i,j)==0){
                LAM(i,j)=0;
            }
        }
    }
    
    
}

/*
void loop_lam(MatrixXd& LAM, MatrixXd& Y,MatrixXd& EX,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI, int s_n, int d_y, int nf, MatrixXd& partL, int j){

        int count_indexALL=0;
        VectorXd indexALL = VectorXd::Constant(nf,0);
        for(int i=0;i<nf;i++){
            if(THETA(j,i)!=0 && PHI(i)!=0){
                indexALL(count_indexALL)=i;
                count_indexALL++;
            }
        }
        
        if(count_indexALL==0){
            LAM.row(j).setZero();
            //continue;
        }else{
            MatrixXd vl = MatrixXd::Constant(count_indexALL, count_indexALL,0);
            
            MatrixXd partV = MatrixXd::Constant(nf,nf,0);
            cal_vl(vl,indexALL,PSI_INV,EXX,Z,THETA,PHI,partV,count_indexALL,j);
            
            MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
            MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
            MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
            cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
            LAMI=partLI*vl;
            for(int i=0;i<count_indexALL;i++){
                if(THETA(j,indexALL(i))==0){
                    LAMI(i)=0;
                }
            }
            
            cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
            LLI=LAMI.transpose()*LAMI;
            
            //#pragma omp parallel for collapse(2)
            for(int i1=0;i1<count_indexALL;i1++){
                for(int i2=0;i2<count_indexALL;i2++){
                    LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j)*(LLI(i1,i2)+vl(i1,i2));
                    //vLXL(j,j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
                }
            }
        }
       
}

void cal_lam_multi_thread2(MatrixXd& LAM, MatrixXd& Y,MatrixXd& EX,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI, int s_n, int d_y, int nf){
 
    MatrixXd PSIY=MatrixXd::Constant(s_n,d_y,0);
    for(int i=0;i<s_n;i++){
        PSIY.row(i) = PSI_INV(i)*Y.row(i);
    }
    
    MatrixXd partL = PSIY*EX.transpose();
    
    LPL.setZero();
    
    //Eigen::initParallel();
    //#pragma omp parallel for schedule(static) reduction(+:LPL)
    for(int j=0;j<s_n;j++){
        loop_lam(LAM, Y, EX,PSI_INV,EXX,Z,LPL,THETA,PHI, s_n, d_y, nf, partL, j);
    }
    
    for(int i=0;i<s_n;i++){
        for(int j=0;j<nf;j++){
            if(THETA(i,j)==0){
                LAM(i,j)=0;
            }
        }
    }
    
    
}
*/

void cal_lam_element_wise(MatrixXd& Y, MatrixXd& LAM, MatrixXd& EX, MatrixXd& THETA, MatrixXd& EXX, MatrixXd& Z, VectorXd& PSI_INV, VectorXd& PHI, int s_n, int nf){
    
    MatrixXd LAM_bak = LAM;
    MatrixXd LAMXX = LAM * EXX;
    
    //MatrixXd TOP = VectorXd::Constant(s_n,0);
    #pragma omp parallel for
    for(int i=0;i<nf;i++){
        VectorXd TOP = VectorXd::Constant(s_n,0);
        //pragma omp parallel for
        for(int j=0;j<s_n;j++){
            TOP(j) = Y.row(j).dot(EX.row(i)) - LAMXX(j,i) + LAM(j,i) * EXX(i,i);
            if(Z(0,i)==0){
                LAM(j,i) = TOP(j)/(EXX(i,i)+double(1)/PSI_INV(j)*Z(1,i)/PHI(i));
            }
            else if(Z(1,i)==0){
                LAM(j,i) = TOP(j)/(EXX(i,i)+double(1)/PSI_INV(j)*Z(0,i)/THETA(j,i));
            }
            else{
                LAM(j,i) = TOP(j)/(EXX(i,i)+double(1)/PSI_INV(j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
            }
        }
        LAMXX = LAMXX  + (LAM.col(i) - LAM_bak.col(i)) * EXX.row(i);
    }
}


// the folling is slow, because LAMXX_tmp didn't update LAMXX immediately
void cal_lam_element_wise_multi_threaded(MatrixXd& Y, MatrixXd& LAM, MatrixXd& EX, MatrixXd& THETA, MatrixXd& EXX, MatrixXd& Z, VectorXd& PSI_INV, VectorXd& PHI, int s_n, int nf){
    
    MatrixXd LAM_bak = LAM;
    MatrixXd LAMXX = LAM * EXX;
    
    MatrixXd TOP = MatrixXd::Constant(s_n,nf,0);

    #pragma omp parallel
    {
        MatrixXd LAMXX_tmp = MatrixXd::Constant(s_n, nf,0);
        #pragma omp for nowait
        for(int i=0;i<nf;i++){
            for(int j=0;j<s_n;j++){
                TOP(j,i) = Y.row(j).dot(EX.row(i)) - LAMXX(j,i) + LAM(j,i) * EXX(i,i);
                if(Z(0,i)==0){
                    LAM(j,i) = TOP(j,i)/(EXX(i,i)+double(1)/PSI_INV(j)*Z(1,i)/PHI(i));
                }
                else if(Z(1,i)==0){
                    LAM(j,i) = TOP(j,i)/(EXX(i,i)+double(1)/PSI_INV(j)*Z(0,i)/THETA(j,i));
                }
                else{
                    LAM(j,i) = TOP(j,i)/(EXX(i,i)+double(1)/PSI_INV(j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
                }
            }
            LAMXX_tmp = LAMXX_tmp  + (LAM.col(i) - LAM_bak.col(i)) * EXX.row(i);
        }
        #pragma omp critical
        {
           LAMXX += LAMXX_tmp;
        }
    }
}

/*
void cal_lam_element_wise2(MatrixXd& Y, MatrixXd& LAM, MatrixXd& EX, MatrixXd& THETA, MatrixXd& EXX, MatrixXd& Z, VectorXd& PSI_INV, VectorXd& PHI, int s_n, int nf){
    // suspect that cal_lam_element_wise is too slow, because it didn't effectively use the updated the parameters. come up with this version, so that the next element is based on the updated previous parameter values, turned out didn't speed up much.
    
    MatrixXd LAM_bak = LAM;
    MatrixXd LAMXX = LAM * EXX;
    //VectorXd TOP = VectorXd::Constant(nf,0);
    double top = 0;
    
    for(int j=0;j<s_n;j++){
        for(int i=0;i<nf;i++){
            double lamji_bak = LAM(j,i);
            top = Y.row(j).dot(EX.row(i)) - LAMXX(j,i) + LAM(j,i) * EXX(i,i);
            
            if(Z(0,i)==0){
                LAM(j,i) = top/(EXX(i,i)+double(1)/PSI_INV(j)*Z(1,i)/PHI(i));
            }
            else if(Z(1,i)==0){
                LAM(j,i) = top/(EXX(i,i)+double(1)/PSI_INV(j)*Z(0,i)/THETA(j,i));
            }
            else{
                LAM(j,i) = top/(EXX(i,i)+double(1)/PSI_INV(j)*(Z(1,i)/PHI(i)+Z(0,i)/THETA(j,i)));
            }
            LAMXX.row(j) = LAMXX.row(j)  + (LAM(j,i) - lamji_bak) * EXX.row(i);
        }
        //YLX=YLX-LAM.col(i)*EX.row(i);
    }
}
*/

void cal_ex_simple(MatrixXd& EX, MatrixXd& EXX, MatrixXd& LAM, VectorXd& PSI_INV, MatrixXd& Y, int s_n, int nf, int d_y){
    
    //Eigen::initParallel();
    
    MatrixXd LP = MatrixXd::Constant(nf,s_n,0);
    MatrixXd ID=MatrixXd::Identity(nf,nf);
    //#pragma omp parallel for
    for(int i=0;i<s_n;i++){
        LP.col(i)=LAM.transpose().col(i)*PSI_INV(i);
    }
    MatrixXd VX=(LP*LAM+ID).lu().solve(ID);
    EX=VX*LP*Y;
    EXX=EX*EX.transpose()+VX*d_y;
    //EXX=EX*EX.transpose();
    //TEX=EX.transpose();
}

void cal_z(MatrixXd& logZ,MatrixXd& LOGV,MatrixXd& LAM, MatrixXd& THETA, MatrixXd& DELTA, VectorXd& PHI, MatrixXd& Z, int nf, int s_n, double a, double b, double alpha, double beta, bool cln){
    
    // logZ
    //Eigen::initParallel();
    //#pragma omp parallel for
    for(int i=0;i<nf;i++){
        logZ(0,i)=LOGV(0,0);
        logZ(1,i)=LOGV(1,0);
        
        //#pragma omp parallel for
        for(int j=0;j<s_n;j++){
            if(cln){
                logZ(0,i)=logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+log_gamma(THETA(j,i),a,DELTA(j,i))+log_gamma(DELTA(j,i),b,PHI(i));
                logZ(1,i)=logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
            }else{
                logZ(0,i)=logZ(0,i)+log_norm(LAM(i,j),0,THETA(i,j))+log_gamma(THETA(i,j),a,DELTA(i,j))+log_gamma(DELTA(i,j),b,PHI(i));
                logZ(1,i)=logZ(1,i)+log_norm(LAM(i,j),0,PHI(i));
            }
        }
    }
    
    // Z
    //#pragma omp parallel for
    for(int i=0;i<nf;i++){
        Z(0,i)=double(1)/(1+exp(logZ(1,i)-logZ(0,i)));
        Z(1,i)=1-Z(0,i);
    }
    
    double ps1=alpha;
    double ps2=beta;
    
    // Probability
    //#pragma omp parallel for
    for(int i=0;i<nf;i++){
        ps1=ps1+Z(0,i);
        ps2=ps2+Z(1,i);
    }
    double dgama = gsl_sf_psi(ps1+ps2);
    LOGV(0,0)=gsl_sf_psi(ps1)-dgama;
    LOGV(1,0)=gsl_sf_psi(ps2)-dgama;
    
    
    
}

void cal_psi(VectorXd& PSI,MatrixXd& LX, MatrixXd& Y,MatrixXd& LAM, MatrixXd& EXX, int s_n, int d_y){
    //Eigen::initParallel();
    //#pragma omp parallel for
    for(int i=0;i<s_n;i++){
        PSI(i)=(0.5*(Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i)))+1)/(double(d_y)/2+1);

    }
    //PSI=PSI/d_y;
}

void cal_lam_all(MatrixXd& LAM, MatrixXd& Y,MatrixXd& EX,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,MatrixXd& THETA,VectorXd& PHI, int s_n, int d_y, int nf, double a, double b, double c, double d, double g, double h, double GAMMA, double ETA, double nu, VectorXd& TAU, MatrixXd& DELTA, double alpha, double beta,string lam_method){
    
    GAMMA=double(g+h)/(ETA+nu);
    ETA=double((d*nf+g))/(GAMMA+TAU.sum());
    // TAU
    cal_tau(TAU, PHI, ETA, nf, c, d);
    
    // phi for loading
    cal_phi(PHI, Z, TAU, DELTA,LAM, nf, s_n, b, c,true);
    
    // calculate DELTA
    cal_delta(DELTA, THETA, PHI, a, b, nf, s_n,true);
    
    // calculate THETA
    cal_theta(THETA, LAM, DELTA, a, s_n, nf,true);
    
    if(lam_method.compare("matrix") == 0){
        cal_lam(LAM, Y, EX, PSI_INV, EXX, Z, LPL, THETA, PHI, s_n,  d_y, nf);
    }
    if(lam_method.compare("element") == 0){
        
        cal_lam_element_wise(Y, LAM, EX, THETA, EXX, Z, PSI_INV, PHI, s_n, nf);
        //cal_lam_element_wise(Y, LAM, EX, THETA, EXX, Z, PSI_INV, PHI, s_n, nf);
        
    }
    
}

void cal_ex_all(MatrixXd& LAM, MatrixXd& Y,MatrixXd& EX,VectorXd& PSI_INV,MatrixXd& EXX,MatrixXd& O,MatrixXd& LPL,MatrixXd& SIGMA,VectorXd& LAMX, int s_n, int d_y, int nf, double a, double b, double c, double d, double g, double h, double VARSIG, double OMEGA, double XI, VectorXd& KAPPA, MatrixXd& RHO, MatrixXd& logO, MatrixXd& LOGVO, double alpha, double beta,string x_method){
    
    if(x_method.compare("sparse") == 0){
        VARSIG=double(g+h)/(OMEGA+XI);
        OMEGA=double((d*nf+g))/(VARSIG+KAPPA.sum());
        
        // KAPPA
        cal_tau(KAPPA, LAMX, OMEGA, nf, c, d);
        
        // PHI for x
        //cal_lamx(LAMX, O, KAPPA, RHO, EX, nf, d_y, b, c);
        cal_phi(LAMX, O, KAPPA, RHO, EX, nf, d_y, b, c,false);
        
        //cal_rho(RHO, SIGMA, LAMX, a, b, d_y, nf);
        cal_delta(RHO, SIGMA, LAMX, a, b, nf, d_y, false);
        
        // calculate SIGMA
        //cal_sigma(SIGMA,EX, RHO, a, d_y, nf);
        cal_theta(SIGMA,EX, RHO, a, d_y, nf,false);
        
        cal_ex(EX, LAM, Y, PSI_INV, EXX, SIGMA, LAMX, O, LPL, d_y, s_n, nf);
        cal_z(logO, LOGVO, EX,  SIGMA,  RHO, LAMX,  O, nf, d_y, a, b, alpha, beta,false);
        //cal_o(logO,LOGVO,EX, SIGMA, RHO, LAMX, O, nf, d_y, a, b, alpha, beta);
    }
    
    if(x_method.compare("dense") == 0){
        cal_ex_simple(EX, EXX, LAM, PSI_INV, Y, s_n, nf, d_y);
    }
}
