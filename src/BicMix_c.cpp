#include <stdlib.h>
#include <ctime>
#include <stdio.h>
#include <math.h>
#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <time.h>

#include <R.h>
#include <Rmath.h>

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/LU>
#include "myHeader.cpp"
#include <Eigen/Core>

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <gsl/gsl_sf_psi.h>
#include <gsl/gsl_sf_gamma.h>
#include <unistd.h>
//#include <unsupported/Eigen/MatrixFunctions>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

//const double PI  =3.141592653589793238462;
/*
Chuan Gao C++
*/

/* 
usage
./bicluster_mixture_simul_up --nf 50 --y ./sim_data/Bicluster/Y_noise1_2.txt --out result --sep space --interval 100 --a 0.5 --b 0.5
./bicluster_mixture_simul_up --nf 50 --y /gpfs/fs0/data/engelhardtlab/cg148/data/CAP/S480RjQ2N_t_noCrossHyb.dat --out result --sep space --interval 50
./bicluster_mixture_simul_up --nf 100 --y /nfs/labs/engelhardtlab/cg148/data/hm3/HM3_gex_raw.t.txt --out result --sep tab
./bicluster_mixture_simul_up --nf 100 --y /nfs/labs/engelhardtlab/cg148/data/genotypes/hapmap2/hapmap.sfa --out result --sep space
*/



extern "C" void BicMix(double *Y_TMP_param ,int *nrow_param, int *ncol_param, double *a_param,double *b_param, int *nf_param, int *itr_param, double *LAM_out, double *EX_out, double *Z_out, double *O_out,double *EXX_out, double *PSI_out, int *nf_out, int *out_itr, char **output_dir,int *rsd, char **x_method, double *tol){
    
    Eigen::initParallel();
    
    double a = *a_param;
    double b = *b_param;
    int nf = *nf_param;
    int s_n = *nrow_param;
    int d_y = *ncol_param;
    int n_itr = *itr_param;
    
    int write_itr = *out_itr;
    
    int rsd_in = *rsd;
    
    double tol_in = *tol;
    
    string out_dir = *output_dir;
    std::replace( out_dir.begin(), out_dir.end(), '%', '/');
    
    stringstream ss;
    
    string x_method_in = *x_method;

    //cout << "a " << a << endl;
    
    
    //cout << "nf " <<  nf << endl;
    //cout << "s_n " << s_n << endl;
    //cout << "a " << a << endl;
    //cout << "a " << a << endl;
    
    
    double c=0.5,d=0.5,g=0.5,h=0.5,alpha=1,beta=1;
    
    int interval = 1000;
    
    MatrixXd Y=MatrixXd::Constant(s_n,d_y,0);
    
    for(int i = 0; i < d_y; i++){
        for(int j = 0; j < s_n; j++){
            Y(j,i) = Y_TMP_param[i*s_n+j];
            //cout << "y_i_j" << Y(j,i) << endl;
        }
    }
 
    // Declare variables independent of factor number to prepare for the EM algorithm
    
    long seed;
    //seed = time (NULL) * getpid();
    //seed = 1000;
    
    seed = (long)rsd;
    
    ss.str("");
    ss.clear();
    ss << out_dir << "/seed";
    ofstream f_seed (ss.str().c_str());
    f_seed << seed << endl;
    f_seed.close();
    

    VectorXd PSI=VectorXd::Constant(s_n,1);
    VectorXd PSI_INV=VectorXd::Constant(s_n,1);
   
    int nt=nf;
    
    MatrixXd LAM=MatrixXd::Constant(s_n,nt,0);
    MatrixXd THETA=MatrixXd::Constant(s_n,nf,1);
    MatrixXd DELTA=MatrixXd::Constant(s_n,nf,1);
    VectorXd PHI = VectorXd::Constant(nf,1);
    VectorXd TAU = VectorXd::Constant(nf,1);
    
    double nu = 1;
    double ETA = 1;
    double GAMMA = 1;
    
    VectorXd count_lam = VectorXd::Constant(nf,0);
    VectorXd index = VectorXd::Constant(nf,0);
    
    double nmix=2;
    double zi = double(1)/nmix;
    MatrixXd Z = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logZ = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGV = MatrixXd::Constant(nmix,1,log(zi));
    
    
    // fill in the lambda matrix
    
    init_lam(LAM, s_n, nf, seed);
    
    //VectorXd lam_count_v = VectorXd::Constant(n_itr,0);
    
    //declare and initialize parameters related to X
    double XI=1,VARSIG=1,OMEGA=1;
    
    MatrixXd EX=MatrixXd::Constant(nt,d_y,0);
    //MatrixXd TEX=MatrixXd::Constant(d_y,nt,0);
    MatrixXd EXX=MatrixXd::Constant(nt,nt,0);
    
    VectorXd KAPPA = VectorXd::Constant(nf,1);
    VectorXd LAMX = VectorXd::Constant(nf,1);
    MatrixXd RHO=MatrixXd::Constant(nf,d_y,1);
    MatrixXd SIGMA=MatrixXd::Constant(nf,d_y,1);
    
    VectorXd count_x = VectorXd::Constant(nf,0);
    
    MatrixXd O = MatrixXd::Constant(nmix,nf,zi);
    MatrixXd logO = MatrixXd::Constant(nmix,nf,log(zi));
    MatrixXd LOGVO = MatrixXd::Constant(nmix,1,log(zi));
    
    MatrixXd LPL = MatrixXd::Constant(nf,nf,0);
    //MatrixXd vLXL = MatrixXd::Constant(s_n,s_n,0);
    //MatrixXd partR = MatrixXd::Constant(nf,d_y,0);
    //MatrixXd partL = MatrixXd::Constant(s_n,nf,0);
    
    // fill in the EX matrix
    init_ex(EX, d_y, nf, seed);
    
    EXX=EX*EX.transpose();
    
    //VectorXd x_count_v = VectorXd::Constant(n_itr,0);
    //MatrixXd LAM_T=LAM.transpose();
    
    LPL.setZero();
    for(int i=0;i<s_n;i++){
        LPL += LAM.transpose().col(i)*PSI_INV(i)*LAM.row(i);
    }
    
    VectorXd det_psi = VectorXd::Constant(n_itr,0);
    
    //cout << "You passed method " << method << endl;
    
    string lam_method = "matrix";
    
    for(int itr=0;itr<n_itr;itr++){
        
        cal_lam_all(LAM, Y,EX,PSI_INV,EXX,Z,LPL,THETA,PHI, s_n, d_y, nf,
                    a, b, c, d, g, h, GAMMA, ETA, nu, TAU, DELTA, alpha, beta, lam_method);
        
        cal_z(logZ, LOGV, LAM,  THETA,  DELTA, PHI,  Z, nf, s_n, a, b, alpha, beta,true);
        
        //if(x_method.compare("bicmix") == 0){
        //cout << "You passed bicmix" << endl;
        cal_ex_all( LAM,  Y, EX, PSI_INV, EXX, O, LPL, SIGMA, LAMX,  s_n,  d_y,  nf,
                   a,  b,  c,  d,  g,  h,  VARSIG,  OMEGA,  XI,  KAPPA,  RHO,  logO,  LOGVO,  alpha,  beta, x_method_in);
        //}
        //if(method.compare("sfamix") == 0){
        //cal_ex_simple(EX, EXX, LAM, PSI_INV, Y, s_n, nf, d_y);
        //}
        
        
        
        // count the number of non-zero values in each row of the x matrix
        int count_nonzero = 0;
        count_nonzero = count_nonzero_lam_ex(count_x, count_lam, index, LAM, EX, PHI, LAMX, nf, s_n, d_y);
        
        // remove factors, loadings that are exclusively zeros, and assign to new matrix
        if(count_nonzero != nf){
            nf=count_nonzero;
            nt=nf;
            //red_dim(EX, EXX, LAM, THETA, DELTA, PHI, TAU, Z, logZ, count_lam, index, KAPPA, LAMX,
            //        RHO, SIGMA, count_x, O, logO, LPL, partR, partL, nf, nt, s_n, d_y, nmix, zi);
            red_dim(EX, EXX, LAM, THETA, DELTA, PHI, TAU, Z, logZ, count_lam, index, KAPPA, LAMX,
                    RHO, SIGMA, count_x, O, logO, LPL, nf, nt, s_n, d_y, nmix, zi,true);
            //red_dim(EX_merge, EXX_merge, LAM_merge, THETA_merge, DELTA_merge, PHI_merge, TAU_merge, Z_merge, logZ_merge, count_lam, index, KAPPA_merge, LAMX_merge,
            //       RHO_merge, SIGMA_merge, count_x, O_merge, logO_merge, LPL_merge, nfcov, nt, s_n, d_y, nmix, zi,false);
            
        }
        
        MatrixXd LX=LAM*EX;
        cal_psi(PSI, LX, Y, LAM, EXX, s_n, d_y);
        inv_psi_vec(PSI,PSI_INV,s_n);
        
        for(int i=0;i<s_n;i++){
            det_psi(itr) = det_psi(itr)+log(PSI(i,i));
        }
        
        if(itr%10==0){
            cout << "itr " << itr << endl;
            
            cout << "number of factors " << nf << endl;
            cout << "count_lam" << endl << count_lam.transpose() << endl;
            //cout << "number of betas " << ncov << endl;
            //cout << "count_beta" << endl << count_lam_cov.transpose() << endl;
            cout << "count_x" << endl << count_x.transpose() << endl;
        }
        if(itr>10){
            if(abs(det_psi(itr) - det_psi(itr-1)) < tol_in){
                //itr_at = itr;
                //write_final_beta(out_dir, LAM_cov, PSI, itr, seed);
                //write_final_hidden(out_dir, LAM, Z, EX, EXX, O, PSI, lam_count_v, LAMX, PHI, itr, seed);
                string sin = out_dir + "/itr";
                write_file <int> (itr,sin);
                
                // write LAM
                sin = out_dir + "/LAM";
                write_file <MatrixXd> (LAM,sin);
                
                sin = out_dir + "/Z";
                write_file <MatrixXd> (Z,sin);
                
                sin = out_dir + "/EX";
                write_file <MatrixXd> (EX,sin);
                
                sin = out_dir + "/EXX";
                write_file <MatrixXd> (EXX,sin);
                
                sin = out_dir + "/PSI";
                write_file <VectorXd> (PSI,sin);
                
                break;
            }
        }
        
    }
    
    
    convert_results(LAM,EX,EXX,Z,O,PSI,nf,LAM_out,EX_out,EXX_out, Z_out,O_out,PSI_out, nf_out,s_n,d_y);
    
    
}

