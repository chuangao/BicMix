
#include <stdlib.h>
#include <gsl/gsl_sf_gamma.h>
#include <math.h>

#include <Eigen/Dense>
#include <limits>
#include <stdio.h>

#include <fstream>
#include <ostream>
#include <string>
#include <sstream>
#include <iostream>

#include "gig_par.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include <float.h>
#include <gsl/gsl_sf_psi.h>

using namespace Eigen;
using namespace std;

using Eigen::MatrixXd;

const double PI_mine  =3.141592653589793238462;

void fInvert(MatrixXd&X,MatrixXd& M1,MatrixXd& M2,int n){
	X.setZero();
	int j=0;
	//MatrixXd SUM=MatrixXd::Constant(k,n,0);
	for(int j=0;j<n;j++){
		for(int i=0;i<n;i++){
			//cout << "i " << i << endl;
			if(i==0){
				X(i,j)=M2(i,j)/M1(i,i);
			}
			else{
				double sum=0;				
				for(int k=0;k<i;k++){
					//cout << "k " << k << endl;
					sum += M1(i,k)*X(k,j);
					X(i,j)=(M2(i,j)-sum)/M1(i,i);
				}
			}
		}
	}
}
void range_colwise(MatrixXd& range,MatrixXd& x,int n,int m){
	int i,j;
	for(i=0;i<m;i++){
		//cout << "i" << i << endl;
		int s;
		double min=1e100,max=1e-100;
		// for(s=0;s<n;s++){
		// 	if(x(s,i)!=0){
		// 		min=abs(x(s,i));
		// 		max=abs(x(s,i));
		// 		continue;
		// 	}
		// 	//cout << "s" << s << endl;
		// }
		for(j=0;j<n;j++){
			if(x(j,i)!=0){
				if(min>abs(x(j,i))){
					min=abs(x(j,i));
				}
				if(max<abs(x(j,i))){
					max=abs(x(j,i));
				}
			}
		}
		range(0,i)=min;
		range(1,i)=max;
	}
	//cout << "range_coleise_finished" << endl;
}

void range_rowwise(MatrixXd& range,MatrixXd& x,int n,int m){
	int i,j;
	for(j=0;j<n;j++){
		int s;
		double min=1e100,max=1e-100;
		// for(s=0;s<m;s++){
		// 	if(x(j,s)!=0){
		// 		min=abs(x(j,s));
		// 		max=abs(x(j,s));
		// 		continue;
		// 	}
			
		// }
		for(i=0;i<m;i++){
			if(x(j,i)!=0){
				if(min>abs(x(j,i))){
					min=abs(x(j,i));
				}
				if(max<abs(x(j,i))){
					max=abs(x(j,i));
				}
			}
		}
		range(0,j)=min;
		range(1,j)=max;
	}
}

// void range_colwise(MatrixXd& range,MatrixXd& x,int n,int m){
// 	int i,j;
// 	for(i=0;i<m;i++){
// 		//cout << "i" << i << endl;
// 		int s;
// 		double min=0,max=0;
// 		// for(s=0;s<n;s++){
// 		// 	if(x(s,i)!=0){
// 		// 		min=abs(x(s,i));
// 		// 		max=abs(x(s,i));
// 		// 		break;
// 		// 	}
// 		// 	cout << "s" << s << endl;
// 		// }
// 		for(j=0;j<n;j++){
// 			if(x(j,i)!=0){
// 				if(min>abs(x(j,i))){
// 					min=abs(x(j,i));
// 				}
// 				if(max<abs(x(j,i))){
// 					max=abs(x(j,i));
// 				}
// 			}
// 		}
// 		range(0,i)=min;
// 		range(1,i)=max;
// 	}
// 	//cout << "range_coleise_finished" << endl;
// }

// void range_rowwise(MatrixXd& range,MatrixXd& x,int n,int m){
// 	int i,j;
// 	for(j=0;j<n;j++){
// 		int s;
// 		double min=0,max=0;
// 		// for(s=0;s<m;s++){
// 		// 	if(x(j,s)!=0){
// 		// 		min=abs(x(j,s));
// 		// 		max=abs(x(j,s));
// 		// 		break;
// 		// 	}
			
// 		// }
// 		for(i=0;i<m;i++){
// 			if(x(j,i)!=0){
// 				if(min>abs(x(j,i))){
// 					min=abs(x(j,i));
// 				}
// 				if(max<abs(x(j,i))){
// 					max=abs(x(j,i));
// 				}
// 			}
// 		}
// 		range(0,j)=min;
// 		range(1,j)=max;
// 	}
// }

void inv_psi(MatrixXd& psi,MatrixXd& psi_inv,int n){
    for(int i=0;i<n;i++){
        psi_inv(i,i)=double(1)/psi(i,i);
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

void set_thresh_matrix(MatrixXd& M,int n,int p){
    for(int i=0;i<n;i++){
		for(int j=0;j<p;j++){
			if(M(i,j)<0.0001){
				M(i,j)=0.0001;
			}
		}
    }
}
void set_thresh_vector(VectorXd& M,int n){
    for(int i=0;i<n;i++){
		if(M(i)<0.0001){
			M(i)=0.0001;
		}
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



/*
void draw_ex(MatrixXd &LAM,MatrixXd &PSI_INV,MatrixXd &Y, MatrixXd &O, MatrixXd &SIGMA, VectorXd &LAMX,MatrixXd &EX, int nf, int s_n, int d_y, gsl_rng *r ){
	MatrixXd LP=MatrixXd::Constant(nf,s_n,0);
	MatrixXd partR=MatrixXd::Constant(nf,d_y,0);
	MatrixXd LPL=MatrixXd::Constant(nf,nf,0);
	MatrixXd IDNF = MatrixXd::Identity(nf, nf);
	MatrixXd partV = MatrixXd::Constant(nf, nf,0);
	MatrixXd P = MatrixXd::Constant(nf, nf,0);
	VectorXd D = VectorXd::Constant(nf,0);
	MatrixXd D_inv = MatrixXd::Constant(nf, nf,0);
	MatrixXd Ll = MatrixXd::Constant(nf, nf,0);
	MatrixXd Ll_inv = MatrixXd::Constant(nf,nf,0);
	MatrixXd L_inv = MatrixXd::Constant(nf, nf,0);
	MatrixXd vx = MatrixXd::Constant(nf, nf,0);
	MatrixXd X2I = MatrixXd::Constant(nf,1,0);
	MatrixXd EXI=MatrixXd::Constant(nf, 1,0);
	
	LP=LAM.transpose()*PSI_INV;
	partR=LP*Y;
	LPL=LP*LAM;

	for(int j=0;j<d_y;j++){
		partV.setZero();
		for(int i1=0;i1<nf;i1++){
			for(int i2=0;i2<nf;i2++){
				partV(i1,i2) = LPL(i1,i2);
			}
			if(O(0,i1)==1){
				partV(i1,i1) += (double)1/SIGMA(i1,j);
			}else{
				partV(i1,i1) += (double)1/LAMX(i1);
			}
		}
			
		LDLT<MatrixXd> ldltOfA(partV);
		P = PermutationMatrix<Dynamic,Dynamic>(ldltOfA.transpositionsP());
		D = ldltOfA.vectorD();
		D = double(1)/D.array().sqrt();
		D_inv = D.asDiagonal();
		Ll = ldltOfA.matrixL();
	
		fInvert(Ll_inv,Ll,IDNF,nf);
		L_inv=P.transpose()*Ll_inv.transpose()*D_inv;
		vx=L_inv*L_inv.transpose();
	
		for(int i1=0;i1<nf;i1++){
			X2I(i1,0)=gsl_ran_gaussian (r, 1);
		}
		EXI=vx*partR.col(j)+L_inv*X2I;				
		EX.col(j)=EXI.col(0);	
	}
			
}

void draw_sigma(MatrixXd &SIGMA,MatrixXd &EX, MatrixXd &RHO,int nf,int d_y,int a,gig *mygig){
	for(int j=0;j<nf;j++){
		//if(O(0,j)==1){
		for(int i=0;i<d_y;i++){
			double sig=0;
			double pg,ag,bg;
			pg=a-0.5;
			ag=EX(j,i)*EX(j,i);
			bg=2*RHO(i);
			mygig->rgig(sig,&pg,&bg,&ag);
			SIGMA(j,i)=sig;
			//if(SIGMA(j,i)<1e-15){SIGMA(j,i)=1e-15;}
		}
		//}
	}
}
void draw_rho(MatrixXd &RHO,MatrixXd &SIGMA, VectorXd &LAMX, double a, double b, int nf, int d_y, gsl_rng *r){
	for(int j=0;j<nf;j++){
		for(int i=0;i<d_y;i++){
			RHO(j,i)=gsl_ran_gamma(r,a+b,double(1)/(SIGMA(j,i)+LAMX(j)));
		}
				
	}
			
}

void draw_lamx(VectorXd &LAMX,MatrixXd &O, MatrixXd &RHO, VectorXd &KAPPA,MatrixXd &EX, MatrixXd &SIGMA, gsl_rng *r, gig *mygig, int nf, int d_y, double b, double c){
	for(int i=0;i<nf;i++){
		if(O(0,i)==1){
			LAMX(i)=gsl_ran_gamma(r,d_y*b+c,double(1)/(RHO.row(i).sum()+KAPPA(i)));
		}
		else{
			double lamx;
			double pg,ag,bg;
			pg=c-0.5*d_y;
			ag=EX.row(i).dot(EX.row(i));
			bg=2*KAPPA(i);
			mygig->rgig(lamx,&pg,&bg,&ag);
			LAMX(i)=lamx;
						
			for(int j=0;j<d_y;j++){
				SIGMA(i,j)=LAMX(i);
			}
					
		}
		//if(LAMX(i)<1e-10){LAMX=1e-10;}
	}
}
void draw_kappa(VectorXd &KAPPA, VectorXd &LAMX, double OMEGA, double c, double d, int nf, gsl_rng *r){
	for(int i=0;i<nf;i++){
		KAPPA(i)=gsl_ran_gamma(r,c+d,double(1)/(LAMX(i)+OMEGA));
	}
}

void draw_o(MatrixXd &O, MatrixXd &LOGVO, MatrixXd &EX, MatrixXd &SIGMA,VectorXd &LAMX, MatrixXd &RHO, int nf, int d_y,double a, double b,gsl_rng *r){
	MatrixXd logO = MatrixXd::Constant(2,nf,0);
	for(int i=0;i<nf;i++){
		logO(0,i)=LOGVO(0,0);
		logO(1,i)=LOGVO(1,0);
		for(int j=0;j<d_y;j++){
			logO(0,i) = logO(0,i)+log_norm(EX(i,j),0,SIGMA(i,j))+lgamma(a+b)-(a+b)*log(LAMX(i)+SIGMA(i,j))+b*log(LAMX(i))+(a-1)*log(SIGMA(i,j))-lgamma(a)-lgamma(b);
			logO(1,i) = logO(1,i)+log_norm(EX(i,j),0,LAMX(i));
		}
	}
	for(int i=0;i<nf;i++){
		O(0,i)= gsl_ran_bernoulli(r,double(1)/(1+exp(logO(1,i)-logO(0,i))));
		O(1,i)=1-O(0,i);
								
		if(O(0,i)==0){
			for(int j=0;j<d_y;j++){
				SIGMA(i,j)=LAMX(i);
			}
		}
		else{
			for(int j=0;j<d_y;j++){
				RHO(i,j)=gsl_ran_gamma(r,a+b,double(1)/(SIGMA(i,j)+LAMX(i)));
			}
		}
				
	}
}
void draw_v(MatrixXd &LOGVO,MatrixXd &O,double alpha, double beta, int nf, gsl_rng *r){
	double ps1=alpha;
	double ps2=beta;
        
	for(int i=0;i<nf;i++){
		ps1=ps1+O(0,i);
		ps2=ps2+O(1,i);
	}
	LOGVO(0,0)=log(gsl_ran_beta(r,ps1,ps2));
	LOGVO(1,0)=log(1-exp(LOGVO(0,0)));
						
}
void draw_lam(MatrixXd &LAM, MatrixXd &PSI_INV, MatrixXd &Y,MatrixXd &EX, MatrixXd &Z, MatrixXd &THETA, VectorXd &PHI, int s_n, int nf, gsl_rng *r){
	MatrixXd partL=PSI_INV*Y*EX.transpose();
	MatrixXd EXX=EX*EX.transpose();
	MatrixXd IDNF = MatrixXd::Identity(nf, nf);
	MatrixXd partV = MatrixXd::Constant(nf,nf,0);
	for(int j=0;j<s_n;j++){		
		partV.setZero();	
		for(int i1=0;i1<nf;i1++){
			for(int i2=0;i2<nf;i2++){
				partV(i1,i2) = PSI_INV(j,j)*EXX(i1,i2);
			}
			if(Z(0,i1)==1){
				partV(i1,i1) += double(1)/THETA(j,i1);
			}
			else if(Z(1,i1)==1){
				partV(i1,i1) += double(1)/PHI(i1);
			}
		}  
								
		LDLT<MatrixXd> ldltOfA(partV);
		
		MatrixXd P = PermutationMatrix<Dynamic,Dynamic>(ldltOfA.transpositionsP());
		VectorXd D = ldltOfA.vectorD();
		D = double(1)/D.array().sqrt();
		MatrixXd D_inv = D.asDiagonal();
				
		MatrixXd Ll = ldltOfA.matrixL();
		MatrixXd Ll_inv = MatrixXd::Constant(nf,nf,0);
		fInvert(Ll_inv,Ll,IDNF,nf);
		//MatrixXd L_inv=P.transpose()*Ll_inv.transpose()*D_inv;
		MatrixXd L_inv=P.transpose()*Ll_inv.transpose()*D_inv;
		MatrixXd vl=L_inv*L_inv.transpose();
		MatrixXd LAM2I = MatrixXd::Constant(1,nf,0);
		for(int i1=0;i1<nf;i1++){
			LAM2I(0,i1)=gsl_ran_gaussian (r, 1);
		}
		MatrixXd LAMI=partL.row(j)*vl+LAM2I*(L_inv.transpose());				
		LAM.row(j)=LAMI.row(0);			

	}

}

void draw_theta(MatrixXd& THETA,MatrixXd& LAM, MatrixXd& DELTA, double a, int nf, int s_n,gig *mygig){
	for(int j=0;j<nf;j++){	
		//if(Z(0,j)==1){
		for(int i=0;i<s_n;i++){					
			double the=0;
			double pg,ag,bg;
			pg=a-0.5;
			ag=LAM(i,j)*LAM(i,j);
			bg=2*DELTA(i,j);
			mygig->rgig(the,&pg,&bg,&ag);
			THETA(i,j)=the;
			//if(THETA(i,j)<1e-15){THETA(i,j)=1e-15;}
		}
		//}
				
	}

}

void draw_delta(MatrixXd &DELTA,MatrixXd &THETA,VectorXd &PHI,double a, double b,int nf, int s_n, gsl_rng *r){
	for(int j=0;j<nf;j++){			
		for(int i=0;i<s_n;i++){	
			DELTA(i,j)=gsl_ran_gamma(r,a+b,double(1)/(THETA(i,j)+PHI(j)));
		}
				
	}


}

void draw_phi(VectorXd &PHI,MatrixXd &Z, MatrixXd &DELTA,VectorXd &TAU,MatrixXd &LAM, MatrixXd &THETA, double b, double c, int nf, int s_n, gsl_rng *r, gig *mygig){
	for(int i=0;i<nf;i++){
		if(Z(0,i)==1){
			PHI(i)=gsl_ran_gamma(r,s_n*b+c,double(1)/(DELTA.col(i).sum()+TAU(i)));
					
		}
		else if(Z(1,i)==1){
			double phi=0;
			double pg,ag,bg;
			pg=c-0.5*s_n;
			ag=LAM.col(i).dot(LAM.col(i));
			bg=2*TAU(i);
			mygig->rgig(phi,&pg,&bg,&ag);
			PHI(i)=phi;
					
			for(int j=0;j<s_n;j++){
				THETA(j,i)=PHI(i);
			}
					
		}
	}
}
void draw_tau(VectorXd &TAU, VectorXd &PHI, double ETA, double c, double d, int nf, gsl_rng *r){
	for(int i=0;i<nf;i++){
		TAU(i)=gsl_ran_gamma(r,c+d,double(1)/(PHI(i)+ETA));
	}
}

void draw_z(MatrixXd &Z, MatrixXd &LOGV, MatrixXd &LAM,MatrixXd &THETA, VectorXd &PHI, MatrixXd &DELTA, double a, double b, int s_n, int nf, gsl_rng *r){
	MatrixXd logZ = MatrixXd::Constant(2,nf,0);
	for(int i=0;i<nf;i++){
		logZ(0,i)=LOGV(0,0);
		logZ(1,i)=LOGV(1,0);
		for(int j=0;j<s_n;j++){
			//lgamma(a+b)-(a+b)*log(PHI[i]+THETA[j,i])+b*log(PHI[i])+(a-1)*log(THETA[j,i])-lgamma(a)-lgamma(b);
			logZ(0,i) = logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+lgamma(a+b)-(a+b)*log(PHI(i)+THETA(j,i))+b*log(PHI(i))+(a-1)*log(THETA(j,i))-lgamma(a)-lgamma(b);
			logZ(1,i) = logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
				
		}
	}
	for(int i=0;i<nf;i++){
		//cout << "pi " << double(1)/(1+exp(logZ(1,i)-logZ(0,i))) << endl;
		Z(0,i)= gsl_ran_bernoulli(r,double(1)/(1+exp(logZ(1,i)-logZ(0,i))));
		Z(1,i)=1-Z(0,i);

	
		if(Z(0,i)==0){
			for(int j=0;j<s_n;j++){
				THETA(j,i)=PHI(i);
			}
		}
		else{
			for(int j=0;j<s_n;j++){
				DELTA(j,i)=gsl_ran_gamma(r,a+b,double(1)/(THETA(j,i)+PHI(i)));
			}
		}
					
	}			
}

void draw_psi(MatrixXd &PSI,MatrixXd &Y,MatrixXd &LAM,MatrixXd &EX, int s_n, int d_y,gsl_rng *r){

	MatrixXd YLX=Y-LAM*EX;
	for(int i=0;i<s_n;i++){
		PSI(i,i)=double(1)/gsl_ran_gamma(r,double(d_y)/2+1,double(1)/(0.5*(YLX.row(i).dot(YLX.row(i)))+1));
	}
	//PSI=PSI/d_y;

			
}

void rescale(MatrixXd &EX, MatrixXd &LAM, MatrixXd &SIGMA, MatrixXd &RHO, VectorXd &LAMX, VectorXd &KAPPA, MatrixXd &THETA, MatrixXd &DELTA, VectorXd &PHI, VectorXd &TAU,int nf, int d_y, int s_n){
	MatrixXd vx_m = MatrixXd::Constant(2,nf,0);
	for(int i=0;i<nf;i++){
			
		double ex2=0,e2x=0,vx=0,nv=0;			
		for(int j=0;j<d_y;j++){
			if(EX(i,j)!=0){
				ex2 += EX(i,j)*EX(i,j);
				e2x += EX(i,j);
				nv++;
			}
		}

		if(nv==0){
			continue;
		}
		if(nv==1){
			vx=ex2;
			vx_m(0,i)=vx;
		}else{
			vx=(ex2/nv-(e2x/nv)*(e2x/nv));
			vx_m(0,i)=vx;
		}
			
		for(int j=0;j<d_y;j++){
			EX(i,j)=EX(i,j)/sqrt(vx_m(0,i));
			SIGMA(i,j)=SIGMA(i,j)/vx_m(0,i);
			RHO(i,j)=RHO(i,j)*vx_m(0,i);
		}
		LAMX(i)=LAMX(i)/vx_m(0,i);
		KAPPA(i)=KAPPA(i)*vx_m(0,i);
		for(int j=0;j<s_n;j++){
			LAM(j,i)=LAM(j,i)*sqrt(vx_m(0,i));
			THETA(j,i)=THETA(j,i)*vx_m(0,i);
			DELTA(j,i)=DELTA(j,i)/vx_m(0,i);
		}
		TAU(i)=TAU(i)*vx_m(0,i);
		PHI(i)=PHI(i)/vx_m(0,i);
	}
		
}

*/


void cal_lamx(VectorXd& LAMX, MatrixXd& O, VectorXd& KAPPA, MatrixXd& RHO,MatrixXd& EX, int nf, int d_y, double b, double c){
	for(int i=0;i<nf;i++){
		double x_sum=0;
		double sum_c=d_y*b*O(0,i)+c-1-0.5*d_y*O(1,i);
		double at = 2*(KAPPA(i)+O(0,i)*(RHO.row(i).sum()));
			
		//x_sum=EXX(i,i);
			
		x_sum=EX.row(i).dot(EX.row(i));
		
		double bt = O(1,i)*x_sum;
		LAMX(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
		
		 if(LAMX(i)<1e-10){
		 	LAMX(i)=1e-10;
		 }
		
	}
}

/*
void cal_lamx_VXM(VectorXd& LAMX, MatrixXd& O, VectorXd& KAPPA, MatrixXd& RHO,MatrixXd& EX, MatrixXd& VXM, int nf, int d_y, double b, double c){
	for(int i=0;i<nf;i++){
		double x_sum=0;
		double sum_c=d_y*b*O(0,i)+c-1-0.5*d_y*O(1,i);
		double at = 2*(KAPPA(i)+O(0,i)*(RHO.row(i).sum()));
			
		//x_sum=EXX(i,i);
			
		x_sum=EX.row(i).dot(EX.row(i))+VXM.col(i).sum();
		
		double bt = O(1,i)*x_sum;
		LAMX(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
		
		 if(LAMX(i)<1e-10){
		 	LAMX(i)=1e-10;
		 }
		
	}
}
*/
void cal_phi(VectorXd& PHI, MatrixXd& Z, VectorXd& TAU, MatrixXd& DELTA,MatrixXd& LAM, int nf, int s_n, double b, double c){
	for(int i=0;i<nf;i++){
		double lam_sum=0;
		double sum_c=s_n*b*Z(0,i)+c-1-0.5*s_n*Z(1,i);
		double at = 2*(TAU(i)+Z(0,i)*(DELTA.col(i).sum()));			
		
		lam_sum=LAM.col(i).dot(LAM.col(i));
			
			
		double bt = Z(1,i)*lam_sum;
		PHI(i)=double(sum_c+sqrt(sum_c*sum_c+at*bt))/at;
		
		 if(PHI(i)<1e-10){
		 	PHI(i)=1e-10;
		 }
		
	}
}

void red_dim(MatrixXd& EX,MatrixXd& EXX, MatrixXd& ELL, MatrixXd& LAM, MatrixXd& THETA, MatrixXd& DELTA, VectorXd& PHI, VectorXd& TAU, MatrixXd& Z, MatrixXd& logZ, VectorXd& count_lam, MatrixXd& count_lam_M, VectorXd& index, VectorXd& KAPPA, VectorXd& LAMX, MatrixXd& RHO, MatrixXd& SIGMA, VectorXd& count_x, MatrixXd& count_x_M, MatrixXd& O, MatrixXd& logO, MatrixXd& LPL, MatrixXd& partR, MatrixXd& partL, int nf, int s_n, int d_y, int n_itr){
	int nt=nf;
	int nmix=2;
	double zi=0.5;
	
	MatrixXd EX2=MatrixXd::Constant(nt,d_y,0);
	//MatrixXd TEX2=MatrixXd::Constant(d_y,nt,0);
	//MatrixXd VX2=MatrixXd::Constant(nt,nt,0);
	MatrixXd EXX2=MatrixXd::Constant(nt,nt,0);
	MatrixXd ELL2=MatrixXd::Constant(nt,nt,0);
            
	MatrixXd LAM2=MatrixXd::Constant(s_n,nt,0);
            
	MatrixXd THETA2=MatrixXd::Constant(s_n,nf,0);
	MatrixXd DELTA2=MatrixXd::Constant(s_n,nf,0);
	VectorXd PHI2 = VectorXd::Constant(nf,0);
	VectorXd TAU2 = VectorXd::Constant(nf,0);
	MatrixXd Z2 = MatrixXd::Constant(nmix,nf,zi);
	MatrixXd logZ2 = MatrixXd::Constant(nmix,nf,log(zi));
	VectorXd count_lam2 = VectorXd::Constant(nf,0);
	MatrixXd count_lam_M2 = MatrixXd::Constant(n_itr,nf,0);
	VectorXd index2 = VectorXd::Constant(nf,0);
	//ID2.diagonal()=id_v2;
            
	// for EX related10
	VectorXd KAPPA2 = VectorXd::Constant(nf,0);
	VectorXd LAMX2 = VectorXd::Constant(nf,0);
	MatrixXd RHO2=MatrixXd::Constant(nf,d_y,0);
	MatrixXd SIGMA2=MatrixXd::Constant(nf,d_y,0);
            
	VectorXd count_x2 = VectorXd::Constant(nf,0);
	MatrixXd count_x_M2 = MatrixXd::Constant(n_itr,nf,0);
	MatrixXd O2 = MatrixXd::Constant(nmix,nf,zi);
	MatrixXd logO2 = MatrixXd::Constant(nmix,nf,log(zi));
	//MatrixXd LOGVO2 = MatrixXd::Constant(nmix,1,log(0.5));
            
	MatrixXd LPL2 = MatrixXd::Constant(nf,nf,0);
	MatrixXd partR2 = MatrixXd::Constant(nf,d_y,0);
	MatrixXd partL2 = MatrixXd::Constant(s_n,nf,0);

	// MatrixXd VXM2=MatrixXd::Constant(d_y,nf,0); 
	// MatrixXd VLM2=MatrixXd::Constant(s_n,nf,0);
	
 
	for(int i=0;i<nf;i++){
		
		//vx_m2.col(i)=vx_m.col(index(i));
		
		// VLM2.col(i)=VLM.col(index(i));
		// VXM2.col(i)=VXM.col(index(i));
		
		EX2.row(i)=EX.row(index(i));
		//TEX2=EX2.transpose();
		for(int j=0;j<nf;j++){
			EXX2(i,j)=EXX(index(i),index(j));
			ELL2(i,j)=ELL(index(i),index(j));
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
		count_lam_M2.col(i)=count_lam_M.col(index(i));
		index2(i)=index(i);
                
		// EX related
		KAPPA2(i)=KAPPA(index(i));
		LAMX2(i)=LAMX(index(i));
		RHO2.row(i)=RHO.row(index(i));
		SIGMA2.row(i)=SIGMA.row(index(i));
		count_x2(i)=count_x(index(i));
		count_x_M2.col(i)=count_x_M.col(index(i));
		O2.col(i)=O.col(index(i));
		logO2.col(i)=logO.col(index(i));
                
		partR2.row(i)=partR.row(index(i));
		partL2.col(i)=partL.col(index(i));
                
	}
         
	// Assign the new parameters back
			
	//vx_m=vx_m2;
	
	  // VLM=VLM2;
	  // VXM=VXM2;
	
	EX=EX2;
	//TEX=TEX2;
	//VX=VX2;
	EXX=EXX2;
	//LP=LP2;
	//ID=ID2;
	LAM=LAM2;
	ELL=ELL2;
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
	count_lam_M=count_lam_M2;
	index=index2;
	// EX related
	KAPPA=KAPPA2;
	LAMX=LAMX2;
	RHO=RHO2;
	SIGMA=SIGMA2;
	count_x=count_x2;
	count_x_M=count_x_M2;
	O=O2;
	logO=logO2;
            
	LPL=LPL2;
	partR=partR2;
	partL=partL2;
	//partV=partV2;
}

void red_dim_cov(MatrixXd& EX,MatrixXd& EXX, MatrixXd& LAM, MatrixXd& THETA, MatrixXd& DELTA, VectorXd& PHI, VectorXd& TAU, MatrixXd& Z, MatrixXd& logZ, VectorXd& count_lam, MatrixXd& count_lam_M, VectorXd& index, VectorXd& index_cov_output, MatrixXd& LPL, MatrixXd& partL, int nf, int s_n, int d_y, int n_itr){
	int nt=nf;
	int nmix=2;
	double zi=0.5;
	
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
	MatrixXd count_lam_M2 = MatrixXd::Constant(n_itr,nf,0);
	VectorXd index2 = VectorXd::Constant(nf,0);
	VectorXd index_cov_output2 = VectorXd::Constant(nf,0);
	//ID2.diagonal()=id_v2;
            

	MatrixXd LPL2 = MatrixXd::Constant(nf,nf,0);
	MatrixXd partL2 = MatrixXd::Constant(s_n,nf,0);
 
	for(int i=0;i<nf;i++){
		
		//vx_m2.col(i)=vx_m.col(index(i));
		/*
		  VLM2.col(i)=VLM.col(index(i));
		  VXM2.col(i)=VXM.col(index(i));
		*/
		EX2.row(i)=EX.row(index(i));
		//TEX2=EX2.transpose();
		for(int j=0;j<nf;j++){
			EXX2(i,j)=EXX(index(i),index(j));
			//ELL2(i,j)=ELL(index(i),index(j));
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
		count_lam_M2.col(i)=count_lam_M.col(index(i));
		index2(i)=index(i);
		index_cov_output2(i)=index_cov_output(index(i));
                
		// EX related
	
		partL2.col(i)=partL.col(index(i));
                
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
	count_lam_M=count_lam_M2;
	index=index2;
	index_cov_output=index_cov_output2;
	// EX related

	LPL=LPL2;
	//partR=partR2;
	partL=partL2;
	//partV=partV2;
}

void cal_rho(MatrixXd& RHO, MatrixXd& SIGMA, VectorXd& LAMX, double a, double b, int d_y, int nf_nc){
	
	for(int i=0;i<d_y;i++){
		for(int j=0;j<nf_nc;j++){
			RHO(j,i)=double((a+b))/(SIGMA(j,i)+LAMX(j));
		}
	}
}

void cal_sigma(MatrixXd& SIGMA,MatrixXd& EX, MatrixXd& RHO, double a, int d_y, int nf_nc){
	
	for(int i=0;i<d_y;i++){
		for(int j=0;j<nf_nc;j++){
			double a23=(2*a-3);
			SIGMA(j,i)=double(a23+sqrt(a23*a23+8*(EX(j,i)*EX(j,i))*RHO(j,i)))/4/RHO(j,i);
		}
	}
}

/*
void cal_sigma_VXM(MatrixXd& SIGMA,MatrixXd& EX, MatrixXd& RHO, MatrixXd& VXM, double a, int d_y, int nf_nc){
	
	for(int i=0;i<d_y;i++){
		for(int j=0;j<nf_nc;j++){
			double a23=(2*a-3);
			SIGMA(j,i)=double(a23+sqrt(a23*a23+8*(EX(j,i)*EX(j,i)+VXM(i,j))*RHO(j,i)))/4/RHO(j,i);
		}
	}
}
*/

void cal_delta(MatrixXd& DELTA,MatrixXd& THETA,VectorXd& PHI,double a, double b,int nf, int s_n){
	for(int i=0;i<s_n;i++){
		for(int j=0;j<nf;j++){
			DELTA(i,j)=double((a+b))/(THETA(i,j)+PHI(j));
		}
	}	
}

void cal_theta(MatrixXd& THETA,MatrixXd& LAM,MatrixXd& DELTA,double a, int s_n, int nf){
	for(int i=0;i<s_n;i++){
		for(int j=0;j<nf;j++){
			double a23=(2*a-3);
			THETA(i,j)=double(a23+sqrt(a23*a23+8*(LAM(i,j)*LAM(i,j))*DELTA(i,j)))/4/DELTA(i,j);				
		}
	}
}
/*
void draw_ex(MatrixXd &LAM,MatrixXd &PSI_INV,MatrixXd &Y, MatrixXd &partV, MatrixXd &O, MatrixXd &SIGMA, VectorXd &LAMX,MatrixXd &Ll_inv,MatrixXd &X2I, MatrixXd &EX, int nf, int d_y, gsl_rng *r ){
	MatrixXd LP=LAM.transpose()*PSI_INV;
	MatrixXd partR=LP*Y;
	MatrixXd LPL=LP*LAM;
	MatrixXd IDNF = MatrixXd::Identity(nf, nf);
	for(int j=0;j<d_y;j++){
		partV.setZero();
		for(int i1=0;i1<nf;i1++){
			for(int i2=0;i2<nf;i2++){
				partV(i1,i2) = LPL(i1,i2);
			}
			if(O(0,i1)==1){
				partV(i1,i1) += (double)1/SIGMA(i1,j);
			}else{
				partV(i1,i1) += (double)1/LAMX(i1);
			}
		}
			
		LDLT<MatrixXd> ldltOfA(partV);
		MatrixXd P = PermutationMatrix<Dynamic,Dynamic>(ldltOfA.transpositionsP());
		VectorXd D = ldltOfA.vectorD();
		D = double(1)/D.array().sqrt();
		MatrixXd D_inv = D.asDiagonal();
		MatrixXd Ll = ldltOfA.matrixL();
		fInvert(Ll_inv,Ll,IDNF,nf);
		MatrixXd L_inv=P.transpose()*Ll_inv.transpose()*D_inv;
		MatrixXd vx=L_inv*L_inv.transpose();
		for(int i1=0;i1<nf;i1++){
			X2I(i1,0)=gsl_ran_gaussian (r, 1);
		}
		MatrixXd EXI=vx*partR.col(j)+L_inv*X2I;				
		EX.col(j)=EXI.col(0);	
	}			
}

void draw_sigma(MatrixXd &SIGMA,MatrixXd &EX, MatrixXd &RHO,int nf,int d_y,int a,gig *mygig){
	for(int j=0;j<nf;j++){
		//if(O(0,j)==1){
		for(int i=0;i<d_y;i++){
			double sig=0;
			double pg,ag,bg;
			pg=a-0.5;
			ag=EX(j,i)*EX(j,i);
			bg=2*RHO(i);
			mygig->rgig(sig,&pg,&bg,&ag);
			SIGMA(j,i)=sig;
			//if(SIGMA(j,i)<1e-15){SIGMA(j,i)=1e-15;}
		}
		//}
	}
}
void draw_lamx(VectorXd &LAMX,MatrixXd &O, MatrixXd &RHO, VectorXd &KAPPA,MatrixXd &EX, gsl_rng *r, gig *mygig, int nf, int d_y, double b, double c){
	for(int i=0;i<nf;i++){
		if(O(0,i)==1){
			LAMX(i)=gsl_ran_gamma(r,d_y*b+c,double(1)/(RHO.row(i).sum()+KAPPA(i)));
		}
		else{
			double lamx;
			double pg,ag,bg;
			pg=c-0.5*d_y;
			ag=EX.row(i).dot(EX.row(i));
			bg=2*KAPPA(i);
			mygig->rgig(lamx,&pg,&bg,&ag);
			LAMX(i)=lamx;
						
			// for(int j=0;j<d_y;j++){
			// SIGMA(i,j)=LAMX(i);
			// }
					
		}
	}
}
*/

/*
void draw_o(MatrixXd &O,MatrixXd &logO, MatrixXd &LOGVO MatrixXd &EX, ){		
	for(int i=0;i<nf;i++){
		logO(0,i)=LOGVO(0,0);
		logO(1,i)=LOGVO(1,0);
		for(int j=0;j<d_y;j++){
			logO(0,i) = logO(0,i)+log_norm(EX(i,j),0,SIGMA(i,j))+lgamma(a+b)-(a+b)*log(LAMX(i)+SIGMA(i,j))+b*log(LAMX(i))+(a-1)*log(SIGMA(i,j))-lgamma(a)-lgamma(b);
			logO(1,i) = logO(1,i)+log_norm(EX(i,j),0,LAMX(i));
		}
	}
			
	for(int i=0;i<nf;i++){
		O(0,i)= gsl_ran_bernoulli(r,double(1)/(1+exp(logO(1,i)-logO(0,i))));
		O(1,i)=1-O(0,i);
								
		// if(O(0,i)==0){
		// 	for(int j=0;j<d_y;j++){
		// 		SIGMA(i,j)=LAMX(i);
		// 	}
		// }
		// else{
		// 	for(int j=0;j<d_y;j++){
		// 		RHO(i,j)=gsl_ran_gamma(r,a+b,double(1)/(SIGMA(i,j)+LAMX(i)));
		// 	}
		// }
				
	}
}
*/
void cal_ex(MatrixXd& EX, MatrixXd& LAM,MatrixXd& Y,MatrixXd& PSI_INV,MatrixXd& partR,MatrixXd& EXX,VectorXd& indexALL,MatrixXd& SIGMA,VectorXd& LAMX,MatrixXd& partV,MatrixXd& O,MatrixXd& LPL,int d_y, int nf_nc){
	
	partR=LAM.transpose()*PSI_INV*Y;
        
	EXX.setZero();
	//EX.setZero();
               
	//cout << "indexLAM " << endl << indexLAM << endl;
	for(int j=0;j<d_y;j++){
		int count_indexALL=0;
            
		for(int i=0;i<nf_nc;i++){
			if(SIGMA(i,j)!=0&&LAMX(i)!=0){
				indexALL(count_indexALL)=i;
				count_indexALL++;
			}
		}
			
		//cout << "indexSIG " << endl << indexSIG << endl;
		//cout << "indexALL " << endl << indexALL.transpose() << endl;
            
		if(count_indexALL==0){
			for(int i=0;i<nf_nc;i++){
				EX(i,j)=0;
			}
			continue;
		}
            
		partV.setZero();
            
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = LPL(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += O(0,indexALL(i1))/SIGMA(indexALL(i1),j)+O(1,indexALL(i1))/LAMX(indexALL(i1));
		}

		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2)=partV(indexALL(i1),indexALL(i2));
			}
		}
			
		MatrixXd partRI=MatrixXd::Constant(count_indexALL,1,0);
		MatrixXd EXI=MatrixXd::Constant(count_indexALL,1,0);
		MatrixXd EXXI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		//MatrixXd EXXI_cov=MatrixXd::Constant(nf-nf_nc,nf-nf_nc,0);
		// if(nf>nf_nc){
		// 	EXXI_cov=MatrixXd::Constant(nf-nf_nc,nf-nf_nc,0);
		// }
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
		// if(nf>nf_nc){
		// 	EXXI_cov=EX.block(nf_nc,j,nf-nf_nc,1)*EX.block(nf_nc,j,nf-nf_nc,1).transpose();
		// }
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2)+vx(i1,i2);
				//EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2);
			}
		}
		// if(nf>nf_nc){
		// 	for(int i1=0;i1<(nf-nf_nc);i1++){
		// 		for(int i2=0;i2<nf-nf_nc;i2++){
		// 			EXX(nf_nc+i1,nf_nc+i2) += EXXI_cov(i1,i2);
		// 			//EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2);
		// 		}
		// 	}
		// }
		
	}

}
/*
void cal_ex_VXM(MatrixXd& EX, MatrixXd& LAM,MatrixXd& Y,MatrixXd& PSI_INV,MatrixXd& partR,MatrixXd& EXX,VectorXd& indexALL,MatrixXd& SIGMA,VectorXd& LAMX,MatrixXd& partV,MatrixXd& O,MatrixXd& LPL, MatrixXd VXM,int d_y, int nf_nc){
	
	partR=LAM.transpose()*PSI_INV*Y;
        
	EXX.setZero();
	//EX.setZero();
               
	//cout << "indexLAM " << endl << indexLAM << endl;
	for(int j=0;j<d_y;j++){
		int count_indexALL=0;
            
		for(int i=0;i<nf_nc;i++){
			if(SIGMA(i,j)!=0&&LAMX(i)!=0){
				indexALL(count_indexALL)=i;
				count_indexALL++;
			}
		}
			
		//cout << "indexSIG " << endl << indexSIG << endl;
		//cout << "indexALL " << endl << indexALL.transpose() << endl;
            
		if(count_indexALL==0){
			for(int i=0;i<nf_nc;i++){
				EX(i,j)=0;
			}
			continue;
		}
            
		partV.setZero();
            
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = LPL(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += O(0,indexALL(i1))/SIGMA(indexALL(i1),j)+O(1,indexALL(i1))/LAMX(indexALL(i1));
		}

		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2)=partV(indexALL(i1),indexALL(i2));
			}
		}
			
		MatrixXd partRI=MatrixXd::Constant(count_indexALL,1,0);
		MatrixXd EXI=MatrixXd::Constant(count_indexALL,1,0);
		MatrixXd EXXI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		//MatrixXd EXXI_cov=MatrixXd::Constant(nf-nf_nc,nf-nf_nc,0);
		// if(nf>nf_nc){
		// 	EXXI_cov=MatrixXd::Constant(nf-nf_nc,nf-nf_nc,0);
		// }
		cpy_col_matrix(partRI,partR,indexALL,count_indexALL,0,j);

		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);
		MatrixXd vx=MatrixXd::Constant(count_indexALL,count_indexALL,0);

		vx = partVI.lu().solve(IDNF);
		EXI=vx*partRI;
		
		for(int i=0;i<count_indexALL;i++){
			VXM(j,indexALL(i))=vx(i,i);
		}
		
	
		for(int i=0;i<count_indexALL;i++){
			if(SIGMA(indexALL(i),j)==0){
				EXI(i)=0;
			}
		}
			
		cpy_col_matrix_bak(EX,EXI,indexALL,count_indexALL,j,0);
			
		EXXI=EXI*EXI.transpose();
		// if(nf>nf_nc){
		// 	EXXI_cov=EX.block(nf_nc,j,nf-nf_nc,1)*EX.block(nf_nc,j,nf-nf_nc,1).transpose();
		// }
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2)+vx(i1,i2);
				//EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2);
			}
		}
		// if(nf>nf_nc){
		// 	for(int i1=0;i1<(nf-nf_nc);i1++){
		// 		for(int i2=0;i2<nf-nf_nc;i2++){
		// 			EXX(nf_nc+i1,nf_nc+i2) += EXXI_cov(i1,i2);
		// 			//EXX(indexALL(i1),indexALL(i2)) += EXXI(i1,i2);
		// 		}
		// 	}
		// }
		
	}

}

*/


void cal_lam(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,VectorXd& vLXL,MatrixXd& ELL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV, int s_n,int nf){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	vLXL.setZero();
	ELL.setZero();
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
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
            
		partV.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}
							
		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			
		MatrixXd vl = partVI.lu().solve(IDNF);
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
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
			
	}

}

/*
void cal_lam_VLM(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,VectorXd& vLXL,MatrixXd& ELL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV,MatrixXd VLM, int s_n,int nf){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	vLXL.setZero();
	ELL.setZero();
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
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
            
		partV.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}
							
		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			
		MatrixXd vl = partVI.lu().solve(IDNF);
		LAMI=partLI*vl;

		for(int i=0;i<count_indexALL;i++){
			VLM(j,indexALL(i))=vl(i,i);
		}
			
		for(int i=0;i<count_indexALL;i++){
			if(THETA(j,indexALL(i))==0){
				LAMI(i)=0;
			}
		}
            
		cpy_row_matrix_bak(LAM,LAMI,indexALL,count_indexALL,j);
			
		LLI=LAMI.transpose()*LAMI;


		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
			
	}

}

void cal_lam_LPL(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,VectorXd& vLXL,MatrixXd& ELL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV,int s_n,int nf){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	vLXL.setZero();
	ELL.setZero();
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
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
            
		partV.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}
							
		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			
		MatrixXd vl = partVI.lu().solve(IDNF);
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
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
			
	}

}


void cal_lam_diag(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,VectorXd& vLXL,MatrixXd& ELL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV,int s_n,int nf){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	vLXL.setZero();
	ELL.setZero();
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
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
            
		partV.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}

		
		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			
		MatrixXd vl = MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i=0;i<count_indexALL;i++){
			vl(i,i)=1/partVI(i,i);
		}
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
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
			
	}

}


void cal_lam_wbr(MatrixXd& LAM,VectorXd& indexALL,MatrixXd& partL,MatrixXd& Y,MatrixXd& EX,MatrixXd& PSI_INV,MatrixXd& EXX,MatrixXd& Z,MatrixXd& LPL,VectorXd& vLXL,MatrixXd& ELL,MatrixXd& THETA,VectorXd& PHI,MatrixXd& partV,MatrixXd& sumPhiTheta,int s_n,int nf, int d_y){
	indexALL.setZero();
	partL=PSI_INV*Y*EX.transpose();
	//cout << "partL " << endl << partL.block(0,0,5,5) << endl;
	//LPL=LAM.transpose()*PSI_INV*LAM;
	LPL.setZero();
	vLXL.setZero();
	ELL.setZero();
	//LAM.setZero();
        
	//MatrixXd partV=MatrixXd::Constant(nf,nf,0);
        
	for(int j=0;j<s_n;j++){
		int count_indexALL=0;
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
		
		partV.setZero();
		sumPhiTheta.setZero();
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partV(indexALL(i1),indexALL(i2)) = PSI_INV(j,j)*EXX(indexALL(i1),indexALL(i2));
			}
			partV(indexALL(i1),indexALL(i1)) += Z(0,indexALL(i1))/THETA(j,indexALL(i1))+Z(1,indexALL(i1))/PHI(indexALL(i1));
			sumPhiTheta(indexALL(i1),indexALL(i1)) += THETA(j,indexALL(i1))/Z(0,indexALL(i1))+PHI(indexALL(i1))/Z(1,indexALL(i1));
		}  
       
		MatrixXd partVI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		for(int i1=0;i1<count_indexALL;i1++){
			for(int i2=0;i2<count_indexALL;i2++){
				partVI(i1,i2) = partV(indexALL(i1),indexALL(i2));
			}
		}

		MatrixXd partLI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LAMI=MatrixXd::Constant(1,count_indexALL,0);
		MatrixXd LLI=MatrixXd::Constant(count_indexALL,count_indexALL,0);
		cpy_row_matrix(partLI,partL,indexALL,count_indexALL,j);
            
		MatrixXd IDNF = MatrixXd::Identity(count_indexALL, count_indexALL);

		//LDLT<MatrixXd> ldltOfA(partVI);
		//MatrixXd vl=ldltOfA.solve(IDNF);
			

		MatrixXd vl = MatrixXd::Constant(count_indexALL,count_indexALL,0);
		if(count_indexALL<d_y){
			vl=partVI.lu().solve(IDNF);
		}
		else{
			
		}
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
				LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2)+vl(i1,i2));
				ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2)+vl(i1,i2);
				//LPL(indexALL(i1),indexALL(i2)) += PSI_INV(j,j)*(LLI(i1,i2));
				//ELL(indexALL(i1),indexALL(i2)) += LLI(i1,i2);
				vLXL(j) += vl(i1,i2)*EXX(indexALL(i1),indexALL(i2));
					
			}
		}
			
	}
		
}
*/

void cal_o(MatrixXd& logO,MatrixXd& LOGVO,MatrixXd& EX, MatrixXd& SIGMA, MatrixXd& RHO, VectorXd& LAMX, MatrixXd& O, int nf, int d_y, double a, double b){
    for(int i=0;i<nf;i++){
		logO(0,i)=LOGVO(0,0);
		logO(1,i)=LOGVO(1,0);
		for(int j=0;j<d_y;j++){
			logO(0,i)=logO(0,i)+log_norm(EX(i,j),0,SIGMA(i,j))+log_gamma(SIGMA(i,j),a,RHO(i,j))+log_gamma(RHO(i,j),b,LAMX(i));
			logO(1,i)=logO(1,i)+log_norm(EX(i,j),0,LAMX(i));
		}
	}
        
	// O
	for(int i=0;i<nf;i++){
		O(0,i)=double(1)/(1+exp(logO(1,i)-logO(0,i)));
		O(1,i)=1-O(0,i);
	}
}

void cal_z(MatrixXd& logZ,MatrixXd& LOGV,MatrixXd& LAM, MatrixXd& THETA, MatrixXd& DELTA, VectorXd& PHI, MatrixXd& Z, int nf, int s_n, double a, double b){
	for(int i=0;i<nf;i++){
		logZ(0,i)=LOGV(0,0);
		logZ(1,i)=LOGV(1,0);
		for(int j=0;j<s_n;j++){
			logZ(0,i)=logZ(0,i)+log_norm(LAM(j,i),0,THETA(j,i))+log_gamma(THETA(j,i),a,DELTA(j,i))+log_gamma(DELTA(j,i),b,PHI(i));
			logZ(1,i)=logZ(1,i)+log_norm(LAM(j,i),0,PHI(i));
		}
	}
        
	// Z 
	for(int i=0;i<nf;i++){
		Z(0,i)=double(1)/(1+exp(logZ(1,i)-logZ(0,i)));
		Z(1,i)=1-Z(0,i);
	}
}

void cal_v(MatrixXd& Z, MatrixXd& LOGV, double alpha, double beta, int nf){
	double ps1=alpha;
	double ps2=beta;
        
	// Probability 
	for(int i=0;i<nf;i++){
		ps1=ps1+Z(0,i);
		ps2=ps2+Z(1,i);
	}
	double dgama = gsl_sf_psi(ps1+ps2);
	LOGV(0,0)=gsl_sf_psi(ps1)-dgama;
	LOGV(1,0)=gsl_sf_psi(ps2)-dgama;
}

void cal_psi(MatrixXd& PSI,MatrixXd& LX, MatrixXd& Y,MatrixXd& LAM,MatrixXd& EX, MatrixXd& EXX, VectorXd& vLXL, int s_n, int d_y){
	
	for(int i=0;i<s_n;i++){
		//PSI(i,i)=Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i))+vLXL(i);
		//PSI(i,i)=Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i));
		PSI(i,i)=(0.5*(Y.row(i).dot(Y.row(i))-2*(LX.row(i)).dot(Y.row(i))+(LAM.row(i)*EXX).dot(LAM.row(i)))+1)/(double(d_y)/2+1);
	}
	 //PSI=PSI/d_y;
}

void write_itr(string dir_out, MatrixXd& LAM,MatrixXd& Z,MatrixXd& EX,MatrixXd& EXX,MatrixXd& O,MatrixXd& PSI,VectorXd& lam_count_v,MatrixXd& count_lam_M, MatrixXd& count_x_M, VectorXd& LAMX,VectorXd& PHI, int itr, long seed){
	stringstream ss;
	
		
	ss.str("");
	ss.clear();
	ss << dir_out << "/itr";
	ofstream f_itr (ss.str().c_str());
	if (f_itr.is_open()){
		f_itr << itr << endl;
	}
	f_itr.close();
				

	ss.str("");
	ss.clear();
	ss << dir_out << "/LAM_" << itr;
	ofstream f_lam (ss.str().c_str());
	if (f_lam.is_open()){
		f_lam << LAM << endl;
	}
	f_lam.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/LAMX_" << itr;
	ofstream f_lamx (ss.str().c_str());
	if (f_lamx.is_open()){
		f_lamx << LAMX << endl;
	}
	f_lamx.close();
                
	ss.str("");
	ss.clear();
	ss << dir_out << "/Z_" << itr;
	ofstream f_Z (ss.str().c_str());
	if (f_Z.is_open()){
		f_Z << Z << endl;
	}
	f_Z.close();
                
                
	ss.str("");
	ss.clear();
	ss << dir_out << "/EX_" << itr;
	ofstream f_EX (ss.str().c_str());
	if (f_EX.is_open()){
		f_EX << EX << endl;
	}
	f_EX.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/PHI_" << itr;
	ofstream f_phi (ss.str().c_str());
	if (f_phi.is_open()){
		f_phi << PHI << endl;
	}
	f_phi.close();
                

	ss.str("");
	ss.clear();
	ss << dir_out << "/EXX_" << itr;
	ofstream f_EXX (ss.str().c_str());
	if (f_EXX.is_open()){
		f_EXX << EXX << endl;
	}
	f_EXX.close();

	/*
	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/vLXL_" << itr;
	  ofstream f_vLXL (ss.str().c_str());
	  if (f_vLXL.is_open()){
	  f_vLXL << vLXL.diagonal() << endl;
	  }
	  f_vLXL.close();


	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/ELL_" << itr;
	  ofstream f_ELL (ss.str().c_str());
	  if (f_ELL.is_open()){
	  f_ELL << ELL << endl;
	  }
	  f_ELL.close();

	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/LPL_" << itr;
	  ofstream f_LPL (ss.str().c_str());
	  if (f_LPL.is_open()){
	  f_LPL << LPL << endl;
	  }
	  f_LPL.close();
	*/

	ss.str("");
	ss.clear();
	ss << dir_out << "/O_" << itr;
	ofstream f_O (ss.str().c_str());
	if (f_O.is_open()){
		f_O << O << endl;
	}
	f_O.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/seed";
	ofstream f_seed (ss.str().c_str());
	if (f_seed.is_open()){
		f_seed << seed << endl;
	}

	f_seed.close();

            
	ss.str("");
	ss.clear();
	ss << dir_out << "/PSI_" << itr;
	ofstream f_PSI (ss.str().c_str());
	if (f_PSI.is_open()){
		f_PSI << PSI.diagonal() << endl;
	}
	f_PSI.close();
                
   

	ss.str("");
	ss.clear();
	ss << dir_out << "/count_lam_v";
	ofstream f_count_lam (ss.str().c_str());
	if (f_count_lam.is_open()){
		f_count_lam << lam_count_v << endl;
	}
	f_count_lam.close();

	/*
	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/log_det_v";
	  ofstream f_det (ss.str().c_str());
	  if (f_det.is_open()){
	  f_det << std::setprecision(10) << log_det << endl;
	  }
	  f_det.close();

	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/bic_v";
	  ofstream f_bic (ss.str().c_str());
	  if (f_bic.is_open()){
	  f_bic << std::setprecision(10) << bic << endl;
	  }

	  f_bic.close();
	*/
    /*
	ss.str("");
	ss.clear();
	ss << dir_out << "/count_lam_M_" << itr;
	ofstream f_count_lam_M (ss.str().c_str());
	if (f_count_lam_M.is_open()){
		f_count_lam_M << std::setprecision(10) << count_lam_M << endl;
	}

	f_count_lam_M.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/count_x_M_" << itr;
	ofstream f_count_x_M (ss.str().c_str());
	if (f_count_x_M.is_open()){
		f_count_x_M << std::setprecision(10) << count_x_M << endl;
	}

	f_count_x_M.close();
     */
}

void write_final(string dir_out, MatrixXd& LAM,MatrixXd& Z,MatrixXd& EX,MatrixXd& EXX,MatrixXd& O,MatrixXd& PSI,VectorXd& lam_count_v,MatrixXd& count_lam_M, MatrixXd& count_x_M, VectorXd& LAMX,VectorXd& PHI, int itr, long seed){
	stringstream ss;
			
	ss.str("");
	ss.clear();
	ss << dir_out << "/itr";
	ofstream f_itr (ss.str().c_str());
	if (f_itr.is_open()){
		f_itr << itr << endl;
	}
	f_itr.close();
				

	ss.str("");
	ss.clear();
	ss << dir_out << "/LAM";
	ofstream f_lam (ss.str().c_str());
	if (f_lam.is_open()){
		f_lam << LAM << endl;
	}
	f_lam.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/LAMX";
	ofstream f_lamx (ss.str().c_str());
	if (f_lamx.is_open()){
		f_lamx << LAMX << endl;
	}
	f_lamx.close();

	
	ss.str("");
	ss.clear();
	ss << dir_out << "/Z";
	ofstream f_Z (ss.str().c_str());
	if (f_Z.is_open()){
		f_Z << Z << endl;
	}
	f_Z.close();
                
                
	ss.str("");
	ss.clear();
	ss << dir_out << "/EX";
	ofstream f_EX (ss.str().c_str());
	if (f_EX.is_open()){
		f_EX << EX << endl;
	}
	f_EX.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/PHI";
	ofstream f_phi (ss.str().c_str());
	if (f_phi.is_open()){
		f_phi << PHI << endl;
	}
	f_phi.close();
                
	ss.str("");
	ss.clear();
	ss << dir_out << "/EXX";
	ofstream f_EXX (ss.str().c_str());
	if (f_EXX.is_open()){
		f_EXX << EXX << endl;
	}
	f_EXX.close();

	/*
	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/vLXL_" << itr;
	  ofstream f_vLXL (ss.str().c_str());
	  if (f_vLXL.is_open()){
	  f_vLXL << vLXL.diagonal() << endl;
	  }
	  f_vLXL.close();


	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/ELL_" << itr;
	  ofstream f_ELL (ss.str().c_str());
	  if (f_ELL.is_open()){
	  f_ELL << ELL << endl;
	  }
	  f_ELL.close();

	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/LPL_" << itr;
	  ofstream f_LPL (ss.str().c_str());
	  if (f_LPL.is_open()){
	  f_LPL << LPL << endl;
	  }
	  f_LPL.close();
	*/

	ss.str("");
	ss.clear();
	ss << dir_out << "/O";
	ofstream f_O (ss.str().c_str());
	if (f_O.is_open()){
		f_O << O << endl;
	}
	f_O.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/seed";
	ofstream f_seed (ss.str().c_str());
	if (f_seed.is_open()){
		f_seed << seed << endl;
	}

	f_seed.close();

            
	ss.str("");
	ss.clear();
	ss << dir_out << "/PSI";
	ofstream f_PSI (ss.str().c_str());
	if (f_PSI.is_open()){
		f_PSI << PSI.diagonal() << endl;
	}
	f_PSI.close();
                
   

	ss.str("");
	ss.clear();
	ss << dir_out << "/count_lam_v";
	ofstream f_count_lam (ss.str().c_str());
	if (f_count_lam.is_open()){
		f_count_lam << lam_count_v << endl;
	}
	f_count_lam.close();

	/*
	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/log_det_v";
	  ofstream f_det (ss.str().c_str());
	  if (f_det.is_open()){
	  f_det << std::setprecision(10) << log_det << endl;
	  }
	  f_det.close();

	  ss.str("");
	  ss.clear();
	  ss << dir_out << "/bic_v";
	  ofstream f_bic (ss.str().c_str());
	  if (f_bic.is_open()){
	  f_bic << std::setprecision(10) << bic << endl;
	  }

	  f_bic.close();
	*/
    /*
	ss.str("");
	ss.clear();
	ss << dir_out << "/count_lam_M";
	ofstream f_count_lam_M (ss.str().c_str());
	if (f_count_lam_M.is_open()){
		f_count_lam_M << std::setprecision(10) << count_lam_M << endl;
	}

	f_count_lam_M.close();

	ss.str("");
	ss.clear();
	ss << dir_out << "/count_x_M";
	ofstream f_count_x_M (ss.str().c_str());
	if (f_count_x_M.is_open()){
		f_count_x_M << std::setprecision(10) << count_x_M << endl;
	}

	f_count_x_M.close();
	
	ss.str("");
	ss.clear();
	ss << dir_out << "/final";
	ofstream f_FINAL (ss.str().c_str());
	if (f_FINAL.is_open()){
		f_FINAL << "done" << endl;
	}
	f_FINAL.close();
     */
}

double like(MatrixXd& PSI,MatrixXd& EXX,MatrixXd& LAM,MatrixXd& THETA,MatrixXd& DELTA,VectorXd& PHI,VectorXd& TAU,MatrixXd& Z,MatrixXd& V,double ETA, double GAMMA,double alpha, double beta, int n,int p, int nf,double a, double b, double c, double d, double e, double f, double nu){
  double det_psi=0;
  
  for(int i=0;i<n;i++){
    det_psi = det_psi+log(PSI(i,i));
  }
  
  double like=(-1)*0.5*n*p*log(2*PI_mine)-0.5*p*det_psi;
  
  double sum_x=0;
  for(int i=0;i<nf;i++){
    sum_x = sum_x + (-1)*0.5*EXX(i,i);
  }
  
  like = like - 0.5*nf*p*log(2*PI_mine) + sum_x;
  
  for(int i=0;i<nf;i++){
    like=like + (Z(0,i)+alpha)*V(0,i) + (Z(1,i)+beta)*V(1,i);
  }
  
  for(int i=0;i<n;i++){
    for(int j=0;j<nf;j++){
      if(THETA(i,j)!=0){
        like=like+Z(0,j)*log_norm(LAM(i,j),0,THETA(i,j));
        //if(DELTA(i,j)!=0){
          like=like+Z(0,j)*log_gamma(THETA(i,j),a,DELTA(i,j));
          //}
      }
      like=like+Z(0,j)*log_gamma(DELTA(i,j),b,PHI(j));
      like=like+Z(1,j)*log_norm(LAM(i,j),0,PHI(j));
    }
  }
 
 for(int i=0;i<nf;i++){
   like=like+log_gamma(PHI(i),c,TAU(i));
   like=like+log_gamma(TAU(i),d,ETA);
 }
 
 like=like+log_gamma(ETA,e,GAMMA);
 like=like+log_gamma(GAMMA,f,nu);
  
  return like;
 
}
/*
#include <iostream>     // cout
#include <math.h>       // acos
#include <float.h>      // DBL_MAX
#include <limits>       // numeric_limits
*/
template<typename T>
bool is_infinite( const T &value )
{
    // Since we're a template, it's wise to use std::numeric_limits<T>
    //
    // Note: std::numeric_limits<T>::min() behaves like DBL_MIN, and is the smallest absolute value possible.
    //
 
    T max_value = std::numeric_limits<T>::max();
    T min_value = - max_value;
 
    return ! ( min_value <= value && value <= max_value );
}
 
template<typename T>
bool is_nan( const T &value )
{
    // True if NAN
    return value != value;
}
 
template<typename T>
bool is_valid( const T &value )
{
    return ! is_infinite(value) && ! is_nan(value);
}

/*
int main()
{
    using std::cout;
 
    double a, b, c, d, e;
 
    a = 1.0;
    b = 0.0;
    c = a / c;          // divide by zero
    d = acos(-1.001);   // domain for acos is [-1, 1], anything else is #IND or inf
    e = b / b;          // zero / zero
 
    cout << "Value of a: " << a << " " << is_valid(a) << " " << (is_nan(a) ? " nan " : "") << (is_infinite(a) ? " infinite " : "") << "n";
    cout << "Value of b: " << b << " " << is_valid(b) << " " << (is_nan(b) ? " nan " : "") << (is_infinite(b) ? " infinite " : "") << "n";
    cout << "Value of c: " << c << " " << is_valid(c) << " " << (is_nan(c) ? " nan " : "") << (is_infinite(c) ? " infinite " : "") << "n";
    cout << "Value of d: " << d << " " << is_valid(d) << " " << (is_nan(d) ? " nan " : "") << (is_infinite(d) ? " infinite " : "") << "n";
    cout << "Value of e: " << e << " " << is_valid(e) << " " << (is_nan(e) ? " nan " : "") << (is_infinite(e) ? " infinite " : "") << "n";
 
    return 0;
}
*/
