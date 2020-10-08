// bfrm.cpp: implementation of the Bfrm class.
//Author: Quanli Wang quanli@stat.duke.edu
//
//////////////////////////////////////////////////////////////////////
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines

#include "newran.h"
#include "tnt.h"

using namespace TNT;

#include <math.h>
#include "stdafx.h"

#include "mdp.h"

#include "Model.h"
#include "bfrmrawdata.h"
#include "bfrmprior.h"

#include "bfrmresult.h"

#include "bfrm.h"
#include "util.h"


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
Bfrm::Bfrm()
{
        Bfrm(-1.0);
}
Bfrm::Bfrm(double seed)
{
	//initialise the newran library
        if (seed < 0){
	      srand( (unsigned)time( NULL ) );
	      Random::Set(((double)rand()) / ((double)RAND_MAX));
	} else {
	      Random::Set(seed);
        }
	//initialise the optimization variables
	mpabj = NULL; mpafj = NULL; mpaibj = NULL; 
	mpaMeanBJ = NULL;
	mpaVarBJ = NULL;
	mnAppend = 0;
}

Bfrm::~Bfrm()
{
	
}



void Bfrm::LoadData(Model& model) {
	if (!InitialiseH(model.mstrHFile,model.mnNObservations,model.mnNDesignControls, model.mnNControlFactors)) {
		std::cout << "Loading H file failed!" << endl;
		exit (1);
	}
	
	if (!InitialiseWeight(model.mstrWeightFile,model.mnNObservations)) {
		std::cout << "Loading Weight file failed!" << endl;
		exit (1);
	}
	
	if (!LoadDataXY(model.mstrDataFile,model.GetDataRowCount(),model.mnNObservations, model.mnNVarBinary,
			model.mnNVarCategorical, model.mnNVarSurvival, model.mnNVarContinuous)) {
		std::cout << "Loading data file failed!" << endl;
		exit (1);
	}
	
	if (!InitialiseYMask(model.mstrResponseMaskFile)) {
		std::cout << "Loading response factor mask file failed!" << endl;
		exit (1);
	}
	
	if (!InitialiseXMask(model.mstrXMaskFile)) {
		std::cout << "Loading X mask file failed!" << endl;
		exit (1);
	}
	if (model.mnNDesignControls == 0) {
		SubtractSampleMeans();
	}
}

void Bfrm::Run(int nStartIter, Model& model) {
	int NonGaussian = model.mnNonGaussianFactors;
	int nExtra = NonGaussian?50:0;
	if (model.mnEvol == 0) { nExtra = 0; }
	model.mnNonGaussianFactors = 0; //force to be gaussian first
	//SaveIter(0);
	for (int it=1 + nStartIter; it <= model.mnMC + model.mnBurnin + nExtra; it++) {
		if (NonGaussian && it > nExtra ) {
			model.mnNonGaussianFactors = NonGaussian;
		}
		if(it% model.mnPrintIteration == 0 || model.mnPrintIteration < 2) { 
			std::cout << "MCMC iteration " << it << endl;
			std::cout << "Sparsity " << mdSparsity << endl;
		}

		Iterate(model.mnNLatentFactors, model.mnNDesignControls, model.mnPriorPiStandard, 0, 
			model.mdDPAlpha, model.mnNonGaussianFactors, it);

		if (it>model.mnBurnin + nExtra) {
			Accumulate(model.mnNLatentFactors, model.mnNDesignControls, model.mnPriorPiStandard);
			if (it > model.mnMC + model.mnBurnin + nExtra - model.mnHistory) {
				SaveHist(model.mnPriorPiStandard, mnAppend);
				mnAppend = 1;
			}
		}
		
	}
	Summarise(model.mnNLatentFactors, model.mnNDesignControls ,model.mnMC, model.mnPriorPiStandard, 0);
}

void Bfrm::Save(int joel) {
	mMeanResult.Save(joel);
}
void Bfrm::SaveHist(int joel, int append) {
	BfrmData::SaveData("mA_Hist.txt",mIter.mB, append);
	BfrmData::SaveData("mF_Hist.txt",mIter.mF, append);
	BfrmData::SaveData("mPsi_Hist.txt", mIter.mPsi, append);
	BfrmData::SaveData("mTau_Hist.txt", mIter.maTau, append);
	BfrmData::SaveData("mAlpha_Hist.txt", mIter.mAlpha,append);
	BfrmData::SaveData("mPib_Hist.txt", mIter.maRho,append);
	BfrmData::SaveData("mPostPib_Hist.txt", mIter.mPostPib,append);
	if (mnYFactors > 0) {
		BfrmData::SaveData("mZ_Hist.txt", mIter.mZ,append);
	}
	if (joel) {
		BfrmData::SaveData("mPib2_Hist.txt", mIter.maPi,append);
	}
	if (mbHasXMask) {
		BfrmData::SaveData("mX_Hist.txt", mIter.mX,append);
	}
}
//accumulate the iteration result
void Bfrm::Accumulate(int factor, int design, int joel) {
	mIter.mAlpha[0] = mDP.mAlpha;
	mResult.Accumulate(mIter, joel);
}

void Bfrm::Initialise(Model& model, int FindFounder)
{
	int i, j;
	int design = model.mnNDesignControls;
	int nTotalFactor = design + model.mnNLatentFactors + mnYFactors;
	
	//mIter components will be initilised seperately
	
	mPrior.Initialise(model, mnYFactors, mnVariables);
	mPrior.InitialisePsi(mIter.mPsi, mnVariables, mnYFactors - mnVarSurvival - mnVarContinuous, mnYFactors - mnVarContinuous );
	mPrior.InitialiseTau(mIter.maTau, nTotalFactor);
	mPrior.InitialiseRho(mIter.maRho, nTotalFactor);
	
	mDP.init(model.mnNObservations,model.mdDPAlpha,model.mdDPAlphaX, model.mdaPriorAlpha[0], model.mdaPriorAlpha[1]);
	
	if (model.mnPriorPiStandard) {
		mPrior.InitialisePibIJ(mIter.maPi, nTotalFactor,mnVariables);
	}

	mIter.mPostPib = Array2D<double>(mnVariables, nTotalFactor); mIter.mPostPib = 0;
	if (model.mnNDesignControls) {			//only do this if there is an intercept factor
		for (i = 0; i < mnVariables; i++) {
			mIter.mPostPib[i][0] = 1.0;			
		}
	}
	
	SetBMasks(model.mnShapeOfB, FindFounder);
	InitiliseX();
	Array2D<double> Mean = RowMeanXY();

	//initilize B matrix
	int ntotalzeros = 0;
	mIter.mB = Array2D<double>(mnVariables, nTotalFactor); mIter.mB = 0;
	for (i = 0; i < mnVariables; i++) {
		for (j = 0; j < nTotalFactor; j++) {	
			if (design && j == 0) {
				mIter.mB[i][0] = Mean[i][0];
			} else {
				int nMask = GetBMask(i,j,model.mnNDesignControls, 0);
				if (nMask == 0) {
					mIter.mB[i][j] = 0;
					ntotalzeros++;
				} else if (nMask == 1) {
					//July 29 2005, use carlos's initialisation instead 
					mIter.mB[i][j] = 1.0; //mPrior.SampleRestrictedNormal(0,1,1,0);
				} else {
					//July 29 2005, use carlos's initialisation instead
					mIter.mB[i][j] = 0.0;
				}
			}
		}
	}
	mdSparsity = (double)(ntotalzeros-mnVariables)/(double)(mnVariables * nTotalFactor);
	
	//initilise  F matrix
	mIter.mF = Array2D<double>(nTotalFactor, mnSampleSize);
	//design part
	for (i = 0; i < design; i++) {
		for (j = 0; j < mnSampleSize; j++) {
			mIter.mF[i][j] = mWeight[j] ? mH[j][i] : 0;			//assume H: mnSampleSize * design
		}
	}

	//set YFactors to zero based on Mike's suggestion, July 28,2005
	for (i = design; i < design + mnYFactors; i++) {
		for (j = 0; j < mnSampleSize; j++) {
			mIter.mF[i][j] = mWeight[j] ? 0 : 0;
		}
	}
	//end Mike's suggestion

	//latent factor part
	for (i = design + mnYFactors; i < nTotalFactor; i++) {
		for (j = 0; j < mnSampleSize; j++) {
			mIter.mF[i][j] = mWeight[j] ? mX[i - design][j] - Mean[i-design][0] : 0;
		}
	}
	
	
	//initilise residules, part 1: X 
	mResX = mX - matmult(mIter.mB,mIter.mF);
	if (mnLeaveout != 0) {			//if there are leave-out samples
		for (j = 0; j < mnSampleSize; j++) {
			if (mWeight[j] == 0) {
				for (i = 0; i < mnVariables; i++) {
					mResX[i][j] = 0;
				}
			}
		}
	}

	//initilize Z if Y factors exist
	//initilise residules, part 2: Z, if Y factors exist
	if (mnYFactors > 0) {	
		mIter.mZ = Array2D<double>(mnYFactors, mnSampleSize);
		SampleZ(model.mnNLatentFactors, model.mnNDesignControls);
	}
	
	//estimate tau
	mIter.maTau = BfrmData::Var(mIter.mF);
	for (i = 0; i < nTotalFactor; i++) {
		if (!mIter.maTau[i]) { mIter.maTau[i] = 1e-9; }
	}
	

	//addjust F and B so F ~ N(0,1)
	Array2D<double> T(nTotalFactor,nTotalFactor); T = 0; T[0][0] = 1.0;
	for (i = 1; i < nTotalFactor; i++) {
		T[i][i] = sqrt(mIter.maTau[i]);
	}
	mIter.mB = matmult(mIter.mB, T);

	
	for (i = 1; i < nTotalFactor; i++) {
		T[i][i] = 1.0 / sqrt(mIter.maTau[i]);
	}
	mIter.mF = matmult(T, mIter.mF);

	mIter.mAlpha = Array1D<double>(1);
	
	
	ManageHelperVariables(true, model.mnNLatentFactors);

	mResult.Init(mnVariables, nTotalFactor, mnSampleSize, mnYFactors, model.mnPriorPiStandard,mbHasXMask, model.mnNDesignControls);
	mMeanResult.Init(mnVariables, nTotalFactor, mnSampleSize, mnYFactors, model.mnPriorPiStandard,mbHasXMask, model.mnNDesignControls);
	if (mbHasXMask) {
		mIter.mX = Array2D<double>(mnVariables - mnYFactors, mnSampleSize); mIter.mX = 0;
	}
}

void Bfrm::RankVariables(int nSkip1, int nSkip2, int VarInCount, bool screener) {
	int i;
	int nVarSkip = mnYFactors + VarInCount;
	int nVar = mnVariables - nVarSkip;

	maOrderedVariables.clear();
	maOrderedScores.clear();

	double *pScores = new double[nVar];
	double *pPostScores = new double[nVar];
	for (i = 0; i < nVar; i++) {
		pScores[i] = 0; pPostScores[i] = 0;
		for (int j = 1; j < mMeanResult.mPostPib.dim2(); j++) {
			if (j < nSkip1 || j >= nSkip2) {		//ignore control factors
				double dTemp = fabs(mMeanResult.mPostPib[i+nVarSkip][j] * mMeanResult.mB[i+nVarSkip][j]);
				if (dTemp > pScores[i]) {
					pScores[i] = dTemp;
					pPostScores[i] = mMeanResult.mPostPib[i+nVarSkip][j];
				}
			}			
		}
	}
	//note: here we rank them based on mB.* mPib but record scores as mPib
	
	util::sort(nVar, pScores, maOrderedVariables, -1); //descending
	
	for (i = 0; i < nVar; i++) {
		maOrderedScores.push_back(pPostScores[maOrderedVariables[i] - 1]);
		if (screener) {
			maOrderedVariables[i] = maOriginalIndex[maOrderedVariables[i] - 1 + VarInCount]; //so make it zero-based
		} else {
			maOrderedVariables[i] = maOrderedVariables[i] - 1;
		}
	}
	delete [] pScores;
	delete [] pPostScores;
}

void Bfrm::SaveIter(int i) {
	
//	BfrmData::SaveData("Z.txt", mIter.mZ);
	BfrmData::SaveData("A"  + util::ToString(i) + ".txt",mIter.mB);
	BfrmData::SaveData("F"  + util::ToString(i) + ".txt",mIter.mF);
	BfrmData::SaveData("Psi" + util::ToString(i) + ".txt", mIter.mPsi);
	BfrmData::SaveData("Tau" + util::ToString(i) + ".txt", mIter.maTau);
	BfrmData::SaveData("Pib.txt", mIter.maRho);
	BfrmData::SaveData("PostPib.txt", mIter.mPostPib);
	BfrmData::SaveData("X.txt", mIter.mX);
	BfrmData::SaveData("T"  + util::ToString(i) +  ".txt", mIter.mT);
}

void Bfrm::Iterate(int factor, int design, int joel, int mode, double alpha, int nongaussian,  int iter) {

	ManageHelperVariables(true, factor);			//usefull when COM mode is used
	if (!mode) {
		if (nongaussian) {
			SampleF(factor, design, &mDP, iter, nongaussian);
		} else {
			SampleF(factor,design);
		}
	}
	
	SampleB(factor, design, joel, mode);

	SamplePsi();
	
	if (!mode) {
		SampleZ(factor, design);
	}
	
	SampleX(factor, design);

}

void Bfrm::ManageHelperVariables(bool allocate, int factor) {
	if (allocate) {
		if (mpaibj == NULL) { //usefull in com mode
			mpabj = new double[mnVariables];
			mpafj = new double[mnSampleSize];
			mpaibj = new int[mnVariables];

			mpaMeanBJ = new double[mnVariables]; 
			mpaVarBJ = new double[mnVariables]; 
		}
	} else {
		delete [] mpabj; mpabj = NULL;
		delete [] mpafj; mpafj = NULL;
		delete [] mpaibj; mpaibj = NULL;
		delete [] mpaMeanBJ; mpaMeanBJ = NULL;
		delete [] mpaVarBJ; mpaVarBJ = NULL;
	}
}

//update mB mpib mPostPib mtau mResX
//this is the main function for sampling, 
//splitting it into small pieces will slow down the code quite a bit
void Bfrm::SampleB(int factor, int design, int joel, int mode) {
	int i,j,k;
	int ntotalBz = 0;

	int nTotalFactor = design + mnYFactors + factor;
	for ( j = 0; j < nTotalFactor; j++) {
		//init fj
		//maybe we need to go back to [] instead of *() 
		double v = 0.0;
		k=0;
		while(k < mnSampleSize) {		//can be optimized here since the first desgin v's are constant!!!
			double dt = mIter.mF[j][k];
			v+= dt * dt;
			*(mpafj+k++) = dt;
		}
		
		int nnz = 0;
		if (v) {	
			//calculate MeanBJ and VarBJ
			v = 1.0 / v;
			CalculateMeanVarBJ(j,v);			//multithreaded

			int m = 0; 
			double bsumsquare = 0;
			
			//used when all loadings in a factor are zero
			double oddratio = (j || !design) ? mIter.maRho[j] / (1 - mIter.maRho[j]):0;
			
			
			bool ConstantTauCondition = ((j == 0) && (design > 0));
			
			//use for the special case where all Bij's are zero
			int minindex = -1; 
			double minv = 1e10;
			double oldminbij = 0;
			int mPi_ij = 0;
			double postpi;

			for ( i = 0; i < mnVariables; i++) {
				double newvij;
				SampleBij(i,j,v,design, oddratio, newvij, mPrior.GetPriorMean(i),postpi,0, 
					joel, j == nTotalFactor - 1, mode);
				int nMask = GetBMask(i,j,design, mode);
				//count non-zero Pi_ij's
				if (joel && mIter.maPi[i][j]) {		
					mPi_ij++;
				}
				if (mIter.mB[i][j]) {
					m++; 
					double mv = mPrior.GetPriorMean(i);
					if (j || !design) {
						bsumsquare += mIter.mB[i][j] * mIter.mB[i][j];
					} else {
						bsumsquare += (mIter.mB[i][j] - mv)  * (mIter.mB[i][j] - mv);
					}
				} else {				
					//only usefull when all Bij == 0
					if (m==0 && nMask != 0) { 
						if (minindex == -1 || newvij < minv) {				
							minindex = i;
							minv = newvij;
							oldminbij = mpabj[i];
						}
					}
				}

				//record changeovers
				double deltaBij = mIter.mB[i][j]-mpabj[i];
				if (deltaBij) {
					mpabj[nnz] = deltaBij;
					mpaibj[nnz++] = i;
				}
			}
		
			//special case, all B[][j] are zeros
			if (m == 0) {   
				if (minindex < 0) {
					minindex = 0;
				}
				mIter.mB[minindex][j] = oldminbij;
				double newvij;
				SampleBij(minindex,j,v,design, oddratio, newvij, 
					mPrior.GetPriorMean(minindex),postpi,1, joel, j == nTotalFactor - 1, mode);
				m++;
				if (joel) { 
					if (mIter.maPi[minindex][j]) {mPi_ij++;}
				}

				//record changeovers
				mpabj[nnz] = mIter.mB[minindex][j] - oldminbij;
				mpaibj[nnz++] = minindex;
			}
			ntotalBz += m;

			//sampling rho/pi and tau
			if (j || !design) {	//if not intercept factor
				int jaddjust = GetBMaskAddjust(j,design);
				if (joel) { 
					mIter.maRho[j] = mPrior.SampleRho(j, mPi_ij, mnVariables - jaddjust - mPi_ij);
				} else {
					mIter.maRho[j] = mPrior.SampleRho(j, m, mnVariables - jaddjust - m);
				}
				
				mIter.maTau[j] = mPrior.SampleTau(j, (double)m,  bsumsquare);
				
				
			} else { // the intercept
				mIter.maRho[j] = 1.0;		
				if (ConstantTauCondition) {
					mIter.maTau[j] = 0.0;	//dummy value, not really necessary since it is never used
				} else {
					mIter.maTau[j] = mPrior.SampleTau(j, (double)mnVariables,  bsumsquare);
				}
			}
		} else {			//variance  == 0
			for ( i = 0; i < mnVariables; i++) {
				if (mIter.mB[i][j]) {
					mpabj[nnz] = -mIter.mB[i][j];
					mpaibj[nnz++] = i;
					mIter.mB[i][j] = 0;
				}
				if (j || !design) {
					mIter.mPostPib[i][j] = 0;
					if (joel) { mIter.maPi[i][j] = 0;}
				}
			}
			mIter.maTau[j] = 0;
			mIter.maRho[j] = 0;	
		}

		//update mResX
		BUpdateResX(nnz); //multithreaded
	}
	mdSparsity = (double)(ntotalBz-mnVariables) / (double)(mnVariables * nTotalFactor);
}

void Bfrm::BUpdateResX(int nStart, int nEnd) {
	for (int i = nStart; i < nEnd; i++) {
		int index = mpaibj[i];
		for (int k = 0; k < mnSampleSize; k++) {
			if (mWeight[k]) {
				mResX[index][k] -=  mpabj[i] * mpafj[k];
			}
		}
	}
}
void Bfrm::BUpdateResX(int nnz) {
	for (int i = 0; i < nnz; i++) {
		int index = mpaibj[i];
		for (int k = 0; k < mnSampleSize; k++) {
			if (mWeight[k]) {
				mResX[index][k] -=  mpabj[i] * mpafj[k];
			}
		}
	}
}
void Bfrm::SamplePsi(int nStart, int nEnd) {
	//NO initilisation here 

	int nMask = mnYFactors - mnVarSurvival - mnVarContinuous;
	for (int i = nStart; i < nEnd; i++) {
		double sumE = 0;
		if (i >= nMask) {
			for (int k = 0; k < mnSampleSize; k++) {
				if (mWeight[k]) {
					double temp = mResX[i][k];
					sumE += temp * temp;
				}
			}
		}
		mIter.mPsi[i] = mPrior.SamplePsi(sumE, i);
	}
}
//update Psi
void Bfrm::SamplePsi() {
	int i,k;	
	mPrior.InitGammaPsi(mnSampleSize - mnLeaveout);
	
	
	//SamplePsi(0, mnVariables);
	int nMask = mnYFactors - mnVarSurvival - mnVarContinuous;
	for (i = 0; i < mnVariables; i++) {
		double sumE = 0;
		if (i >= nMask) {
			for (k = 0; k < mnSampleSize; k++) {
				if (mWeight[k]) {
					double temp = mResX[i][k];
					sumE += temp * temp;
				}
			}
		}
		mIter.mPsi[i] = mPrior.SamplePsi(sumE, i);
	}
	
}

void Bfrm::CalculateMeanVarBJ(int j, double v) {
	
	//CalculateMeanVarBJ(j,v,0,mnVariables);
	for (int i = 0; i < mnVariables; i++) {
		mpabj[i] = mIter.mB[i][j];
		double oldBij = mpabj[i];
		double Bij = 0;
		//here comes the most time consuming line
		for (int k = 0; k < mnSampleSize; k++) {
			Bij += (mResX[i][k] + oldBij * mpafj[k]) * mpafj[k];
		}
		mpaMeanBJ[i] = Bij * v;
		mpaVarBJ[i] = mIter.mPsi[i] * v;
	}
	
}
void Bfrm::CalculateMeanVarBJ(int j, double v, int nStart, int nEnd) {
	for (int i = nStart; i < nEnd; i++) {
		mpabj[i] = mIter.mB[i][j];
		double oldBij = mpabj[i];
		double Bij = 0;
		//here comes the most time consuming line
		for (int k = 0; k < mnSampleSize; k++) {
			Bij += (mResX[i][k] + oldBij * mpafj[k]) * mpafj[k];
		}
		mpaMeanBJ[i] = Bij * v;
		mpaVarBJ[i] = mIter.mPsi[i] * v;
	}
}
void Bfrm::SampleBij(int i, int j, double v, int design, double oddratio, double& newv, double meanvalue, 
	double &postpi, int specialcase, int joel, int lastfactor, int mode) {
	
	bool ConstantTauCondition = ((j == 0) && (design > 0));
	double tauj = ConstantTauCondition?mPrior.GetPriorVar(i) : mIter.maTau[j];	
	int nMask = GetBMask(i,j,design, mode);
	newv = 1e9;
	
	if (j || !design) { 
		if (nMask == -1) { //no changes tp be made
			postpi = mIter.mPostPib[i][j];
			//no changes to be made to mB[i][j] & mPostPib{i][j]
		} else if (nMask == 0) {
			postpi = 0;
			mIter.mPostPib[i][j] = 0;
			mIter.mB[i][j] = 0;
			if (joel) { mIter.maPi[i][j] = 0; }
		} else if (nMask == 1) {
			postpi = 1;
			mIter.mPostPib[i][j] = 1;
			double e = 1 / (tauj + mpaVarBJ[i]);
			double dMu =  e * mpaMeanBJ[i] * tauj;
			double sqrtEtau = sqrt(tauj * mpaVarBJ[i] * e);
			if (sqrtEtau <= 1e-9) { sqrtEtau = 1e-9;}
			//mIter.mB[i][j] = dMu + sqrtEtau * mNormal.Next();
			
			if (lastfactor && mnBCarlosMask) {
				mIter.mB[i][j] = dMu + sqrtEtau * mNormal.Next();
			} else {
				mIter.mB[i][j] = mPrior.SampleRestrictedNormal(dMu, sqrtEtau * sqrtEtau, 1, 0);
			}			
			
			if (joel) { mIter.maPi[i][j] = 1.0; }
		} else {
			double e, sqrtE,sqrtEtau;
			try {//might work for some platforms
				double temp = (mpaMeanBJ[i] * mpaMeanBJ[i]*tauj)/(2 *(mpaVarBJ[i] *(tauj + mpaVarBJ[i])));
				if (joel) { //need to recalculate oddratio for joeL
					oddratio = mIter.maPi[i][j] / (1 - mIter.maPi[i][j]);
				}
				e = 1 / (tauj + mpaVarBJ[i]);
				sqrtE = sqrt(mpaVarBJ[i] * e);
				sqrtEtau = sqrtE * sqrt(tauj);
				newv = oddratio * sqrtE * exp(temp);
				#if defined(BFRM_LINUX)
					if (!finite(newv)) {
						newv  = 1e9;
					}
				#else
					if (!_finite(newv)) {
						newv  = 1e9;
					}
				#endif
			} catch (...) {
				newv = 1e9;
			}
			if (newv > 1e9) {newv = 1e9;}
			if (newv < 1e-9) {newv = 1e-9;}
			
			postpi = newv / ( 1 + newv);
			mIter.mPostPib[i][j] = postpi;

			if (specialcase) {
				mIter.mB[i][j] = e * mpaMeanBJ[i] * tauj + 0.01 * sqrtEtau * mNormal.Next();
			} else {
				mIter.mB[i][j] =  ((mRandom.Next() <= postpi)?(e * mpaMeanBJ[i] * tauj + sqrtEtau * mNormal.Next()):0);
			}
			if (joel) { mIter.maPi[i][j] = mPrior.SamplePiIJ(j, mIter.mB[i][j], mIter.maRho[j]); }

		}
	} else { //the intercept
		if (nMask >=0) {
			double e = 1 / (tauj + mpaVarBJ[i]);
			double sqrtEtau = sqrt(tauj * mpaVarBJ[i] * e);
			mIter.mB[i][j] =  meanvalue + e * (mpaMeanBJ[i] - meanvalue) * tauj + sqrtEtau * mNormal.Next();
		}
		postpi = 1;
		mIter.mPostPib[i][j] = 1;
		if (joel) {mIter.maPi[i][j] = 1; }
	}
}

//update mF mResX
void Bfrm::SampleF(int factor, int design) {
	int i,j,k;
	int N = factor + mnYFactors;
	if (N <= 0) return;

	//get B' * inv(diag(psi)) 
	double* BtPsiInv = new double[N * mnVariables];
	int nCount = 0;
	for ( j = design; j < design + N; j++) {
		for (i = 0; i < mnVariables; i++) {
			BtPsiInv[nCount++] = mIter.mB[i][j] / mIter.mPsi[i];
		}
	}
	SymmetricMatrix QInv(N);
	double *p = QInv.Store(); nCount = 0;
	for ( i = 0; i < N; i++) {
		for ( j = 0; j <= i; j++) {
			double sum = 0;
			int column = j + design;
			int anchor = i * mnVariables;
			for (k = 0; k < mnVariables; k++) {
				sum += BtPsiInv[anchor + k] * mIter.mB[k][column];
			}
			if (i == j) {
				sum += 1;
			}
			p[nCount++] = sum;
		}
	}

	LowerTriangularMatrix LInv = Cholesky(QInv);
	LowerTriangularMatrix L = LInv.i();
	
	double * xp = new double[mnVariables];
	double* ui = new double[N];
	double* Ltui = new double[N];
	double* lambda_i = new double[N];
	for (i = 0; i < mnSampleSize; i++) { 

		//now calculate x_i
		for (j = 0; j < mnVariables; j++) {
			double sum = mResX[j][i];
			for (k = design; k < design + N; k++) {
				sum += mIter.mB[j][k] * mIter.mF[k][i];		
			}
			xp[j] = sum;
		}
		
		for (j = 0; j < N;j++) {
			double sum = 0;
			int anchor = j * mnVariables;
			for (k = 0; k < mnVariables; k++) {
				sum += BtPsiInv[anchor + k] * xp[k];
			}
			ui[j] = sum;
		}
	
		for (j = 0; j < N; j++) {
			Ltui[j] = 0;
		}
		p = L.Store(); nCount = 0;			
		for (j = 0; j < N; j++) {
			Ltui[j] = mNormal.Next();
			for (k = 0; k <= j; k++) {
				Ltui[j] += p[nCount++] * ui[k];
			}
		}
		//now times by L
		//calculate L * (Ltui + nr);
		for (j = 0; j < N; j++) {
			lambda_i[j] = 0;
		}
		p = L.Store(); nCount = 0;			
		for (j = 0; j < N; j++) {
			for (k = 0; k <= j; k++) {
				lambda_i[k] += p[nCount++] * Ltui[j];
			}
		}

		for (k = design; k < design + N; k++) {
			double temp = mIter.mF[k][i];
			int ind = k - design;
			mIter.mF[k][i] = lambda_i[ind];
			lambda_i[ind] =temp - lambda_i[ind]; 
		}


		//now updating mResX
		for (j = 0; j < mnVariables; j++) {
			for (k = 0; k < N; k++) {
				mResX[j][i] += mIter.mB[j][k + design] * lambda_i[k];
			}
		}
	}
	delete [] xp;
	delete [] ui;
	delete [] Ltui;
	delete [] lambda_i;
	delete [] BtPsiInv;
}


void Bfrm::SampleF(int factor, int design, mdp* pdp, int iter, int nongaussian) {
	int i,j,k,m;
	int N = factor + mnYFactors;
	if (N <= 0) return;

	//get B' * inv(diag(psi)) 
	double* BtPsiInv = new double[N * mnVariables];
	int nCount = 0;
	for ( j = design; j < design + N; j++) {
		for (i = 0; i < mnVariables; i++) {
			BtPsiInv[nCount++] = mIter.mB[i][j] / mIter.mPsi[i];
		}
	}

	SymmetricMatrix QInv(N);
	SymmetricMatrix QInv_j(N);
	IdentityMatrix Id(N);
	double *p = NULL;
	p = QInv.Store(); nCount = 0;
	double * p_i = QInv_j.Store();
	for ( i = 0; i < N; i++) {
		for ( j = 0; j <= i; j++) {
			double sum = 0;
			int column = j + design;
			int anchor = i * mnVariables;
			for (k = 0; k < mnVariables; k++) {
				sum += BtPsiInv[anchor + k] * mIter.mB[k][column];
			}
			p_i[nCount] = sum;
			if (i == j) {
				sum += 1;
			}
			p[nCount++] = sum;
		}
	}
	
	SymmetricMatrix Q = QInv.i();
	LowerTriangularMatrix L = Cholesky(Q);
	LowerTriangularMatrix LInv = L.i();

	double QDetRoot = 1.0;
	for (i = 0; i < N; i++) {
		QDetRoot *= L[i][i];
	}
	
	//first copy mT into this symmetric matrix
	
	double* BtStarPsiInv = NULL;
	LowerTriangularMatrix LStar;
	LowerTriangularMatrix LM;
	SymmetricMatrix TInv;

	SymmetricMatrix T(N);
	double* Ltui = NULL;
	double* LtStarui = NULL;
	
	Ltui= new double[N];
	//-------------------------------------------------------------------

	double* xp= new double[mnVariables];
	double* theta_j = new double[N];
	double* deltaF_j = new double[N]; 
	double* alpha_ij = new double[mnSampleSize+1]; //alpha_ij
	double* ui = new double[N];
	double* uiStar = new double[N];
	double* nr = new double[N];

	double sum;
	//--------------------------------------Start of the configuration step
	for (i = 0; i < mnSampleSize; i++) { 	
		//now calculate x_i
		for (j = 0; j < mnVariables; j++) {
			double sum = mResX[j][i];
			for (k = design; k < design + N; k++) {
				sum += mIter.mB[j][k] * mIter.mF[k][i];		
			}
			xp[j] = sum;
		}
		
		//calculate alpha_ij's
		for (j = 0; j < N;j++) {
			double sum = 0;
			int anchor = j * mnVariables;
			for (k = 0; k < mnVariables; k++) {
				sum += BtPsiInv[anchor + k] * xp[k];
			}
			ui[j] = sum;
		}
		if (nongaussian) {
			
			//step 1: calculate alpha0
			//alpha_ij[0] = pdp->mAlpha * QDetRoot * exp(0.5 * DotProduct(ui,Q * ui));
			//transpose(L) * ui
			
			for (j = 0; j < N; j++) {
				Ltui[j] = 0;
			}
			p = L.Store(); nCount = 0;			
			for (j = 0; j < N; j++) {
				for (k = 0; k <= j; k++) {
					Ltui[k] +=  p[nCount++] * ui[j];
				}
			}

			double sum2 = 0;
			for (j = 0; j < N; j++) {
				sum2 += Ltui[j] * Ltui[j];
			}
			alpha_ij[0] = pdp->mAlpha * QDetRoot * exp(0.5 * sum2);
			

			//now the other alpha's
			int group = pdp->mGroupIndicators[i];
			sum = alpha_ij[0];
			for ( j = 0; j < pdp->mnGroups; j++) {					
				int nj = (int)pdp->mpGroups[j].size();
				if (j == group) {
					nj--;
				} 
				if (nj) {
					for ( k = design; k < design + N; k++) {
						theta_j[k - design] = mIter.mF[k][j]; 
					}

					//ColumnVector temp  = ui + 0.5 * theta_j - 0.5 * QInv * theta_j;
					//alpha_ij[j + 1] = nj * exp(DotProduct(theta_j, temp));


					//calculate LInv * theta_j
					for (m = 0; m < N; m++) {
						nr[m] = 0;
					}
					
					p = LInv.Store(); nCount = 0;			//LInv * theta_j
					for (m = 0; m < N; m++) {
						for (k = 0; k <= m; k++) {
							nr[m] += p[nCount++] * theta_j[k];
						}
					}
					double temp2 = 0;
					for (m = 0; m < N; m++) {
						temp2 += (2 * ui[m]  + theta_j[m]) * theta_j[m] - nr[m] * nr[m];
					}		
					alpha_ij[j + 1] = nj * exp(temp2 / 2);
					sum += alpha_ij[j + 1];
				} else {
					alpha_ij[j + 1] = 0;
				}
			}
		}
		//figure out which theta to pick
		int index = 0;
		if (nongaussian) {
			if (pdp != NULL) {
				double prob = mRandom.Next() * sum;
				double cumsum = 0;
				for ( j = 0; j <= pdp->mnGroups; j++) {
					cumsum += alpha_ij[j];
					if (prob < cumsum) {
						index = j;
						break;
					}
				}
			}
		}
	
		//index = 0;
		//generate new theta here only if index == 0;
		if (index == 0) {		//sample from a normal model
			pdp->AddGroup(i);
			
			//draw a new value for this new lambda/theta
			//theta_j =  M * ui + LM * nr;  theta_j = LM * transpose(LM) * ui + LM * nr = LM * (Ltui + nr);

			for ( k = 0; k < N; k++) {
				nr[k] = Ltui[k] + mNormal.Next();
			}
			//calculate L * (Ltui + nr);
			for (j = 0; j < N; j++) {
				theta_j[j] = 0;
			}
			
			
			p = L.Store();
			
			
			nCount = 0;			
			for (j = 0; j < N; j++) {
				for (k = 0; k <= j; k++) {
					theta_j[j] += p[nCount++] * Ltui[k];
				}
			}

			for (k = design; k < design + N; k++) {
				double temp = mIter.mF[k][i];
				int ind = k - design;
				mIter.mF[k][i] = theta_j[ind];
				theta_j[ind] =temp - theta_j[ind]; 
			}
			//now updating mResX
			for (j = 0; j < mnVariables; j++) {
				for (k = 0; k < N; k++) {
					mResX[j][i] += mIter.mB[j][k + design] * theta_j[k];
				}
			}

		} else {
			pdp->SwitchObservationToGroup(i,index-1); //minus one because of the one is located to alpha_0
		}
	}
	//--------------------------------------End of the configuration step

	//--------------------------------------Start of the updating theta step
	SymmetricMatrix ThetaSquare(N); ThetaSquare = 0;		//only for wishart
	for (j = 0; j < pdp->mnGroups;j++) {
		//calculate the group mean xi
		for ( i = 0; i < mnVariables; i++) { xp[i] = 0; } 
		set<int>::iterator it;
		for (it = pdp->mpGroups[j].begin(); it != pdp->mpGroups[j].end(); it++) {
			int column = *it;
			for (i = 0; i < mnVariables; i++) {
				xp[i] += mResX[i][column];
			}
		}
		int nGroupSize = (int)pdp->mpGroups[j].size();
		int nFirstSampleIndexInGroup = *(pdp->mpGroups[j].begin());
		for (i = 0; i < mnVariables; i++) {
			xp[i] /= nGroupSize;
			double sum = 0;
			for (k = design; k < design + N; k++) {
				sum += mIter.mB[i][k] * mIter.mF[k][nFirstSampleIndexInGroup];		
			}
			xp[i] += sum;
		}

		//draw a new theta_j
		for (i = 0; i < N;i++) {
			double sum = 0;
			int anchor = i * mnVariables;
			for (k = 0; k < mnVariables; k++) {
				sum += BtPsiInv[anchor + k] * xp[k];
			}
			ui[i] = nGroupSize * sum;
		}

		
		//theta_j =  Q * ui + L * nr / sqrt((double)nGroupSize); 

		//Get the updated Ltui + nr / sqrt((double)nGroupSize)
		for (i = 0; i < N; i++) {
			Ltui[i] = 0;
		}
		
		SymmetricMatrix QTemp_j;  
		
		QTemp_j << Id + QInv_j * nGroupSize;
		
	
		SymmetricMatrix QTemp_j_inv = QTemp_j.i();
		LowerTriangularMatrix L_j = Cholesky(QTemp_j_inv);

		//LowerTriangularMatrix LInv_j = Cholesky(QTemp_j);
		//LowerTriangularMatrix L_j = LInv_j.i();
		p = L_j.Store();
		
		nCount = 0;			
		for (i = 0; i < N; i++) {
			Ltui[i] = mNormal.Next();
			for (k = 0; k <= i; k++) {
				Ltui[k] += p[nCount++] * ui[i];
			}
		}
		//now times by L
		//calculate L * (Ltui + nr);
		for (i = 0; i < N; i++) {
			theta_j[i] = 0;
		}
		p = L_j.Store();
		nCount = 0;			
		for (i = 0; i < N; i++) {
			for (k = 0; k <= i; k++) {
				theta_j[i] += p[nCount++] * Ltui[k];
			}
		}
		
		//update mF and mResX
		for (it = pdp->mpGroups[j].begin(); it != pdp->mpGroups[j].end(); it++) {
			int obs = *it;
			for (k = design; k < design + N; k++) {
				double temp = mIter.mF[k][obs];
				mIter.mF[k][obs] = theta_j[k - design];
				deltaF_j[k - design] =temp - mIter.mF[k][obs]; 
			}
			//now updating mResX
			for (i = 0; i < mnVariables; i++) {
				for (k = 0; k < N; k++) {
					mResX[i][obs] += mIter.mB[i][k + design] * deltaF_j[k];
				}
			}
		}

	}
	//--------------------------------------End of the updating theta step

	//now update alpha
	pdp->UpdateAlpha();
	delete [] BtPsiInv; BtPsiInv = NULL; 
	delete [] xp; xp = NULL;
	delete [] theta_j; theta_j = NULL;
	delete [] deltaF_j; deltaF_j = NULL;
	delete [] alpha_ij; alpha_ij = NULL;
	delete [] ui; ui = NULL;
	delete [] Ltui; Ltui = NULL;
}

void Bfrm::SampleX(int factor, int design) {
	if (mbHasXMask) { 
		int nTotalFactor = design + mnYFactors + factor;	//total number of columns	
		//for all X variables
		int i,j;
		for (i = mnYFactors; i < mnVariables; i++) {
			int index = i - mnYFactors;
			for (j = 0; j < mnSampleSize; j++) {
				if (mWeight[j]) {
					double mu_zij = 0;
					for (int k = 0; k < nTotalFactor; k++) {
						mu_zij += mIter.mB[i][k] * mIter.mF[k][j];
					}
					if (mXMask[index][j]>0) { 	
						mIter.mX[index][j] = mu_zij + sqrt(mIter.mPsi[i]) * mNormal.Next(); 
					} else {
						mIter.mX[index][j] = mX[i][j];		
					}
					mResX[i][j] = mIter.mX[index][j] - mu_zij;
				} else {
					mIter.mX[index][j] =0;
				} 
			}
		}
	}
}

//update Z mResX(Z)
void Bfrm::SampleZ(int factor, int design) {
	int i,j;
	int nTotalFactor = design + mnYFactors + factor;	//total number of columns	

	//binary variables
	for (i = 0; i < mnVarBinary; i++) {
		for (j = 0; j < mnSampleSize; j++) {

			if (mWeight[j]) {
				//calculate Z_mu_i_j first
				double mu_zij = 0;
				for (int k = 0; k < nTotalFactor; k++) {
					mu_zij += mIter.mB[i][k] * mIter.mF[k][j];
				}

				if (mYMask[i][j]) { //missing Y, so impute it
					mIter.mZ[i][j] = mu_zij + mNormal.Next();	
				} else {
					//sample Z[i][j] according to Y (the first part of mX)
					mIter.mZ[i][j] = mPrior.SampleRestrictedNormal(mu_zij, 1, mX[i][j] - 0.5, 0);
				}

				//initilise residules, part 2: Z
				mResX[i][j] = mIter.mZ[i][j] - mu_zij;
			} else {
				mIter.mZ[i][j] = 0;
			}
		}
	}
	
	//categorical variables
	int nCount = mnVarBinary;
	for (i = mnVarBinary; i < mnVarBinary +  mnVarCategorical; i++) {
		int nCurrentCV = i - mnVarBinary;
		int nVar = mCategoryIndex[nCurrentCV][2] - mCategoryIndex[nCurrentCV][1];
		for (int nCat = 0; nCat < nVar; nCat++) {
			int nCurrent = nCount + nCat;
			for (j = 0; j < mnSampleSize; j++) {
				if (mWeight[j]) {
					//calculate Z_mu_i_j first
					double mu_zij = 0;
					for (int k = 0; k < nTotalFactor; k++) {
						mu_zij += mIter.mB[nCurrent][k] * mIter.mF[k][j];
					}
					
					//sample Z[i][j] according to Y (the first part of mX)
					//example
					//			0		1		2
					//0			-1		?		?
					//1			1		-1		?
					//2			1		1		-1
					//3			1		1		1
					if (mYMask[nCurrent][j]) {
						mIter.mZ[nCurrent][j] = mu_zij + mNormal.Next();
					} else {
						double nY = (int)mCategorical[nCurrentCV][j];
						if (nY == nCat) {				//-1		
							mIter.mZ[nCurrent][j] = mPrior.SampleRestrictedNormal(mu_zij, 1, -1, 0);  
						} else if (nY > nCat)	 {		//1
							mIter.mZ[nCurrent][j] = mPrior.SampleRestrictedNormal(mu_zij, 1,  1, 0);
						} else {						//?	
							mIter.mZ[nCurrent][j] = mu_zij + mNormal.Next();			//normal(nu_zij, 1);
						}
					}
					
					//initilise residules, part 2: Z
					mResX[nCurrent][j] = mIter.mZ[nCurrent][j] - mu_zij;
				} else {			//can be optimized
					mIter.mZ[nCurrent][j] = 0;
				}
			}
		}
		nCount += nVar;
	}

	//survival variables
	for (i = mnYFactors - mnVarContinuous - mnVarSurvival; i < mnYFactors - mnVarContinuous;  i++) {
		for (j = 0; j < mnSampleSize; j++) {
			if (mWeight[j]) {
				//calculate Z_mu_i_j first
				double mu_zij = 0;
				for (int k = 0; k < nTotalFactor; k++) {
					mu_zij += mIter.mB[i][k] * mIter.mF[k][j];
				}
				if (mYMask[i][j]) { 
					if (mYMask[i][j] != 1) {			//censored data
						mIter.mZ[i][j] = mPrior.SampleRestrictedNormal(mu_zij, mIter.mPsi[i], mX[i][j] + 1, mX[i][j]);
					} else {		//missing value, inpute it here
						mIter.mZ[i][j] = mPrior.SampleRestrictedNormal(mu_zij, mIter.mPsi[i], 1, 0);
					}
				} else {			//observated survival time
					mIter.mZ[i][j] = mX[i][j];
					
				}
				//initilise residules, part 2: Z
				mResX[i][j] = mIter.mZ[i][j] - mu_zij;
			} else {
				mIter.mZ[i][j] = 0;
			}
		}
	}
	
	//now all continuous Y variables
	for (i = mnYFactors - mnVarContinuous; i < mnYFactors; i++) {
		for (j = 0; j < mnSampleSize; j++) {
			if (mWeight[j]) {
				double mu_zij = 0;
				for (int k = 0; k < nTotalFactor; k++) {
					mu_zij += mIter.mB[i][k] * mIter.mF[k][j];
				}
				if (mYMask[i][j]) { 	
					mIter.mZ[i][j] = mu_zij + sqrt(mIter.mPsi[i]) * mNormal.Next(); 
				} else {
					mIter.mZ[i][j] = mX[i][j];		
				}
				mResX[i][j] = mIter.mZ[i][j] - mu_zij;
			} else {
				mIter.mZ[i][j] =0;
			}
		}
	}

}

//output: mA, mBz, mF
void Bfrm::Summarise(int factor, int design, int nmc, int joel, int mode) {
	int i,j;
	mMeanResult = mResult.Summarise(nmc, joel);
	if (!mode) {
		//transform H back
		if (design > 1) {
			//recover H
			for (i = 0; i < design; i++) {
				for (j = 0; j < mnSampleSize; j++) {
					mMeanResult.mF[i][j] = mH[j][i];			//assume H: mnSampleSize * design
				}
			}

			mIter.maTau = BfrmData::Var(mMeanResult.mF);
			for (i = 1; i < design; i++) {
				if (!mIter.maTau[i]) { mIter.maTau[i] = 1e-9; }
				for (j = 0; j < mnVariables; j++) {
					mMeanResult.mB[j][i] /= sqrt(mIter.maTau[i]);
				}
				mMeanResult.maTau[i] /= mIter.maTau[i];
			}
		}
	}
	ManageHelperVariables(false, factor);
}

