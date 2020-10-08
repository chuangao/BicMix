// bfrmprior.cpp: implementation of the BfrmPrior class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////
#define WANT_STREAM                  // include.h will get stream fns
#define WANT_MATH                    // include.h will get math fns
#include "newmatap.h"                // need matrix applications
#include "newmatio.h"                // need matrix output routines

#include "newran.h"
#include "tnt.h"
using namespace TNT;

#include "Model.h"
#include "stdafx.h"
#include "util.h"

#include "bfrmprior.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BfrmPrior::BfrmPrior()
{
	mpGammaPsi = NULL;
	mpGammaPsiSurvival = NULL;
	
	mpGammaPia0 = NULL;
	mpGammaPia1 = NULL;
	mpGammaPib0 = NULL;
	mpGammaPib1 = NULL;

	//srand( (unsigned)time( NULL ) );
	//Random::Set(((double)rand()) / ((double)RAND_MAX));
}

void BfrmPrior::Destroy() {
	if (mpGammaPsi != NULL) {
		delete mpGammaPsi; mpGammaPsi = NULL;
		delete mpGammaPsiSurvival; mpGammaPsiSurvival = NULL;

		delete mpGammaPia0; mpGammaPia0 = NULL;
		delete mpGammaPia1; mpGammaPia1 = NULL;
		delete mpGammaPib0; mpGammaPib0 = NULL;
		delete mpGammaPib1; mpGammaPib1 = NULL;
	}
}
BfrmPrior::~BfrmPrior()
{
	Destroy();
}

void BfrmPrior::Initialise(Model& model, int yfactor, int TotalNewVariables) {
	int i;

	Destroy();

	int design = model.mnNDesignControls;
	int factor = model.mnNLatentFactors;
	int TotalFactors = design + yfactor + factor;

	//pi_j or r_j
	maRho = Array2D<double>(TotalFactors,3);
	double dp = model.mdPriorRhoMean;
	double dM = model.mdPriorRhoN;
	double A = dM * dp;
	double B = dM * ( 1- dp); 
	for (i = 0; i < TotalFactors; i++) {
		maRho[i][0]= A;
		maRho[i][1] = B;
		maRho[i][2] = dp;
	}
	if (design) {
		maRho[0][0] = 1;
		maRho[0][1] = 1;
		maRho[0][2] = 1;
	}
	maPi[0] = model.mdPriorPiMean * model.mdPriorPiN;
	maPi[1] = model.mdPriorPiN * (1 - model.mdPriorPiMean);
	maPi[2] = model.mdPriorPiMean;
	mpGammaPia0 = new Gamma(maPi[0]);
	mpGammaPia1 = new Gamma(maPi[0] + 1);
	mpGammaPib0 = new Gamma(maPi[1] + 1);
	mpGammaPib1 = new Gamma(maPi[1]);

	maMeanValue = Array1D<double>(TotalNewVariables);
	maVarValue = Array1D<double>(TotalNewVariables);
	
	for (i = 0; i < model.mnNVarBinary; i++) {
		maMeanValue[i] = 0.0;
		maVarValue[i] = 1.0;
	}

	//default for categorical is 0.5, need to be changed
	for (i = model.mnNVarBinary; i < yfactor - model.mnNVarContinuous - model.mnNVarSurvival; i++) {
		maMeanValue[i] = 0;
		maVarValue[i] = 1.0;
	}
	
	for (i = yfactor - model.mnNVarContinuous - model.mnNVarSurvival; i < yfactor - model.mnNVarContinuous; i++) {
		maMeanValue[i] = model.mdMeanSurvival;
		maVarValue[i] = model.mdVarSurvival;
	}

	for (i = yfactor - model.mnNVarContinuous; i < yfactor; i++) {
		maMeanValue[i] = model.mdMeanContinuous;
		maVarValue[i] = model.mdVarContinuous;
	}
	
	for (i = yfactor; i < TotalNewVariables; i++) {
		maMeanValue[i] = model.mdMeanValue;
		maVarValue[i] = model.mdVarValue;
	}

	for (i = 0; i < 2; i++) {
		//prior params for Gamma(apsi/2,bpsi/2) on ideo precs 
		maPsi[i] = model.mdaPriorPsi[i];
		maPsiSurvival[i] = model.mdaPriorSurvivalPsi[i];
		
		//beta proir for design factor
		mTauDesignControl[i] = model.mdaPriorTauDesign[i];
		//beta prior for Y factor
		for (int j = 0; j < 4; j++) {
			maTauResponse[j][i] = model.mdaPriorTauResponse[j][i];
		}
		//beta prior for latent factor
		maTauLatent[i] = model.mdaPriorTauLatent[i];
	}
	
	//now deal with gammaB and gammaA
	maGammaB = Array1D<double>(TotalFactors);
	maGammaA = Array1D<double>(TotalFactors);	maGammaA = 1.0;
	for (i = design + yfactor; i < TotalFactors; i++) {
		maGammaA[i] = maTauLatent[0];
		maGammaB[i] = maTauLatent[1];
	}
	
	//Contniuous Y factors if any
	for (i = design + yfactor - model.mnNVarContinuous; i < design + yfactor;  i++) {
		maGammaA[i] = maTauResponse[3][0];
		maGammaB[i] = maTauResponse[3][1];
	}
	
	//Survival Y factors if any
	for (i = design + yfactor - model.mnNVarContinuous - model.mnNVarSurvival; 
		i < design + yfactor - model.mnNVarContinuous;  i++) {
		maGammaA[i] = maTauResponse[2][0];
		maGammaB[i] = maTauResponse[2][1];
	}

	//Categorical Y factors if any
	for (i = design + model.mnNVarBinary; 
		i < design + yfactor - model.mnNVarContinuous - model.mnNVarSurvival;  i++) {
		maGammaA[i] = maTauResponse[1][0];
		maGammaB[i] = maTauResponse[1][1];
	}
	
	//Binary Y factors if any
	for (i = design; i < design + model.mnNVarBinary;  i++) {
		maGammaA[i] = maTauResponse[0][0];
		maGammaB[i] = maTauResponse[0][1];
	}
	
	//design factor if any
	for (i = 0; i < design; i++) {
		maGammaA[i] = mTauDesignControl[0];
		maGammaB[i] = mTauDesignControl[1];
	}

	
}

void BfrmPrior::InitialiseRho(Array1D<double>& rho,int nFactors) {
	rho = Array1D<double>(nFactors); 
	for (int i = 0; i < nFactors; i++) {
		rho[i] = maRho[i][2];
	}
}

void BfrmPrior::InitialisePibIJ(Array2D<double>& pibIJ,int nFactors, int nVariables) {
	pibIJ = Array2D<double>(nVariables, nFactors);
	for (int j = 0; j < nVariables; j++) {
		for (int i = 0; i < nFactors; i++) {
			pibIJ[j][i] = maRho[i][2];
		}
	}
}

void BfrmPrior::InitialiseTau(Array1D<double>& tau, int nFactors) {
	tau = Array1D<double>(nFactors);
	for (int i = 0; i < nFactors; i++) {	
		tau[i] = maGammaB[i] / maGammaA[i];	
	}
}

void BfrmPrior::InitialisePsi(Array1D<double>& psi,int nVariables, int PsiMask, int SurvivalPsiMask) {
	psi = Array1D<double>(nVariables);
	int i;
	for (i = 0; i < PsiMask; i++) {
		psi[i] = 1.0;
	}
	for (i = PsiMask; i < SurvivalPsiMask; i++) {
		psi[i] = maPsiSurvival[1] / maPsiSurvival[0];
	}
	for (i = SurvivalPsiMask; i < nVariables; i++) {
		psi[i] = maPsi[1] / maPsi[0];
	}
	mnPsiMask = PsiMask;
	mnPsiMaskSurvival = SurvivalPsiMask;
}

double BfrmPrior::SampleTau(int factor, double m, double sumsquare) {
	Gamma gtau(maGammaA[factor] + m / 2.0);
	return (1.0 / (gtau.Next() / (maGammaB[factor] + sumsquare / 2))); //gamma(a,b) = gamma(a,1) / b
}

double BfrmPrior::SamplePsi(double sumsquare, int index) {
	
	if (index < mnPsiMask) {
		return 1.0;
	} else if (index < mnPsiMaskSurvival) {
		return (1.0 / (mpGammaPsiSurvival->Next() * 2.0 / (maPsiSurvival[1]+sumsquare))); //gamma(a,b) = gamma(a,1) / b
	} else {
		return (1.0 / (mpGammaPsi->Next() * 2.0 / (maPsi[1]+sumsquare))); //gamma(a,b) = gamma(a,1) / b
	}
}

void BfrmPrior::InitGammaPsi(int nSampleCount) {
	if (mpGammaPsi != NULL) {
		delete mpGammaPsi; mpGammaPsi = NULL;
		delete mpGammaPsiSurvival; mpGammaPsiSurvival = NULL;
	} 
	mpGammaPsi = new Gamma((maPsi[0]+nSampleCount)/2.0);
	mpGammaPsiSurvival = new Gamma((maPsiSurvival[0]+nSampleCount)/2.0);
}

double BfrmPrior::SampleRho(int factor, int m, int n) { //sample from a beta distribution
	Gamma ga(m + maRho[factor][0]);
	Gamma gb(n + maRho[factor][1]);
	double ra = ga.Next();
	double rb = gb.Next();
	return (ra / (ra  + rb));	
}

double BfrmPrior::SamplePiIJ(int factor, double bij, double rho) {
	if (bij != 0) {
		double ra = mpGammaPia1->Next();
		double rb = mpGammaPib1->Next();
		return (ra / (ra  + rb));	
	} else {
		double dp =  mRandom.Next();
		rho = rho * ( 1 - maPi[2]) / ( 1 - rho * maPi[2]);
		if (dp < rho) {
			double ra = mpGammaPia0->Next();
			double rb = mpGammaPib0->Next();
			return (ra / (ra  + rb));	
		} else {
			return 0;
		}
	}
}

double BfrmPrior::GetPriorVar(int index) {
	return maVarValue[index];
}

double BfrmPrior::GetPriorMean(int index) {
	return maMeanValue[index];
}

double BfrmPrior::SampleRestrictedNormal(double mu, double sigma2, double obs, double cutoff) {
	double sigma = sqrt(sigma2);
	double dX = (cutoff - mu) / sigma;
	double dC;
	
	if (dX < -1e9)  { 
		dC = 0; 
	} else if (dX > 1e9) { 
		dC = 1; 
	} else {
		dC = util::normcdf(dX);
	}

	double dA, dB;
	if (obs > cutoff) {
		dA = dC;
		dB = 1.0;
	} else {
		dA = 0.0;
		dB = dC;
	}
	
	double dP = dA + (dB - dA) * mRandom.Next(); 
	if (dP < 0.000001) {
		dP = 0.000001;
	} else if (dP > 0.999999) {
		dP = 0.999999;
	}
	
	double retmu = mu + sigma * util::norminv(dP);
	if (obs > cutoff && retmu <=0) {
		retmu = 1e-9;
	} else if (obs < cutoff && retmu >= 0) {
		retmu = -1e-9;
	}
	return retmu;
}
