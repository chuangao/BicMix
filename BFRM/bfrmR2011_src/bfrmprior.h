// bfrmprior.h: interface for the BfrmPrior class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BFRMPRIOR_H)
#define AFX_BFRMPRIOR_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class BfrmPrior  
{
public:
	BfrmPrior();
	virtual ~BfrmPrior();
	void Destroy();
	void Initialise(Model& model, int yfactor, int TotalNewVariables);
	void InitialiseTau(Array1D<double>& tau,int nFactors); 
	void InitialiseRho(Array1D<double>& rho,int nFactors);
	//void InitialiseG(Array1D<double>& g);
	void InitialisePibIJ(Array2D<double>& pibIJ,int nFactors, int nVariables);
	void InitialisePsi(Array1D<double>& psi,int nVariables, int PsiMask, int SurvivalPsiMask); 
	void InitGammaPsi(int nSampleCount);
	double SampleTau(int factor, double m, double sumsquare);
	double SampleG(int sk, double sumsquare);
	double SamplePsi(double sumsquare, int index);
	double SampleRho(int factor, int m, int n);
	double SamplePiIJ(int factor, double bij, double rho);
	double SampleRestrictedNormal(double mu, double sigma2, double obs, double cutoff);
	double GetPriorMean(int index);
	double GetPriorVar(int index);

private:
	//prior info
	double mTauDesignControl[2];
	double maTauResponse[4][2];
	double maTauLatent[2];

	Gamma *mpGammaPsi;
	Gamma *mpGammaPsiSurvival;

	//optimization, Joe's model
	Gamma *mpGammaPia0;
	Gamma *mpGammaPia1;
	Gamma *mpGammaPib0;
	Gamma *mpGammaPib1;

	int mnPsiMask;
	int mnPsiMaskSurvival;
	Array1D<double> maMeanValue;
	Array1D<double> maVarValue;
	double maPi[3];

	Uniform mRandom;
	
protected:
	Array1D<double> maGammaB;
	Array1D<double> maGammaA;
	double maPsi[2];
	double maPsiSurvival[2];
	Array2D<double> maRho;
public:
	SymmetricMatrix mT;
};

#endif // !defined(AFX_BFRMPRIOR_H)
