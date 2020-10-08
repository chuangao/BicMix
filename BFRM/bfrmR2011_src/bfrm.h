// bfrm.h: interface for the Bfrm class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#if !defined(_BFRM_H)
#define _BFRM_H

#if _MSC_VER > 1000
#pragma once
#endif 

class Model;
class Bfrm  : public BfrmRawData
{
public:
	Bfrm();
        Bfrm(double seed);
	virtual ~Bfrm();
	virtual void Initialise(Model& model, int FindFounder);
	virtual void LoadData(Model& model);
	void Run(int nStartIter, Model& model);
	virtual void Iterate(int factor, int design, int joel, int mode, double alpha, int nongaussian, 
		int iter);
	void SampleF(int factor, int design, mdp* pdp, int iter, int nongaussian);
	void SampleF(int factor, int design);
	virtual void SampleB(int factor, int design, int joel, int mode);
	void SamplePsi();
	void SampleZ(int factor, int design);
	void SampleX(int factor, int design);
	void Accumulate(int factor, int design, int joel);
	void Summarise(int factor, int design, int nmc, int joel, int mode);
	void Save(int joel);
	void SaveHist(int joel, int append);
	void RankVariables(int nSkip1, int nSkip2, int VarInCount, bool screener);

	void SaveIter(int i);
	BfrmResult mMeanResult;
	BfrmPrior mPrior;
	virtual void SampleBij(int i, int j, double v, int design, double oddratio, 
		double& newv, double meanvalue, double &postpi, int specialcase, int joel,
		int lastfactor, int mode);
	void CalculateMeanVarBJ(int j, double v, int nStart, int nEnd);
	void BUpdateResX(int nnz);
	void BUpdateResX(int nStart, int nEnd);
	void SamplePsi(int nStart, int nEnd);

	double* mpaMeanBJ;
	double* mpaVarBJ;

protected:
	//working variables for optimization
	int* mpaibj;
	double* mpafj;
	double* mpabj;
 
	Uniform mRandom;
	Normal mNormal;
	Array2D<double> mResX;

protected:
	void ManageHelperVariables(bool allocate, int factor);
	void CalculateMeanVarBJ(int j, double v);
	
	BfrmResult mResult;
	BfrmResult mIter;
	vector<int> maOriginalIndex;
	mdp mDP;
	int mnAppend; 
	
public:
	double mdSparsity;
	vector<int> maOrderedVariables;
	vector<double> maOrderedScores;

};

#endif 
