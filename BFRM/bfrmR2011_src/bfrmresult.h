// BfrmResult.h: interface for the BfrmResult class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BFRMRESULT_H)
#define AFX_BFRMRESULT_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class BfrmResult  
{
public:
	BfrmResult();
	virtual ~BfrmResult();
	void Init(int nVariavles, int nFactors, int nSamples, int nYFactors, int joel, bool bHasXMask, int nDesignControls);
	double Get(int i, int j, int which);
	void Save(int joel);
	void Accumulate(BfrmResult& Iter, int joel);
	int GetSignificantLoadingCount(int nfactor, double dthreshold);
	BfrmResult Summarise(int nmc, int joel);
	int GetMaxBjIndex(int j, int design);

	Array2D<double> mB;			
	Array2D<double> mBz;
	Array2D<double> mF;
	Array1D<double> maTau;			
	Array1D<double> mAlpha;	
	//Array1D<double> mG;
	Array1D<double> mPsi;	
	Array1D<double> maRho;
	Array2D<double> maPi;
	Array2D<double> mPostPib;	
	Array2D<double> mT;
	
	Array2D<double> mZ;
	Array2D<double> mX;

	int mnVariables;
	int mnFactors;
	int mnSamples;
	int mnYFactors;
	bool mbHasXMask;
	int mnDesignControls;
};
#endif // !defined(AFX_BFRMRESULT_H)
