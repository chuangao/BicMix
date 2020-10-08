// bfrmrawdata.h: interface for the BfrmRawData class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BFRMRAWDATA_H)
#define AFX_BFRMRAWDATA_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "bfrmdata.h"

class BfrmRawData : public BfrmData  
{
public:
//	virtual void Copy(BfrmRawData & other, int VariablesIn, Array1D<int> & Index);
	BfrmRawData();
	virtual ~BfrmRawData();
	bool LoadDataXY(string FileName, int nVariables, int nSampleSize, int nBinary,
			int nCategorical, int nSurvival, int nContinuous);
	bool LoadDataXY(double *data, int rows, int columns, int nBinary,
			int nCategorical, int nSurvival, int nContinuous);
	bool SaveData(string FileName);
	void SubtractSampleMeans();
	bool InitialiseH(string FileName, int nSampleSize, int design, int control);
	bool InitialiseH(double* data, int nSampleSize, int design);
	bool InitialiseWeight(string FileName, int nSampleSize);
	bool InitialiseWeight(double* data, int nSampleSize);
	bool InitialiseYMask(string FileName);
	bool InitialiseXMask(string FileName);
	void InitiliseX();
	bool InitialiseYMask(double* data);

	void SetBMasks(int nBMask, int nCarlosBMask);

	int GetBMask(int i, int j, int design, int mode); 
	int GetBMaskAddjust(int j, int design);
	Array2D<double> RowMeanXY();
protected:
	//input
	Array2D<double> mX;				//Y and X matrix
	Array2D<double> mCategorical;				//Y and X matrix
	Array2D<double> mYMask;			//mask matrix for Y variables
	Array2D<double> mXMask;			//mask matrix for X variables
	bool mbHasXMask;				//a convinence variable
	int mnSampleSize;	//total number of samples
	int mnVariables;	//number of variables in X
	int mnYFactors;	//number of Y factors

	//the next variables are defined within ReformatData
	int mnVarBinary;		//number of binary variables in Y
	int mnVarCategorical;	//number of categorical variables in Y
	int mnVarSurvival;		//number of Survival variables in Y
	int mnVarContinuous;	//number of Continuous variables in X
	Array2D<int> mCategoryIndex;
	int mnExtraVariables;
	
	Array1D<double> mWeight;		//weight matrix
	int mnLeaveout;	//number of leave out samples

	//H related
	Array2D<double> mH;				//design matrix
	int mnControls;  //number of control factors

	void ReformatYMask(Array2D<double>& RawYMask);
	bool ReformatData(int nbinary, int nCategorical, int nSurvival, int nContinuous);  
	int mnBMask;	//indicator for B mask
	int mnBCarlosMask; // Carlos' 'include factor' mask

};

#endif // !defined(AFX_BFRMRAWDATA_H)
