// bfrmevo.h: interface for the BfrmEvo class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BFRMEVO_H__8EB32FFC_629F_4AA1_A52C_0E5204D68845__INCLUDED_)
#define AFX_BFRMEVO_H__8EB32FFC_629F_4AA1_A52C_0E5204D68845__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "bfrm.h"

class BfrmEvo : public Bfrm  
{
public:
	void SwitchXVariables(int i, int j);
	void GetFounder(int k, int design);
	void LoadDataEvol(Model& model);
	void RemoveDesignEffect(Model& model);
	int AddVariables(Model& model);
	int AddFactors(Model& model);
	void Evol(Model& model);
	void CalculateExternalProbs(Model& model);
	void SelectData(int factor);
	BfrmEvo(double seed);
	void GetInitialModel(Model& model);
	int FindIncludeVariables(vector<int> & include, vector<int> & highscorefactors, Model& model);
	bool SaveVariablesInIndex(string FileName);
	virtual ~BfrmEvo();

protected:
	//track the variables that are present in the model
	int mnTotalVariables;
	int mnVariablesInCount;
	Array1D<int> maVariablesInIndex;
	Array2D<double> mXAll;
	Array2D<double> mXMaskAll;
	Array2D<double> mXAllCorrected;
	Array2D<double> mExternalProb;
	Array1D<int> maVariablesOutIndicator;

	double GetPijHat(double v, int Vin, int VTotal, double r, double Rhoj);

};

#endif // !defined(AFX_BFRMEVO_H__8EB32FFC_629F_4AA1_A52C_0E5204D68845__INCLUDED_)
