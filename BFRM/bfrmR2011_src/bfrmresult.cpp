// bfrmresult.cpp: implementation of the BfrmResult class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#include "tnt.h"
using namespace TNT;
#include "stdafx.h"
#include "bfrmdata.h"
#include "bfrmresult.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
BfrmResult::BfrmResult()
{
	mnVariables = 0;
	mnFactors = 0;
	mnSamples = 0;
	mnYFactors = 0;
	mbHasXMask = 0;
}

BfrmResult::~BfrmResult()
{

}
void BfrmResult::Init(int nVariavles, int nFactors, int nSamples, int nYFactors, int joel, bool bHasXMask, int nDesignControls) {
	mnVariables = nVariavles;
	mnFactors = nFactors;
	mnSamples = nSamples;
	mnYFactors = nYFactors;
	mbHasXMask = bHasXMask;
	mnDesignControls = nDesignControls;
	mB = Array2D<double>(mnVariables,mnFactors); mB = 0; 
	mBz = Array2D<double>(mnVariables,mnFactors); mBz = 0;
	mF = Array2D<double>(mnFactors, mnSamples); mF = 0;
	maTau = Array1D<double>(mnFactors); maTau = 0;			
	mAlpha = Array1D<double>(1); mAlpha = 0;
//	mG = Array1D<double>(1); mG = 0;
	mPsi = Array1D<double>(mnVariables);	mPsi = 0;
	maRho = Array1D<double>(mnFactors); maRho = 0;
	mPostPib = Array2D<double>(mnVariables,mnFactors); mPostPib = 0;
	if (mnYFactors > 0) {
		mZ = Array2D<double>(mnYFactors, mnSamples); mZ = 0; 
	}
	if (joel) {
		maPi = Array2D<double>(mnVariables,mnFactors); maPi = 0;
	}
	if (mbHasXMask) {
		mX = Array2D<double>(mnVariables - mnYFactors, mnSamples); mX = 0;
	}
	int nResponseLatentFactors = mnFactors - mnDesignControls;
	if (nResponseLatentFactors > 0) {
		mT = Array2D<double>(nResponseLatentFactors, nResponseLatentFactors); mT= 0;
	}
}

//used in Carlos' initialModel function to identify founders
int BfrmResult::GetMaxBjIndex(int j, int design) {
	int index = j - design;
	double d = fabs(mB[j - design][j]);
	for (int i = j - design + 1; i < mnVariables; i++) {
		if (fabs(mB[i][j]) > d) {
			d = fabs(mB[i][j]);
			index = i;
		}
	}
	return index - mnYFactors;			//index in genes
}
int BfrmResult::GetSignificantLoadingCount(int nfactor, double dthreshold) {
	int nCount = 0; //don't exclude the response loadings for now
	for (int i = 0; i < mnVariables; i++) {
		if (mPostPib[i][nfactor] >= dthreshold) {
			nCount++;
		}
	}
	return nCount;
}

double BfrmResult::Get(int i, int j, int which)
{
	if (which == 0) { //mA,mB
		return mB[i][j];
	} else if (which == 1) { //mBz
		return mBz[i][j];
	} else if (which == 2) { //mF
		return mF[i][j];
	} else if (which == 3) { //mPsi
		return mPsi[i];
	} else if (which == 4) { //mtau
		return maTau[i];
	} else if (which == 5) { //mrho
		return maRho[i];
	} else if (which == 6) {//mPostPib
		return mPostPib[i][j];
	} else if (which == 7) {
		return maPi[i][j];
	} else {//mZ
		return mZ[i][j];
	}
}

void BfrmResult::Save(int joel) {
	BfrmData::SaveData("mA.txt",mB);
	BfrmData::SaveData("mF.txt",mF);
	BfrmData::SaveData("mPsi.txt", mPsi);
	BfrmData::SaveData("mTau.txt", maTau);
	BfrmData::SaveData("mAlpha.txt", mAlpha);
	BfrmData::SaveData("mPib.txt", maRho);
	BfrmData::SaveData("mPostPib.txt", mPostPib);
	if (mnYFactors > 0) {
		BfrmData::SaveData("mZ.txt", mZ);
	}
	
	if (mbHasXMask) {
		BfrmData::SaveData("mX.txt", mX);
	}
}

BfrmResult BfrmResult::Summarise(int nmc, int joel) {
	BfrmResult Mean;
	Mean.Init(mnVariables,mnFactors,mnSamples, mnYFactors, joel, mbHasXMask, mnDesignControls);
	for (int j = 0; j < mnFactors; j++) { 
		for (int i = 0; i < mnVariables; i++) {
			if (mBz[i][j] > 0) {
				Mean.mB[i][j] = mB[i][j] / mBz[i][j];
			} else {
				Mean.mB[i][j] = 0;
			}
		}
	}
	
	double dmc = nmc;
	Mean.mF += mF; Mean.mF /= dmc; 
	Mean.mPsi += mPsi; Mean.mPsi /=  dmc; 
	Mean.maTau += maTau; Mean.maTau /= dmc; 
	Mean.mAlpha += mAlpha; Mean.mAlpha /= dmc;
//	Mean.mG += mG; Mean.mG /= dmc;
	Mean.mBz += mBz; Mean.mBz /= dmc;
	Mean.maRho += maRho; Mean.maRho /= dmc;   
	Mean.mPostPib += mPostPib; Mean.mPostPib /= dmc;
	if (joel) {
		Mean.maPi += maPi; Mean.maPi /= dmc;
	}
	if (mnYFactors > 0) {
		Mean.mZ += mZ; Mean.mZ /= dmc;
	}
	if (mbHasXMask) {
		Mean.mX += mX; Mean.mX /= dmc;
	}
	if (mnFactors - mnDesignControls > 0) {
		Mean.mT += mT; Mean.mT /= dmc;
	}
	return Mean;
}

void BfrmResult::Accumulate(BfrmResult& Iter, int joel) {
	mF += Iter.mF; 
	mPsi += Iter.mPsi; 
	maTau += Iter.maTau;
	mAlpha+= Iter.mAlpha;
//	mG+= Iter.mG;
	maRho += Iter.maRho; 
	mPostPib += Iter.mPostPib;
	if (mnYFactors > 0) {
		mZ += Iter.mZ;
	}
	if (joel) { maPi += Iter.maPi; }
	if (mbHasXMask) {
		mX += Iter.mX;
	}
	
	for (int i = 0; i < mnVariables; i++) { 
		for (int j = 0; j < mnFactors ; j++) { 
			if (Iter.mB[i][j]) {
				mBz[i][j] += 1.0;
				mB[i][j] += Iter.mB[i][j];
			}
		}
	}
	if (mnFactors - mnDesignControls > 0) {
		mT += Iter.mT;
	}
}
