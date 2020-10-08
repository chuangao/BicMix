// bfrmrawdata.cpp: implementation of the BfrmRawData class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////
#include "tnt.h"
using namespace TNT;
#include "stdafx.h"
#include "bfrmrawdata.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BfrmRawData::BfrmRawData()
{
	mnBMask = 0;
	mnBCarlosMask = 0;
	mnYFactors = 0;
	mnControls = 0;
	mbHasXMask = false;
}

BfrmRawData::~BfrmRawData()
{

}
//now extended to allow a strip
//2 2 2 1 0 0 0
//2	2 2 2 1 0 0
//2 2 2 2 2 1 0
//2 2 2 2 2 2 1
//2 2 2 2 2 2 2 
//2 2 2 2 2 2 2
int BfrmRawData::GetBMask(int i, int j, int design, int mode) {
	if (mode) {				//variable selection mode
		if (i < mnYFactors) {
			return -1;			//-1 means keep it unchanged
		}
	}
	if (mnBMask) {
		int jaddjust = j - design;
		if (i==jaddjust)  {
			return 1;
		} else if (i > jaddjust) {
			//mask out block coresponding to control--phenotype interaction to zero
			//Mike 08/21/2005
			if (mnControls > 0 && mnBMask == 2) {
				if (i < mnYFactors && j >= design - mnControls && j < design) {
					return 0;
				}
			} 
			//End of Mike 08/21/2005
			return 2;
			
		} else {		//mnBMask == 2
			if (mnBMask == 2) {
				if (i < mnYFactors && jaddjust >= mnYFactors) {
					return 2;
				}
			}
		}
		return 0;
	}
	return 2;
}

int BfrmRawData::GetBMaskAddjust(int j, int design) {
	int jaddjust = j - design;
	if (!mnBMask || jaddjust < 0) { jaddjust = 0;} 
	
	if (mnBMask == 2) {
		if ( jaddjust >= mnYFactors) {
			jaddjust = jaddjust - mnYFactors;
		}
	}

	//mask out block coresponding to control--phenotype interaction to zero
	//Mike 08/21/2005
	if (mnControls > 0 && mnBMask == 2) {
		if (j >= design - mnControls && j < design) {
			jaddjust = mnYFactors;
		}
	} 
	//End of Mike 08/21/2005

	return jaddjust;
}

bool BfrmRawData::InitialiseH(double* data, int nSampleSize, int design) {
	try {
		mH = Array2D<double>(nSampleSize, design);
		if (design <= 1) {
			mH = 1.0;
		} else {
			int index = 0;
			for (int j = 0; j < design; j++) {
				for (int i = 0; i < nSampleSize; i++) {
					mH[i][j] = data[index++];
				}
			}
		}
		return true;
	} catch (...) {
		return false;
	}
}

bool BfrmRawData::InitialiseH(string FileName, int nSampleSize, int design, int control) {
	mnControls = control;
	if (design <= 1) {
		mH = Array2D<double>(nSampleSize, design);
		mH = 1.0;
	} else {
		return LoadData(FileName, mH, nSampleSize, design);
	}
	return true;
}

bool BfrmRawData::ReformatData(int nbinary, int nCategorical, int nSurvival, int nContinuous) {
	
	//if we have survival variables
	int i;
	mnExtraVariables = 0;
	for (i = nbinary + nCategorical; i <  nbinary + nCategorical + nSurvival; i++) {
		for (int j  = 0; j < mnSampleSize; j++) {
			if (mWeight[j] > 0 ) {
				if (mX[i][j] > 0) {
					mX[i][j] = log(mX[i][j]);
				} else {
					std::cout << "Survival time must be positive!" << endl;
					return false;
				}
			}
		}
	}
	

	//figure out how many extra categarical variables we will need
	if (nCategorical > 0) {
		int nExtraVariables = 0;
		int* nmax= new int[nCategorical];
		for (i = nbinary; i < nbinary +  nCategorical; i++) {
			nmax[i-nbinary] = 0;
			for (int j = 0; j < mnSampleSize; j++) {
				if (mX[i][j] > nmax[i-nbinary]) {
					nmax[i-nbinary] = (int)mX[i][j];
				}
			}
			nExtraVariables += nmax[i-nbinary]; 
		}
		nExtraVariables -= nCategorical;
		Array2D<double> newX = Array2D<double>(mnVariables + nExtraVariables, mnSampleSize);
		mCategorical = Array2D<double>(nCategorical, mnSampleSize);
		mCategoryIndex = Array2D<int>(nCategorical, 3);
		int nCount = 0;
		for (i = 0; i < nbinary; i++) {		//binary variables
			for (int j = 0; j < mnSampleSize; j++) {
				newX[nCount][j] = mX[i][j];
			}
			nCount++;
		}

		for (i = nbinary; i < nbinary +  nCategorical; i++) {
			int nCurrentCV = i - nbinary;
			for (int j = 0; j < mnSampleSize; j++) {
				mCategorical[nCurrentCV][j] = mX[i][j];
				for (int k = 0; k < nmax[nCurrentCV]; k++) {
					newX[nCount + k][j] = (((int)mX[i][j]) == (k)) ? 1 : 0;
				}
			}
			mCategoryIndex[nCurrentCV][0] = i;
			mCategoryIndex[nCurrentCV][1] = nCount;		// [)
			nCount += nmax[nCurrentCV];
			mCategoryIndex[nCurrentCV][2] = nCount;
		}	
		//now all other variables
		for (i = nbinary +  nCategorical; i < mnVariables; i++) {
			for (int j = 0; j < mnSampleSize; j++) {
				newX[nCount][j] = mX[i][j];
			}
			nCount++;
		}
		mX = Array2D<double>(mnVariables + nExtraVariables, mnSampleSize);
		for (i = 0; i < mnVariables + nExtraVariables; i++) {
			for (int j = 0; j < mnSampleSize; j++) {
				mX[i][j] = newX[i][j];
			}
		}
		mnVariables += nExtraVariables;
		mnExtraVariables = nExtraVariables;
		delete [] nmax; 
	}
	mnVarBinary = nbinary;
	mnVarCategorical = nCategorical;
	mnVarSurvival = nSurvival;
	mnVarContinuous = nContinuous;
	mnYFactors = mnVarBinary + mnVarCategorical + mnExtraVariables + mnVarSurvival + mnVarContinuous;
	return true;
}

void BfrmRawData::ReformatYMask(Array2D<double>& RawYMask) {
	int i;
	int nCount = 0;
	mYMask = Array2D<double>(mnYFactors, mnSampleSize);
	for (i = 0; i < mnVarBinary; i++) {		//binary variables
		for (int j = 0; j < mnSampleSize; j++) {
			mYMask[nCount][j] = RawYMask[i][j];
		}
		nCount++;
	}
	
	for (i = mnVarBinary; i < mnVarBinary +  mnVarCategorical; i++) {
		int nCurrentCV = i - mnVarBinary;
		int nVar = mCategoryIndex[nCurrentCV][2] - mCategoryIndex[nCurrentCV][1];
		for (int j = 0; j < mnSampleSize; j++) {
			for (int k = 0; k < nVar; k++) {
				mYMask[nCount + k][j] = RawYMask[i][j];
			}
		}
		nCount += nVar;
	}	

	//now all other Y variables
	for (i = mnVarBinary +  mnVarCategorical; i < mnVarBinary +  mnVarCategorical + mnVarSurvival + mnVarContinuous; i++) {
		for (int j = 0; j < mnSampleSize; j++) {
			mYMask[nCount][j] = RawYMask[i][j];
		}
		nCount++;
	}
}

bool BfrmRawData::SaveData(string FileName) {
	return BfrmData::SaveData(FileName,mX);
}

bool BfrmRawData::InitialiseWeight(string FileName, int nSampleSize) {
	if (FileName == "") {
		mWeight = Array1D<double>(nSampleSize);
		mWeight = 1.0;
		mnLeaveout = 0;
	} else {
		if (LoadData(FileName,mWeight,nSampleSize)) {
			mnLeaveout = 0;
			for (int i = 0; i < nSampleSize; i++) {
				if (mWeight[i] == 0) {
					mnLeaveout++;
				}
			}
			std::cout << "Number of samples used " << (nSampleSize - mnLeaveout) << endl; 
			return true;
		}
		return false;
	}
	return true;
}

bool BfrmRawData::InitialiseWeight(double* data, int nSampleSize) {
	try {
		mWeight = Array1D<double>(nSampleSize);
		mnLeaveout = 0;
		for (int i = 0; i < nSampleSize; i++) {
			mWeight[i] = data[i];
			if (mWeight[i] == 0) {
				mnLeaveout++;
			}
		}
		std::cout << "Number of samples used " << (nSampleSize - mnLeaveout) << endl; 
		return true;
	} catch (...) {
		return false;
	}
	
}

bool BfrmRawData::InitialiseYMask(double* data) {
	//Only if Y exists
	if (mnYFactors > 0) {
		try {
			int i , j;
			mYMask = Array2D<double>(mnYFactors, mnSampleSize);
			if (data == NULL) { // no masks
				mYMask = 0;
			} else {
				if (!mnExtraVariables) {
					int index = 0;
					for (i = 0; i < mnYFactors; i++) {
						for (j = 0; j < mnSampleSize; j++) {
							mYMask[i][j] = data[index++];
						}
					}
				} else {
					Array2D<double> RawYMask(mnYFactors - mnExtraVariables,mnSampleSize);
					int index = 0;
					for (i = 0; i < mnYFactors; i++) {
						for (j = 0; j < mnSampleSize; j++) {
							RawYMask[i][j] = data[index++];
						}
					}
					ReformatYMask(RawYMask);
				}
				return true;
			}
		} catch (...) {
			return false;
		}
	}
	return true;
}

bool BfrmRawData::InitialiseXMask(string FileName) {
	if (FileName == "") { // no mask file
		mbHasXMask = false;
		return true;
	} else {
		mXMask = Array2D<double>(mnVariables - mnYFactors, mnSampleSize);
		mbHasXMask = LoadData(FileName, mXMask, mnVariables - mnYFactors, mnSampleSize);
		return mbHasXMask;
	}
	return true;
}

void BfrmRawData::InitiliseX() {
	if (mbHasXMask) {			//fill the missings with mean values
		int i;
		//need to take missing values into consideration
		for ( i = mnYFactors; i < mnVariables; i++) { 
			double sum =0;
			int nCount = 0;
			int j;
			for (j = 0; j < mnSampleSize; j++) { //if observed
				if ((mWeight[j] > 0) && (mXMask[i-mnYFactors][j] == 0)) {
					nCount ++ ;
					sum += mX[i][j];
				}
			}
			double dmean = sum / nCount;
			for (j = 0; j < mnSampleSize; j++) { //if missing
				if ((mWeight[j] > 0) && (mXMask[i-mnYFactors][j] > 0)) {
					mX[i][j] = dmean;
				}
			}
		}
	}
}

bool BfrmRawData::InitialiseYMask(string FileName) {
	//Only if Y exists
	if (mnYFactors > 0) {
		mYMask = Array2D<double>(mnYFactors, mnSampleSize);
		if (FileName == "") { // no mask file
			mYMask = 0;
		} else {
			if (!mnExtraVariables) {
				return LoadData(FileName, mYMask, mnYFactors, mnSampleSize);
			} else {
				Array2D<double> RawYMask;
				if (LoadData(FileName, RawYMask, mnYFactors - mnExtraVariables, mnSampleSize)) {
					ReformatYMask(RawYMask);
					return true;
				}
			}
			return false;
		}
	}
	return true;
}

bool BfrmRawData::LoadDataXY(double *data, int rows, int columns, int nBinary,
			int nCategorical, int nSurvival, int nContinuous) {
	try {
		mnVariables = rows;
		mnSampleSize = columns;
		mX = Array2D<double>(rows,columns);
		int index = 0;
		for (int j = 0; j < columns; j++) {
			for (int i = 0; i < rows; i++) {
				mX[i][j] = data[index++];
			}
		}
		//reformat survival and categorical data if necessary, mnVariables could be changed here
		return ReformatData(nBinary, nCategorical, nSurvival, nContinuous); 
	} catch (...) {
		return false;
	}
}

bool BfrmRawData::LoadDataXY(string FileName, int nVariables, int nSampleSize, int nBinary,
			int nCategorical, int nSurvival, int nContinuous) {
	mnVariables = nVariables;
	mnSampleSize = nSampleSize;
	if (!LoadData(FileName,mX,nVariables,nSampleSize)) { return false; }
	//reformat survival and categorical data if necessary, mnVariables could be changed here
	return ReformatData(nBinary, nCategorical, nSurvival, nContinuous); 
}

void BfrmRawData::SetBMasks(int nBMask, int nCarlosBMask) {
	mnBMask = nBMask;
	mnBCarlosMask = nCarlosBMask;
}

Array2D<double> BfrmRawData::RowMeanXY() {
	//calculate row means
	Array2D<double> Ones(mnSampleSize,1); 
	try {
		for (int i = 0; i < mnSampleSize; i++) {
			Ones[0][i] = mWeight[i];
		}
	} catch (...) {
		Ones = 1;
		mnLeaveout = 0;
	}
	Array2D<double> Mean = matmult(mX,Ones); Mean /= (double)(mnSampleSize-mnLeaveout);
	return Mean;
}

void BfrmRawData::SubtractSampleMeans()
{
	int i;
	if (!mbHasXMask) {
		Array2D<double> Ones(mnSampleSize,1);
		for (i = 0; i < mnSampleSize; i++) {
			Ones[0][i] = mWeight[i];
		}
		Array2D<double> Mean = matmult(mX,Ones); 
		Mean /=  (double)(mnSampleSize-mnLeaveout);
		mX -= matmult(Mean,transpose(Ones));
	} else {
		//need to take missing values into consideration
		for ( i = 0; i < mnVariables; i++) { //for each variable including response variables
			double sum =0;
			int nCount = 0;
			int j;
			for (j = 0; j < mnSampleSize; j++) {
				if ((mWeight[j] > 0) && (i < mnYFactors || ((i >= mnYFactors) && mXMask[i-mnYFactors][j] == 0))) {
					nCount ++ ;
					sum += mX[i][j];
				}
			}
			if (nCount <= 2) {
				std::cout << "At least two observed values are required for each variable" << endl;
				exit(1);
			}
			double dmean = sum / nCount;
			for (j = 0; j < mnSampleSize; j++) {
				if ((mWeight[j] > 0) && (i < mnYFactors || ((i >= mnYFactors) && mXMask[i-mnYFactors][j] == 0))) {
					mX[i][j] = mX[i][j] - dmean;
				} else {
					mX[i][j] = dmean;		//does not matter here, would be overwritten anyway
				}
			}
		}
	}
}
