// bfrmevo.cpp: implementation of the BfrmEvo class.
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

#include <math.h>
#include "stdafx.h"
#include "mdp.h"
#include "Model.h"
#include "bfrmrawdata.h"
#include "bfrmprior.h"

#include "bfrmresult.h"
#include "bfrm.h"
#include "util.h"

#include "bfrmevo.h"



//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BfrmEvo::BfrmEvo(double seed):Bfrm(seed)
{

}

BfrmEvo::~BfrmEvo()
{

}

void BfrmEvo::LoadDataEvol(Model& model) {
	if (model.mnEvol) {
		std::cout << "Evolving model is activated" << endl;
	}

	int i;

	//keep a copy of all data
	mnTotalVariables = mnVariables;
	mXAll = mX.copy();
	if (mbHasXMask) {
		mXMaskAll = mXMask.copy();
	}
	
	//record variables initially in
	mnVariablesInCount = model.mnEvolVarIn;
	BfrmData::LoadData(model.mstrEvolVarInFile, maVariablesInIndex, mnVariablesInCount, mnTotalVariables);
	maVariablesOutIndicator = Array1D<int>(mnTotalVariables - mnYFactors); maVariablesOutIndicator = 1;
	for (i = 0; i < mnVariablesInCount; i++) {
		maVariablesInIndex[i] -= 1;
		maVariablesOutIndicator[maVariablesInIndex[i]] = 0;
	}
	
	SelectData(model.mnNLatentFactors);
}

void BfrmEvo::SelectData(int factor) {
	int i,j;
	mnVariables = mnYFactors + mnVariablesInCount;
	mX = Array2D<double>(mnVariables, mnSampleSize);
	if (mbHasXMask) {
		mXMask = Array2D<double>(mnVariablesInCount, mnSampleSize);
	}
	for (j = 0; j < mnSampleSize; j++) {
		for (i = 0; i < mnYFactors; i++) {
			mX[i][j] = mXAll[i][j];
		}
		for (i = 0; i < mnVariablesInCount; i++) {
			mX[i+mnYFactors][j] = mXAll[mnYFactors + maVariablesInIndex[i]][j];
		}
		if (mbHasXMask) {
			for (i = 0; i < mnVariablesInCount; i++) {
				mXMask[i][j] = mXMaskAll[maVariablesInIndex[i]][j];
			}
		}
	}
	std::cout << "Selecting data......Done: p = "<< mnVariablesInCount << " variables, k =  " << factor
		<< " factors" << endl;

}

//basic duplicates FindIncludeVariables
void BfrmEvo::CalculateExternalProbs(Model& model) {
	
	int i,j;
	ColumnVector Pi_j_hat_ratio(model.mnNDesignControls + mnYFactors + model.mnNLatentFactors);
	//for each factor
	for (i = model.mnNDesignControls; i < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; i++) {
		Pi_j_hat_ratio(i+1) = GetPijHat(model.mdPriorRhoN,mnVariablesInCount + mnYFactors,
			mnTotalVariables,model.mdPriorPiMean,mMeanResult.maRho[i]);
		Pi_j_hat_ratio(i+1) = Pi_j_hat_ratio(i+1) / (1 - Pi_j_hat_ratio(i+1));
	}
	mExternalProb = Array2D<double>(mnTotalVariables - mnYFactors - mnVariablesInCount, mnYFactors + model.mnNLatentFactors + 1);
	int nVariablesCount = 0;
	if (mbHasXMask) { 
		//get the scores for all genes
		for (i = mnYFactors; i < mnTotalVariables; i++) {
			int g = i - mnYFactors;  //gene i
			if (maVariablesOutIndicator[g]) {
				mExternalProb[nVariablesCount][0] = g + 1;
				//for each factor
				for (j = model.mnNDesignControls; j < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; j++) {
					//need to calculate fdot for each gene and factor
					 //(1)caculate beta_g_hat
					double fdot = 0;
					double betag_hat = 0;
					int k;
					for (k = 0; k < mnSampleSize; k++) {
						if (mXMaskAll[g][k] == 0) {
							fdot += mMeanResult.mF[j][k] * mMeanResult.mF[j][k]; 
							betag_hat += mMeanResult.mF[j][k] * mXAllCorrected[i][k];
						}
					}
					double fdot_inv = (fdot > 0)? (1.0 / fdot):0;
					betag_hat *= fdot_inv;

					//(2) caculate psi_g_hat
					double psig_hat = 0;
					double psig_null_hat = 0;
					int nValidSamples = 0;
					for (k = 0; k < mnSampleSize; k++) {
						if (mXMaskAll[g][k] == 0) {
							double temp = mXAllCorrected[i][k]  - betag_hat * mMeanResult.mF[j][k];
							psig_hat += temp * temp;
							temp = mXAllCorrected[i][k];
							psig_null_hat += temp * temp;
							nValidSamples++;
						}
					}

					psig_hat = psig_hat  / (nValidSamples - 1);
					psig_null_hat = psig_null_hat  / nValidSamples;


					double w = mMeanResult.maTau[j] + psig_hat / fdot;
					double v = psig_null_hat / fdot;


					//log score
					double log_score = 0.5 * log(v) - 0.5 *log(w) + log(Pi_j_hat_ratio(j+1));
					log_score += -0.5 * betag_hat * betag_hat * ((v-w) / (v * w));

					double dScore = exp(log_score); dScore = dScore / ( 1 + dScore);
					mExternalProb[nVariablesCount][j - model.mnNDesignControls + 1] = dScore;

				}
				nVariablesCount++;
				
			} 

		}
	} else {
		ColumnVector fdot(model.mnNDesignControls + mnYFactors + model.mnNLatentFactors);
		ColumnVector fdot_inv(model.mnNDesignControls + mnYFactors + model.mnNLatentFactors);
		//for each factor
		for (i = model.mnNDesignControls; i < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; i++) {
			fdot(i+1) = 0;
			for (j = 0; j < mnSampleSize; j++) {
				fdot(i+1) += mMeanResult.mF[i][j] * mMeanResult.mF[i][j]; 
			}
			fdot_inv(i+1) = 0;
			if (fdot(i+1) > 0) {
				fdot_inv(i+1) = 1.0 / fdot(i+1);
			}
		}

		//get the scores for all genes
		for (i = mnYFactors; i < mnTotalVariables; i++) {
			int g = i - mnYFactors;  //gene i
			if (maVariablesOutIndicator[g]) {
				mExternalProb[nVariablesCount][0] = g + 1;
				for (j = model.mnNDesignControls; j < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; j++) {

				     //(1)caculate beta_g_hat
					double betag_hat = 0;
					int k;
					for (k = 0; k < mnSampleSize; k++) {
						betag_hat += mMeanResult.mF[j][k] * mXAllCorrected[i][k];
					}
					betag_hat *= fdot_inv(j + 1);

					//(2) caculate psi_g_hat
					double psig_hat = 0;
					for (k = 0; k < mnSampleSize; k++) {
						double temp = mXAllCorrected[i][k]  - betag_hat * mMeanResult.mF[j][k];
						psig_hat += temp * temp;
					}
					psig_hat = psig_hat  / (mnSampleSize - 1);
					
					double psig_null_hat = 0;
					for (k = 0; k < mnSampleSize; k++) {
						double temp = mXAllCorrected[i][k];
						psig_null_hat += temp * temp;
					}

					psig_null_hat = psig_null_hat  / mnSampleSize;

					double w = mMeanResult.maTau[j] + psig_hat / fdot(j+ 1);
					double v = psig_null_hat / fdot(j+1);


					//log score
					double log_score = 0.5 * log(v) - 0.5 *log(w) + log(Pi_j_hat_ratio(j+1));
					log_score += -0.5 * betag_hat * betag_hat * ((v-w) / (v * w));
				
					double dScore = exp(log_score); dScore = dScore / ( 1 + dScore);
					mExternalProb[nVariablesCount][j - model.mnNDesignControls + 1] = dScore;
				}
				nVariablesCount++;
			} 
		}
	}
	BfrmData::SaveData("mExternalProb.txt",mExternalProb);
}
int BfrmEvo::FindIncludeVariables(vector<int> & include, vector<int> & highscorefactors, Model& model) {
	include.clear();
	highscorefactors.clear();
	if (mnVariablesInCount == maVariablesOutIndicator.dim()) { return 0;}
	int i,j;

	double *pScores = new double[mnTotalVariables - mnYFactors];
	int *pHighScorsFactors = new int[mnTotalVariables - mnYFactors];
	

	//create a factor mask to decide which factors will be used for variable selection
	int  *pFactorIndicator = new int[mnYFactors + model.mnNLatentFactors];
	for (i = 0; i < mnYFactors + model.mnNLatentFactors; i++) {
		pFactorIndicator[i] = 1;
	}
	if (model.mnEvolMaximumVariablesPerFactor > 0) {
		double dthreshold = model.mdEvolVariableThreshold;
		int nStartIndex = model.mnNDesignControls+mnYFactors;
		//for each latent factor
		for (i = nStartIndex; i < nStartIndex + model.mnNLatentFactors; i++) {			
			int nCount = mMeanResult.GetSignificantLoadingCount(i,dthreshold);
			if (nCount >= model.mnEvolMaximumVariablesPerFactor) {
				pFactorIndicator[i - nStartIndex + mnYFactors] = 0;
			}
		}
	}
	//------------------------------------------------------------------------------------


	ColumnVector Pi_j_hat_ratio(model.mnNDesignControls + mnYFactors + model.mnNLatentFactors);
	//for each factor
	for (i = model.mnNDesignControls; i < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; i++) {
		Pi_j_hat_ratio(i+1) = GetPijHat(model.mdPriorRhoN,mnVariablesInCount + mnYFactors,
			mnTotalVariables,model.mdPriorPiMean,mMeanResult.maRho[i]);
		Pi_j_hat_ratio(i+1) = Pi_j_hat_ratio(i+1) / (1 - Pi_j_hat_ratio(i+1));
	}
	if (mbHasXMask) { 
		//get the scores for all genes
		for (i = mnYFactors; i < mnTotalVariables; i++) {
			int g = i - mnYFactors;  //gene i
			pScores[g] = -1e9;
			pHighScorsFactors[g] = 0;
			if (maVariablesOutIndicator[g]) {
				//for each factor
				for (j = model.mnNDesignControls; j < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; j++) {
					if (pFactorIndicator[j - model.mnNDesignControls] > 0) {
						//need to calculate fdot for each gene and factor
						 //(1)caculate beta_g_hat
						double fdot = 0;
						double betag_hat = 0;
						int k;
						for (k = 0; k < mnSampleSize; k++) {
							if (mXMaskAll[g][k] == 0) {
								fdot += mMeanResult.mF[j][k] * mMeanResult.mF[j][k]; 
								betag_hat += mMeanResult.mF[j][k] * mXAllCorrected[i][k];
							}
						}
						double fdot_inv = (fdot > 0)? (1.0 / fdot):0;
						betag_hat *= fdot_inv;

						//(2) caculate psi_g_hat
						double psig_hat = 0;
						double psig_null_hat = 0;
						int nValidSamples = 0;
						for (k = 0; k < mnSampleSize; k++) {
							if (mXMaskAll[g][k] == 0) {
								double temp = mXAllCorrected[i][k]  - betag_hat * mMeanResult.mF[j][k];
								psig_hat += temp * temp;
								temp = mXAllCorrected[i][k];
								psig_null_hat += temp * temp;
								nValidSamples++;
							}
						}

						psig_hat = psig_hat  / (nValidSamples - 1);
						psig_null_hat = psig_null_hat  / nValidSamples;


						double w = mMeanResult.maTau[j] + psig_hat / fdot;
						double v = psig_null_hat / fdot;


						//log score
						double log_score = 0.5 * log(v) - 0.5 *log(w) + log(Pi_j_hat_ratio(j+1));
						log_score += -0.5 * betag_hat * betag_hat * ((v-w) / (v * w));

						if (log_score > pScores[g]) {
							pScores[g] = log_score;
							pHighScorsFactors[g] = j - model.mnNDesignControls + 1;
						}
					}

				}
			} 
		}
	} else {
		ColumnVector fdot(model.mnNDesignControls + mnYFactors + model.mnNLatentFactors);
		ColumnVector fdot_inv(model.mnNDesignControls + mnYFactors + model.mnNLatentFactors);
		//for each factor
		for (i = model.mnNDesignControls; i < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; i++) {
			fdot(i+1) = 0;
			for (j = 0; j < mnSampleSize; j++) {
				fdot(i+1) += mMeanResult.mF[i][j] * mMeanResult.mF[i][j]; 
			}
			fdot_inv(i+1) = 0;
			if (fdot(i+1) > 0) {
				fdot_inv(i+1) = 1.0 / fdot(i+1);
			}
		}

		//get the scores for all genes
		for (i = mnYFactors; i < mnTotalVariables; i++) {
			int g = i - mnYFactors;  //gene i
			pScores[g] = -1e9;
			pHighScorsFactors[g] = 0;
			if (maVariablesOutIndicator[g]) {
				for (j = model.mnNDesignControls; j < model.mnNDesignControls + mnYFactors + model.mnNLatentFactors; j++) {
					if (pFactorIndicator[j - model.mnNDesignControls] > 0) {
						 //(1)caculate beta_g_hat
						double betag_hat = 0;
						int k;
						for (k = 0; k < mnSampleSize; k++) {
							betag_hat += mMeanResult.mF[j][k] * mXAllCorrected[i][k];
						}
						betag_hat *= fdot_inv(j + 1);

						//(2) caculate psi_g_hat
						double psig_hat = 0;
						for (k = 0; k < mnSampleSize; k++) {
							double temp = mXAllCorrected[i][k]  - betag_hat * mMeanResult.mF[j][k];
							psig_hat += temp * temp;
						}
						psig_hat = psig_hat  / (mnSampleSize - 1);
						
						double psig_null_hat = 0;
						for (k = 0; k < mnSampleSize; k++) {
							double temp = mXAllCorrected[i][k];
							psig_null_hat += temp * temp;
						}

						psig_null_hat = psig_null_hat  / mnSampleSize;

						double w = mMeanResult.maTau[j] + psig_hat / fdot(j+ 1);
						double v = psig_null_hat / fdot(j+1);


						//log score
						double log_score = 0.5 * log(v) - 0.5 *log(w) + log(Pi_j_hat_ratio(j+1));
						log_score += -0.5 * betag_hat * betag_hat * ((v-w) / (v * w));
					
						if (log_score > pScores[g]) {
							pScores[g] = log_score;
							pHighScorsFactors[g] = j - model.mnNDesignControls+1;
						}
					}

				}
			} 
		}
	}
	//order genes by scores
	maOrderedVariables.clear();
	maOrderedScores.clear();

	util::sort(mnTotalVariables - mnYFactors, pScores, maOrderedVariables, -1); //descending
	
	for (i = 0; i < mnTotalVariables - mnYFactors; i++) {
		maOrderedScores.push_back(pScores[i]);
		maOrderedVariables[i] = maOrderedVariables[i] - 1;
	}
	delete [] pScores;
	delete [] pFactorIndicator;
	

	double log_threshold = log(model.mdEvolVariableThreshold / (1 - model.mdEvolVariableThreshold));
	int nCount = 0; 
	for (i = 0; i < (int)maOrderedScores.size(); i++) {
		//if ((maOrderedScores[i]) > model.mdEvolVariableThreshold) {
		if ((maOrderedScores[i]) > log_threshold) {
			include.push_back(maOrderedVariables[i]);
			highscorefactors.push_back(pHighScorsFactors[maOrderedVariables[i]]);
			nCount++;
			std::cout << (maOrderedVariables[i] + 1)  << '\t' << maOrderedScores[i]  <<  '\t' << highscorefactors[i] << endl;
		}
		if (nCount >= model.mnEvolMaximumVariablesPerIteration) break;
	}
	delete [] pHighScorsFactors;
	return (int)include.size();
	
}

void BfrmEvo::GetInitialModel(Model& model) {
	std::cout << "Searching for an initial model..." << endl;
	int K = model.mnNLatentFactors;
	for (int k = 0; k < K; k++) {
		std::cout << "Searching for founder variable #" << (k+1) << endl;
		model.mnNLatentFactors = k + 1;
		Initialise(model,1);
		Run(0, model);
		GetFounder(k,model.mnNDesignControls);
	}
	std::cout << "Fitting the initial model..." << endl;
	model.mnNLatentFactors = K;
	Initialise(model, 0);
	Run(0, model);
	std::cout << "Fitting the initial model......Done" << endl;
}

//update mnVariablesInIndex and mX , mXMask matrices
void BfrmEvo::SwitchXVariables(int i, int j) {
	int IndexTemp = maVariablesInIndex[i];
	maVariablesInIndex[i] = maVariablesInIndex[j];
	maVariablesInIndex[j] = IndexTemp;
	
	int IndexI = i + mnYFactors;
	int IndexJ = j + mnYFactors;
	for (int k = 0; k < mnSampleSize; k++) {
		double temp = mX[IndexI][k];
		mX[IndexI][k] = mX[IndexJ][k];
		mX[IndexJ][k] = temp;
		if (mbHasXMask) {
			temp = mXMask[i][k];
			mXMask[i][k] = mXMask[j][k];
			mXMask[j][k] = temp;
		}
	}
}

void BfrmEvo::GetFounder(int k, int design) {
	int j = mMeanResult.GetMaxBjIndex(design + mnYFactors + k, design);
	SwitchXVariables(k,j);
}

double BfrmEvo::GetPijHat(double v, int Vin, int VTotal, double r, double Rhoj) {
	return r * Rhoj * (v + Vin) / (v + VTotal);
}

void BfrmEvo::RemoveDesignEffect(Model& model) {
	
	mXAllCorrected = mXAll.copy();
	//now reomve the design effects
	if (model.mnNDesignControls> 0) {
		int i,j,k;
		if (mbHasXMask) {
			for (i = 0; i < mnTotalVariables; i++) {
				//in this case, need to calculate HAT matrix for each variable
				//may not be necessary for response variables, but just to be consistent

				//figure out the number of observed data for this variable
				int nObserved = 0;
				for (j = 0; j < mnSampleSize; j++) {
					if ( i < mnYFactors) {
						if (mYMask[i][j] == 0) { nObserved++; }
					} else {
						if (mXMaskAll[i-mnYFactors][j] == 0) { nObserved++; }
					}
				}
				if (nObserved < 2) {
					std::cout << "Two few observed data!" << endl;
					exit(1);
				}

				//now get the design matrix for observed data only
				bool bObserved;
				Matrix FT(nObserved, model.mnNDesignControls);
				int nCount = 0;
				for (j = 0; j < mnSampleSize; j++) {
					bObserved = (i < mnYFactors) ? (mYMask[i][j] == 0) : (mXMaskAll[i-mnYFactors][j] == 0);
					if (bObserved) {
						for (int k = 0; k < model.mnNDesignControls; k++) {
							FT[nCount][k] = mIter.mF[k][j];
						}
						nCount++;
					}
				}

				//calculate the HAT for this specific variable
				Matrix HAT = FT * (FT.t() * FT).i() * FT.t(); //nObserved by nObserved
				for (j = 0; j < nObserved; j++) {
					HAT[j][j] = HAT[j][j] - 1;
				}

				//do the correction
				nCount = 0;
				for (j = 0; j < mnSampleSize; j++) {
					bObserved = (i < mnYFactors) ? (mYMask[i][j] == 0) : (mXMaskAll[i-mnYFactors][j] == 0);
					if (bObserved) {
						double sum = 0;
						int nCount1 = 0;
						for (k=0; k< mnSampleSize; k++) {
							bObserved = (i < mnYFactors) ? (mYMask[i][k] == 0) : (mXMaskAll[i-mnYFactors][k] == 0);
							if (bObserved) {
								sum += HAT[nCount][nCount1] * mXAll[i][k];
								nCount1++;
							}
						}
						mXAllCorrected[i][j] = -sum;
						nCount++;
					}
				}
			}
			
			
		} else {
			Matrix FT(mnSampleSize, model.mnNDesignControls);
			for (i = 0; i < mnSampleSize; i++) {
				for (j = 0; j < model.mnNDesignControls; j++) {
					FT(i+1, j+1) = mIter.mF[j][i];
				}
			}
			Matrix HAT = FT * (FT.t() * FT).i() * FT.t(); 
			for (i = 1; i <= mnSampleSize; i++) {
				HAT(i,i) = HAT(i,i) - 1;
			}
			for (i = 0; i < mnTotalVariables; i++) {
				for (j = 1; j <= mnSampleSize; j++) {
					double sum = 0;
					for (k=1; k<= mnSampleSize; k++) {
						sum += HAT(j,k) * mXAll[i][k-1];
					}
					mXAllCorrected[i][j-1] = -sum;
				}
			}
		}
	}
		
}

int BfrmEvo::AddVariables(Model& model) {
	//selecting new variables
	int nNewVariables = 0;
	vector<int> IncludeVariables;
	vector<int> HighScoreFactors;
	if (mnVariablesInCount < model.mnEvolMaximumVariables) {
		std::cout << "Selecting new variables..." << endl << endl;
		nNewVariables = FindIncludeVariables(IncludeVariables, HighScoreFactors,model);
		if (nNewVariables > 0) {
			//a bit detour here according to MW's May 26 request
			ofstream theFile;
			theFile.open( "mEvolInfo.txt", ios_base::app);
			if  (theFile.fail()) {
				std::cout << "Failed to create file" << "mEvolInfo.txt" << endl;
				exit(1);
			}
			for (int i = 0; i < nNewVariables; i++) {
				theFile << model.mnNLatentFactors << "\t" << (IncludeVariables[i] + 1) << "\t" << HighScoreFactors[i] << endl;
			}
			theFile.close();
		}
		
	} else {
		nNewVariables = 0;
		IncludeVariables.clear();
	}

	if (nNewVariables) {
		std::cout << "Number of newly selected variables = " << nNewVariables <<endl;
		int i;
		for (i = 0; i < (int)IncludeVariables.size(); i++) {
			maVariablesInIndex[mnVariablesInCount++] = IncludeVariables[i];
			maVariablesOutIndicator[IncludeVariables[i]] = 0;
			if (mnVariablesInCount >= model.mnEvolMaximumVariables) break;
		}
		if (IncludeVariables.size()) {
			SelectData(model.mnNLatentFactors);
		}
	}
	return nNewVariables;
}

int BfrmEvo::AddFactors(Model& model) {
	int nNewLatentFactors = 0;
	if (model.mnNLatentFactors < model.mnEvolMaximumFactors) {
		std::cout << "Trying to include a new factor..." << endl; 
		int K = model.mnNLatentFactors;
		model.mnNLatentFactors = K + 1;
		Initialise(model,1);
		Run(0, model);
		GetFounder(K,model.mnNDesignControls);

		int nCount = 0;
		int Index = K + model.mnNDesignControls + mnYFactors;
		for (int i = mnYFactors;  i < mnVariables; i++) {
			if (mMeanResult.mPostPib[i][Index] > model.mdEvolFactorThreshold) {
				nCount++;
			}
		}

		if (nCount <= model.mnEvolMiniumVariablesInFactor) {
			model.mnNLatentFactors = K;
			std::cout << "No factor can be added this time" << endl;
			nNewLatentFactors = 0;
		} else {
			std::cout << endl << "Factor " << (K+1) << " has been added" << endl << endl;
			nNewLatentFactors = 1;
		}
	}
	return nNewLatentFactors;
}

void BfrmEvo::Evol(Model& model) {
	int nNewVariables = 0;
	int nNewLatentFactors = 0;
	do {

		nNewVariables = AddVariables(model); //selecting new variables
		bool bTriedNewFactors = (model.mnNLatentFactors < model.mnEvolMaximumFactors);
		nNewLatentFactors = AddFactors(model);

		if (nNewVariables || nNewLatentFactors || bTriedNewFactors) {
			std::cout << "Refitting the model..." << endl;
			Initialise(model, 0);
			Run(0, model);
			std::cout << "Refitting the model......done" << endl;
		} else {
			std::cout << "No more variables or LatentFactors can be added" << endl;
		}

		
		if (mnVariablesInCount >= model.mnEvolMaximumVariables) {
			std::cout << "Maximum number of variables allowed in the model is reached" << endl;
			if (model.mnNLatentFactors < model.mnEvolMaximumFactors) {
				std::cout << "Trying to add latent factors only from now on." << endl;
			}
			nNewVariables = 0;
		} 

		if (model.mnNLatentFactors >= model.mnEvolMaximumFactors) {
			std::cout << "Maximum number of latent factors allowed in the model is reached" << endl;
			if (mnVariablesInCount < model.mnEvolMaximumVariables) {
				std::cout << "Trying to add variables only from now on." << endl;
			}
			nNewLatentFactors = 0;
		} 
		SaveVariablesInIndex("mVariablesIn.txt");
	} while (nNewVariables || nNewLatentFactors);
//	BfrmData::SaveData("aaa.txt",maVariablesOutIndicator);
	CalculateExternalProbs(model);
}

bool BfrmEvo::SaveVariablesInIndex(string FileName) {
	ofstream theFile(FileName.c_str());
	if  (theFile.fail()) {
		std::cout << "Failed to create file" << FileName.c_str() << endl;
		exit(1);
	}
	for (int i = 0; i < mnVariablesInCount; i++) {
		theFile << (maVariablesInIndex[i] + 1) << "\t";
	}
	theFile << endl;
	theFile.close();
	return true;
}


int main(int argc,char* argv[]) {
	
	Model model;
		
	string pfile = "parameters.txt";
	if (argc > 1) {
		pfile = string(argv[1]);
		if (pfile == "-DEFAULT" || pfile == "-Default" || pfile == "-default") {
			model.Save("default.parameters.txt");
			std::cout << "Default Parameters are saved in 'default.parameters.txt'" << endl;
			exit(0);
		}
	}
	if (!model.Load(pfile)) {
		std::cout << "Loading parameters file failed" << endl;
		exit (1);
	}
	BfrmEvo T(model.mdRandomSeed);
	//read data
	T.LoadData(model);
	
	if (model.mnEvol) {
		T.LoadDataEvol(model);
		T.GetInitialModel(model);
		T.RemoveDesignEffect(model);
		T.Evol(model);
		T.Save(model.mnPriorPiStandard);
		T.SaveVariablesInIndex("mVariablesIn.txt");
	} else {
		//Initialise
		T.Initialise(model, 0);
		std::cout<<"Done init" << endl;

		int nDoneIter = 0;
		T.Run(nDoneIter, model);
		T.Save(model.mnPriorPiStandard);
	}
	std::cout << "done" << endl;
	return 0;
}
