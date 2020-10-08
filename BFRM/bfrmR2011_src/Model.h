// Model.h: interface for the Model class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////
#include <map>
#include <cstring>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>

using namespace std;

class Model  
{
public:
	Model();
	bool Load(string FileName);
	bool Save(string FileName);
public:
	string mstrDataFile;
	string mstrHFile;
	string mstrWeightFile;
	string mstrResponseMaskFile;
	string mstrEvolVarInFile;
	string mstrXMaskFile;

	int mnNObservations;
	int mnNLatentFactors;;
	int mnNDesignControls;

	int mnNVarBinary;
	int mnNVarCategorical;
	int mnNVarSurvival;
	int mnNVarContinuous;

	int mnBurnin;
	int mnMC;
	int mnHistory;

	int mnBurninSelect;
	int mnMCSelect;

	int mnPrintIteration;
	int mnShapeOfB;
	int mnEvol;
	int mnEvolMiniumVariablesInFactor;
	double mdEvolFactorThreshold;
	int mnEvolMaximumVariables;
	int mnEvolMaximumVariablesPerIteration;
	int mnEvolMaximumVariablesPerFactor;
	int mnEvolMaximumFactors;
	double mdEvolVariableThreshold;
	int mnEvolVarIn;

	int mnPriorPiStandard;						
	int mnNControlFactors;

	//Rho prior
	double mdPriorRhoN;
	double mdPriorRhoMean;

	//Pi prior
	double mdPriorPiMean;
	double mdPriorPiN;
	
	
	//Psi
	double mdaPriorPsi[2];
	double mdaPriorSurvivalPsi[2];
	
	//Tau
	double mdaPriorTauDesign[2];
	double mdaPriorTauResponse[4][2];
	double mdaPriorTauLatent[2];
	double mdaPriorAlpha[2];

	//new intercept prriors
	double mdMeanValue;
	double mdMeanSurvival;
	double mdMeanContinuous;
	double mdVarValue;
	double mdVarSurvival;
	double mdVarContinuous;

	double mdDPAlpha;
	double mdDPAlphaX;
	int mnNonGaussianFactors;

        double mdRandomSeed;
	
	int GetDataRowCount() { return mnNVariables; }
	static double ToDouble(string str);
	static int ToInt(string str);
	static string ToLower(string str);
	static string ToString(int value);
	static string ToString(double value);
	static string trim(string s);
private:
	int mnNVariables;

};

