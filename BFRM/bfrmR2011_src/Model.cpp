// Model.cpp: implementation of the Model class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#include "Model.h"
#pragma warning(disable:4786)
#pragma warning(disable : 4996)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Model::Model()
{
	mnNObservations = 0;		//total number of observations within the data
	mnNVariables = 0;			//total number of variables within the model (include Z and X)
	mnNLatentFactors = 0;			//number of latent factors
	mnNDesignControls = 1;			//for the constant or design(control) factors
	mnNControlFactors = 0;

	//Z (Y) variables, by default, these terms are not included
	mnNVarBinary = 0;		//number of binary variables in the model
	mnNVarCategorical = 0;	//number of categorical variabes in the model
	mnNVarSurvival = 0;		//number of survival variables in the model
	mnNVarContinuous = 0;	//number of continuous variables in the model


	mstrDataFile = "";			//X (gene expression matrix)
	mstrHFile = "";				//H matrix (constant, design, control etc.)
	mstrWeightFile = "";		//Sample inclusion indicators
	mstrResponseMaskFile = "";		//Z (Y) mask variables, used to indicate missing values, censorships etc
	mstrXMaskFile = "";
	mstrEvolVarInFile = "";		//variables in in evolving model
	
	mnBurnin = 2000;
	mnMC = 5000;
	mnHistory = 0;

	mnBurninSelect = 1000;
	mnMCSelect = 2000;
	
	mnPrintIteration = 100;
	mnShapeOfB = 2;			//0 = no mask
	mnEvol = 0;
	mnEvolMiniumVariablesInFactor = 5;
	mnEvolMaximumVariables = 1000;
	mnEvolMaximumVariablesPerIteration = 25;
	mnEvolMaximumVariablesPerFactor = -1;
	mnEvolMaximumFactors = 50;

	mdMeanValue = 8.0;			//grand mean for X matrix
	mdMeanSurvival = 4.0;
	mdMeanContinuous = 0.0;
	mdVarValue = 100.0;
	mdVarSurvival = 1.0;
	mdVarContinuous = 1.0;
	mdEvolFactorThreshold = 0.7;

	mdEvolVariableThreshold = 0.95;
	mdaPriorPsi[0] = 10.0;			
	mdaPriorPsi[1] = 2.0;

	mdaPriorSurvivalPsi[0] = 2.0;
	mdaPriorSurvivalPsi[1] = 0.5;
	mdPriorRhoMean = 0.001;
	mdPriorRhoN = 200;

	mnPriorPiStandard = 1;		
	mdPriorPiMean = 0.9; 
	mdPriorPiN = 10.0;


	mdaPriorTauDesign[0] = 5.0;
	mdaPriorTauDesign[1] = 1.0;
	for (int i = 0; i < 4; i++) {
		mdaPriorTauResponse[i][0] = 5.0;
		mdaPriorTauResponse[i][1] = 1.0;
	}
	mdaPriorTauLatent[0] = 5.0;
	mdaPriorTauLatent[1] = 1.0;

	mnEvolVarIn = 0;
	mdDPAlpha  = 1.0;
	mdDPAlphaX = 0.5;
	mdaPriorAlpha[0] = 1.0;
	mdaPriorAlpha[1] = 1.0;

	mnNonGaussianFactors = 1;
        mdRandomSeed = -1.0;
}

string Model::ToLower(string str)
{
	//new implementation for GNU
	char *newstr = strdup(str.c_str());
	int i = 0;
	while (newstr[i] != '\0') {
		newstr[i] = tolower(newstr[i]);
		i++;
	}
	return newstr;
	//return strlwr(strdup(str.c_str())); 
}

int Model::ToInt(string str)
{
	return atoi(str.c_str());
}

double Model::ToDouble(string str)
{
	return atof(str.c_str());
}

string Model::ToString(int value) {
	char  buffer[10]; 
	sprintf(buffer, "%d", value);
	string result(buffer);
	return result;
	cout << result.c_str() << endl;
	
}

string Model::ToString(double value) {
	char  buffer[10]; 
	sprintf(buffer, "%4.1f", value);
	string result(buffer);
	return result;
	
}

string Model::trim(string s)
{
	if (s.empty()) return s;
	string ret;
	for (int i = 0; i < (int)s.length();  i++) {
		if (!isspace(s.at(i)) && !iscntrl(s.at(i)))
			ret.append(s.substr(i,1));
	}
	return ret;
}


bool Model::Load(string FileName){
	int BufferSize = 4096;
	char* theLine = new char[BufferSize];
	ifstream theFile(FileName.c_str());
	if (theFile.fail()) {
		std::cout <<  "Failed to open the model file!";
		return false;
	}

	int nLineCount = 0;
	while (!theFile.eof()) {
		theFile.getline(theLine, BufferSize);
		nLineCount++;
		string theline(theLine);
		string Name(""), Value("");
		theline = trim(theline);			//so space is not allowed, should be improved later
		if (theline.length() && (theline.c_str()[0] != '#'))
		{
			int pos = 0;
			if ((pos = (int)theline.find("=")) != -1) {
				Name = theline.substr(0, pos);
				Value = theline.substr(pos + 1);
			}
			if (Name == "" && Value == "") {
			} else if (Name == "" || Value == "") {
				cout << "Invalid Parameter!" << endl;
				cout << theLine << endl;
				return false;
			} else {
				string name = ToLower(Name);
				string value = ToLower(Value);
				if (name == "datafile") {
					mstrDataFile = Value;
				} else if (name == "hfile") {
					mstrHFile = Value;
				} else if (name == "weightfile") {
					mstrWeightFile = Value;
				} else if (name == "responsemaskfile") {
					mstrResponseMaskFile = Value;
				} else if (name == "xmaskfile") {
					mstrXMaskFile = Value;
				} else if (name == "evolvarinfile") {
					mstrEvolVarInFile = Value;
				} else if (name == "nobservations") {
					mnNObservations = ToInt(value);
				} else if (name == "shapeofb") {
					mnShapeOfB = ToInt(value);
				} else if (name == "nvariables") {
					mnNVariables = ToInt(value);
				} else if (name == "nongaussianfactors") {
					mnNonGaussianFactors = ToInt(value);
				} else if (name == "nbinaryresponses") {
					mnNVarBinary = ToInt(value);
				} else if (name == "ncategoricalresponses") {
					mnNVarCategorical = ToInt(value);
				} else if (name == "nsurvivalresponses") {
					mnNVarSurvival = ToInt(value);
				} else if (name == "ncontinuousresponses") {
					mnNVarContinuous = ToInt(value);
				} else if (name == "evolvarin") {
					mnEvolVarIn = ToInt(value);
				} else if ((name == "evolminiumvariablesinfactor") || (name == "evolminimumvariablesinfactor"))  {
					mnEvolMiniumVariablesInFactor = ToInt(value);
				} else if (name == "evolmaximumvariables") {
					mnEvolMaximumVariables = ToInt(value);
				} else if (name == "evolmaximumvariablesperiteration") {
					mnEvolMaximumVariablesPerIteration = ToInt(value);
				} else if (name == "evolmaximumvariablesperfactor") {
					mnEvolMaximumVariablesPerFactor = ToInt(value);
				} else if (name == "evolmaximumfactors") {
					mnEvolMaximumFactors = ToInt(value);
				} else if (name == "evol") {
					mnEvol = ToInt(value);
				} else if (name == "printiteration") {
					mnPrintIteration = ToInt(value);
				} else if (name == "nlatentfactors") {
					mnNLatentFactors = ToInt(value);
				} else if (name == "ndesignvariables") {
					mnNDesignControls = ToInt(value);
				} else if (name == "ncontrolvariables") {
					mnNControlFactors = ToInt(value);
				} else if (name == "priorpistandard") {
					mnPriorPiStandard = ToInt(value);
				} else if (name == "burnin") {
					mnBurnin = ToInt(value);
				} else if (name == "burnin_select") {
					mnBurninSelect = ToInt(value);
				} else if (name == "history") {
					mnHistory = ToInt(value);
				} else if (name == "priorrhon") {
					mdPriorRhoN = ToDouble(value);
				} else if (name == "alpha") {
					mdDPAlpha = ToDouble(value);
				} else if (name == "alphax") {
					mdDPAlphaX = ToDouble(value);
				} else if (name == "priorrhomean") {
					mdPriorRhoMean = ToDouble(value);
				} else if (name == "priorpin") {
					mdPriorPiN = ToDouble(value);
				} else if (name == "priorpimean") {
					mdPriorPiMean = ToDouble(value);
				} else if (name == "priorpsia") {
					mdaPriorPsi[0] = ToDouble(value);
				} else if (name == "priorpsib") {
					mdaPriorPsi[1] = ToDouble(value);
				} else if (name == "priorsurvivalpsia") {
					mdaPriorSurvivalPsi[0] = ToDouble(value);
				} else if (name == "priorsurvivalpsib") {
					mdaPriorSurvivalPsi[1] = ToDouble(value);
				} else if (name == "priortaudesigna") {
					mdaPriorTauDesign[0] = ToDouble(value);
				} else if (name == "priortauresponsebinarya") {
					mdaPriorTauResponse[0][0] = ToDouble(value);
				} else if (name == "priortauresponsecategoricala") {
					mdaPriorTauResponse[1][0] = ToDouble(value);
				} else if (name == "priortauresponsesurvivala") {
					mdaPriorTauResponse[2][0] = ToDouble(value);
				} else if (name == "priortauresponsecontinuousa") {
					mdaPriorTauResponse[3][0] = ToDouble(value);
				} else if (name == "priortaulatenta") {
					mdaPriorTauLatent[0] = ToDouble(value);
				} else if (name == "prioralphaa") {
					mdaPriorAlpha[0] = ToDouble(value);
				} else if (name == "prioralphab") {
					mdaPriorAlpha[1] = ToDouble(value);
				} else if (name == "priortaudesignb") {
					mdaPriorTauDesign[1] = ToDouble(value);
				} else if (name == "priortauresponsebinaryb") {
					mdaPriorTauResponse[0][1] = ToDouble(value);
				} else if (name == "priortauresponsecategoricalb") {
					mdaPriorTauResponse[1][1] = ToDouble(value);
				} else if (name == "priortauresponsesurvivalb") {
					mdaPriorTauResponse[2][1] = ToDouble(value);
				} else if (name == "priortauresponsecontinuousb") {
					mdaPriorTauResponse[3][1] = ToDouble(value);
				} else if (name == "priortaulatentb") {
					mdaPriorTauLatent[1] = ToDouble(value);
				} else if (name == "evolincludevariablethreshold") {
					mdEvolVariableThreshold = ToDouble(value);
					if (mdEvolVariableThreshold <= 0) {
						mdEvolVariableThreshold  = 0.00000001; 
					}
				} else if (name == "evolincludefactorthreshold") {
					mdEvolFactorThreshold = ToDouble(value); 
					if (mdEvolFactorThreshold <=0) {
						mdEvolFactorThreshold = 0.00000001;
					}
				} else if (name == "priorinterceptmean") {
					mdMeanValue = ToDouble(value);
				} else if (name == "priorinterceptvar") {
					mdVarValue = ToDouble(value);
				} else if (name == "priorsurvivalmean") {
					mdMeanSurvival = ToDouble(value);
				} else if (name == "priorcontinuousmean") {
					mdMeanContinuous = ToDouble(value);
				} else if (name == "priorsurvivalvar") {
					mdVarSurvival = ToDouble(value);
				} else if (name == "priorcontinuousvar") {
					mdVarContinuous = ToDouble(value);
				} else if (name == "nmcsamples") {
					mnMC = ToInt(value);
				} else if (name == "nmcsamples_select") {
					mnMCSelect = ToInt(value);
				} else if (name == "randomseed") {
					mdRandomSeed = ToDouble(value);
				} else {
					cout << "Unknown Parameter" << endl; //to be refined later
					cout << theLine << endl;
					return false;
				}
			}
		}
	}
	mnNDesignControls += mnNControlFactors;
	delete[] theLine;
	if (mstrDataFile == "") {
		std::cout << endl << "Warning: Text file "  << FileName.c_str() << " might come from a different platform" << endl;
		std::cout << "Suggestion: Use dos2unix/unis2dos or similar tools to convert the file first" << endl;
		return false;
	}
	return true;
}

bool Model::Save(string FileName){
	ofstream theFile(FileName.c_str());
	if (theFile.fail()) {
		cout << "Failed to create file!" << endl;
		return false;
	}
	theFile << "#Version 2.0" << endl << endl;

	theFile << "#data section" << endl;
	theFile << "NObservations = " << mnNObservations << endl;
	theFile << "NVariables = " << mnNVariables << endl;
	theFile << "NBinaryResponses = " << mnNVarBinary << endl;
	theFile << "NCategoricalResponses = " << mnNVarCategorical << endl;
	theFile << "NSurvivalResponses = " << mnNVarSurvival << endl;
	theFile << "NContinuousResponses = " << mnNVarContinuous << endl;
	theFile << "NDesignVariables = " << mnNDesignControls - mnNControlFactors<< endl;
	theFile << "NControlVariables = " << mnNControlFactors << endl;
	theFile << "NLatentFactors = " << mnNLatentFactors << endl;
	if (mstrDataFile != "") {
		theFile << "DataFile = " << mstrDataFile.c_str() << endl;
	} else {
		theFile << "#DataFile = " << endl;
	}
	if (mstrHFile != "") {
		theFile << "HFile = " << mstrHFile.c_str() << endl;
	} else {
		theFile << "#HFile = " << endl;
	}
	if (mstrResponseMaskFile != "") {
		theFile << "ResponseMaskFile = " << mstrResponseMaskFile.c_str() << endl;
	} else {
		theFile << "#ResponseMaskFile = " << endl;
	}
	if (mstrXMaskFile != "") {
		theFile << "XMaskFile = " << mstrXMaskFile.c_str() << endl;
	} else {
		theFile << "#XMaskFile = " << endl;
	}

	theFile << endl;

	theFile << "#prior section" << endl;

	theFile << "#model specification" << endl;
	theFile << "ShapeOfB = " << mnShapeOfB << endl << endl;
	theFile << "NonGaussianFactors = " << mnNonGaussianFactors << endl;
	theFile << "PriorPiStandard = " << mnPriorPiStandard << endl;
	

	theFile << "#prior Psi" << endl;
	theFile << "PriorPsia = " << mdaPriorPsi[0] << endl;
	theFile << "PriorPsib = " << mdaPriorPsi[1] << endl;
	theFile << "PriorSurvivalPsia = " << mdaPriorSurvivalPsi[0] << endl;
	theFile << "PriorSurvivalPsib = " << mdaPriorSurvivalPsi[1] << endl;
	theFile << endl;

	theFile << "#prior Rho" << endl;
	theFile << "PriorRhoMean = " << mdPriorRhoMean << endl;
	theFile << "PriorRhoN = " << mdPriorRhoN << endl;
	theFile << endl;
	
	theFile << "#prior Pi" << endl;
	theFile << "PriorPiMean = " << mdPriorPiMean << endl;
	theFile << "PriorPiN = " << mdPriorPiN << endl << endl;
	
	theFile << "#prior Tau" << endl;
	theFile << "PriorTauDesigna = " << mdaPriorTauDesign[0] << endl;
	theFile << "PriorTauDesignb = " << mdaPriorTauDesign[1] << endl;
	theFile << endl;
	
	theFile << "PriorTauResponseBinarya = " << mdaPriorTauResponse[0][0] << endl;
	theFile << "PriorTauResponseBinaryb = " << mdaPriorTauResponse[0][1] << endl;
	theFile << endl;
	
	theFile << "PriorTauResponseCategoricala = " << mdaPriorTauResponse[1][0] << endl;
	theFile << "PriorTauResponseCategoricalb = " << mdaPriorTauResponse[1][1] << endl;
	theFile << endl;

	theFile << "PriorTauResponseSurvivala = " << mdaPriorTauResponse[2][0] << endl;
	theFile << "PriorTauResponseSurvivalb = " << mdaPriorTauResponse[2][1] << endl;
	theFile << endl;

	theFile << "PriorTauResponseContinuousa = " << mdaPriorTauResponse[3][0] << endl;
	theFile << "PriorTauResponseContinuousb = " << mdaPriorTauResponse[3][1] << endl;
	theFile << endl;

	theFile << "PriorTauLatenta = " << mdaPriorTauLatent[0] << endl;
	theFile << "PriorTauLatentb = " << mdaPriorTauLatent[1] << endl;
	theFile << endl;

	theFile << "#priors on Intercept" << endl;
	theFile << "PriorInterceptMean = " << mdMeanValue << endl;
	theFile << "PriorSurvivalMean = " << mdMeanSurvival << endl;
	theFile << "PriorContinuousMean = " << mdMeanContinuous << endl;
	theFile << "PriorInterceptVar = " << mdVarValue << endl;
	theFile << "PriorSurvivalVar = " << mdVarSurvival << endl;
	theFile << "PriorContinuousVar = " << mdVarContinuous << endl << endl;

	theFile << "#evolving mode section" << endl;
	theFile << "Evol = " << mnEvol << endl;
	theFile << "EvolVarIn = " << mnEvolVarIn << endl;
	if (mstrEvolVarInFile != "") {
		theFile << "EvolVarInFile = " << mstrEvolVarInFile.c_str() << endl;
	} else {
		theFile << "#EvolVarInFile = " << endl;
	}
	theFile << "EvolIncludeVariableThreshold = " << mdEvolVariableThreshold << endl;
	theFile << "EvolIncludeFactorThreshold = " << mdEvolFactorThreshold << endl;
	theFile << "EvolMinimumVariablesInFactor = " << mnEvolMiniumVariablesInFactor << endl;
	theFile << "EvolMaximumFactors = " << mnEvolMaximumFactors << endl;
	theFile << "EvolMaximumVariables = " << mnEvolMaximumVariables << endl;
	theFile << "EvolMaximumVariablesPerIteration = " << mnEvolMaximumVariablesPerIteration << endl;
	theFile << endl;

	theFile << "#mcmc section" << endl;
	theFile << "Burnin = " << mnBurnin << endl;
	theFile << "Burnin_Select = " << mnBurninSelect << endl;
	theFile << "nMCSamples = " << mnMC << endl;
	theFile << "nMCSamples_Select = " << mnMCSelect << endl;	
	theFile << "History = " << mnHistory << endl;	
	theFile << endl;

	theFile << "#monitoring section" << endl;
	theFile << "PrintIteration = " << mnPrintIteration << endl;
	theFile << endl;

	theFile << "#DP parameters" << endl;
	theFile << "alpha = " << mdDPAlpha << endl;
	theFile << "alphaX = " << mdDPAlphaX << endl;
	theFile << "PriorAlphaa = " << mdaPriorAlpha[0] << endl;
	theFile << "PriorAlphab = " << mdaPriorAlpha[1] << endl;

	theFile.flush();
	theFile.close();
	return true;
}

