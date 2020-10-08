// util.cpp: implementation of the util class.
//
//////////////////////////////////////////////////////////////////////
#include "stdafx.h"
#include "util.h"
#include <math.h>
#pragma warning(disable : 4996)

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

util::util()
{

}

util::~util()
{

}

string util::Trim(string s)
{
	if (s.empty()) return s;
	string ret;
	for (int i = 0; i < (int)s.length();  i++) {
		if (!isspace(s.at(i)) && !iscntrl(s.at(i)))
			ret.append(s.substr(i,1));
	}
	return ret;
}

string util::ToLower(string str)
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

int util::ToInt(string str)
{
	return atoi(str.c_str());
}

double util::ToDouble(string str)
{
	return atof(str.c_str());
}

string util::ToString(int value) {
	char  buffer[10]; 
	sprintf(buffer, "%d", value);
	string result(buffer);
	return result;
	cout << result.c_str() << endl;
	
}

string util::ToString(double value) {
	char  buffer[10]; 
	sprintf(buffer, "%4.1f", value);
	string result(buffer);
	return result;
	
}

//Multinomial sampling from cells 0,1,....,k-1 with probabilities p0,p1,....,p(k-1)
//Result is (k-1) -vector of counts of #draws in cells 0,1,2,...,k-1
void util::RandSamples(int nMC, vector<double> P, int nSize, vector<int>& result)
{
	result.clear();
	//int nSize = P.size();
	if ((nMC > 0) && (nSize > 0)) {
		double dTotal = 0;
		int i;

		for (i =0 ; i < nSize; i++) {
			dTotal += P[i];
		}
		for (i = 0; i < nSize; i++) {
			P[i] = P[i] /dTotal;
		}
		vector<double> Q;
		cumsum(P,Q);
		double *pRand = new double[nMC];
		cmRand(nMC, pRand);
		vector<int> temp = vector<int>(nMC, 0);
		
		
		for (i = 0; i < nSize; i++) {
			for (int j = 0; j < nMC; j++) {
				if (pRand[j] >= Q[i]) temp[j]++;		
			}
		}

		result = vector<int>(nSize, 0);
		for (i = 0; i < nMC; i++) {
			result[temp[i]]++;
		}
		delete [] pRand;
	}
}

//change to integer result later ????
void util::cmRandPerm(const int nSize, double* pResult ) {
	vector<int> index;
	cmRand(nSize,pResult);
	sort(nSize, pResult, index, 1);
	for (int i = 0; i < nSize; i++) {
		pResult[i] = index[i];
	}
}

void util::SelectCVGroup(int Total, int GroupSize, vector<int>& subject, vector<int>& current, vector<int>& train, int GroupID)
{
	int nSize = (int)subject.size();
	if (GroupSize == 0) {		//no cross validation
		current.clear();
		train = subject;
		subject.clear();
		return;
	}
	
	int i;

	double *pRandPerm = new double[nSize];
	if (GroupID == 0) {
		cmRandPerm(nSize, pRandPerm);
	} else {
		for (int i = 0; i < nSize; i++) {
			pRandPerm[i] = i + 1;
		}
	}
	
	list<int> TempCurrent;
	list<int> TempSubject;
	int nCVSize = (nSize > GroupSize) ? GroupSize : nSize;
	
	for (i = 0; i < nSize; i++) {
		if (i < nCVSize) {
			TempCurrent.push_back(subject[(int)pRandPerm[i] - 1]);
		} else {
			TempSubject.push_back(subject[(int)pRandPerm[i] - 1]);
		}
	}

	subject.clear();
	list<int>::iterator it;
	if (!TempSubject.empty()) {
		TempSubject.sort();
		for (it = TempSubject.begin(); it != TempSubject.end(); it++) {
			subject.push_back(*it);
		}
	}
	
	train.clear();
	current.clear();
	if (!TempCurrent.empty()) {
		i = 1;
		TempCurrent.sort();
		while (!TempCurrent.empty()) {
			int nCurrent = TempCurrent.front();
			TempCurrent.pop_front();
			current.push_back(nCurrent);
			while (i < nCurrent) {
				train.push_back(i++);
			}
			i++;
		}
		while (i <= Total) {
			train.push_back(i++);
		}
	}
	delete [] pRandPerm;
}

vector<int>& util::colon(vector<int>& result, int start, int end, int step /*= 1 */) {
	result.clear();
	for (int i = start; i <= end; i+= step) {
		result.push_back(i);
	}	return result;
}

/*
vector<int>& util::sort2(int size, double* pData, vector<int>& index, int nDirection) {
	index.clear();
	double *pIndex = new double[size];
	cmSort(size, pData, pIndex, nDirection);
	for (int i = 0; i <size; i++) {
		index.push_back((int)pIndex[i]);
	}
	return index;
	delete [] pIndex;
}
*/

vector<double>& util::cumsum(vector<double>& source, vector<double>& result)
{
	result.clear();
	if (!source.empty()) {
		int nSize = (int)source.size();
		result.push_back(source[0]);
		for (int i = 1; i < nSize; i++) {
			result.push_back(result[i - 1] + source[i]);
		}
	}
	return result;
}

//nDirection = 1 :: ascending
vector<int>& util::sort(int size, double* pData, vector<int>& index, int nDirection/* = 1*/) {
	index.clear();
	int *pfrom = new int[size];
	int *pto   = new int[size];
	double *pDataCopy = new double[size];
	int i;
	for (i = 0; i < size; i++) {
		pfrom[i] = i;
	}
	memcpy(pDataCopy, pData, sizeof(double) * size);
	memcpy(pto, pfrom, sizeof(int) * size);
	int ascending = 1;		
	if (nDirection != 1) ascending = 0;  
	shuttlesort(pData, pfrom, pto, 0, size, ascending);
	
	for (i = 0; i <size; i++) {
		index.push_back(pto[i] + 1);   //convert to one based as Matlab does
		pData[i] = pDataCopy[pto[i]];
	}
	
	delete [] pfrom;
	delete [] pto;
	delete [] pDataCopy;
	return index;
}

void util::shuttlesort(double* pData, int* from, int* to, int low, int high, int asceding) {
	if (high - low < 2) return;
    int middle = (low + high ) / 2;
    shuttlesort(pData, to, from, low, middle, asceding);
    shuttlesort(pData, to, from, middle, high, asceding);
    int p = low;
    int q = middle;

    if (high - low >= 4 && compare(pData[from[middle-1]], pData[from[middle]], asceding) <= 0) {
        for (int i = low; i < high; i++) {
            to[i] = from[i];
        }
        return;
    }
    // A normal merge. 
    for (int i = low; i < high; i++) {
        if (q >= high || (p < middle && compare(pData[from[p]], pData[from[q]], asceding) <= 0)) {
            to[i] = from[p++];
        }
        else {
            to[i] = from[q++];
        }
    }
}

int util::compare(double d1, double d2, int ascending) {
	int result;
    if (d1 < d2) {
        result = -1;
    } else if (d1 > d2) {
        result = 1;
    } else {
        result = 0;
    }
    if (result != 0) {
        return ascending ? result : -result;
    }
	return 0;
}

void util::cmPower2(int nSize, double *px, double* py, double* pResult) {
	for (int i = 0; i < nSize; i++) {
		pResult[i] = pow(px[i], py[i]);
	}
}
