// bfrmdata.cpp: implementation of the BfrmData class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////
#include "tnt.h"
using namespace TNT;
#include "stdafx.h"
#include "bfrmdata.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BfrmData::BfrmData()
{

}

BfrmData::~BfrmData()
{

}

bool BfrmData::LoadData(string FileName, Array2D<double>& A, int rows, int columns) {
	std::cout << "Start loading file " << FileName.c_str() << " ...";
	ifstream theFile(FileName.c_str());
	if (theFile.fail()) {
		std::cout << "Failed to open the file " << FileName.c_str() << endl;
		exit(1);
	}
	std::cout << " done" << endl;
	A = Array2D<double>(rows,columns);
	double dCurrent = 0;
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < columns; j++) {
			if (!theFile.eof()) {
				theFile >> dCurrent;
				A[i][j] = dCurrent;
			} else {
				std::cout << "Not enough numbers to be read" << endl;
				exit(1);
			}
		}
	}
	theFile.close();
	return true;
}

bool BfrmData::SaveData(string FileName, Array2D<double>& A, int append) {
	int M = A.dim1();
    int N = A.dim2();
	//int mode = (append > 0) ? ios::app : ios::out;
	ofstream theFile(FileName.c_str(), (append > 0) ? ios::app : ios::out);
	if (theFile.fail()) {
		std::cout << "Failed to create file " << FileName.c_str()  << endl;
		exit(1);
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			theFile << A[i][j] << "\t";
		}
		theFile << endl;
	}
	theFile << endl;
	theFile.close();
	return true;
}
bool BfrmData::SaveData(string FileName, Array1D<double>& A, int append) {
	int M = A.dim();
	//int mode = (append > 0) ? ios::app : ios::out;
	ofstream theFile(FileName.c_str(), (append > 0) ? ios::app : ios::out);
	//ofstream theFile(FileName.c_str());
	if  (theFile.fail()) {
		std::cout << "Failed to create file" << FileName.c_str() << endl;
		exit(1);
	}
	for (int i = 0; i < M; i++) {
		theFile << A[i] << "\t";
	}
	theFile << endl;
	theFile.close();
	return true;
}

bool BfrmData::SaveData(string FileName, Array2D<double>& A) {
	int M = A.dim1();
    int N = A.dim2();

	ofstream theFile(FileName.c_str());
	if (theFile.fail()) {
		std::cout << "Failed to create file " << FileName.c_str()  << endl;
		exit(1);
	}
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			theFile << A[i][j] << "\t";
		}
		theFile << endl;
	}
	theFile << endl;
	theFile.close();
	return true;
}

bool BfrmData::SaveData(string FileName, Array1D<double>& A) {
	int M = A.dim();
	ofstream theFile(FileName.c_str());
	if  (theFile.fail()) {
		std::cout << "Failed to create file" << FileName.c_str() << endl;
		exit(1);
	}
	for (int i = 0; i < M; i++) {
		theFile << A[i] << "\t";
	}
	theFile << endl;
	theFile.close();
	return true;
}

bool BfrmData::LoadData(string FileName, Array1D<int>& A, int size, int capacity) {
	A = Array1D<int>(capacity); A = 0;
	if (size == 0) { return true;}
	if (FileName == "") {
		for (int i = 0; i < size; i++) {
			A[i] = i + 1;
		}
		return true;
	}
	ifstream theFile(FileName.c_str());
	if (theFile.fail()) {
		std::cout << "Failed to open the file " << FileName.c_str() << endl;
		exit(1);
	}

	double dCurrent = 0;
	for (int i = 0; i < size; i++) {
		if (!theFile.eof()) {
			theFile >> dCurrent;
			A[i] = (int)dCurrent;
		} else {
			std::cout << "Failed to read the file " << FileName.c_str() << endl;
			exit(1);
		}
	}
	theFile.close();
	return true;
}


bool BfrmData::LoadData(string FileName, Array1D<double>& A, int size) {
	ifstream theFile(FileName.c_str());
	if (theFile.fail()) {
		std::cout << "Failed to open the file" << FileName.c_str() << endl;
		exit(1);
	}
	A = Array1D<double>(size);
	double dCurrent = 0;
	for (int i = 0; i < size; i++) {
		if (!theFile.eof()) {
			theFile >> dCurrent;
			A[i] = dCurrent;
		} else {
			break;
		}
	}
	theFile.close();
	return true;
}

//calculate variance for row
Array1D<double> BfrmData::Var(Array2D<double>& A) {
	
	int M = A.dim1();
	int N = A.dim2();
	
	//subtract means  -- subtractsamplemeans
	Array2D<double> Ones(N,1); Ones = 1.0;
	Array2D<double> Mean = matmult(A,Ones);
	Ones = 1.0 / N;
	Array2D<double> B = A  -  matmult(Mean,transpose(Ones));

	Array1D<double> V(M); V = 0;
	for (int i = 0; i < M; i++) {
		for (int j = 0; j < N; j++) {
			V[i] += B[i][j] * B[i][j]; 	
		}
		V[i] = V[i] / (N -1);
	}
	return V;
}
