// bfrmdata.h: interface for the BfrmData class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_BFRMDATA_H)
#define AFX_BFRMDATA_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class BfrmData  
{
public:
	BfrmData();
	virtual ~BfrmData();
	static bool LoadData(string FileName, Array2D<double>& A, int rows, int columns);
	static bool LoadData(string FileName, Array1D<double>& A, int size);
	static bool LoadData(string FileName, Array1D<int>& A, int size, int capacity);
	static bool SaveData(string FileName, Array2D<double>& A);
	static bool SaveData(string FileName, Array2D<double>& A, int append);
	static bool SaveData(string FileName, Array1D<double>& A, int append);
	static bool SaveData(string FileName, Array1D<double>& A);
	static Array1D<double> Var(Array2D<double>& A);
};

#endif // !defined(AFX_BFRMDATA_H)
