// util.h: interface for the util class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_UTIL_H__256A235D_9557_4858_B315_88FB13A9044B__INCLUDED_)
#define AFX_UTIL_H__256A235D_9557_4858_B315_88FB13A9044B__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class util  
{
public:
	util();
	virtual ~util();
	
	static double ToDouble(string str);
	static int ToInt(string str);
	static string ToLower(string str);
	static string Trim(string s);
	static string ToString(int value);
	static string ToString(double value);

	//c++ base power function
	static void cmPower2(int nSize, double *px, double* py, double* pResult);
	static void SelectCVGroup(int Total, int GroupSize, vector<int>& subject, vector<int>& current, vector<int>& train, int GroupID);
	static void cmRandPerm(const int nSize, double* pResult);
	static vector<int>& colon(vector<int>& result, int start, int end, int step = 1); 

	//based on home made shuttle sort
	static vector<int>& sort(int size, double* pData, vector<int>& index, int nDirection = 1); 
	static vector<double>& cumsum(vector<double>& source, vector<double>& result);
	static void RandSamples(int nMC, vector<double> P, int nSize, vector<int>& result);
	static bool GammaRand(double a, double b, int nSize, vector<double>& result);
	static bool BetaRand(double a, double b, int nSize, vector<double>& result);
	static vector<double> Percentile(vector<double>& X, int nCount, double * pP); 
	static double Gammaln(double x);
	static double Betaln(double x, double y);
	static void cmRand(int nSize, double* pResult);
	static void cmInterp2(int nSize, double *px, double *py, int nCount, double *pP, double *pResult);

	static double norminv(double p);			//inverse normal cdf
	static double normcdf(double u);
	//static double cmMean(vector<double>& data);
private:
	static void shuttlesort(double* pData, int* from, int* to, int low, int high, int asceding);
	static int compare(double d1, double d2, int ascending); 
	
};

#endif // !defined(AFX_UTIL_H__256A235D_9557_4858_B315_88FB13A9044B__INCLUDED_)
