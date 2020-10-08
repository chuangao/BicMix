// mdp.h: interface for the mdp class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_MDP_H__589CB87D_E356_47F0_800B_958B01547E81__INCLUDED_)
#define AFX_MDP_H__589CB87D_E356_47F0_800B_958B01547E81__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class mdp  
{
public:
	void init(int nobs, double alpha, double x, double priora, double priorb); 
	mdp();
	virtual ~mdp();

	int AddGroup(int obs);
	void SwitchObservationToGroup(int obs, int group);
	void UpdateAlpha();
	int mnObservations;
	int mnGroups;
	set<int>* mpGroups;
	vector<int> mGroupIndicators;
	double mAlpha;
	double mX;

	double mdPriora;
	double mdPriorb;

	Uniform mRandom;
};

#endif // !defined(AFX_MDP_H__589CB87D_E356_47F0_800B_958B01547E81__INCLUDED_)
