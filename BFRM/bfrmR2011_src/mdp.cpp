// mdp.cpp: implementation of the mdp class.
//Author: Quanli Wang quanli@stat.duke.edu
//Department of Statistical Science, Duke University
//////////////////////////////////////////////////////////////////////
#include "newran.h"

#include "stdafx.h"
#include "mdp.h"
#include "math.h"
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

void mdp::init(int nobs, double alpha, double x, double priora, double priorb)
{
	mnObservations = nobs;
	mAlpha = alpha;
	mX = x;

	mdPriora = priora;
	mdPriorb = priorb;
	if (mpGroups) {
		delete [] mpGroups; mpGroups = NULL;
	}

	//orginally each observation belongs to its own group
	mpGroups = new set<int>[nobs];
	mnGroups = nobs;
	mGroupIndicators.clear();
	for (int i = 0; i < mnGroups; i++) {
		mpGroups[i].insert(i);
		mGroupIndicators.push_back(i);
	}
}

mdp::mdp()
{
	srand( (unsigned)time( NULL ) );
	Random::Set(((double)rand()) / ((double)RAND_MAX));
	mpGroups = NULL;
}

int mdp::AddGroup(int obs) {
	int from = mGroupIndicators[obs];
	if (mpGroups[from].size() == 1) { return from; }
	mpGroups[from].erase(obs);
	
	mpGroups[mnGroups].clear();
	mpGroups[mnGroups].insert(obs);
	mGroupIndicators[obs] = mnGroups;
	mnGroups++;
	return mnGroups - 1;
}

void mdp::SwitchObservationToGroup(int obs, int group) {
	int from = mGroupIndicators[obs];
	//if already there, do nothing
	if (from == group) { return;}

	//first assign the new group ownership
	mpGroups[group].insert(obs);
	mGroupIndicators[obs] = group;


	//now remove from the old group
	mpGroups[from].erase(obs);
	//if the old group is empty them get ride of it
	if (mpGroups[from].empty()) {
		//move group only if the 'from' group is not the last group
		if (mnGroups - 1 > from) {
			//swap the last group into this empty group
			mpGroups[from].swap(mpGroups[mnGroups-1]);
			//then change the indicators for that last group to the empty group
			for (set<int>::iterator it = mpGroups[from].begin(); it != mpGroups[from].end(); it++ ) {
				mGroupIndicators[*it] = from;
			}	
		}
		//and remove the empty group(now the last one)
		mnGroups--;
	}
}

void mdp::UpdateAlpha() {
	//update x
	Gamma ga(mAlpha + 1);
	Gamma gb(mnObservations);
	double ra = ga.Next();
	double rb = gb.Next();
	mX =  (ra / (ra  + rb));

	//calculate Pi_x
	double logx = log(mX);
	double ratio = (mdPriora + mnGroups - 1) / (mnObservations * (mdPriorb - logx));
	double pi_x = ratio / ( 1 + ratio);

	double a = mdPriora + mnGroups;
	double prob = mRandom.Next();
	if (prob >pi_x) {
		a = a - 1;
	} 
	Gamma g(a);
	mAlpha = g.Next() / (mdPriorb - logx); //gamma(a,b) = gamma(a,1) / b; ***************
}


mdp::~mdp()
{
	if (mpGroups) {
		delete [] mpGroups;
		mpGroups = NULL;
	}
}
