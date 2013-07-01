// $Id: bblEMSeperate.h 962 2006-11-07 15:13:34Z privmane $
#ifndef ___BBL_EM_SEPERATE_H
#define ___BBL_EM_SEPERATE_H

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "sequenceContainer.h"

#include <vector>
using namespace std;


class bblEMSeperate {
public:
	explicit bblEMSeperate(vector<tree>& et,
									const vector<sequenceContainer>& sc,
									const vector<stochasticProcess> &sp,
									const vector<Vdouble *> * weights,
									const int maxIterations=50,
									const MDOUBLE epsilon=0.05,
									const MDOUBLE tollForPairwiseDist=0.0001);
	MDOUBLE getTreeLikelihood() const {return _treeLikelihood;}

private:
	MDOUBLE _treeLikelihood;
	
};

#endif
