#ifndef ___BB_COMPUTE_UP_ALG__
#define ___BB_COMPUTE_UP_ALG__

#include "computePijComponent.h"
#include "suffStatComponent.h"

// the only different from computeUpAlg is that here char assignments to 
// internal nodes are taken into account while calculating compute up.

#include "tree.h"
#include "sequenceContainer.h"
#include "computePijComponent.h" 
#include "suffStatComponent.h"
#include "sequence.h"
#include <vector>
using namespace std;

void BBfillComputeUp(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const computePijHom& pi,
				   suffStatGlobalHomPos& ssc,
				   const vector<sequence>& ancS);

#endif

