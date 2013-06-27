#ifndef ___BB_COMPUTE_DOWN_ALG__
#define ___BB_COMPUTE_DOWN_ALG__

#include "tree.h"
#include "sequenceContainer.h"
#include "computePijComponent.h" 
#include "suffStatComponent.h"
#include "sequence.h"
#include <vector>
using namespace std;

void BBfillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const vector<sequence>& ancS);



#endif

