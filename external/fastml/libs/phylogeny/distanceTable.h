// $Id: distanceTable.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___DISTANCE_TABLE
#define ___DISTANCE_TABLE

#include "definitions.h"
#include "distanceMethod.h"
#include "sequenceContainer.h"

void giveDistanceTable(const distanceMethod* dis,
					   const sequenceContainer& sc,
					   VVdouble& res,
					   vector<string>& names,
					   const vector<MDOUBLE> * weights = NULL);


#endif
