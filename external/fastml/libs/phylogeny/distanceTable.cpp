// $Id: distanceTable.cpp 1740 2007-02-26 13:53:10Z itaymay $

#include "definitions.h"
#include "distanceTable.h"

void giveDistanceTable(const distanceMethod* dis,
		       const sequenceContainer& sc,
		       VVdouble& res,
		       vector<string>& names,
		       const vector<MDOUBLE> * weights){
	res.resize(sc.numberOfSeqs());
	for (int z=0; z< sc.numberOfSeqs();++z) res[z].resize(sc.numberOfSeqs(),0.0);

	for (int i=0; i < sc.numberOfSeqs();++i) {
		for (int j=i+1; j < sc.numberOfSeqs();++j) {
			res[i][j] = dis->giveDistance(sc[sc.placeToId(i)],sc[sc.placeToId(j)],weights,NULL);
			//LOG(5,<<"res["<<i<<"]["<<j<<"] ="<<res[i][j]<<endl);
		}
		names.push_back(sc[sc.placeToId(i)].name());
	}
}
