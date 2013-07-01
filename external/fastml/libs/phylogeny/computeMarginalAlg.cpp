// $Id: computeMarginalAlg.cpp 1735 2007-02-26 13:46:37Z itaymay $

#include "definitions.h"
#include "treeIt.h"
#include "computeMarginalAlg.h"
#include <iostream>
#include <cassert>
using namespace std;


void computeMarginalAlg::fillComputeMarginal(const tree& et,
					   const sequenceContainer& sc,
					   const stochasticProcess& sp,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const suffStatGlobalHomPos& cdown,
					   doubleRep & posProb){

	// filling the exact probs.
	tree::nodeP mynode = NULL;
	ssc.allocatePlace(et.getNodesNum(),pi.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		assert (mynode != NULL);
		int letter;
		if (mynode->isLeaf()) {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				doubleRep val=convert(cup.get(mynode->id(),letter))?1.0:0.0;
				ssc.set(mynode->id(),letter,val);
			}
			continue;
		}
		doubleRep sumProb =0;
		for(letter=0; letter<pi.alphabetSize();letter++) {
			doubleRep prob=0.0;
			if (mynode->father()==NULL) prob=1.0; // special case of the root.
			else {
				for(int letter_in_f=0; letter_in_f<pi.alphabetSize();letter_in_f++) {
					prob +=cdown.get(mynode->id(),letter_in_f)*
					pi.getPij(mynode->id(),letter,letter_in_f);
				}
			}
			
			prob = prob*sp.freq(letter)*
				cup.get(mynode->id(),letter);
			ssc.set(mynode->id(),letter,prob);
			sumProb += prob;
		}
		for(letter=0; letter<pi.alphabetSize();letter++) {
			doubleRep getV = ssc.get(mynode->id(),letter);
			ssc.set(mynode->id(),letter,getV/sumProb);
		}

	

		// CHECKING:
/*		LOG(5,<<" checking marginal of node: "<<mynode->name()<<endl);
		MDOUBLE SSum =0;
		for (int u=0; u < pi.alphabetSize(); ++u) {
			LOG(5,<<ssc.get(mynode->id(),u)<<" ");
			SSum +=ssc.get(mynode->id(),u);
		}
		LOG(5,<<"\nsum of marginals = "<<SSum<<endl);
*/		
	if (mynode->isRoot()) posProb = convert(sumProb);
	}
}




/*
if (val>1) {
					LOG(5,<<"x val = " << val<<endl);
					LOG(5,<<" my node = " << mynode->name()<<endl);
					LOG(5,<<" let = " << let << endl);
					LOG(5,<<" up = " << cup.get(mynode->id(),let));
					LOG(5,<< "pos prob = " << posProb<<endl);
					LOG(5,<<" root of tree = " << et.getRoot()->name()<<endl);
					errorMsg::reportError(" error in compute marginal >1 ");
				}
if (val>1) {
					LOG(5,<<" val = " << val<<endl);
					LOG(5,<<" pos = " << pos<<endl);
					LOG(5,<<" my node = " << mynode->name()<<endl);
					LOG(5,<<" let = " << let << endl);
					LOG(5,<<" up = " << cup.get(mynode->id(),let)<<endl);
					LOG(5,<<" down[sameLetter] = " << cdown.get(mynode->id(),let)<<endl);
					LOG(5,<<" pij[sameLetter] = " << pi.getPij(mynode->id(),let,let)<<endl);
					LOG(5,<< "pos prob = " << posProb<<endl);
					LOG(5,<<" root of tree = " << et.getRoot()->name()<<endl);
					LOG(5,<<"sp.freq(letter) = "<<sp.freq(let)<<endl);
					errorMsg::reportError(" error in compute marginal >1 ");
				}


  */

