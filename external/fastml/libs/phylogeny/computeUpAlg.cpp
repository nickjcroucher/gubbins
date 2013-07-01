// $Id: computeUpAlg.cpp 5988 2009-03-18 18:20:05Z itaymay $

#include "definitions.h"
#include "computeUpAlg.h"
#include "treeIt.h"
#include "seqContainerTreeMap.h"
#include "logFile.h"
#include <iostream>
#include <cassert>
using namespace std;

void computeUpAlg::fillComputeUp(const tree& et,
								 const sequenceContainer & sc,
								 const computePijGam& pi,
								 suffStatGlobalGam& ssc) {
	computeUpAlg cupAlg;
	ssc.allocatePlace(sc.seqLen(),pi.categories(),et.getNodesNum(),pi.alphabetSize());
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		for (int categor = 0; categor < pi.categories(); ++categor) {
			cupAlg.fillComputeUp(et,sc,pos,pi[categor],ssc[pos][categor]);
		}
	}
}

void computeUpAlg::fillComputeUp(const tree& et,
								 const sequenceContainer& sc,
								 const int pos,
								 const computePijHom& pi,
								 suffStatGlobalHomPos& ssc) {

	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),pi.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter;
		if (mynode->isLeaf()) {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
		}
		else {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				doubleRep total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons();++i){				
					doubleRep prob=0.0;
					for(int letInSon=0; letInSon<pi.alphabetSize();letInSon++) {
						prob += ssc.get(mynode->getSon(i)->id(), letInSon)*
							pi.getPij(mynode->getSon(i)->id(),letter,letInSon);
					}
				total_prob*=prob;
				}
				ssc.set(mynode->id(),letter,total_prob);
			}
		}
	}
}
/*
void computeUpAlg::fillComputeUp(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc) {

	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),sp.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter;
		if (mynode->isLeaf()) {// leaf
			for(letter=0; letter<sp.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				MDOUBLE val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
		}
		else {
			for(letter=0; letter<sp.alphabetSize();letter++) {
				MDOUBLE total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons();++i){				
					MDOUBLE prob=0.0;
					for(int letInSon=0; letInSon<sp.alphabetSize();letInSon++) {
						prob += ssc.get(mynode->getSon(i)->id(),letInSon)*
							sp.Pij_t(letter,letInSon,mynode->getSon(i)->dis2father()*sp.getGlobalRate());// taking care of the glubal is new.
					}
					assert(prob>=0.0);
				total_prob*=prob;
				}
				ssc.set(mynode->id(),letter,total_prob);
			}
		}
	}
}
*/
void computeUpAlg::fillComputeUpSpecificGlobalRate(const tree& et,
												   const sequenceContainer& sc,
												   const int pos,
												   const stochasticProcess& sp,
												   suffStatGlobalHomPos& ssc,
												   const MDOUBLE gRate) {
	if (sp.categories() >1) {// because we do not multiply all branch lengths by the rate[categories])
		errorMsg::reportError("the function fillComputeUpSpecificGlobalRate should not be used with a gamma model");
	}

	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),sp.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
#ifdef VERBOS
		LOG(15,<<endl<<endl<<"doing node: "<<mynode->name()<<endl);
#endif
		int letter;
		if (mynode->isLeaf()) {
			for(letter=0; letter<sp.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
		}
		else {
			int letterWithTotalProbEqZero =0;
			for(letter=0; letter<sp.alphabetSize();letter++) {
				doubleRep total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons();++i){				
					doubleRep prob=0.0;
					for(int letInSon=0; letInSon<sp.alphabetSize();letInSon++) {
						assert(ssc.get(mynode->getSon(i)->id(),letInSon)>=0);
						assert(sp.Pij_t(letter,letInSon,mynode->getSon(i)->dis2father()*gRate)>=0);
						prob += ssc.get(mynode->getSon(i)->id(),letInSon)*
							sp.Pij_t(letter,letInSon,mynode->getSon(i)->dis2father()*gRate);
					}
				assert(prob>=0.0);
				total_prob*=prob;
				}
				if (total_prob==0.0) ++letterWithTotalProbEqZero;
				
				ssc.set(mynode->id(),letter,total_prob);
			} // end of else
			if (letterWithTotalProbEqZero == sp.alphabetSize() && (mynode->getNumberOfSons() > 0)) {
				LOG(5,<<" total prob =0");
				for (int z=0; z <mynode->getNumberOfSons(); ++z) {
					LOG(5,<<"son "<<z<<" is "<<mynode->getSon(z)->name()<<endl);
					LOG(5,<<"dis2father is "<<mynode->getSon(z)->dis2father()<<endl);
					for(int letInSon=0; letInSon<sp.alphabetSize();letInSon++) {
						LOG(5,<<"let = "<<letInSon<<endl);
						LOG(5,<<"ssc.get(mynode->getSon(z)->id(),letInSon) = "<<convert(ssc.get(mynode->getSon(z)->id(),letInSon))<<endl);
					}
				}
				return;
			}
		}
	}
}
