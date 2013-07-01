// $Id: computeUpAlgFactors.cpp 1738 2007-02-26 13:49:16Z itaymay $

#include "definitions.h"
#include "computeUpAlg.h"
#include "seqContainerTreeMap.h"
#include "logFile.h"
#include <iostream>
#include <cassert>
#include <cmath>
#include <cstdlib>
using namespace std;

void computeNodeFactorAndSetSsc(MDOUBLE & minFactor,suffStatGlobalHomPos& ssc, int nodeId, const int alphSize){
	// given a number = probability (val), it is changed to a new number which is 10 to the power of factor + val.
	// for example if val = 0.001, it is changed to 0.1 and factor 2.
	minFactor=100000;
	for (int i=0; i < alphSize; ++i) {
		MDOUBLE tmpfactor=0;
		doubleRep val = ssc.get(nodeId,i);
		if (val >0) {
			while (val < 0.1) {
				val *=10;
				tmpfactor++;
			}
		}
		else tmpfactor=minFactor;
		if (tmpfactor<minFactor) minFactor=tmpfactor;
	}
	for (int j=0; j < alphSize; ++j) {
		doubleRep tmp = ssc.get(nodeId,j);
		tmp = tmp * pow(static_cast<MDOUBLE>(10.0),minFactor);
		ssc.set(nodeId,j,tmp);
	}
}

void computeUpAlg::fillComputeUpWithFactors(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const computePijHom& pi,
				   suffStatGlobalHomPos& ssc,
				   vector<MDOUBLE>& factors) {
	factors.resize(et.getNodesNum(),0.0);
	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),pi.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter;
		if (mynode->getNumberOfSons() == 0) {// leaf
			for(letter=0; letter<pi.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
			computeNodeFactorAndSetSsc(factors[mynode->id()],ssc,mynode->id(),pi.alphabetSize());
		}
		else {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				doubleRep total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons(); ++i){				
					doubleRep prob=0.0;
					for(int letInSon=0; letInSon<pi.alphabetSize();letInSon++) {
						prob += ssc.get(mynode->getSon(i)->id(),letInSon)*
							pi.getPij(mynode->getSon(i)->id(),letter,letInSon);
					}
					total_prob*=prob;
				}
				ssc.set(mynode->id(),letter,total_prob);
			}
			computeNodeFactorAndSetSsc(factors[mynode->id()],ssc,mynode->id(),pi.alphabetSize());
			for(int k=0; k < mynode->getNumberOfSons();++k) {
				factors[mynode->id()]+=factors[mynode->getSon(k)->id()];
			}
		}
	}
}

void computeUpAlg::fillComputeUpWithFactors(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   vector<MDOUBLE>& factors) {
	factors.resize(et.getNodesNum(),0.0);
	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),sp.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter;
		if (mynode->getNumberOfSons() == 0) {// leaf
			for(letter=0; letter<sp.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
			computeNodeFactorAndSetSsc(factors[mynode->id()],ssc,mynode->id(),sp.alphabetSize());
		}
		else {
			for(letter=0; letter<sp.alphabetSize();letter++) {
				doubleRep total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons();++i){				
					doubleRep prob=0.0;
					for(int letInSon=0; letInSon<sp.alphabetSize();letInSon++) {
						prob += ssc.get(mynode->getSon(i)->id(),letInSon)*
							sp.Pij_t(letter,letInSon,mynode->getSon(i)->dis2father()*sp.getGlobalRate());// taking care of the glubal is new.
					}
					assert(prob>=0);
				total_prob*=prob;
				}
				ssc.set(mynode->id(),letter,total_prob);
			}
			computeNodeFactorAndSetSsc(factors[mynode->id()],ssc,mynode->id(),sp.alphabetSize());
			for(int k=0; k < mynode->getNumberOfSons();++k) {
				factors[mynode->id()]+=factors[mynode->getSon(k)->id()];
			}
		}
	}
}

void computeUpAlg::fillComputeUpSpecificGlobalRateFactors(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const stochasticProcess& sp,
				   suffStatGlobalHomPos& ssc,
				   const MDOUBLE gRate,
				   vector<MDOUBLE>& factors) {
	factors.resize(et.getNodesNum(),0.0);
	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),sp.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
#ifdef VERBOS
		LOG(5,<<endl<<endl<<"doing node: "<<mynode->name()<<endl);
#endif
		int letter;
		if (mynode->getNumberOfSons() == 0) {// leaf
			for(letter=0; letter<sp.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				doubleRep val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
			computeNodeFactorAndSetSsc(factors[mynode->id()],ssc,mynode->id(),sp.alphabetSize());
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
				assert(prob>=0);
				total_prob*=prob;
				}
				if (total_prob ==0) ++letterWithTotalProbEqZero;
				
				ssc.set(mynode->id(),letter,total_prob);
			} // end of else
			computeNodeFactorAndSetSsc(factors[mynode->id()],ssc,mynode->id(),sp.alphabetSize());
			for(int k=0; k < mynode->getNumberOfSons();++k) {
				factors[mynode->id()]+=factors[mynode->getSon(k)->id()];
			}
			if (letterWithTotalProbEqZero == sp.alphabetSize() && (mynode->getNumberOfSons() > 0)) {
				LOG(5,<<" total prob =0");
				for (int z=0; z <mynode->getNumberOfSons(); ++z) {
					LOG(5,<<"son "<<z<<" is "<<mynode->getSon(z)->name()<<endl);
					LOG(5,<<"dis2father is "<<mynode->getSon(z)->dis2father()<<endl);
					for(int letInSon=0; letInSon<sp.alphabetSize();letInSon++) {
						LOG(5,<<"let = "<<letInSon<<endl);
						LOG(5,<<"ssc.get(mynode->sons[z]->id(),letInSon) = "<<convert(ssc.get(mynode->getSon(z)->id(),letInSon))<<endl);
//						LOG(5,<<"sp.Pij_t(letter,letInSon,mynode->getSon(i)->dis2father()*gRate) = "<<sp.Pij_t(letter,letInSon,mynode->sons[i]->dis2father()*gRate)<<endl);
//						LOG(5,<<"mynode->getSon(i)->dis2father() = "<<mynode->getSon(i)->dis2father()<<endl);

					
					
					
					
					}
				}
				exit(3);
			}
		}
	}
}
