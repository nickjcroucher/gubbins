// $Id: computeCounts.cpp 4583 2008-08-05 15:02:26Z cohenofi $

#include "computeCounts.h"
void computeCounts::computeCountsNodeFatherNodeSonHomPos(const sequenceContainer& sc,
														 const computePijHom& pi,
														 const stochasticProcess& sp,
														 const suffStatGlobalHomPos& cup,
														 const suffStatGlobalHomPos& cdown,
														 const MDOUBLE weight,
														 const doubleRep posProb,
														 const tree::nodeP nodeSon,
														 countTableComponentHom& _ctc,
														 const MDOUBLE rateCategorProb 
												   ) 
{
	assert(posProb>0.0);
	if (weight == 0) return;
	int alph1,alph2;
	for (alph1 =0; alph1< pi.alphabetSize(); ++alph1) {
		for (alph2 =0; alph2< pi.alphabetSize(); ++alph2) {
			doubleRep tmp = cup.get(nodeSon->id(),alph1) *
				cdown.get(nodeSon->id(),alph2) *
				pi.getPij(nodeSon->id(),alph1,alph2)*
				sp.freq(alph1)
				* rateCategorProb
				/
				posProb;
			_ctc.addToCounts(alph1,alph2,convert(tmp)*weight);
		}
	}
}

void computeCounts::computeCountsNodeFatherNodeSonHomPos(const sequenceContainer& sc,
												   const computePijHom& pi,
												   const stochasticProcess& sp,
												   const suffStatGlobalHomPos& cup,
												   const suffStatGlobalHomPos& cdown,
												   const MDOUBLE weight,
												   const doubleRep posProb,
												   const tree::nodeP nodeSon,
												   countTableComponentHom& _ctc,
												   const MDOUBLE rateCategorProb,
												   const int letterInRoot
												   ) 
{
	assert(posProb>0.0);
	if (weight == 0) return;
	int alph1,alph2;
	for (alph1 =0; alph1< pi.alphabetSize(); ++alph1) {
		for (alph2 =0; alph2< pi.alphabetSize(); ++alph2) {
			doubleRep tmp = cup.get(nodeSon->id(),alph1) *
			cdown.get(nodeSon->id(),alph2) *				// down was given with specific root
			pi.getPij(nodeSon->id(),alph1,alph2)*
			sp.freq(letterInRoot)							// fixed root
			* rateCategorProb 
			* sp.freq(letterInRoot)						// to account for the additional letterAtRoot loop 
			/posProb;
			_ctc.addToCounts(alph1,alph2,convert(tmp)*weight);
		}
	}
}



void computeCounts::fillCountTableComponentGam(countTableComponentGam& ctcGam,
								const stochasticProcess& sp,
								const sequenceContainer& sc,
								const computePijGam& pij0,
								const suffStatGlobalGam& cup,
								const suffStatGlobalGam& cdown,
								const Vdouble * weights,
								tree::nodeP nodeSon,
								const VdoubleRep& posProbVec) {
	ctcGam.countTableComponentAllocatePlace(sp.alphabetSize(),sp.categories());
	for (int rateCat =0; rateCat< sp.categories(); ++ rateCat) {
		fillCountTableComponentGamSpecRateCategor(rateCat,ctcGam[rateCat],sp,
													sc,pij0[rateCat],
													cup,cdown,weights,posProbVec,nodeSon);
	}
}

void computeCounts::fillCountTableComponentGamSpecRateCategor(const int rateCategor,
											   countTableComponentHom& ctcHom,
											   const stochasticProcess& sp,
											   const sequenceContainer& sc,
											   const computePijHom& pi,
											   const suffStatGlobalGam& cup,
												const suffStatGlobalGam& cdown,
												const Vdouble * weights,
												const VdoubleRep& posProbVec, //prob of the position with gamma
												tree::nodeP nodeSon) {
	computeCounts cc;
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		MDOUBLE weig = (weights ? (*weights)[pos] : 1.0);
		cc.computeCountsNodeFatherNodeSonHomPos(sc,pi,sp,cup[pos][rateCategor],
												cdown[pos][rateCategor],
												weig,posProbVec[pos],nodeSon,
												ctcHom,sp.ratesProb(rateCategor)); 
	}
}
/*
void computeCounts::computeCountsNodeXNodeYHomPos(
						const tree::nodeP nodeX,
						const tree::nodeP nodeY) {

	const tree::nodeP nodeFather = nodeSon->father();
	_ctc.zero();
	if (_weight!=NULL) {	// this is one of the MAIN LOOPS.  no "if"s deep inside it!
		for (int pos=0; pos< _pi.seqLen(); ++pos) {
			if ((*_weight)[pos] == 0) continue;
			for (int alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
				for (int alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
					for (int rate =0; rate< _pi.categories(); ++rate) {
						MDOUBLE tmp = _cup.get(nodeSon->id(),pos,rate,alph1) *
						_cdown.get(nodeSon->id(),pos,rate,alph2) *
						_pi.pij(pos)->getPij(nodeSon->id(),alph1,alph2,rate)*
						_pi.stocProcessFromPos(pos)->freq(alph1)/
						_cprobAtEachPos.getProb(pos);
						_ctc.addToCounts(alph1,alph2,rate,tmp*(*_weight)[pos]);
					}
				}
			}
		}
	}
	else {
		for (int pos=0; pos< _pi.seqLen(); ++pos) {
			for (int alph1 =0; alph1< _pi.alphabetSize(); ++alph1) {
				for (int alph2 =0; alph2< _pi.alphabetSize(); ++alph2) {
					for (int rate =0; rate< _pi.categories(); ++rate) {
						MDOUBLE tmp = _cup.get(nodeSon->id(),pos,rate,alph1) *
						_cdown.get(nodeSon->id(),pos,rate,alph2) *
						_pi.pij(pos)->getPij(nodeSon->id(),alph1,alph2,rate)*
						_pi.stocProcessFromPos(pos)->freq(alph1)/
						_cprobAtEachPos.getProb(pos);
						_ctc.addToCounts(alph1,alph2,rate,tmp);
					}
				}
			}
		}
	}
	*/

