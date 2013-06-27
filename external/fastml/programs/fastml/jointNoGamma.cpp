#include "jointNoGamma.h"
#include "treeIt.h"
#include "seqContainerTreeMap.h"
#include <fstream>
#include <cmath>
using namespace std;

jointNoGamma::jointNoGamma(const tree& et,
						   const stochasticProcess& sp,
						   const sequenceContainer& sc) 
						   : _et(et), _sp(sp), _sc(sc) {
	_cpih.fillPij(_et,_sp);
}

void jointNoGamma::compute() {
	
	suffStatGlobalHomPos ssc;
	suffStatGlobalHomPosJointNoGamma sscJointNoGam;
	ssc.allocatePlace(_et.getNodesNum(),_sc.alphabetSize());
	sscJointNoGam.allocatePlace(_et.getNodesNum(),_sc.alphabetSize());

	vector<string> ancestralSequences(_et.getNodesNum());
	MDOUBLE totalLikelihoodOfReconstruction = 0;
	cout<<"doing position (joint): ";
	for (int pos=0; pos<_sc.seqLen(); ++pos) {
		cout<<pos+1<<" ";
		fillComputeUp(pos,ssc,sscJointNoGam);
		doubleRep likelihoodOfPos = 0;

		vector<int> res =computeJointAncestralFromSSC(pos,ssc,sscJointNoGam,likelihoodOfPos);
		treeIterDownTopConst tIt(_et);
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
			if (mynode->isInternal()) { 
				ancestralSequences[mynode->id()]+=_sc.getAlphabet()->fromInt(res[mynode->id()]);
			}
		}
		_jointLikelihoodOfPositions.push_back(likelihoodOfPos);
	}
	cout<<endl;
	fromJointReconstructionToSequenceContainer(ancestralSequences);
}

void jointNoGamma::fillComputeUp(const int pos,
				   suffStatGlobalHomPos& ssc,
				   suffStatGlobalHomPosJointNoGamma& sscJointNoGam) {
	seqContainerTreeMap sctm(_sc,_et);
	ssc.allocatePlace(_et.getNodesNum(),_cpih.alphabetSize());
	treeIterDownTopConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf()) {// leaf
			for(int letterInFather=0; letterInFather<_cpih.alphabetSize();letterInFather++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				MDOUBLE totalVal = 0.0;
				for (int let=0; let<_cpih.alphabetSize();let++) {
					MDOUBLE val = _sc.getAlphabet()->relations(_sc[seqID][pos],let);
					if (val>0) {
						val*=_cpih.getPij(mynode->id(),letterInFather,let);
						totalVal +=val;
					}
				}
				//cerr<<"val =" << val <<" "; // REMOVE!
				//cerr<<"_pi->data(mynode->id(),pos)= "<<_pi->data(mynode->id(),pos)<<" ";//REMOVE
				ssc.set(mynode->id(),letterInFather,totalVal);
				sscJointNoGam.set(mynode->id(),letterInFather,_sc[seqID][pos]);
			}
		}
		else {
			for(int letterInFather=0; letterInFather<_cpih.alphabetSize();letterInFather++) {
				doubleRep maxProb=0.0;
				int bestLet = -1;
				for (int let=0; let<_cpih.alphabetSize();++let) {
					doubleRep tmpProb = 1;
					if (mynode->isRoot() == false) {
						tmpProb *= _cpih.getPij(mynode->id(),letterInFather,let);
					}
					for(int i=0; i < mynode->getNumberOfSons();++i){				
						tmpProb *= ssc.get(mynode->getSon(i)->id(),let);
					}
					if (tmpProb>maxProb) {
						maxProb = tmpProb;
						bestLet = let;
					}
				}
				ssc.set(mynode->id(),letterInFather,maxProb);
				assert(bestLet>=0);
				assert(bestLet<_cpih.alphabetSize());

				sscJointNoGam.set(mynode->id(),letterInFather,bestLet);
				if (mynode->isRoot()) break; // there's no meening to letterInFather in case of root.
			}
		}
	}
}

vector<int> jointNoGamma::computeJointAncestralFromSSC(
				   const int pos,
				   const suffStatGlobalHomPos& ssc,
				   const suffStatGlobalHomPosJointNoGamma& sscFASTML,
				   doubleRep & likelihoodOfReconstruction) {
	treeIterTopDownConst tIt(_et);
	vector<int> res(_et.getNodesNum());
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isRoot() == false) {
			int letterInFather = res[mynode->father()->id()];
			int tmp = sscFASTML.get(mynode->id(),letterInFather);
			res[mynode->id()] = tmp;
		} else {//special case of the root
			MDOUBLE maxL = VERYSMALL;
			int bestCharInRoot = sscFASTML.get(mynode->id(),0);
			likelihoodOfReconstruction = ssc.get(mynode->id(),0)*_sp.freq(bestCharInRoot);;
			res[mynode->id()] = bestCharInRoot;
		}
	}
	return res;
}

void jointNoGamma::fromJointReconstructionToSequenceContainer(const vector<string> & ancestralSequences){
	_resultSec = _sc;
	treeIterDownTopConst tIt2(_et);
	for (tree::nodeP mynode = tIt2.first(); mynode != tIt2.end(); mynode = tIt2.next()) {
		if (mynode->isInternal()) { 
			sequence tmp(ancestralSequences[mynode->id()],mynode->name(),"joint reconstruction",_resultSec.numberOfSeqs(),_sc.getAlphabet());
			_resultSec.add(tmp);
		}
	}
}

void jointNoGamma::outputTheJointProbAtEachSite(const string & outputFileProbJoint) {
	ofstream jointProbOutput(outputFileProbJoint.c_str());
	MDOUBLE totalLogLikelihood =0;
	for (int j=0; j < _jointLikelihoodOfPositions.size(); ++j) {
		totalLogLikelihood+=log(_jointLikelihoodOfPositions[j]);
		jointProbOutput<<"Joint log likelihood of position "<<j+1;// j+1 so that positions start from 1, and not from 0.
		jointProbOutput<<": "<<log(_jointLikelihoodOfPositions[j])<<endl;
	}
	jointProbOutput<<"total log likelihood of joint reconstruction: "<<totalLogLikelihood<<endl;
	jointProbOutput.close();
}


