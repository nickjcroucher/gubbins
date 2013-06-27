#include "computeMarginalReconstruction.h"
#include "computeUpAlg.h"
#include "computePijComponent.h"

#include "computeDownAlg.h"
#include "computeMarginalAlg.h"
#include "treeIt.h"
#include <algorithm>
#include <iostream>
#include <fstream>

using namespace std;


computeMarginalReconstruction::computeMarginalReconstruction(const tree& et, 
															 vector<stochasticProcess>& spVec, 
															 const sequenceContainer& sc) : _et(et), _spVec(spVec), _sc(sc) {
	_resultProb.resize(_sc.seqLen());
	_bestProb.resize(_sc.seqLen());
	for (int i=0; i < _sc.seqLen(); ++i) {
		_resultProb[i].resize(et.getNodesNum());
		_bestProb[i].resize(et.getNodesNum());
		for (int j=0; j < et.getNodesNum(); ++j) {
			_resultProb[i][j].resize(_spVec[0].alphabetSize(),0.0);
		}
	}
}



void computeMarginalReconstruction::compute(const distribution * forceDistr){
	computePijGam pi;
	if (_spVec.size()>1) {//w codon model + gamma special case
		pi._V.resize(forceDistr->categories());
		for (int i=0; i < _spVec.size(); ++i) 
			pi._V[i].fillPij(_et,_spVec[i]);
		_spVec[0].setDistribution(forceDistr);//update the first process with gamma distr 
												//for all the functions that needs no catregor and categor probabilty
	}
	else{
		pi.fillPij(_et,_spVec[0]); 
	}

	//pi.fillPij(_et,_sp);
	MDOUBLE totalLikelihoodOfReconstruction = 0;
	cout<<"doing position (marginal): ";
	for (int pos=0; pos<_sc.seqLen(); ++pos) {
		suffStatGlobalGamPos sscUp;// this is for a specific position.
		suffStatGlobalGamPos sscDown;// this is for a specific position.
		suffStatGlobalGamPos sscMarginal; // this is for a specific position.
		sscUp.allocatePlace(_spVec[0].categories(),_et.getNodesNum(),_sc.alphabetSize());
		sscDown.allocatePlace(_spVec[0].categories(),_et.getNodesNum(),_sc.alphabetSize());
		sscMarginal.allocatePlace(_spVec[0].categories(),_et.getNodesNum(),_sc.alphabetSize());

		cout<<pos+1<<" ";
		computeUpAlg computeUpAlg1;
		computeDownAlg computeDownAlg1;
		computeMarginalAlg computeMarginalAlg1;
	
		for (int cat = 0; cat < _spVec[0].categories(); ++cat) {
            computeUpAlg1.fillComputeUp(_et,_sc,pos,pi[cat],sscUp[cat]);
			computeDownAlg1.fillComputeDown(_et,_sc,pos,pi[cat],sscDown[cat],sscUp[cat]);
			doubleRep posProb =0;
			computeMarginalAlg1.fillComputeMarginal(_et,_sc,_spVec[0],pos,pi[cat],sscMarginal[cat],sscUp[cat],sscDown[cat],posProb);
		}

		MDOUBLE likelihoodOfPos = 0;

		fillResultProb(sscMarginal,_spVec[0],_et,pos);
		fillMarginalReconstruction();
	}
	cout<<endl;
}

void computeMarginalReconstruction::fillResultProb(
	const suffStatGlobalGamPos& ssc,
	const stochasticProcess & sp,
	const tree& et,
	const int pos){
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		for (int i=0; i < sp.alphabetSize(); ++i) {
			doubleRep tmp=0; // the value for this letter in this node.
			for (int j=0; j < sp.categories(); ++j) {
				tmp += ssc.get(j,mynode->id(),i)*sp.ratesProb(j);
			}
			_resultProb[pos][mynode->id()][i] = convert(tmp);
		}
	}
}

void computeMarginalReconstruction::fillMarginalReconstruction() {
	_resultSec = _sc;
	treeIterTopDownConst tIt(_et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		if (mynode->isLeaf()) continue;
		// creating the place for this sequence in the resulting sequence container
		sequence tmp("",mynode->name(),"",_resultSec.numberOfSeqs(),_sc.getAlphabet());
		_resultSec.add(tmp);
		fillMarginalReconstructionSpecificNode(mynode);
	}
}

void computeMarginalReconstruction::fillMarginalReconstructionSpecificNode(tree::nodeP mynode) {
	for (int pos=0; pos < _sc.seqLen(); ++pos) {
		MDOUBLE bestP =-1.0;
		int bestChar = -1;
		for (int letter=0; letter < _spVec[0].alphabetSize(); ++letter) {
			if (_resultProb[pos][mynode->id()][letter] > bestP) {
				bestP = _resultProb[pos][mynode->id()][letter];
				bestChar = letter;
			}
		}
		_bestProb[pos][mynode->id()] = bestP;

		// adding bestChar to the resulting sequence container.
		string res = _sc.getAlphabet()->fromInt(bestChar);
		int id = _resultSec.getId(mynode->name());
		_resultSec[id].addFromString(res);
	}
}

void computeMarginalReconstruction::outputTheMarginalProbForEachCharForEachNode(const string& outputFileName) {
	ofstream out(outputFileName.c_str());
	for (int pos=0; pos<_sc.seqLen(); ++pos) {
        outputTheMarginalProbForEachCharForEachNodePos(out,pos);
	}
	out.close();
}

void computeMarginalReconstruction::outputTheMarginalProbForEachCharForEachNodePos(ostream& out,const int pos){//(DEFAULT = JPF, same file as above).
	treeIterDownTopConst tIt(_et);
	out<<"marginal probabilities at position: "<<pos+1<<endl;
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		//if (mynode->isLeaf()) continue;
		out<<"of node: "<<mynode->name()<<": ";
		vector<pair< MDOUBLE,string> > pres;
		int c=0;
		for (c=0; c < _spVec[0].alphabetSize(); ++c) {
			pres.push_back(pair<MDOUBLE,string>(_resultProb[pos][mynode->id()][c],_sc.getAlphabet()->fromInt(c)));
		}
		sort(pres.begin(),pres.end());
		for (c=pres.size()-1; c >=0 ; --c) {
			if (pres[c].first<0.0001) continue;
			out<<"p("<<pres[c].second;
			out<<")="<<pres[c].first<<" ";
		}
		out<<endl;
	}
	out<<endl;
}

