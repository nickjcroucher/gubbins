#include "bbNodeOrderAlg.h"
#include "bbComputeUpAlg.h"
#include "bbComputeDownAlg.h"
#include "computeMarginalAlg.h"
#include <algorithm>
using namespace std;

bbNodeOrderAlg::bbNodeOrderAlg(const tree& et,
					const stochasticProcess &sp,
					const sequenceContainer& sc,
					const computePijGam& cpij,
					const MDOUBLE computeAgainExactTreshold) :_et(et),_sp(sp),_sc(sc),_cpij(cpij){
	_alphabetSize=_sp.alphabetSize();
	_computeAgainExactTreshold = computeAgainExactTreshold;
	cupbb.allocatePlace(sp.categories(),et.getNodesNum(),sp.alphabetSize());
	cdownbb.allocatePlace(sp.categories(),et.getNodesNum(),sp.alphabetSize());
	cmarginalbb.allocatePlace(sp.categories(),et.getNodesNum(),sp.alphabetSize());
}

bbNodeOrderAlg::~bbNodeOrderAlg(){}

// note: there is a way to dynamically correct exact.
// it is not implemented here.
void bbNodeOrderAlg::getNextNodeAndCharOrder(tree::nodeP &nextNode,
							 vector<int> &charOrder,
							 vector<sequence> &ancestralSequences,
							 const int pos,
							 const bool firstTime,
							 doubleRep& exactVal){
	doubleRep highestProb=0;
	if (firstTime) {
		_et.getAllHTUs(_nodesLeft,_et.getRoot());
		recalculateExact(ancestralSequences,pos);
		rankRemainingNodesAccordingToTheirMarginalProb(pos);
	}
	assert(_nodesLeftExact.size()>=1);
	assert(_nodesLeftExact.size()==_nodesLeft.size());
	highestProb = _nodesLeftExact[_nodesLeftExact.size()-1];
	if (highestProb<_computeAgainExactTreshold) {
		recalculateExact(ancestralSequences,pos);
		rankRemainingNodesAccordingToTheirMarginalProb(pos);
		highestProb = _nodesLeftExact[_nodesLeftExact.size()-1];
	}
	_nodesLeftExact.pop_back();
	nextNode = _nodesLeft[_nodesLeft.size()-1];
	_nodesLeft.pop_back();
	charOrder = findBestOrderInNode(nextNode,pos);
	exactVal = highestProb;
}

void bbNodeOrderAlg::putBack(tree::nodeP& node2check,const doubleRep & exactVal) {
	_nodesLeft.push_back(node2check);
	_nodesLeftExact.push_back(exactVal);
}


void bbNodeOrderAlg::rankRemainingNodesAccordingToTheirMarginalProb(
			const int pos) {

	typedef pair<doubleRep,tree::nodeP> sortedElement;
	vector<sortedElement> sortVec;
	int i;
	doubleRep tmpVal;
	for ( i = 0 ; i < _nodesLeft.size() ; ++i) {
		tmpVal = getNodeHighestMarginal(_nodesLeft[i]);
		sortedElement elem(tmpVal,_nodesLeft[i]);
		sortVec.push_back(elem);
	}

	sort(sortVec.begin(), sortVec.end());
	_nodesLeft.clear();
	_nodesLeftExact.clear();
	_nodesLeft.resize(sortVec.size());
	_nodesLeftExact.resize(sortVec.size());
	for ( i = 0 ; i < _nodesLeft.size() ; ++i ) {
		_nodesLeft[i] = sortVec[i].second;
		_nodesLeftExact[i]=sortVec[i].first;
	}
}

// this function gets as input the "exact" sufficient statistic for a given node
// for a given position. It goes over all the alphabet, and computes
// the marginal at each position. Then he returns the highest marginal.
doubleRep bbNodeOrderAlg::getNodeHighestMarginal(const tree::nodeP& inNodeP) {
	doubleRep highestProb =0.0;

	int j,s;
	for (j=0;j<_alphabetSize;++j) {
		doubleRep tmpVal = 0;
		for (s=0; s< _sp.categories();++s ) {
			tmpVal += cmarginalbb.get(s,inNodeP->id(),j)*_sp.ratesProb(s);
		}
		if (highestProb<tmpVal) {
			highestProb=tmpVal;
		}
	}
	return highestProb;
}

void bbNodeOrderAlg::recalculateExact(vector<sequence> &ancestralSequences,
									  const int pos) {
	for (int i=0; i < _sp.categories(); ++i) {
		BBfillComputeUp(_et,_sc,pos,_cpij[i],cupbb[i],ancestralSequences);
		BBfillComputeDown(_et,_sc,pos,_cpij[i],cdownbb[i],cupbb[i],ancestralSequences);
		doubleRep posProb = 0.0;
		computeMarginalAlg cmalg;
		cmalg.fillComputeMarginal(_et,_sc,_sp,pos,_cpij[i],cmarginalbb[i],cupbb[i],cdownbb[i],posProb);
	}
}

vector<int> bbNodeOrderAlg::findBestOrderInNode(const tree::nodeP node2check,
												const int pos) const {
	assert (node2check != NULL);
	typedef pair<doubleRep,int> sortedElement; // (marginal, letter)
	vector<sortedElement> sortVec;
	int i,s;
	for ( i = 0 ; i < _alphabetSize ; i++ ) {
		doubleRep tmpVal = 0;
		for (s=0; s< _sp.categories();++s ) {
			tmpVal += cmarginalbb.get(s,node2check->id(),i)*_sp.ratesProb(s);
		}
		sortedElement elem(tmpVal,i);
		sortVec.push_back(elem);
	}

	sort(sortVec.begin(), sortVec.end());
	reverse(sortVec.begin(), sortVec.end());
	vector<int> bestCharOrder(_alphabetSize);
	for ( i = 0 ; i < _alphabetSize ; i++ ) {
		bestCharOrder[i] = sortVec[i].second;
	}
	return bestCharOrder;
}

