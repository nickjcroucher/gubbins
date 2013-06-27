// $Id: bootstrap.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"
#include "someUtil.h"
#include "bootstrap.h"
#include "splitTreeUtil.h"
#include <algorithm>
#include <set>
using namespace std;

// -----------------------------------------------------------------------------------------
// ----------------------------- The constructor and its related functions -----------------
// -----------------------------------------------------------------------------------------

bootstrap::bootstrap(const treeVec& treevect):_numTrees(0), _nTaxa(0){
	fillFromTreeVec(treevect);
}
bootstrap::bootstrap (const string& filename):_numTrees(0), _nTaxa(0){
	fillFromTreeVec(getStartingTreeVecFromFile(filename));
}

void bootstrap::fillFromTreeVec(const treeVec& treevect) {
// for each tree, we compute the set of all splits.
// we update for each split in each tree the split-map.
// so we have the frequency of each split.
	for (treeVec::const_iterator i=treevect.begin();i!=treevect.end();++i)
		splitTree(*i);
}

// takes a tree, computes all splits and
// enter them into the Splits map
void bootstrap::splitTree(const tree& T){
	_numTrees++;
	updateNtaxaAndNameMapAndValidateConsistency(T);
	splitSubTreeRecursivly(T.getRoot(), true); // the true because we call the recursion with the root. Otherwise it is false;
}

void bootstrap::updateNtaxaAndNameMapAndValidateConsistency(const tree& T) {
	if (!_nTaxa) { // only for the first tree, this part intializes the _nameMap and the _nTaxa
		_sequenceNames = getSequencesNames(T);
		for (_nTaxa=0;_nTaxa<_sequenceNames.size();++_nTaxa) {
		  _nameMap[_sequenceNames[_nTaxa]] =_nTaxa;
		}
	}
	else {
		vector<string> namesInT1 = getSequencesNames(T);
		if (namesInT1.size() <  _nameMap.size()) {
			string errMs1 = "Not all trees have the same number of sequences. ";
			errMs1 += "tree number 1 has: ";
			errMs1 += int2string(_nameMap.size());
			errMs1 += " while tree number: ";
			errMs1 += int2string(_numTrees);
			errMs1 += " has ";
			errMs1 += int2string(namesInT1.size());
			errMs1 += "\nError in function bootstrap::splitTree";
			errorMsg::reportError(errMs1);
		}
		for (int i=0; i < namesInT1.size(); ++i) {
			if (_nameMap.count(namesInT1[i])==0) {
				string errMs = "The taxa ";
				errMs += namesInT1[i];
				errMs += " found in tree number ";
				errMs += int2string(_numTrees);
				errMs += " is not present in the first tree. Error in function bootstrap::splitTree";
				errorMsg::reportError(errMs);
			}
		}
	}
}

set<int> bootstrap::splitSubTreeRecursivly(const tree::nodeP &n,
										   const bool isRoot) {//false
// this function assumes that the root of the tree is not a leaf
	set<int> s; // the id of all leaves of the subtree of the nodeP n.
	for(int i=0; i<n->getNumberOfSons() ;++i) {
		set<int> sonSet(splitSubTreeRecursivly(n->getSon(i)));
		set<int>::iterator it = sonSet.begin();
		for (; it != sonSet.end(); ++it) s.insert(*it); 
	}
	if(isRoot) return s;
	if (n->isLeaf()) {
		s.insert(idFromName(n->name()));
	} else { // this avoids keeping track of trivial splits.
		set<int>::const_iterator sBeg(s.begin());
		set<int>::const_iterator sEnd(s.end());
		split sp(sBeg,sEnd,_nTaxa);
		_Splits.add(sp);
	}
  return(s);
}

// -----------------------------------------------------------------------------------------
// ----------------------------- getWeightsForTree -----------------------------------------
// -----------------------------------------------------------------------------------------

map<int, MDOUBLE>  bootstrap::getWeightsForTree(const tree& inTree) const {
  map<int, MDOUBLE> v;
  recursivelyBuiltBPMap(inTree.getRoot(), v);
  return (v);
}

// the function returns the ids of the leaves in the subtree defined by rootOfSubtree.
set<int> bootstrap::recursivelyBuiltBPMap(const tree::nodeP &rootOfSubtree, map<int, MDOUBLE> &v) const {
	set<int> s;
	for(int i=0;i<rootOfSubtree->getNumberOfSons();++i) {
		set<int> sonSet(recursivelyBuiltBPMap(rootOfSubtree->getSon(i),v));
		set<int>::iterator it = sonSet.begin();
		for (; it != sonSet.end(); ++it) s.insert(*it); 
	}
	if (rootOfSubtree->isLeaf()) {
		s.insert(idFromName(rootOfSubtree->name()));
	}
	set<int>::const_iterator sBeg(s.begin());
	set<int>::const_iterator sEnd(s.end());
	split sp(sBeg,sEnd,_nTaxa);
	v[rootOfSubtree->id()]=(static_cast<MDOUBLE>(_Splits.counts(sp)))/_numTrees;
  return(s);
}

// We get different trees, and the id's are not consistent among different trees.
// here, we map a name to a single id.
int bootstrap::idFromName(const string & name) const {
	NameMap_t::const_iterator i(_nameMap.find(name));
	if (i==_nameMap.end()) {
		string s="Can not find an Id for the taxa name:";
		s+=name;
		s+="\n error in function bootstrap::idFromName\n";
      errorMsg::reportError(s);
    }
	return(i->second);
}

// -----------------------------------------------------------------------------------------
// ----------------------------- Printing the bp  ------------------------------------------
// -----------------------------------------------------------------------------------------

void bootstrap::print(ostream& sout){//  = cout
  _Splits.print(sout);
}

void bootstrap::printTreeWithBPvalues(ostream &out, const tree &t, const map<int, MDOUBLE> & v, const bool printBranchLenght) const{
  recursivlyPrintTreeWithBPvalues(out,t.getRoot(),v, printBranchLenght);
	out<<";";
}

void bootstrap::recursivlyPrintTreeWithBPvalues(ostream &out, 
									  const tree::nodeP &myNode,
												const map<int, MDOUBLE> &v,
												const bool printBranchLenght) const {
	if (myNode->isLeaf()) {
		out << myNode->name();
		if (printBranchLenght) out << ":"<<myNode->dis2father();
		return;
	} else {
		out <<"(";
		for (int i=0;i<myNode->getNumberOfSons();++i) {
			if (i>0) out <<",";
			recursivlyPrintTreeWithBPvalues(out, myNode->getSon(i),v, printBranchLenght);
		}
		out <<")";
		if (myNode->isRoot()==false) {
			if (printBranchLenght) out<<":"<<myNode->dis2father(); 
			map<int,MDOUBLE>::const_iterator val=v.find(myNode->id());
			if ((val!=v.end()) &&  val->second>0.0) {
				out << "["<<val->second<<"]";
			} 
		}
	}
}

// for DEBUGGING ONLY:
void bootstrap::print_names(ostream &out) const {
  NameMap_t::const_iterator i(_nameMap.begin());
  for (;i!=_nameMap.end();++i)
    out << "{"<<i->first<<" = "<<i->second<<"}"<<endl;
}

// -----------------------------------------------------------------------------------------
// ----------------------------- Building consensus tree  ----------------------------------
// -----------------------------------------------------------------------------------------
// returns the bp values of the consensus tree.
// the idea is to start from the split map, extract a split at a time.
// first, the splits with the highest bp (i.e., in a sorted way).
// Each splits is checked for compatibility with the consensus tree constructed so far.
// if it is compatible, it is added to the consensus.
// Otherwise - it is discarded.
// returns the consensus tree
tree bootstrap::consensusTree(const MDOUBLE threshold) const {// =0.5
// 1. get the names of the sequences
	vector<string> names;
	for (NameMap_t::const_iterator i(_nameMap.begin());i!=_nameMap.end();++i)
		names.push_back(i->first);

// 2. create a star tree
	tree res = starTree(names);

// 3. get the sorted vector of the splits from which the consensus is to be built.
	vector<pair<split,int> > sortedSplits = _Splits.sortSplits();
// 4. get a list of compatible splits
	MDOUBLE thresholdForNumTrees = threshold * _numTrees;

	vector<split> consensus;
	for (int k=0; k < sortedSplits.size(); ++k) { 
		bool compatible = true;
        if (sortedSplits[k].second < thresholdForNumTrees) break;

		for (vector<split>::const_iterator j=consensus.begin(); j != consensus.end(); ++j) {
			if (!(sortedSplits[k].first.compatible(*j))) {
				compatible=false;
				 break;
			}
		}
        if (compatible) {
			consensus.push_back(sortedSplits[k].first);
		}
	}

// 5. Now we build a tree from all the compatible splits

	for (vector<split>::iterator i1 = consensus.begin();i1!=consensus.end();++i1) {
	  applySplit(res,*i1,_nameMap);
	}
	res.create_names_to_internal_nodes();
	res.makeSureAllBranchesArePositive();

	return (res);
}
