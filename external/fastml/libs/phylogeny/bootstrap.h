// $Id: bootstrap.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___BOOTSTRAP
#define ___BOOTSTRAP

#include "definitions.h"
#include "split.h"
#include "splitMap.h"
#include "tree.h"
#include "treeUtil.h"
#include <sstream>
using namespace std;

// this class gets as input many trees and can answer questions such as
// 1. the bootstrap value (bp) of a tree
// 2. the bp of a split
// 3. can reconstruct a multifurcating consensus trees.
// We note that 3 can always be done if done only on those splits with bp > 50%
// In this case there is only one tree.
// If the treshold value is <= 50% there might be more than one tree for which
// all splits on this tree have bp>= treshold. 
// In this case we want to give the tree with the highest sum of bp.
// This is probably NP hard, and we use a greedy search to chose
// this tree.

class bootstrap {
public:
  typedef  vector<tree> treeVec;
  explicit bootstrap(const treeVec& treevect);  // constructor 

  // this construction is the same as above, but it reads the trees from
  // an input file.
  explicit bootstrap (const string& filename);
  
  // give a tree and return a map from each edge to a bp value.
  // edge 5 is the edge between node id 5 and its father.
  map<int, MDOUBLE> getWeightsForTree(const tree& inTree) const;

  
  // give a threshold >= 0.5 and get a concensus tree with all splits
  // that are more confident then the threshold.
	tree consensusTree(const MDOUBLE threshold = 0.5) const;

  void print(ostream& sout = cout);
  void printTreeWithBPvalues(ostream &os, const tree &t, const map<int, MDOUBLE> & v, const bool printBranchLenght=true) const;
  
  void print_names(ostream &os) const;
  

private:




  void fillFromTreeVec(const treeVec& treevect);
  int idFromName (const string & name) const;


  set<int> recursivelyBuiltBPMap(const tree::nodeP &rootOfSubtree, map<int, MDOUBLE> &v) const;
  set<int> splitSubTreeRecursivly(const tree::nodeP &n, const bool isRoot=false);	// this function assumes that the tree is rooted not in a leaf
  // take tree, compute all splits and enter them into the Splits map
  void splitTree(const tree& T);
	void recursivlyPrintTreeWithBPvalues(ostream &os, 
									  const tree::nodeP &nP,
									  const map<int, MDOUBLE> &v,
									  const bool printBranchLenght) const;
  void getTreeNodes(const tree& t)  const ; // note that _allTree_nodes is mutable
  void updateNtaxaAndNameMapAndValidateConsistency(const tree& T);

  int _numTrees;			// total number of trees
  splitMap _Splits;
  typedef map<string,int> NameMap_t; 
  NameMap_t _nameMap; // this is a map from the names of the sequences to integers.
  int _nTaxa;
  mutable  vector<int> _id2TreeId, _treeId2Id;
  vector<string> _sequenceNames; // the names of the sequences.
};



#endif  // ___BOOTSTRAP

