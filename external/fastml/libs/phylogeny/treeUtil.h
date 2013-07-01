// $Id: treeUtil.h 6091 2009-04-20 08:31:23Z rubi $

#ifndef ___TREE_UTIL
#define ___TREE_UTIL
#include "definitions.h"
#include "tree.h"

vector<tree> getStartingTreeVecFromFile(string fileName);

tree starTree(const vector<string>& names);

void getStartingTreeVecFromFile(string fileName,
											vector<tree>& vecT,
											vector<char>& constraintsOfT0);


bool sameTreeTolopogy(tree t1, tree t2);

bool cutTreeToTwo(tree bigTree,
			  const string& nameOfNodeToCut,
			  tree &small1,
			  tree &small2);

tree::nodeP makeNodeBetweenTwoNodes(	tree& et,
										tree::nodeP nodePTR1,
										tree::nodeP nodePTR2,
										const string &interName);

void cutTreeToTwoSpecial(const tree& source,
						tree::nodeP intermediateNode,
						tree &resultT1PTR,
						tree &resultT2PTR);

vector<string> getSequencesNames(const tree& t);

MDOUBLE getSumOfBranchLengths(const tree &t);

void printDataOnTreeAsBPValues(ostream &out, Vstring &data, tree &tr) ;
void printDataOnTreeAsBPValues(ostream &out, Vstring &data, const tree::nodeP &myNode) ;

MDOUBLE getDistanceFromNode2ROOT(const tree::nodeP &myNode);
void fillAllNodesNames(Vstring& Vnames,const tree& tr);

void printTreeWithValuesAsBP(ostream &out, const tree &tr, Vstring values, VVVdouble *probs, int from, int to);
void printTreeWithValuesAsBP(ostream &out, const tree::nodeP &myNode, Vstring values,  VVVdouble *probs, int from, int to);


#endif

