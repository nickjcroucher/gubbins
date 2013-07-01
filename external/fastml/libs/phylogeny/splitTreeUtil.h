// $Id: splitTreeUtil.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___SPLIT_TREE_UTIL
#define ___SPLIT_TREE_UTIL
#include "tree.h"
#include "split.h"

#include <vector>
#include <map>
using namespace std;


tree::nodeP findNodeToSplit(const tree& et,const split& mySplit,const map<string, int> & nameIdMap);
void applySplit(tree& et, const split& mySplit,const map<string, int> & nameIdMap);
void splitSonsFromNode(tree & et, tree::nodeP fatherNode, vector<tree::nodeP> & son2split);
void applySplitToRoot(tree& et, const split& mySplit,const map<string, int> & nameIdMap);
vector<tree::nodeP> findSonsThatHaveToBeSplit(const tree& et,const split& mySplit,const map<string, int> & nameIdMap);
bool childIsInTheSplit(const tree::nodeP & myNode, const split& mySplit,const map<string, int> & nameIdMap);



#endif



