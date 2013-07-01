// $Id: splitTreeUtil.cpp 962 2006-11-07 15:13:34Z privmane $

#include "splitTreeUtil.h"
#include "someUtil.h"

static int idFromName(const string name, const map<string, int> & nameIdMap)
{
	map<string, int>::const_iterator i=nameIdMap.find(name);
	if (i==nameIdMap.end()) errorMsg::reportError(" error in splitTreeUtil. Name not found in nameIdMap");
	return (i->second);
}

// returns true if all the sons of myNode are in the split.
// return false if all the sons of myNode are NOT in the split
// if some of the sons are in and some are not - set foundTheNodeAlready to true.
// and set splitNode to be that node.
static bool findNodeToSplitRecursive(		const tree::nodeP myNode,
												const split& mySplit,
												tree::nodeP& splitNode,
												bool & foundTheNodeAlready,
												const map<string, int> & nameIdMap) {
if (myNode->isLeaf()) return (mySplit.isMember(idFromName(myNode->name(),nameIdMap)));	
bool inSplit = findNodeToSplitRecursive(myNode->getSon(0),mySplit,splitNode,foundTheNodeAlready,nameIdMap);
	if (foundTheNodeAlready) return true;
	for (int i=1; i < myNode->getNumberOfSons(); ++i) {
		bool tmp = findNodeToSplitRecursive(myNode->getSon(i),mySplit,splitNode,foundTheNodeAlready,nameIdMap);
		if (foundTheNodeAlready) return true;
		if (tmp != inSplit) {
			foundTheNodeAlready = true;
			splitNode = myNode;
			return true;
		}
	}
	return inSplit;
}



tree::nodeP findNodeToSplit(const tree& et,
							const split& mySplit,
							const map<string, int> & nameIdMap) {
	tree::nodeP res;
	bool foundTheNodeAlready = false;
	findNodeToSplitRecursive(et.getRoot(),mySplit,res,foundTheNodeAlready,nameIdMap);
	return res;
}

void applySplit(tree& et,
				const split& mySplit,
				const map<string, int> & nameIdMap) {
	tree::nodeP node2split = findNodeToSplit(et,mySplit,nameIdMap);
	et.rootAt(node2split);
	applySplitToRoot(et,mySplit,nameIdMap);
}

void splitSonsFromNode(tree & et, tree::nodeP fatherNode, vector<tree::nodeP> & son2split)
{
	for (int k=0; k < son2split.size(); ++k) {
		if (son2split[k]->father() != fatherNode ) 
		errorMsg::reportError(" error in function bootstrap::splitSonsFromNode - nodes don't have the same father");
	}
	// if the split allready exists, we do not need to do anything.
	if (son2split.size()==fatherNode->getNumberOfSons() // the branch above us is the required split
	    || son2split.size() <=1 // the branch below us is it
	    || (fatherNode->father()==NULL &&  son2split.size()==fatherNode->getNumberOfSons()-1) 
		// the branch above us is the required split
	    )
	  return;

	tree::nodeP theNewNode = et.createNode(fatherNode,et.getNodesNum());
	theNewNode->setName("N"+int2string(theNewNode->id()));
	for (int i=0; i < son2split.size(); ++i) {
		son2split[i]->setFather(theNewNode);
		theNewNode->setSon(son2split[i]);
		// remove from son list of father node.
		fatherNode->removeSon(son2split[i]); 
	}
}

void applySplitToRoot(tree& et,
					  const split& mySplit,
					  const map<string, int> & nameIdMap) {
	vector<tree::nodeP> sonsThatHaveToBeSplit = findSonsThatHaveToBeSplit(et,mySplit,nameIdMap);
	splitSonsFromNode(et, et.getRoot(), sonsThatHaveToBeSplit);
}

vector<tree::nodeP> findSonsThatHaveToBeSplit(const tree& et,
											  const split& mySplit,
											  const map<string, int> & nameIdMap){
// we assume that split is compatible with the tree and that the split is a subset of the children of the root.
// i.e., the node that has to be splitted is the root.
	vector<tree::nodeP> res;
	for (int i=0; i < et.getRoot()->getNumberOfSons(); ++i) {
		if (childIsInTheSplit(et.getRoot()->getSon(i),mySplit,nameIdMap)) {
			res.push_back(et.getRoot()->getSon(i));
		}
	}
	return res;
}

bool childIsInTheSplit(const tree::nodeP & myNode,
					   const split& mySplit,
					   const map<string, int> & nameIdMap) {
	if (myNode->isInternal()) return childIsInTheSplit(myNode->getSon(0),mySplit,nameIdMap);
	else {// we are in a leaf
		return (mySplit.isMember(idFromName(myNode->name(),nameIdMap)));
	}
}

