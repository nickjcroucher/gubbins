// $Id: tree.h 3996 2008-05-13 13:06:16Z adido $

#ifndef ___TREE
#define ___TREE

#include "definitions.h"
#include "readTree.h"
#include "errorMsg.h"
#include "logFile.h"


//***********************************************************************************
// class tree represents only the topology. It has no MSA and assumes no model of evolution.
//***********************************************************************************


class tree {
public:
	static const MDOUBLE FLAT_LENGTH_VALUE;// = 0.3;
	static const int TREE_NULL;// = -1;
	static const MDOUBLE SHORT_LENGTH_VALUE;// = 0.000001f;

//---------------------------- TREE NODE ----------------------
public:
	class TreeNode {
	public:
	  explicit TreeNode(const int id) :_sons(0),_father(NULL),_id(id),_name( (string)"" ),_dis2father(TREE_NULL),_comment((string)"") {}
		const int id() const {return _id;}
		const string name() const {return _name;}
		const MDOUBLE dis2father() const {return _dis2father;}
		TreeNode* father() {return _father;}
		void setName(const string &inS) {_name = inS;}
		void setID(const int inID) {_id = inID;}
		void setDisToFather(const MDOUBLE dis) {_dis2father = dis;}
		void setFather(TreeNode* tn){_father=tn;}
		int getNumberOfSons() const {return _sons.size();}
		TreeNode* getSon (int i) {return _sons[i];}
		TreeNode* getLastSon () {return _sons.back();}
		void removeLastSon() {_sons.pop_back();}
		void removeSon(TreeNode* pSon);
		//setSon: updates only the father pointer to the son!
		void setSon(TreeNode* pSon) {_sons.push_back(pSon);}
		void setSon(TreeNode* pSon, int i) {_sons[i]=pSon;} // this will overwrite previous pointer!
		bool isRoot() const {return (_father == NULL);}
		bool isLeaf() const {
			return (
				(getNumberOfSons() ==0) || 
				(isRoot() && (getNumberOfSons() ==1))
				) ;
		}
		bool isInternal() const {return (!isLeaf());}
		//claimSons: sets the _father pointer of all sons to (this)
		//this function is used after setSon has been called without updating the son pointer.
		void claimSons();
		void removeAllSons() {_sons.clear();}
		void copySons(TreeNode* other) {//copy the vector of nodeP only from one node to the other 
			_sons=other->_sons;
		}
	  void setComment(string comment) {_comment = comment;
		if (comment.length())
		  LOG(16,<<"comment for "<<_name<<" set to "<<comment<<endl );}
	  const string getComment(void) const {return _comment;}
	private:
		vector<TreeNode*> _sons;
		TreeNode* _father;
		int _id;
		string _name;
		MDOUBLE _dis2father;
	    string _comment;
		friend class tree;
	};
//------------------------------------------------------------


public:
	//NEWICK is the standard format
	//ANCESTOR/ANCESTORID are for debugging purposes: output a list of nodes one for each line.
	//for each node print the name, dist2father and its sons. id are printed only in ANCESTORID.
	//PAML is like Newick format but with extra line: #of leaves space and #of trees
	typedef enum { PHYLIP, ANCESTOR, ANCESTORID, PAML } TREEformats;
	typedef TreeNode* nodeP;

public:
//*******************************************************************************
// constructors
//*******************************************************************************
	tree();
	tree(const string& treeFileName); 
	tree(istream &treeFile);
	tree(const vector<char>& tree_contents);

	tree(const string& treeFileName,vector<char>& isFixed); 
	tree(const vector<char>& tree_contents, vector<char>& isFixed);
	tree(istream &in, vector<char>& isFixed);

	tree(const tree &otherTree);
	tree& operator=(const tree &otherTree);

	virtual ~tree() {clear();};
	
//*******************************************************************************
// questions on the tree topology
//*******************************************************************************

	nodeP getRoot() const {return _root;};
	inline int getLeavesNum() const; 
	inline int getNodesNum() const;
	inline int getInternalNodesNum() const;
	//findNodeByName: searches the subtree of myNode for a node with a specified name.
	//if myNode==NULL: the search starts from the root
	nodeP findNodeByName(const string inName, nodeP myNode=NULL) const; 
	nodeP findNodeById(const int inId, nodeP myNode=NULL) const;
	bool withBranchLength() const;
	//getNeigboursOfNode: stores into neighbourVec the father and sons of myNode
	void getNeigboursOfNode(vector<nodeP> &neighbourVec, const nodeP myNode) const;
	void getTreeDistanceTableAndNames(VVdouble& disTab, vector <string>& vNames) const;
	MDOUBLE findLengthBetweenAnyTwoNodes(const nodeP node1,const nodeP node2) const;
	//lengthBetweenNodes: find length between two neighbouring nodes only
	MDOUBLE lengthBetweenNodes(const nodeP i, const nodeP j) const; 

	void getPathBetweenAnyTwoNodes(vector<nodeP> &path,const nodeP node1, const nodeP node2) const;
	void getFromLeavesToRoot(vector<nodeP> &vNeighbourVector) const;
	void getFromRootToLeaves(vector<nodeP> &vec) const;
	void getFromNodeToLeaves(vector<nodeP> &vec, const nodeP fromHereDown) const;

	void getAllHTUs(vector<nodeP> &vec,const nodeP fromHereDown) const ;
	void getAllNodes(vector<nodeP> &vec,const nodeP fromHereDown) const ;
	void getAllLeaves(vector<nodeP> &vec,const nodeP fromHereDown) const;

//*******************************************************************************
// change tree topoplogy parameters - should be applied carefully
//*******************************************************************************
	//rootAt: sets newRoot as the root. updates the iterator order lists.
	void rootAt(const nodeP newRoot);   
	void rootToUnrootedTree();
	void multipleAllBranchesByFactor(const MDOUBLE InFactor);
	void create_names_to_internal_nodes();
	void makeSureAllBranchesArePositive();
	void makeSureAllBranchesAreLargerThanEpsilon(MDOUBLE epsilon);

	// removeNodeFromSonListOfItsFather:
	// removes sonNode from its father according to the name of sonNode
	// this function should ONLY be used when sonNode is to be recycled soon!
	// because this function does not change the number of leaves nor the number of nodes!
	// nor does it change the father of sonNode.
	void removeNodeFromSonListOfItsFather(nodeP sonNode);

	void shrinkNode(nodeP nodePTR); 
	//removeLeaf: removes nodePTR from tree. also deletes nodePTR
	void removeLeaf(nodeP nodePTR); 
	//getAllBranches: returns two vectors such that nodesUp[i] is the father of nodesDown[i]
	void getAllBranches(vector<nodeP> &nodesUP, vector<nodeP> & nodesDown);
	//createRootNode: erase the current tree and create a tree with one node. 
	void createRootNode();
	nodeP createNode(nodeP fatherNode, const int id);
	void updateNumberofNodesANDleaves();

// **********************************************************
//  initialization
// **********************************************************
	
	//createFlatLengthMatrix: sets the distance of all branches to newFlatDistance
	void createFlatLengthMatrix(const MDOUBLE newFlatDistance = FLAT_LENGTH_VALUE);  
	//recursiveBuildTree: copy the information from other_nodePTR to a new node, and set the father to father_nodePTR
	//used by treeUtil
	nodeP recursiveBuildTree(tree::nodeP father_nodePTR,const tree::nodeP other_nodePTR);

//*******************************************************************************
// Input-Output
//*******************************************************************************
	void output(string treeOutFile, TREEformats fmt= PHYLIP,bool withHTU=false) const;
	void output(ostream& os, TREEformats fmt= PHYLIP,bool withHTU=false) const;

private:
	void clear();

	void outputInAncestorTreeFormat(ostream& treeOutStream, bool withDist = false) const;
	void outputInPhylipTreeFormat(ostream& treeOutStream,bool withHTU=false) const;
	void outputInAncestorIdTreeFormat(ostream& treeOutStream, bool withDist = false) const;
	void outputInPamlTreeFormat(ostream& treeOutStream, bool withHTU = false) const;
	int print_from(nodeP from_node, ostream& os, bool withHTU) const;
	int print_from(nodeP from_node, ostream& os, bool withHTU);

	bool readPhylipTreeTopology(istream& in,vector<char>& isFixed); //same as the constructor with file name
	bool readPhylipTreeTopology(const vector<char>& tree_contents,vector<char>& isFixed);
	bool readPhylipTreeTopology(istream& in); //same as the constructor with file name
	bool readPhylipTreeTopology(const vector<char>& tree_contents);
	nodeP readPart(vector<char>::const_iterator& p_itCurrent, int& nextFreeID, vector<char> & isFixed);

	void getAllHTUsPrivate(vector<nodeP> &vec,nodeP fromHereDown) const ;
	void getAllNodesPrivate(vector<nodeP> &vec,nodeP fromHereDown) const ;
	void getAllLeavesPrivate(vector<nodeP> &vec,nodeP fromHereDown) const;


protected:
	TreeNode *_root;
	int _leaves;
	int _nodes;
};

inline int tree::getLeavesNum() const {return _leaves;}
inline int tree::getNodesNum() const {return _nodes;}
inline int tree::getInternalNodesNum() const {return getNodesNum() - getLeavesNum();}

ostream &operator<<(ostream &out, const tree &tr);

#endif 

