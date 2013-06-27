// $Id: nj.h 962 2006-11-07 15:13:34Z privmane $

// version 1.00
// last modified 3 Nov 2002

#ifndef ___NJ
#define ___NJ
#include "definitions.h"
#include "tree.h"
#include "sequenceContainer.h"
#include "njConstrain.h"
#include "distances2Tree.h"
using namespace std;

class NJalg : public distances2Tree {
public:
  virtual   NJalg* clone() const {return new NJalg(*this);}
  // changed from computeNJtree to computeTree for competability to "distances2Tree"
	virtual tree computeTree(VVdouble distances, const vector<string>& names, const tree * const constriantTree = NULL);
	tree startingTree(const vector<string>& names);
	tree startingTree(const tree& inTree);
	void NJiterate(tree& et,vector<tree::nodeP>& currentNodes,
							 VVdouble& distanceTable);
	void NJiterate(tree& et,vector<tree::nodeP>& currentNodes,
		       VVdouble& distanceTable, njConstraint& njc);
	void calc_M_matrix(vector<tree::nodeP>& currentNodes,
						  const VVdouble& distanceTable,
						  const Vdouble & r_values,
						  int& minRaw,int& minCol);
	void calc_M_matrix(vector<tree::nodeP>& currentNodes,
						  const VVdouble& distanceTable,
						  const Vdouble & r_values,
			   int& minRaw,int& minCol, const njConstraint& njc);
	Vdouble calc_r_values(vector<tree::nodeP>& currentNodes,const VVdouble& distanceTable);
	tree::nodeP SeparateNodes(tree& et,tree::nodeP node1,tree::nodeP node2);
	void update3taxaLevel(VVdouble& distanceTable,Vdouble & r_values,vector<tree::nodeP>& currentNodes);
	void updateBranchDistance(const VVdouble& disT,
								const Vdouble& rValues,
								tree::nodeP nodeNew,
							  tree::nodeP nodeI,
							  tree::nodeP nodeJ,
							  int Iplace, int Jplace);

	void UpdateDistanceTableAndCurrentNodes(vector<tree::nodeP>& currentNodes,
											VVdouble& distanceTable,
											tree::nodeP nodeI,
											tree::nodeP nodeJ,
											tree::nodeP theNewNode,
											int Iplace, int Jplace);

};

/*
	//explicit NJalg(const tree& inTree, const computeDistance* cd);
	explicit NJalg();
	tree getNJtree() const {return *_myET;}// return a copy...
	void computeTree(const sequenceContainer& sd,const computeDistance* cd,const vector<MDOUBLE> * weights = NULL);
	VVdouble getDistanceTable(vector<string>& names) {
		names.erase(names.begin(),names.end());
		names = _nodeNames;
		return _startingDistanceTable;}
	VVdouble getLTable(vector<string>& names) {
		names.erase(names.begin(),names.end());
		names = _nodeNames;
		return LTable;}
private:
	//void starTreeFromInputTree(const tree& inTree);
	void starTreeFromInputsequenceContainer(const sequenceContainer& sd);
	void GetDisTable(const sequenceContainer& sd,const vector<MDOUBLE>  * weights);
	MDOUBLE dis(const int i, const int j) const{
		return (i<j) ? distanceTable[i][j] : distanceTable[j][i];
	}
	void findMinM(int& minRaw,int& minCol);


	tree* _myET;
	VVdouble distanceTable;
	VVdouble Mmatrix;
	Vdouble r_values;
	vector<tree::nodeP> currentNodes;
	const computeDistance* _cd;

	VVdouble _startingDistanceTable; // for printing etc... not used by the algorithm.
	vector<string> _nodeNames; // for printing etc... not used by the algorithm.
	VVdouble LTable;// for printing etc... not used by the algorithm.

*/
#endif


