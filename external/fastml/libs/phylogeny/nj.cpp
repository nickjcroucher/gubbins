// $Id: nj.cpp 962 2006-11-07 15:13:34Z privmane $

// version 1.00
// last modified 3 Nov 2002

#include "nj.h"
#include "errorMsg.h"
#include "logFile.h"
#include "treeUtil.h"
#include <cassert>
#include <algorithm>
#include <map>
using namespace std;


//------------------------------------------
// general outline:
// we follow Swofford's book, "Molecular Systematics" pg489.
// currentNodes is the vector of the nodes that are "in process".
// in the beggining, these are all the leaves. Once, 2 leaves are separeted, 
// they are excluded from currentNodes, and their father is added to currentNodes.
// we (almost) finish the algorithm when currentNodes's size is 3. (i.e., we know the topology).
// thus when we start from an evolutionary tree, all we do, is to construct a star (start) tree
//------------------------------------------




//------------------------------------------
// constructor and start
//------------------------------------------
tree NJalg::computeTree(VVdouble distances,const vector<string>& names, const tree * const constriantTree /*= NULL*/){
	assert(distances.size() == names.size());
	tree resTree = startingTree(names);
	if (distances.size()<3) return resTree;
	vector<tree::nodeP> currentNodes;
	resTree.getAllLeaves(currentNodes,resTree.getRoot());
	if (constriantTree) {
		njConstraint njc(resTree, *constriantTree);
		while (currentNodes.size() >= 3) NJiterate(resTree,currentNodes,distances, njc);
	} else {
		while (currentNodes.size() >= 3) NJiterate(resTree,currentNodes,distances);
	}
	resTree.create_names_to_internal_nodes();
	LOGDO(5,resTree.output(myLog::LogFile()));
	return resTree;
}

tree NJalg::startingTree(const vector<string>& names) {
	return starTree(names);
}

tree NJalg::startingTree(const tree& inTree) {
	tree et;
	et.createRootNode();
	vector<tree::nodeP> allLeaves;
	inTree.getAllLeaves(allLeaves,inTree.getRoot());
	
	vector<string> names(allLeaves.size());
	for (int k = 0 ; k < allLeaves.size(); ++k)
	  names[k]=allLeaves[k]->name();
	  
	return startingTree(names);
}

void NJalg::updateBranchDistance(const VVdouble& distanceTable,
								 const Vdouble& rValues,
								 tree::nodeP nodeNew,
								 tree::nodeP nodeI,
								 tree::nodeP nodeJ,
								 int Iplace,
								 int Jplace) {
	MDOUBLE dis= (Iplace<Jplace) ? distanceTable[Iplace][Jplace] : distanceTable[Jplace][Iplace];
	MDOUBLE DisI_new = dis/2.0;
	MDOUBLE tmp = rValues[Iplace] - rValues[Jplace];
	tmp/= (  2.0*(distanceTable.size()-2)  );
	DisI_new = DisI_new+ tmp;
	MDOUBLE DisJ_new = dis - DisI_new;
	if (DisI_new<tree::SHORT_LENGTH_VALUE) DisI_new=tree::SHORT_LENGTH_VALUE; // no negative..
	if (DisJ_new<tree::SHORT_LENGTH_VALUE) DisJ_new=tree::SHORT_LENGTH_VALUE; // no negative..
	nodeI->setDisToFather(DisI_new);
	nodeJ->setDisToFather(DisJ_new);
}

void NJalg::NJiterate(tree& et,
					  vector<tree::nodeP>& currentNodes,
							 VVdouble& distanceTable) {
	Vdouble	rVector = calc_r_values(currentNodes,distanceTable);//CHECK2
	
	if (currentNodes.size() == 3) {
		update3taxaLevel(distanceTable,rVector,currentNodes);
		currentNodes.clear();
		return;
	}
	
	int minRaw,minCol;
	calc_M_matrix(currentNodes,distanceTable,rVector,minRaw,minCol);//CHECK3
	tree::nodeP nodeI = currentNodes[minRaw];
	tree::nodeP nodeJ = currentNodes[minCol];
	tree::nodeP theNewNode;
	theNewNode= SeparateNodes(et,nodeI,nodeJ);
	//CHECK4
	updateBranchDistance(distanceTable,rVector,theNewNode,nodeI,nodeJ,minRaw,minCol);
	//CHECK6
		et.create_names_to_internal_nodes();
	UpdateDistanceTableAndCurrentNodes(currentNodes,distanceTable,nodeI,nodeJ,theNewNode,minRaw,minCol);
}

void NJalg::NJiterate(tree& et,
		      vector<tree::nodeP>& currentNodes,
		      VVdouble& distanceTable,
		      njConstraint& njc) {
	Vdouble	rMatrix = calc_r_values(currentNodes,distanceTable);//CHECK2
	
	if (currentNodes.size() == 3) {
		update3taxaLevel(distanceTable,rMatrix,currentNodes);
		currentNodes.clear();
		return;
	}
	
	int minRaw,minCol;
	calc_M_matrix(currentNodes,distanceTable,rMatrix,minRaw,minCol, njc);//CHECK3
	tree::nodeP nodeI = currentNodes[minRaw];
	tree::nodeP nodeJ = currentNodes[minCol];
	tree::nodeP theNewNode;
	theNewNode= SeparateNodes(et,nodeI,nodeJ);
	njc.join(nodeI, nodeJ, theNewNode);
	//CHECK4
	updateBranchDistance(distanceTable,rMatrix,theNewNode,nodeI,nodeJ,minRaw,minCol);
	//CHECK6
		et.create_names_to_internal_nodes();
	UpdateDistanceTableAndCurrentNodes(currentNodes,distanceTable,nodeI,nodeJ,theNewNode,minRaw,minCol);
	LOGDO(15,et.output(myLog::LogFile(),tree::ANCESTORID));

}



Vdouble NJalg::calc_r_values(vector<tree::nodeP>& currentNodes,
							 const VVdouble& distanceTable) {
	Vdouble r_values(currentNodes.size(),0.0);
	for (int i=0; i <r_values.size();++i) {
		for (int j =0; j < r_values.size();++j) {
			MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
			r_values[i] += dis;
		}
	}
	return r_values;
}

void NJalg::calc_M_matrix(vector<tree::nodeP>& currentNodes,
						  const VVdouble& distanceTable,
						  const Vdouble & r_values,
						  int& minRaw,int& minCol){
	MDOUBLE min = VERYBIG;
	for (int i=0; i < currentNodes.size();++i){
		for (int j =i+1; j < currentNodes.size();++j) {
			MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
			MDOUBLE tmp = dis-(r_values[i]+r_values[j])/(currentNodes.size()-2);
			if (tmp<min) {minRaw = i;minCol=j;min=tmp;}
			
		}
	}
}

void NJalg::calc_M_matrix(vector<tree::nodeP>& currentNodes,
			  const VVdouble& distanceTable,
			  const Vdouble & r_values,
			  int& minRaw,int& minCol, 
			  const njConstraint& njc){
	MDOUBLE min = VERYBIG;
	MDOUBLE min_noc =  VERYBIG;
	int minRaw_noc=-1,minCol_noc=-1;
	for (int i=0; i < currentNodes.size();++i){
	  for (int j =i+1; j < currentNodes.size();++j) {
	    if (njc.isCompatible(currentNodes[i],currentNodes[j])) {
	      MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
	      MDOUBLE tmp = dis-(r_values[i]+r_values[j])/(currentNodes.size()-2);
	      if (tmp<min) {minRaw = i;minCol=j;min=tmp;}
	    }
	    LOGDO(10,{
		    MDOUBLE dis= (i<j) ? distanceTable[i][j] : distanceTable[j][i];
		    MDOUBLE tmp = dis-(r_values[i]+r_values[j])/(currentNodes.size()-2);
		    if (tmp<min_noc) {minRaw_noc = i;minCol_noc=j;min_noc=tmp;}
		  });
	      
	      }
	}
	LOGDO(10, {if (min_noc != min) {myLog::LogFile() 
		    << "NJ-constratin changes outcome " <<
		    currentNodes[minRaw_noc]->name()<<","<<currentNodes[minCol_noc]->name() <<"-> " <<
		    currentNodes[minRaw]    ->name()<<","<<currentNodes[minCol]    ->name()<< 
		    "  ("<<min-min_noc<<")"<<endl; 
		  njc.isCompatible(currentNodes[minRaw_noc], currentNodes[minCol_noc], true);
		  myLog::LogFile() << njc <<endl;
		}
	      });
}

tree::nodeP NJalg::SeparateNodes(tree& et, tree::nodeP node1,
											 tree::nodeP node2) {
	if (node1->father() != node2->father()) 
	 errorMsg::reportError(" error in function NJalg::SeparateNodes - nodes don't have the same father");

	tree::nodeP fatherNode = node1->father();

	tree::nodeP theNewNode = et.createNode(fatherNode,et.getNodesNum());
	node1->setFather(theNewNode);
	theNewNode->setSon(node1);
	node2->setFather(theNewNode);
	theNewNode->setSon(node2);

	// remove from son list of father node.
	fatherNode->removeSon(node1); 

	fatherNode->removeSon(node2); 
	return theNewNode;
}

void NJalg::update3taxaLevel(VVdouble& distanceTable,Vdouble & r_values,
							 vector<tree::nodeP>& currentNodes) {
	// update the distance of the 3 taxa that are left in the end, to the root.
	
	MDOUBLE dis0root = distanceTable[0][1]/2+0.5*(r_values[0]-r_values[1]);
	MDOUBLE dis1root = distanceTable[0][1]/2+0.5*(r_values[1]-r_values[0]);
	MDOUBLE dis2root = distanceTable[0][2]/2+0.5*(r_values[2]-r_values[0]);
	if (dis0root<tree::SHORT_LENGTH_VALUE) dis0root=tree::SHORT_LENGTH_VALUE; // no negative..
	if (dis1root<tree::SHORT_LENGTH_VALUE) dis1root=tree::SHORT_LENGTH_VALUE; // no negative..
	if (dis2root<tree::SHORT_LENGTH_VALUE) dis2root=tree::SHORT_LENGTH_VALUE; // no negative..
	currentNodes[0]->setDisToFather(dis0root);
	currentNodes[1]->setDisToFather(dis1root);
	currentNodes[2]->setDisToFather(dis2root);
}

void NJalg::UpdateDistanceTableAndCurrentNodes(vector<tree::nodeP>& currentNodes,
											   VVdouble& distanceTable,
											   tree::nodeP nodeI,
											   tree::nodeP nodeJ,
											   tree::nodeP theNewNode,
											   int Iplace,
											   int Jplace) {
	//	Iplace is the place of i in the "old" currentNodes vector
	int i,j;
	//	updating currentNodes
	vector<tree::nodeP> newCurrentNode= currentNodes;

	vector<tree::nodeP>::iterator vec_iter1=remove(
		newCurrentNode.begin(),newCurrentNode.end(),nodeI );
	newCurrentNode.erase(vec_iter1,newCurrentNode.end());

	vector<tree::nodeP>::iterator vec_iter2=remove(
	newCurrentNode.begin(),newCurrentNode.end(),nodeJ );
	newCurrentNode.erase(vec_iter2,newCurrentNode.end());
	
	newCurrentNode.push_back(theNewNode);

	map<tree::nodeP,int> nodeIntMap1;
	for (int z=0; z<currentNodes.size();++z) {
		nodeIntMap1.insert(map<tree::nodeP,int>::value_type(currentNodes[z],z));
	}

	VVdouble newDisTable;
	newDisTable.resize(newCurrentNode.size());
	for (int z1=0;z1<newDisTable.size();++z1) newDisTable[z1].resize(newCurrentNode.size(),0.0);

// updatine the table
	for (i=0; i < newCurrentNode.size(); i++) {
		for (j=i+1; j < newCurrentNode.size() ; j++) {
			if ((i!=newCurrentNode.size()-1) && (j!=newCurrentNode.size()-1)) {// both old nodes
				int oldI = nodeIntMap1[newCurrentNode[i]];
				int oldJ = nodeIntMap1[newCurrentNode[j]];
				MDOUBLE dis= (oldI<oldJ) ? distanceTable[oldI][oldJ] : distanceTable[oldJ][oldI];
				newDisTable[i][j] = dis;
			} //else if (i==newCurrentNode.size()-1) { // i is new
			//	newDisTable[i][j] = (dis(Iplace,NewOldPlaces[j])+dis(Jplace,NewOldPlaces[j])-dis(Iplace,Jplace))/2.0;
			//}
			else if (j==newCurrentNode.size()-1) { // j is new
				int oldI = Iplace;
				int oldJ = Jplace;
				int oldK = nodeIntMap1[newCurrentNode[i]];
				MDOUBLE disIK= (oldI<oldK) ? distanceTable[oldI][oldK] : distanceTable[oldK][oldI];
				MDOUBLE disIJ= (oldI<oldJ) ? distanceTable[oldI][oldJ] : distanceTable[oldJ][oldI];
				MDOUBLE disJK= (oldJ<oldK) ? distanceTable[oldJ][oldK] : distanceTable[oldK][oldJ];
				newDisTable[i][j] = 0.5*(disIK+disJK-disIJ); //EQ. 43 SWOFFORD PAGE 489.
			}
		}
	}

	currentNodes=newCurrentNode;
	distanceTable=newDisTable;
}

/*
NJalg::NJalg(){
	_myET = NULL;
}



//-----------------------------
// The algorithm
//-----------------------------

void NJalg::GetDisTable(const sequenceContainer& sd,const vector<MDOUBLE>  * weights) {
	
	VVresize(_startingDistanceTable,distanceTable.size(),distanceTable.size());// for printing stuff later.
	VVresize(LTable,distanceTable.size(),distanceTable.size());// for printing stuff later.

	int i,j;
	_nodeNames.resize(currentNodes.size());
	for ( i=0; i < currentNodes.size(); i++) {
		_nodeNames[i] =(currentNodes[i]->name()); 
		for ( j=i+1; j < currentNodes.size(); j++) {
			MDOUBLE tempDis = -2000.0;
			MDOUBLE resLikelihood;
			int seqnodeI_ID = sd.getId(currentNodes[i]->name());
			int seqnodeJ_ID = sd.getId(currentNodes[j]->name());
			const sequence& snodeI = *sd.getSeqPtr(seqnodeI_ID,true);
			const sequence& snodeJ = *sd.getSeqPtr(seqnodeJ_ID,true);
			tempDis = _cd->giveDistance(snodeI,snodeJ,weights,&resLikelihood);
			distanceTable[i][j] = tempDis;
			LTable[i][j] = resLikelihood;
		}
	}
	if (myLog::LogLevel()>4) {
		for (i=0; i < currentNodes.size(); i++) {
			for (j=i+1; j < currentNodes.size(); j++) {
				LOG(100,<<"nj distance ["<<i<<"]["<<j<<"] ="<<distanceTable[i][j]<<endl);
			}
		}
	}
	//if (myLog::LogLevel()>4) {
	//	for (i=0; i < currentNodes.size(); i++) {
	//		for (j=i+1; j < currentNodes.size(); j++) {
	//			LOG(4,<<"nj likelihood for distance["<<i<<"]["<<j<<"] ="<<LTable[i][j]<<endl);
	//		}
	//	}
	//}
	// for printing stuff later.
	for (int tmp1=0; tmp1<distanceTable.size();++tmp1)
	for (int tmp2=0; tmp2<distanceTable.size();++tmp2) 
	_startingDistanceTable[tmp1][tmp2] = distanceTable[tmp1][tmp2];
}






void NJalg::NJiterate() {
	getMmatrixFromDistanceTable();
	int minRaw,minCol;
	findMinM(minRaw,minCol);
	
	tree::nodeP nodeI = currentNodes[minRaw];
	tree::nodeP nodeJ = currentNodes[minCol];
	tree::nodeP theNewNode;
	theNewNode= SeparateNodes(nodeI,nodeJ);

	//CHECK4

	updateBranchDistance(theNewNode,nodeI,nodeJ,minRaw,minCol);
	//CHECK6

	UpdateDistanceTableAndCurrentNodes(nodeI,nodeJ,theNewNode,minRaw,minCol);
}

		
		












//CHECK1
//cout<<"\n-----------------------------------------------"<<endl;
//for (int h=0; h < currentNodes.size(); h++) cout<<currentNodes[h]->name()<<" = "<<h<<endl;

//CHECK2
//	for (int i =0; i < r_values.size();++i) cout<<"r["<<i<<"] = "<<r_values[i]<<endl;

//CHECK3
//	for (i =0; i < currentNodes.size();++i) 
//		for (int j =i+1; j <currentNodes.size();++j) 
//			cout<<"M["<<i<<"]["<<j<<"] = "<<Mmatrix[i][j]<<endl;

//CHECK4
//	string htuname = "HTU";
//	char k = 'a'+currentNodes.size();
//	htuname+=k;
//	theNewNode->SetName(htuname);		
	
//CHECK5
//_myET->getRoot()->SetName("RootOfStar");		
	
//CHECK6
//	et.output(cout,et.getRoot(),tree::ANCESTOR);

	



*/
