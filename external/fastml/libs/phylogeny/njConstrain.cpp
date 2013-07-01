// $Id: njConstrain.cpp 962 2006-11-07 15:13:34Z privmane $

#include "definitions.h"  
#include <cassert>
#include "njConstrain.h"
#include "logFile.h"



njConstraint::njConstraint(const tree& starttree, const tree& constraintTree):_cTree(constraintTree), _interTreeMap(){
  vector<tree::nodeP> currentNodes;
  starttree.getAllLeaves(currentNodes,starttree.getRoot());
  vector<tree::nodeP> constraintNodes;
  _cTree.getAllLeaves(constraintNodes,_cTree.getRoot());
  assert(currentNodes.size()==constraintNodes.size());
  
  map<string,tree::nodeP> name2Node;
  for (vector<tree::nodeP>::iterator vec_iter=constraintNodes.begin();vec_iter!=constraintNodes.end();++vec_iter){
    //    name2Node[test];//=*vec_iter;
    name2Node[(*vec_iter)->name()]=*vec_iter;
  }
    
  for (vector<tree::nodeP>::iterator vec_iter2=currentNodes.begin();vec_iter2!=currentNodes.end();++vec_iter2){
    assert(name2Node.find((*vec_iter2)->name()) !=  name2Node.end()); // cant find the taxa in the constratin tree!
    _interTreeMap[*vec_iter2]=name2Node[(*vec_iter2)->name()];
  }
}


bool njConstraint::isCompatible(const tree::nodeP& n1, const tree::nodeP& n2, const bool verbose) const
{
  bool compatible;
  assert( _interTreeMap.find(n1) !=  _interTreeMap.end()); // cant find the taxa in the map!
  assert( _interTreeMap.find(n2) !=  _interTreeMap.end()); // cant find the taxa in the map!
 
  tree::nodeP s1=_interTreeMap.find(n1)->second;
  tree::nodeP s2=_interTreeMap.find(n2)->second;
  
  if (s1==_cTree.getRoot()) {	// we are asking undirected questions from a directed tree
    compatible =  (s2 != _cTree.getRoot()) && (s2->father() != _cTree.getRoot()) && (s2->father()->father() == _cTree.getRoot());
    if (verbose) LOG(11,<<"isCompatible - s1 is root"<<endl);
  } else if (s2==_cTree.getRoot()) {	// we are asking undirected questions from a directed tree
    compatible =  (s1 != _cTree.getRoot()) && (s1->father() != _cTree.getRoot()) && (s1->father()->father() == _cTree.getRoot());
    if (verbose) LOG(11,<<"isCompatible - s2 is root"<<endl);
  } else { 
    compatible = (s1->father()==s2->father());
  }

  if (verbose) LOG(11,<<"isCompatible:" <<s1->name()<<" + "<<s2->name()<<"-->"	  <<compatible<< endl);
 return (compatible);
}

tree::nodeP joinNodesToSubtree(tree& t,tree::nodeP& s1, tree::nodeP& s2)
{
  assert (s1->father()==s2->father()); // we can only do this if both nodes have same father

  LOG(10,<<endl<<s1->name()<<" and "<<s2->name()<<endl);

  tree::nodeP fatherNode=s1->father();

  if (fatherNode->getNumberOfSons()==2) {
    //    fatherNode->sons.clear();
    return (fatherNode); // no splitting needed
  }
  
  if (s1->father()==t.getRoot() && t.getRoot()->getNumberOfSons()==3) { // no split needed, but the root needs to change
  
    LOG(10,<<"************************* spacial case of constratin join"<<endl);
    LOGDO(10,t.output(myLog::LogFile(),tree::ANCESTORID));
    LOG(10,<<endl<<s1->name()<<" and "<<s2->name()<<endl);
    LOG(10,<<endl<<s1->father()->name()<<" and father "<<s2->father()->name()<<endl);

    tree::nodeP newFatherNode = s1->father();
    for (int i=0; i<3; ++i)
      if (t.getRoot()->getSon(i)!= s1 && t.getRoot()->getSon(i)!= s2){
	t.rootAt(t.getRoot()->getSon(i));
	LOGDO(10,t.output(myLog::LogFile(),tree::ANCESTORID));
	LOG(10,<<endl<<endl);
	return (newFatherNode);	//  this is the new root;
      }
  }

  tree::nodeP newNode = t.createNode(fatherNode, t.getNodesNum());
  newNode->setSon(s1);
  newNode->setSon(s2);
  newNode->claimSons();
  
  
  int k = fatherNode->getNumberOfSons();
  fatherNode->removeSon(s1);
  fatherNode->removeSon(s2);
  assert (k=fatherNode->getNumberOfSons()+2); // both s1 and s2 should have been skiped
  //  fatherNode->sons.resize(k);
  
  t.updateNumberofNodesANDleaves();
  t.create_names_to_internal_nodes();
  return(newNode);
} 

void njConstraint::join(const tree::nodeP& n1, const tree::nodeP& n2, const tree::nodeP& newFather)
{
  assert(_interTreeMap.find(n1) !=  _interTreeMap.end()); // cant find the taxa in the map!
  assert(_interTreeMap.find(n2) !=  _interTreeMap.end()); // cant find the taxa in the map!
  assert(_interTreeMap.find(newFather) ==  _interTreeMap.end()); // should not find the new father in the map!
  assert(isCompatible(n1,n2));

  //  tree::nodeP origFather=_interTreeMap.find(n1)->father();

  // do tree things
  LOG(10,<<endl<<n1->name()<<" AND "<<n2->name()<<endl);
  tree::nodeP newNode=joinNodesToSubtree(_cTree, _interTreeMap[n1], _interTreeMap[n2]);
  
  
  _interTreeMap.erase(n1);
  _interTreeMap.erase(n2);
  _interTreeMap[newFather]=newNode;


  LOGDO(17,_cTree.output(myLog::LogFile()));

}
void njConstraint::output(ostream &out) const{
  _cTree.output(out,tree::ANCESTORID);
  out <<endl;
}

ostream &operator<<(ostream &out, const njConstraint &c){
  c.output(out);
  return(out);
}
