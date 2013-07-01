#include "bbComputeDownAlg.h"
#include "seqContainerTreeMap.h"

void BBfillComputeDown(const tree& et,
					   const sequenceContainer& sc,
					   const int pos,
					   const computePijHom& pi,
					   suffStatGlobalHomPos& ssc,
					   const suffStatGlobalHomPos& cup,
					   const vector<sequence>& ancS){
	ssc.allocatePlace(et.getNodesNum(), pi.alphabetSize());
	treeIterTopDownConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter,letterInFather,bro,letterInSon;
		if (mynode->father()==NULL) {// if root
			for(letter=0; letter<pi.alphabetSize();letter++) {
				ssc.set(mynode->id(),letter,1.0);
			}
			mynode = tIt.next(); //continue
		}
		tree::nodeP fatherNode=mynode->father();
		const int n_bro=fatherNode->getNumberOfSons();
		for(letter=0; letter<pi.alphabetSize();letter++) {
			if ((ancS[mynode->father()->id()][pos]!=-2)&&(ancS[mynode->father()->id()][pos]!=letter)){
				ssc.set(mynode->id(),letter,0);
				continue;
			} // this if takes care of internal node assignments...

			doubleRep totalProb=1.0;
			doubleRep fatherTerm=0;
			if (fatherNode->father()!=NULL) {
				for(letterInFather=0; letterInFather<pi.alphabetSize();letterInFather++)
					fatherTerm += pi.getPij(fatherNode->id(),letter,letterInFather)*
					ssc.get(fatherNode->id(),letterInFather);
			}
			else {
				fatherTerm=1.0;
			}
			doubleRep brotherTerm=1.0;
			for(bro = 0; bro < n_bro; bro++) {
				tree::nodeP brother = fatherNode->getSon(bro);
				if (brother != mynode) {
					doubleRep tmp_bro=0.0;
					for(letterInSon=0; letterInSon<pi.alphabetSize();letterInSon++) {
						tmp_bro+=pi.getPij(fatherNode->getSon(bro)->id(),letter,letterInSon)*
						cup.get(brother->id(),letterInSon);
					}
					brotherTerm *=tmp_bro;
				}
			}
			totalProb = fatherTerm * brotherTerm;
			ssc.set(mynode->id(),letter,totalProb);
		}
	}
}
/*
const evolTree* bbComputeDownAlg::_et=NULL;
const stochasticProcess* bbComputeDownAlg::_sp=NULL;
const suffStatComponent* bbComputeDownAlg::_cup=NULL;
const computePij* bbComputeDownAlg::_cpij=NULL;
suffStatComponent* bbComputeDownAlg::_ssc=NULL;
const vector<sequence>* bbComputeDownAlg::_ancS = NULL;

void bbComputeDownAlg::bbFillComputeDown(const evolTree* et,
				  const stochasticProcess* sp,
				  const suffStatComponent* cup,
				  const computePij* cpij,
				  suffStatComponent* ssc,
				  vector<sequence>* ancS) {


	_et=et;_sp=sp;_cup=cup;_cpij=cpij, _ssc=ssc;_ancS=ancS;
	_ssc->resize(et->iNodes());
	if (_ssc->size()>0) 
	if ((*_ssc)[0].isEmpty()==true) {// alocating memory for the pij(t)...
		for (vector<suffStatComponent::suffStatComponentCell>::iterator it=ssc->_suffCellVec.begin();
				it !=ssc->_suffCellVec.end();++it) {
			it->allocatePlace(_et->seqLen(),
				_sp->categories(),_et->alphabetSize());
		}
	}
	recursiveFillDown(_et->iRoot());
}

void bbComputeDownAlg::bbFillComputeDownForOnePos(const evolTree* et,
				  const stochasticProcess* sp,
				  const suffStatComponent* cup,
				  const computePij* cpij,
				  suffStatComponent* ssc,
				  vector<sequence>* ancS,
				  const int pos) {


	_et=et;_sp=sp;_cup=cup;_cpij=cpij, _ssc=ssc;_ancS=ancS;
	_ssc->resize(et->iNodes());
	if (_ssc->size()>0) 
	if ((*_ssc)[0].isEmpty()==true) {// alocating memory for the pij(t)...
		for (vector<suffStatComponent::suffStatComponentCell>::iterator it=ssc->_suffCellVec.begin();
				it !=ssc->_suffCellVec.end();++it) {
			it->allocatePlace(_et->seqLen(),
				_sp->categories(),_et->alphabetSize());
		}
	}
	recursiveFillDownPos(_et->iRoot(),pos);
}

void bbComputeDownAlg::recursiveFillDownPos(const evolTree::NodeP& mynode,
											const int pos) {
	fillDownNodePos(mynode,pos);
	for (vector<evolTree::nodeP>::iterator i=mynode->sons.begin(); i != mynode->sons.end();++i) {
		recursiveFillDownPos(*i,pos);
	}
}

void bbComputeDownAlg::recursiveFillDown(const evolTree::NodeP& mynode) {
	fillDownNode(mynode);
	for (vector<evolTree::nodeP>::iterator i=mynode->sons.begin(); i != mynode->sons.end();++i) {
		recursiveFillDown(*i);
	}
}

void bbComputeDownAlg::fillDownNode(
				const evolTree::NodeP& mynode) {
	for(int pos=0; pos<_et->seqLen();pos++) fillDownNodePos(mynode,pos);			
}

void bbComputeDownAlg::fillDownNodePos(
				const evolTree::NodeP& mynode,
				const int pos) {

	int rateCategor,letter,letter_in_father,bro,letter_in_son;
	if (mynode->father==NULL) {// if root
		for (rateCategor = 0; rateCategor<_sp->categories(); ++rateCategor) {
			for(letter=0; letter<_et->alphabetSize();letter++) {
				(*_ssc)[mynode->id()].set(pos,rateCategor,letter,1.0);
			}
		}	
		return;
	}
	for (rateCategor = 0; rateCategor<_sp->categories(); ++rateCategor) {
		evolTree::NodeP father_node=mynode->father;
		const int n_bro=father_node->sons.size();
		for(letter=0; letter<_et->alphabetSize();letter++) {//alpha
			assert(_ancS != NULL);
			//------------------------------------------------------
			if (((*_ancS)[mynode->father->id()][pos]!=letter) && 
				((*_ancS)[mynode->father->id()][pos]!=-2)) 			{
				(*_ssc)[mynode->id()].set(pos,rateCategor,letter,0);
				continue;
			} // this if takes care of internal node assignments...
			//------------------------------------------------------
		
			MDOUBLE total_prob=1.0;
			MDOUBLE father_term=0;
			if (father_node->father!=NULL) {
				for(letter_in_father=0; letter_in_father<_et->alphabetSize();letter_in_father++)
					father_term += _cpij->getPij(father_node->id(),letter,letter_in_father,rateCategor)*
					(*_ssc)[father_node->id()].get(pos,rateCategor,letter_in_father);
			}
			else {
				father_term=1.0;
			}
				MDOUBLE brother_term=1.0;
			for(bro=0;bro<n_bro;bro++) {
				evolTree::NodeP brother=father_node->sons[bro];
				if (brother != mynode) {
					MDOUBLE tmp_bro=0.0;
					for(letter_in_son=0; letter_in_son<_et->alphabetSize();letter_in_son++) {
						tmp_bro+=_cpij->getPij(
							father_node->sons[bro]->id(),
							letter,
							letter_in_son,rateCategor)*
						_cup->get(brother->id(),
						pos,
						rateCategor,
						letter_in_son);
					}
					brother_term *=tmp_bro;
				}
			}
			total_prob = father_term * brother_term;
			(*_ssc)[mynode->id()].set(pos,rateCategor,letter,total_prob);
		}
	}
}
*/





