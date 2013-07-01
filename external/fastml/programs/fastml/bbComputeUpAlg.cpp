#include "bbComputeUpAlg.h"
#include "seqContainerTreeMap.h"

void BBfillComputeUp(const tree& et,
				   const sequenceContainer& sc,
				   const int pos,
				   const computePijHom& pi,
				   suffStatGlobalHomPos& ssc,
				   const vector<sequence>& ancS) {

	seqContainerTreeMap sctm(sc,et);

	ssc.allocatePlace(et.getNodesNum(),pi.alphabetSize());
	treeIterDownTopConst tIt(et);
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		int letter;
		if (mynode->isLeaf()) {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				const int seqID = sctm.seqIdOfNodeI(mynode->id());
				MDOUBLE val = sc.getAlphabet()->relations(sc[seqID][pos],letter);
				ssc.set(mynode->id(),letter,val);
			}
		}
		else {
			for(letter=0; letter<pi.alphabetSize();letter++) {
				if ((ancS[mynode->id()][pos]!=-2) && // if there is already assignments for this node
					(ancS[mynode->id()][pos]!=letter)) {
					ssc.set(mynode->id(),letter,0);
					continue;
				} // this if takes care of internal node assignments...


				doubleRep total_prob=1.0;
				for(int i=0; i < mynode->getNumberOfSons();++i){				
					doubleRep prob=0.0;
					for(int letInSon=0; letInSon<pi.alphabetSize();letInSon++) {
						prob += ssc.get(mynode->getSon(i)->id(), letInSon)*
							pi.getPij(mynode->getSon(i)->id(),letter,letInSon);
					}
				total_prob*=prob;
				}
				ssc.set(mynode->id(),letter,total_prob);
			}
		}
	}
}
