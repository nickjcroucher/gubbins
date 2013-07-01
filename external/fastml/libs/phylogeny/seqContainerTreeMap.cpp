// $Id: seqContainerTreeMap.cpp 5106 2008-10-31 02:17:49Z itaymay $

#include "seqContainerTreeMap.h"
#include "logFile.h"


//if bLeavesOnly == true then checks only leaves, otherwise the sequence container includes also internal nodes (as may be the result of simlations
void checkThatNamesInTreeAreSameAsNamesInSequenceContainer(const tree& et,const sequenceContainer & sc, bool bLeavesOnly){
	treeIterDownTopConst tIt(et);
	//cout<<"tree names:"<<endl;	
	for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
		bool bFound = false;
		if (bLeavesOnly) {
            if (mynode->isInternal()) 
                continue;
		}
		sequenceContainer::constTaxaIterator it=sc.constTaxaBegin();
		for (;it != sc.constTaxaEnd(); ++it) 
		{
			string scName = it->name();
			string treeNodeName = mynode->name();

			if (it->name() == mynode->name()) 
			{
				bFound = true;
				break;
			}
		}
		if (bFound == false) 
		{
			string errMsg = "The sequence name: ";
			errMsg += mynode->name();
			errMsg += " was found in the tree file but not found in the sequence file.\n";
			LOG(4,<<errMsg<<endl);
			errorMsg::reportError(errMsg);
		}
	}
	
	sequenceContainer::constTaxaIterator it=sc.constTaxaBegin();
	for (;it != sc.constTaxaEnd(); ++it){
		bool bFound = false;	
		for (tree::nodeP mynode = tIt.first(); mynode != tIt.end(); mynode = tIt.next()) {
			if (bLeavesOnly)
			{
                if (mynode->isInternal()) 
                    continue;
			}
			if (it->name() == mynode->name()) 
			{
				bFound = true;
				break;
			}
		}
		if (bFound == false) 
		{
			string errMsg = "The sequence name: ";
			errMsg += it->name();
			errMsg += " was found in the sequence file but not found in the tree file.\n";
			errorMsg::reportError(errMsg);
		}
	}
}

