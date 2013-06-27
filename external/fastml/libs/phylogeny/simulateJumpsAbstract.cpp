#include "simulateJumpsAbstract.h"


simulateJumpsAbstract::simulateJumpsAbstract(const tree& inTree, const stochasticProcess& sp, const int alphabetSize)
: _tree(inTree), _sp(sp), _alphabetSize(alphabetSize)	
{
}



//runSimulation: do the actual simulation. iterNum specifies the number of iterations starting from each state
void simulateJumpsAbstract::runSimulation(int iterNum)
{
	init();
	for (int state = 0; state < _alphabetSize; ++state)
	{
		for (int iter = 0; iter < iterNum; ++iter)
		{
			runOneIter(state);
		}
	}
	
	computeExpectationsAndPosterior();	
}

//////////////////////////////////////////////////////////
//combined two characters into a combined state.
//For example. if the alphabet is {0,1,2} then the combined alphabet will be {0,1...8}.
//The states (terminalStart, terminalEnd) = (0,2) then combinedId = 2.
//The states (terminalStart, terminalEnd) = (1,2) then combinedId = 5. etc.
int simulateJumpsAbstract::getCombinedState(int terminalStart, int terminalEnd) const
{
	return (terminalStart * _alphabetSize + terminalEnd);
}
int simulateJumpsAbstract::getStartId(int combinedState) const
{
	return combinedState / _alphabetSize;
}
int simulateJumpsAbstract::getEndId(int combinedState) const
{
	return combinedState % _alphabetSize;
}
//////////////////////////////////////////////////////////

