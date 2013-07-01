// 	$Id: stochasticProcessSSRV.cpp 962 2006-11-07 15:13:34Z privmane $	

#include "stochasticProcessSSRV.h"
#include "replacementModelSSRV.h"

// it's important to call static_cast<replacementModelSSRV*>(_pijAccelerator->getReplacementModel())->updateQ(), after changing
// this returned pointer. (when changing alpha)
distribution* stochasticProcessSSRV::distr() const
{
	return ( static_cast<replacementModelSSRV*>(_pijAccelerator->getReplacementModel())->getDistribution() );
}


void stochasticProcessSSRV::setDistribution(const distribution* in_distr)
{
	static_cast<replacementModelSSRV*>(_pijAccelerator->getReplacementModel())->setDistribution(in_distr);
}


