// $Id: uniformDistribution.cpp 2712 2007-11-19 14:50:12Z itaymay $

#include "uniformDistribution.h"


uniformDistribution::uniformDistribution(const int numOfCategories, MDOUBLE lowerBound, 
										 MDOUBLE upperBound) :distribution() {
	_globalRate=1.0;
	setUniformParameters(numOfCategories, lowerBound, upperBound);
}


//copy constructor
uniformDistribution::uniformDistribution(const uniformDistribution& other) : 
	_rates(other._rates),
	_ratesProb(other._ratesProb),
	_globalRate(other._globalRate),
	_interval(other._interval),
	_upperBound(other._upperBound),
	_lowerBound(other._lowerBound)
{
}



void uniformDistribution::setUniformParameters(const int number_of_categories, 
											   MDOUBLE lowerBound, MDOUBLE upperBound){
	_upperBound = upperBound;
	_lowerBound = lowerBound;
	
	_interval = ((upperBound - lowerBound) / (number_of_categories+0.0));
	_rates.clear();
	_rates.resize(number_of_categories);
	_ratesProb.erase(_ratesProb.begin(),_ratesProb.end());
	_ratesProb.resize(number_of_categories, 1.0/number_of_categories);
	//setting _rates[i] as the middle value of each category	
	for (int i = 0; i < number_of_categories; ++i) { 
		_rates[i] = _lowerBound + (_interval * (i + 0.5));
	}
}

//returns the ith border between categories
//getBorder(0) = _lowerBound, getBorder(categories()) = _upperBound
MDOUBLE uniformDistribution::getBorder(int i) const {
	return (i == categories()) ?  _upperBound : (_rates[i] - (_interval/2));
}

const MDOUBLE uniformDistribution::getCumulativeProb(const MDOUBLE x) const
{
	if (x<_lowerBound)
		return 0;
	else if (x>= _upperBound)
		return 1;
	else
		return ((x-_lowerBound) / (_upperBound - _lowerBound));
}

void uniformDistribution::change_number_of_categories(int in_number_of_categories)
{
	if (in_number_of_categories == categories())
		return;
	setUniformParameters(in_number_of_categories, _lowerBound, _upperBound);
}

