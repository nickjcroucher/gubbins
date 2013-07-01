// $Id: distribution.h 2709 2007-11-19 14:49:21Z itaymay $

// version 2.00
// last modified 21 Mar 2004

/************************************************************
 This is a virtual class from which all types of distribution classes inherit from.
************************************************************/

#ifndef ___DISTRIBUTION
#define ___DISTRIBUTION

#include "definitions.h"

class distribution {
public:
	virtual distribution* clone() const = 0;
	virtual ~distribution() = 0;

	virtual const int categories() const=0;  // @@@@ there is no need to return a const int.
	virtual void change_number_of_categories(int in_number_of_categories);
	virtual const MDOUBLE rates(const int i) const=0; // @@@@ there is no need to return a const MDOUBLE.
	virtual const MDOUBLE ratesProb(const int i) const=0; // @@@@ there is no need to return a const MDOUBLE.
 	virtual void setGlobalRate(const MDOUBLE x)=0;
 	virtual MDOUBLE getGlobalRate()const=0; // @@@@ there is no need to return a const MDOUBLE.
	virtual const MDOUBLE getCumulativeProb(const MDOUBLE x) const = 0; // @@@@ there is no need to return a const MDOUBLE.

};
#endif


