// 	$Id: replacementModelSSRV.cpp 4165 2008-06-04 09:19:48Z osnatz $	

#include "replacementModelSSRV.h"
#include "logFile.h"
#include <iomanip>
#include <iostream>


replacementModelSSRV::replacementModelSSRV(const distribution* dist, const replacementModel* baseRM, MDOUBLE rateOfRate /*= 1 */) :
_dist(dist->clone()),
_baseRM(baseRM->clone()),
_rateOfRate(rateOfRate)
{
	if (_dist->categories() == 0)
		errorMsg::reportError("replacementModelSSRV::replacementModelSSRV : number of categories == 0");
	
	updateFreq();
	updateQ();

	
}

//// similar to goldmanYangModel.cpp 
//replacementModelSSRV::replacementModelSSRV(const replacementModelSSRV& other) :
//_dist(other._dist->clone()),
//_baseRM(other._baseRM->clone()),
//_rateOfRate(other._rateOfRate)
//{
//	int size = alphabetSize();
//	_Q.resize(size);
//	for (int z=0; z < _Q.size();++z) 
//		_Q[z].resize(size,0);
//	updateFreq();
//	updateQ();
//}

// Instead of calling updateQ here, like in goldmanYangModel.cpp,
// this method uses the copy constructor of q2pt and also copies _freq and _Q
replacementModelSSRV::replacementModelSSRV(const replacementModelSSRV& other) :
_dist(other._dist->clone()),
_baseRM(other._baseRM->clone()),	
_rateOfRate(other._rateOfRate),
_q2pt(other._q2pt),
_freq(other._freq),
_Q(other._Q)
{
}

replacementModelSSRV::~replacementModelSSRV()
{
	if (_dist) delete (_dist);
	if (_baseRM) delete (_baseRM);
}


replacementModelSSRV& replacementModelSSRV::operator=(const replacementModelSSRV &other)
{
	if (_dist) delete (_dist);
	if (_baseRM) delete (_baseRM);
	
	_dist = other._dist->clone();
	_baseRM = other._baseRM->clone();	
	_rateOfRate = other._rateOfRate;
	_q2pt = other._q2pt; //@@@@ why doesn't this work ? explicit ?
//	_q2pt.fillFromRateMatrix(other._freq,other._Q);
	_freq = other._freq;
	_Q = other._Q;

	return (*this);
}

const int replacementModelSSRV::alphabetSize() const
{
	return (_baseRM->alphabetSize() * _dist->categories());
}



// The freq of each mulCharacter is its freq in the _baseRM * the freq of the rate-category
void replacementModelSSRV::updateFreq()	
{
	_freq.clear();
	int size = alphabetSize();
	int numCategories = _dist->categories();
	_freq.resize(size);
	int idInCategory;
		
	for(idInCategory=0; idInCategory < _baseRM->alphabetSize() ; ++idInCategory)
	{
		for (int categoryNumber=0; categoryNumber < numCategories; ++categoryNumber)
			_freq[categoryNumber*_baseRM->alphabetSize() + idInCategory] = 
				_baseRM->freq(idInCategory) * _dist->ratesProb(categoryNumber);
	}
}


void replacementModelSSRV::updateQ()
{
	if (_rateOfRate < EPSILON) _rateOfRate = EPSILON; // Temporary - to overcome a bug in QL algorithm, when _rateOfRate == 0

	_Q.clear();
	int size = alphabetSize();
	_Q.resize(size);
	for (int z=0; z < _Q.size();++z) 
		_Q[z].resize(size,0.0);
	
	// fill Q
	int _BaseRM_alphabetSize = _baseRM->alphabetSize();
	int numCategories = _dist->categories();
	// i,j : go over all the base-alphabet.
	// z,w : go over all the categories.
	for (int i=0; i < _BaseRM_alphabetSize; ++i)
	{
		for (int j=0; j < _BaseRM_alphabetSize; ++j)
		{
			for (int z=0; z < numCategories; ++z)
			{
				for (int w=0; w < numCategories; ++w)
				{
					if (i!=j)
					{
						// different alphabet, same rate category
						if (z==w)
							_Q[z*_BaseRM_alphabetSize + i][z*_BaseRM_alphabetSize+j] 
							= _dist->rates(z) * _baseRM->dPij_dt(i,j,0);
					}
					else
					{
						// same alphabet, different rate category
						if (z!=w)
						{
							_Q[z*_BaseRM_alphabetSize+i][w*_BaseRM_alphabetSize+i] =  _rateOfRate * _dist->ratesProb(w);
						}
						// same alphabet, same rate category
						else
							_Q[z*_BaseRM_alphabetSize+i][z*_BaseRM_alphabetSize+i] = 
								_dist->rates(z) * _baseRM->dPij_dt(i,j,0) 
								- ( _rateOfRate * (1.0 - _dist->ratesProb(z)));
					}

				}
			}
		}
	}
	
//	// check OZ
//	LOG(4, <<"THE Q MATRIX IS: "<<endl ) ;
//	VVdouble::iterator itr1 = _Q.begin();
//	Vdouble::iterator itr2;
//	for (; itr1 != _Q.end(); ++itr1)
//	{
//		for (itr2 = itr1->begin(); itr2 != itr1->end(); ++itr2)
//			LOG(4,<< setprecision(3) <<  setw(5) << *itr2 <<'\t');
//		LOG(4,<<endl);
//	}
//	LOG (4,<<endl);
////	 end of check

	_q2pt.fillFromRateMatrix(_freq,_Q); 
	
}

void replacementModelSSRV::setDistribution(const distribution* dist)
 {
	 if (dist->categories() == 0)
		errorMsg::reportError("replacementModelSSRV::setDistribution : number of categories == 0");
	 if (_dist) delete (_dist);
		_dist=dist->clone();
	updateQ();
 }

MDOUBLE replacementModelSSRV::sumPijQij() const{
	MDOUBLE sum=0.0;
	for (int i=0; i < _Q.size(); ++i) {
		sum -= _Q[i][i]*_freq[i];
	}
	return sum;
}


//void replacementModelSSRV::norm(MDOUBLE scale){
//	
//	for (int i=0; i < _Q.size(); ++i) {
//		for (int j=0; j < _Q.size(); ++j) {
//			_Q[i][j]*=scale;
//		}
//	}
//	
//	_q2pt.fillFromRateMatrix(_freq,_Q);
//}








