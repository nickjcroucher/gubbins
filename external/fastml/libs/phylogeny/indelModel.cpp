// 	$Id: indelModel.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "indelModel.h"


void indelModel::setFreqX(const MDOUBLE freq_x)
{
	_freq[0] =freq_x ;
	_alpha = 1/(2*_freq[0]*_freq[1]) ;
}

void indelModel::setFreqG(const MDOUBLE freq_g)
{
	_freq[0] =freq_g ;
	_alpha = 1/(2*_freq[0]*_freq[1]) ;
}
