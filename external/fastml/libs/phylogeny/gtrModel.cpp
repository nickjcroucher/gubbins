#include "gtrModel.h"
#include "readDatMatrix.h" // for the normalizeQ function.
#include "matrixUtils.h"

gtrModel::gtrModel(const Vdouble& freq,
				      const MDOUBLE a2c,
					  const MDOUBLE a2g,
					  const MDOUBLE a2t,
					  const MDOUBLE c2g,
					  const MDOUBLE c2t,
					  const MDOUBLE g2t)
					  :_a2c(a2c),_a2g(a2g),_a2t(a2t),_c2g(c2g),_c2t(c2t),_g2t(g2t),_freq(freq)
	{
	_Q.resize(alphabetSize());
    for (int z=0; z < _Q.size();++z) _Q[z].resize(alphabetSize(),0.0);
	updateQ(a2c,a2g,a2t,c2g,c2t,g2t);
}
 

gtrModel& gtrModel::operator=(const gtrModel &other) 
{
	_Q = other._Q;
	_freq = other._freq;
	_q2pt = other._q2pt;
	_a2c = other._a2c;
	_a2g = other._a2g;
	_a2t = other._a2t;
	_c2g = other._c2g;
	_c2t = other._c2t;
	_g2t = other._g2t;
	return *this;
}

gtrModel::gtrModel(const gtrModel &other)
{
	_Q = other._Q;
	_freq = other._freq;
	_q2pt = other._q2pt;
	_a2c = other._a2c;
	_a2g = other._a2g;
	_a2t = other._a2t;
	_c2g = other._c2g;
	_c2t = other._c2t;
	_g2t = other._g2t;
}

void gtrModel::norm(const MDOUBLE scale)
{
	for (int i=0; i < _Q.size(); ++i) {
		for (int j=0; j < _Q.size(); ++j) {
			_Q[i][j] *= scale; 		
		}
	}
}

MDOUBLE gtrModel::sumPijQij(){
	MDOUBLE sum=0.0;
	for (int i=0; i < _Q.size(); ++i) {
		sum -= (_Q[i][i])*_freq[i];
	}
	return sum;
}

void gtrModel::updateQ(const MDOUBLE a2c,const MDOUBLE a2g,const MDOUBLE a2t,const MDOUBLE c2g,const MDOUBLE c2t,const MDOUBLE g2t) 
{
	_a2c = a2c;	
	_Q[a][c] = (_a2c);
	_Q[c][a] = (_freq[a]*_a2c/_freq[c]);
	_a2g = a2g;
	_Q[a][g] = (_a2g);
	_Q[g][a] = (_freq[a]*_a2g/_freq[g]);
	_a2t = a2t;
	_Q[a][t] = (_a2t);
	_Q[t][a] = (_freq[a]*_a2t/_freq[t]);
	_c2g = c2g;
	_Q[c][g] = (_c2g);
	_Q[g][c] = (_freq[c]*_c2g/_freq[g]);
	_c2t = c2t;
	_Q[c][t] = (_c2t);
	_Q[t][c] = (_freq[c]*_c2t/_freq[t]);
	_g2t = g2t;
	_Q[g][t] = (_g2t);
	_Q[t][g] = (_freq[g]*_g2t/_freq[t]);
	_Q[a][a] = -1.0*(_Q[a][c]+_Q[a][g]+_Q[a][t]);
	_Q[c][c] = -1.0*(_Q[c][a]+_Q[c][g]+_Q[c][t]);
	_Q[g][g] = -1.0*(_Q[g][a]+_Q[g][c]+_Q[g][t]);
	_Q[t][t] = -1.0*(_Q[t][a]+_Q[t][c]+_Q[t][g]);
	norm(1.0/sumPijQij());
	_q2pt.fillFromRateMatrix(_freq,_Q);
}

void gtrModel::set_a2c(const MDOUBLE a2c)
{
	_a2c = a2c;
	updateQ(_a2c,_a2g,_a2t,_c2g,_c2t,_g2t);
}

void gtrModel::set_a2g(const MDOUBLE a2g)
{
	_a2g = a2g;
	updateQ(_a2c,_a2g,_a2t,_c2g,_c2t,_g2t);
}

void gtrModel::set_a2t(const MDOUBLE a2t)
{
	_a2t = a2t;
	updateQ(_a2c,_a2g,_a2t,_c2g,_c2t,_g2t);
}

void gtrModel::set_c2g(const MDOUBLE c2g)
{
	_c2g = c2g;
	updateQ(_a2c,_a2g,_a2t,_c2g,_c2t,_g2t);
}

void gtrModel::set_c2t(const MDOUBLE c2t)
{
	_c2t = c2t;
	updateQ(_a2c,_a2g,_a2t,_c2g,_c2t,_g2t);
}

void gtrModel::set_g2t(const MDOUBLE g2t)
{
	_g2t = g2t;
	updateQ(_a2c,_a2g,_a2t,_c2g,_c2t,_g2t);
}

MDOUBLE gtrModel::get_a2c() const
{
	MDOUBLE result;
	if(_Q.size() < alphabetSize())
		errorMsg::reportError("Attempting to reach an uninitiallized Q matrix in gtrModel::get_a2c");
	else{
		if((_Q[a].size() < alphabetSize())||(_Q[c].size() < alphabetSize()))
			errorMsg::reportError("Attempting to reach an uninitiallzed Q matrix element in Model::get_a2c");
		else
			result = _a2c;
	}
	return result;
}

MDOUBLE gtrModel::get_a2g() const
{
	MDOUBLE result;
	if(_Q.size() < alphabetSize())
		errorMsg::reportError("Attempting to reach an uninitiallized Q matrix in gtrModel::get_a2g");
	else{
		if((_Q[a].size() < alphabetSize())||(_Q[g].size() < alphabetSize()))
			errorMsg::reportError("Attempting to reach an uninitiallzed Q matrix element in Model::get_a2g");
		else
			result = _a2g;
	}
	return result;
}

MDOUBLE gtrModel::get_a2t() const
{
	MDOUBLE result;
	if(_Q.size() < alphabetSize())
		errorMsg::reportError("Attempting to reach an uninitiallized Q matrix in gtrModel::get_a2t");
	else{
		if((_Q[a].size() < alphabetSize())||(_Q[t].size() < alphabetSize()))
			errorMsg::reportError("Attempting to reach an uninitiallzed Q matrix element in Model::get_a2t");
		else
			result = _a2t;
	}
	return result;
}

MDOUBLE gtrModel::get_c2g() const
{
	MDOUBLE result;
	if(_Q.size() < alphabetSize())
		errorMsg::reportError("Attempting to reach an uninitiallized Q matrix in gtrModel::get_c2g");
	else{
		if((_Q[c].size() < alphabetSize())||(_Q[g].size() < alphabetSize()))
			errorMsg::reportError("Attempting to reach an uninitiallzed Q matrix element in Model::get_c2g");
		else
			result = _c2g;
	}
	return result;
}

MDOUBLE gtrModel::get_c2t() const
{
	MDOUBLE result;
	if(_Q.size() < alphabetSize())
		errorMsg::reportError("Attempting to reach an uninitiallized Q matrix in gtrModel::get_c2t");
	else{
		if((_Q[c].size() < alphabetSize())||(_Q[t].size() < alphabetSize()))
			errorMsg::reportError("Attempting to reach an uninitiallzed Q matrix element in Model::get_c2t");
		else
			result = _c2t;
	}
	return result;
}

MDOUBLE gtrModel::get_g2t() const
{
	MDOUBLE result;
	if(_Q.size() < alphabetSize())
		errorMsg::reportError("Attempting to reach an uninitiallized Q matrix in gtrModel::get_g2t");
	else{
		if((_Q[g].size() < alphabetSize())||(_Q[t].size() < alphabetSize()))
			errorMsg::reportError("Attempting to reach an uninitiallzed Q matrix element in Model::get_g2t");
		else
			result = _g2t;
	}
	return result;
}
