// $Id: computePijComponent.h 3064 2007-12-23 10:23:57Z cohenofi $

#ifndef ___COMPUTE_PIJ_COMPONENT
#define ___COMPUTE_PIJ_COMPONENT

#include "definitions.h"
#include "tree.h"
#include "stochasticProcess.h"

class computePijHomSpec {//specific node, no rate variation
public:
	virtual ~computePijHomSpec(){};
	void fillPij(const MDOUBLE dis, const stochasticProcess& sp, int derivationOrder = 0, bool isReversible =true);
	void resize(const int alphabetSize) {
		_V.resize(alphabetSize);
		for (int z=0;z<alphabetSize;++z) _V[z].resize(alphabetSize);
	}

	int alphabetSize() const {return _V.size();}
	MDOUBLE getPij(const int let1,const int let2)const{
		return _V[let1][let2];
	}
	VVdouble _V; // let, let
};

class computePijHom {//all nodes, no rate variation
public:
	virtual ~computePijHom(){};
	void fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder = 0, bool isReversible =true);
	int alphabetSize() const {return _V[0].alphabetSize();}
	int getNodesNum() const {return _V.size();}
	MDOUBLE getPij(const int nodeId,const int let1,const int let2)const{
		return _V[nodeId].getPij(let1,let2);
	}
	vector<computePijHomSpec> _V; // let, let
};

class computePijGam {//
public:
	virtual ~computePijGam(){};
	void fillPij(const tree& et, const stochasticProcess& sp, int derivationOrder = 0, bool isReversible =true);
	int categories() const {return _V.size();}
	int alphabetSize() const {return _V[0].alphabetSize();}
	int getNodesNum() const {return _V[0].getNodesNum();}

	MDOUBLE getPij(const int rateCategor,const int nodeId,const int let1,const int let2)const{
		return _V[rateCategor].getPij(nodeId,let1,let2);
	}
	computePijHom& operator[] (int i) {return _V[i];}
	const computePijHom& operator[] (int i) const {return _V[i];}
	vector<computePijHom> _V; // each rate category
};

#endif
