// 	$Id: C_evalParamUSSRV.h 1915 2007-04-04 15:56:24Z privmane $	
#ifndef ___C_EVAL_PARAM_USSRV
#define ___C_EVAL_PARAM_USSRV

#include "definitions.h"

#include "likelihoodComputation.h"
#include "likelihoodComputation2USSRV.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "gammaDistribution.h"
#include "tree.h"
#include "replacementModelSSRV.h"
#include "tamura92.h"
#include "stochasticProcessSSRV.h"
#include "ussrvModel.h"
#include "logFile.h"

// *********************
// *       USSRV       *
// *********************

class C_evalParamUSSRV {
public:
	C_evalParamUSSRV(const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble* weights = NULL)
    : _et(et),_sc(sc),_baseSc(baseSc),_pModel(pModel),_weights(weights){}

	MDOUBLE operator() (MDOUBLE param) ;
    virtual ~C_evalParamUSSRV(){}

protected:
	const tree& _et;
	const sequenceContainer& _sc;
	const sequenceContainer& _baseSc;
	ussrvModel* _pModel;
	const Vdouble * _weights;


protected:
	virtual void setParam(MDOUBLE param) = 0;
	virtual void print(MDOUBLE param,MDOUBLE res) =0;
};


class C_evalAlphaUSSRV : public C_evalParamUSSRV {
public:
  C_evalAlphaUSSRV(const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble *weights = NULL) 
				 : C_evalParamUSSRV(et,sc,baseSc,pModel,weights) 
				{}
  
protected:
	virtual void setParam(MDOUBLE alpha);
	virtual void print(MDOUBLE alpha,MDOUBLE res);
};



class C_evalNuUSSRV : public C_evalParamUSSRV{
public:
  C_evalNuUSSRV(	const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble * weights = NULL)
    : C_evalParamUSSRV(et,sc,baseSc,pModel,weights){}

protected:
	virtual void setParam(MDOUBLE Nu);
	virtual void print(MDOUBLE nu,MDOUBLE res);
};

class C_evalFUSSRV : public C_evalParamUSSRV{
public:
  C_evalFUSSRV(	const tree& et,
				const sequenceContainer& sc,
				const sequenceContainer& baseSc,
				ussrvModel* pModel,
				const Vdouble * weights = NULL)
    : C_evalParamUSSRV(et,sc,baseSc,pModel,weights){}

protected:
	virtual void setParam(MDOUBLE F);
	virtual void print(MDOUBLE f,MDOUBLE res);
};

// *********************
// *       SSRV        *
// *********************

class C_evalParamSSRV {
public:
	C_evalParamSSRV(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble* weights = NULL)
		: _et(et),_sc(sc),_ssrvSp(ssrvSp),_weights(weights){}

		MDOUBLE operator() (MDOUBLE param) ;
		virtual ~C_evalParamSSRV(){}

protected:
	const tree& _et;
	const sequenceContainer& _sc;
	stochasticProcessSSRV& _ssrvSp;
	const Vdouble * _weights;


protected:
	virtual void setParam(MDOUBLE param) = 0;
	virtual void print(MDOUBLE param,MDOUBLE res) =0;
};


class C_evalAlphaSSRV : public C_evalParamSSRV {
public:
	C_evalAlphaSSRV(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble *weights = NULL) 
		: C_evalParamSSRV(et,sc,ssrvSp,weights) 
	{}

protected:
	virtual void setParam(MDOUBLE alpha);
	virtual void print(MDOUBLE alpha,MDOUBLE res);
};



class C_evalNuSSRV : public C_evalParamSSRV{
public:
	C_evalNuSSRV(	const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights = NULL)
		: C_evalParamSSRV(et,sc,ssrvSp,weights){}

protected:
	virtual void setParam(MDOUBLE Nu);
	virtual void print(MDOUBLE nu,MDOUBLE res);
};

class C_evalTrTvSSRV : public C_evalParamSSRV{
public:
	C_evalTrTvSSRV(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights = NULL)
		: C_evalParamSSRV(et,sc,ssrvSp,weights){}

protected:
	virtual void setParam(MDOUBLE TrTv);
	virtual void print(MDOUBLE TrTv,MDOUBLE res);
};

class C_evalThetaSSRV : public C_evalParamSSRV{
public:
	C_evalThetaSSRV(const tree& et,
		const sequenceContainer& sc,
		stochasticProcessSSRV& ssrvSp,
		const Vdouble * weights = NULL)
		: C_evalParamSSRV(et,sc,ssrvSp,weights){}

protected:
	virtual void setParam(MDOUBLE Theta);
	virtual void print(MDOUBLE Theta,MDOUBLE res);
};

#endif
