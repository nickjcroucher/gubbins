// 	$Id: C_evalParamUSSRV.cpp 1915 2007-04-04 15:56:24Z privmane $	
#include "C_evalParamUSSRV.h"

// *********************
// *       USSRV       *
// *********************

MDOUBLE C_evalParamUSSRV::operator() (MDOUBLE param) {
		
		setParam(param);
		MDOUBLE res = likelihoodComputation2USSRV::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_baseSc,*_pModel,_weights);
		print(param,res);
		return -res;
}

void C_evalAlphaUSSRV::setParam(MDOUBLE alpha)
{
	if (_pModel->noOfCategor() == 1)
		errorMsg::reportError(" one category when trying to optimize alpha");
	_pModel->updateAlpha(alpha);	
}

void C_evalAlphaUSSRV::print(MDOUBLE alpha,MDOUBLE res) {
			LOG(5,<<" with Alpha = "<<alpha<<" logL = " <<res<<endl);
}


void C_evalNuUSSRV::setParam(MDOUBLE Nu)
{
	_pModel->updateNu(Nu);
}

void C_evalNuUSSRV::print(MDOUBLE nu,MDOUBLE res) {
	LOG(5,<<" with Nu = "<<nu<<" logL = " <<res<<endl);
}

void C_evalFUSSRV::setParam(MDOUBLE f)
{
	_pModel->updateF(f);
}

void C_evalFUSSRV::print(MDOUBLE f,MDOUBLE res) {
	LOG(5,<<" with F = "<<f<<" logL = " <<res<<endl);
}


// *********************
// *       SSRV       *
// *********************

MDOUBLE C_evalParamSSRV::operator() (MDOUBLE param) {

	setParam(param);
	MDOUBLE res = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_ssrvSp,_weights);
	print(param,res);
	return -res;
}

void C_evalAlphaSSRV::setParam(MDOUBLE alpha)
{
	if (alpha<0) 
		errorMsg::reportError("ERROR in C_evalAlphaSSRV::setParam, alpha is < 0 ");
		
	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(_ssrvSp.getPijAccelerator()->getReplacementModel());
	gammaDistribution* gammaDist = static_cast<gammaDistribution*>(pMulRM->getDistribution()); 
	gammaDist->setAlpha(alpha);
	pMulRM->updateQ();
}

void C_evalAlphaSSRV::print(MDOUBLE alpha,MDOUBLE res) {
	LOG(5,<<" with Alpha = "<<alpha<<" logL = " <<res<<endl);
}


void C_evalNuSSRV::setParam(MDOUBLE Nu)
{
	if (Nu<0) 
		errorMsg::reportError("C_evalNuSSRV::setParam, nu is < 0 ");
		
	static_cast<replacementModelSSRV*>(_ssrvSp.getPijAccelerator()->getReplacementModel())->setRateOfRate(Nu);
}

void C_evalNuSSRV::print(MDOUBLE nu,MDOUBLE res) {
	LOG(5,<<" with Nu = "<<nu<<" logL = " <<res<<endl);
}

void C_evalTrTvSSRV::setParam(MDOUBLE TrTv)
{
	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(_ssrvSp.getPijAccelerator()->getReplacementModel());
	static_cast<tamura92*>(pMulRM->getBaseRM())->changeTrTv(TrTv);
	pMulRM->updateQ();
}

void C_evalTrTvSSRV::print(MDOUBLE TrTv,MDOUBLE res) {
	LOG(5,<<" with TrTv = "<<TrTv<<" logL = " <<res<<endl);
}

void C_evalThetaSSRV::setParam(MDOUBLE Theta)
{
	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(_ssrvSp.getPijAccelerator()->getReplacementModel());
	static_cast<tamura92*>(pMulRM->getBaseRM())->changeTheta(Theta);
	pMulRM->updateFreq();
	pMulRM->updateQ();
}

void C_evalThetaSSRV::print(MDOUBLE Theta,MDOUBLE res) {
	LOG(5,<<" with Theta = "<<Theta<<" logL = " <<res<<endl);
}




