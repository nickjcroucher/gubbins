#include "computeJumps.h"
#include "talRandom.h"
#include "someUtil.h"
#include "matrixUtils.h"
#include <algorithm>


computeJumps::computeJumps(const MDOUBLE Lambda1, const MDOUBLE Lambda2 , const MDOUBLE r)
: _Lambda1(Lambda1), _Lambda2(Lambda2)
{
	_gFuncStart0		= gFunc(Lambda1, Lambda2, r);
	_gFuncStart0MinusR	= gFunc(Lambda1, Lambda2, -r);
	_gFuncStart1		= gFunc(Lambda2, Lambda1, r);
	_gFuncStart1MinusR	= gFunc(Lambda2, Lambda1, -r);
}
computeJumps::~computeJumps()
{
}


/********************************************************************************************
getExpectation
*********************************************************************************************/
MDOUBLE computeJumps::getExpectation(const MDOUBLE BranchLength, int terminalStart, int terminalEnd, int fromId, int toId)
{
	if(fromId==0 && toId==1 && BranchLength>=0){ // Gain
		if(terminalStart==0 && terminalEnd==1)
			return gainExpGiven01(BranchLength);
		if(terminalStart==0 && terminalEnd==0)
			return gainExpGiven00(BranchLength);
		if(terminalStart==1 && terminalEnd==1)
			return gainExpGiven11(BranchLength);
		else //(terminalStart==1 && terminalEnd==0)
			return gainExpGiven10(BranchLength);
	}
	else
		return 0;

}
//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::gainExpGiven01(MDOUBLE BranchLength){
	return 0.5*(m01(BranchLength) +Pij_t(0,1,BranchLength));
}
MDOUBLE computeJumps::gainExpGiven00(MDOUBLE BranchLength){
	return 0.5*(m00(BranchLength));
}
MDOUBLE computeJumps::gainExpGiven11(MDOUBLE BranchLength){
	return 0.5*(m11(BranchLength) ); //???
}
MDOUBLE computeJumps::gainExpGiven10(MDOUBLE BranchLength){
	return 0.5*(m10(BranchLength) ); //???
}


MDOUBLE computeJumps::lossExpGiven01(MDOUBLE BranchLength){
	return 0.5*(m01(BranchLength) ); //???
}
MDOUBLE computeJumps::lossExpGiven00(MDOUBLE BranchLength){
	return 0.5*(m11(BranchLength)  ); //???
}
MDOUBLE computeJumps::lossExpGiven11(MDOUBLE BranchLength){
	return 0.5*(m11(BranchLength) ); //???
}
MDOUBLE computeJumps::lossExpGiven10(MDOUBLE BranchLength){
	return 0.5*(m10(BranchLength) + Pij_t(1,0,BranchLength) ); //???
}



//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::m01(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart0.gFunc_dr(BranchLength) - _gFuncStart0MinusR.gFunc_dr(BranchLength));
}
MDOUBLE computeJumps::m00(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart0.gFunc_dr(BranchLength) + _gFuncStart0MinusR.gFunc_dr(BranchLength));
}
MDOUBLE computeJumps::m11(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart1.gFunc_dr(BranchLength) - _gFuncStart1MinusR.gFunc_dr(BranchLength));
}
MDOUBLE computeJumps::m10(MDOUBLE BranchLength){
	return 0.5 *( _gFuncStart1.gFunc_dr(BranchLength) + _gFuncStart1MinusR.gFunc_dr(BranchLength));
}


//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::gFunc_dr(MDOUBLE BranchLength){
	return _gFuncStart0.g1Func_dr(BranchLength) + _gFuncStart0.g2Func_dr(BranchLength);
}
MDOUBLE computeJumps::gFunc::gFunc_dr(MDOUBLE BranchLength){
	return g1Func_dr(BranchLength) + g2Func_dr(BranchLength);
}


MDOUBLE computeJumps::gFunc::g1Func_dr(MDOUBLE BranchLength){
	return _g1Part_dr*g1Exp(BranchLength) + _g1Part*g1Exp(BranchLength)*BranchLength*_Alpha1_dr;
}
MDOUBLE computeJumps::gFunc::g2Func_dr(MDOUBLE BranchLength){
	return _g2Part_dr*g2Exp(BranchLength) +  _g2Part*g2Exp(BranchLength)*BranchLength*_Alpha2_dr;
}
//////////////////////////////////////////////////////////////////////////
MDOUBLE computeJumps::gFunc::g1Exp(MDOUBLE BranchLength){
	return exp(_Alpha1*BranchLength);
}
MDOUBLE computeJumps::gFunc::g2Exp(MDOUBLE BranchLength){
	return exp(_Alpha2*BranchLength);
}


//MDOUBLE computeJumps::gainExp(MDOUBLE BranchLength,MDOUBLE prob01,MDOUBLE prob11){
//	return gainExpGiven01(BranchLength)*prob01 + gainExpGiven00(BranchLength)*prob11;
//}


/********************************************************************************************
Pij_t - Based on Analytic solution
*********************************************************************************************/
MDOUBLE computeJumps::Pij_t(const int i,const int j, const MDOUBLE d)  {
	MDOUBLE gain = _Lambda1;
	MDOUBLE loss = _Lambda2;
	MDOUBLE eigenvalue =  -(gain + loss);
	bool withHGT = true;

	MDOUBLE noHGTfactor = 0.0001;

	VVdouble Pt;
	int AlphaSize = 2;
	resizeMatrix(Pt,AlphaSize,AlphaSize);
	int caseNum = i + j*2;
	switch (caseNum) {
		case 0 : Pt[0][0] =  loss/(-eigenvalue) + exp(eigenvalue*d)*(1 - loss/(-eigenvalue)); break;
		case 1 : Pt[1][0] =  loss/(-eigenvalue) - exp(eigenvalue*d)*(1 - gain/(-eigenvalue)); break;
		case 2 : if(withHGT)
				 { Pt[0][1] =  gain/(-eigenvalue) - exp(eigenvalue*d)*(1 - loss/(-eigenvalue));}
				 else 
				 { Pt[0][1] =  (gain/(-eigenvalue) - exp(eigenvalue*d)*(1 - loss/(-eigenvalue)))*noHGTfactor;} break;
		case 3 : Pt[1][1] =  gain/(-eigenvalue) + exp(eigenvalue*d)*(1 - gain/(-eigenvalue));  break;
	}
	MDOUBLE val = (Pt[i][j]);
	return val; 
}



/********************************************************************************************
*********************************************************************************************/
computeJumps::gFunc::gFunc(const MDOUBLE Lambda1, const MDOUBLE Lambda2 , const MDOUBLE r)
: _Lambda1(Lambda1), _Lambda2(Lambda2), _r(r)
{
	_delta = sqrt( pow((_Lambda1+_Lambda2),2) + 4*(_r*_r - 1)*_Lambda1*_Lambda2 );
	_delta_dr = (4*_r*_Lambda1*_Lambda2)/_delta;

	_Alpha1 = 0.5*(-_Lambda1-_Lambda2 +_delta);
	_Alpha2 = 0.5*(-_Lambda1-_Lambda2 -_delta);

	_Alpha1_dr =  0.5*_delta_dr;
	_Alpha2_dr = -0.5*_delta_dr;

	_Alpha1_2 = _delta;			//= _Alpha1-_Alpha2;
	_Alpha1_2_dr = _delta_dr;	//= _Alpha1_dr - _Alpha2_dr;

	_g1Part = ( (_r-1)*_Lambda1 - _Alpha2)/_Alpha1_2;
	_g2Part = (-(_r-1)*_Lambda1 + _Alpha1)/_Alpha1_2;

	_g1Part_dr = ( _Alpha1_2*( _Lambda1-_Alpha2_dr) - ( (_r-1)*_Lambda1 - _Alpha2)*_Alpha1_2_dr )/(_Alpha1_2*_Alpha1_2);
	_g2Part_dr = ( _Alpha1_2*(-_Lambda1+_Alpha1_dr) - (-(_r-1)*_Lambda1 + _Alpha1)*_Alpha1_2_dr )/(_Alpha1_2*_Alpha1_2);
}
