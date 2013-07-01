#include "bestAlphaAndK.h"
#include "computePijComponent.h"
#include "betaOmegaDistribution.h"
#include "codonUtils.h"


optimizeSelectonParameters::optimizeSelectonParameters(tree& et, //find Best params and best BBL
					   const sequenceContainer& sc,
					   vector<stochasticProcess>& spVec,
					   distribution * distr,
					   bool bblFlag,
					   bool isGamma, bool isBetaProbSet,bool isOmegaSet,
					   bool isKappaSet, bool isAlphaSet, bool isBetaSet,
					   const MDOUBLE upperBoundOnAlpha,
					   const MDOUBLE upperBoundOnBeta,
					   const MDOUBLE epsilonAlphaOptimization,
					   const MDOUBLE epsilonKOptimization,
					   const MDOUBLE epsilonLikelihoodImprovment,
					   const int maxBBLIterations,
					   const int maxTotalIterations){
   //initialization	
	MDOUBLE lowerValueOfParamK = 0;
	MDOUBLE lowerValueOfParamAlpha = 0.1;
	MDOUBLE lowerValueOfParamBeta = 0.1;
	MDOUBLE omegaLowerBoundary = 0.99; // this is to allow brent to reach the exact lower bound value
	MDOUBLE omegaUpperBoundary = 5.0; 
	MDOUBLE upperValueOfParamK = 5; // changed from 50, Adi S. 2/1/07 
	
	MDOUBLE initialGuessValueOfParamTr;
	initialGuessValueOfParamTr = _bestK = static_cast<wYangModel*>(spVec[0].getPijAccelerator()->getReplacementModel())->getK();

	MDOUBLE initialGuessValueOfParamAlpha;
	if (isGamma) initialGuessValueOfParamAlpha = _bestAlpha = static_cast<generalGammaDistribution*>(distr)->getAlpha();
	else initialGuessValueOfParamAlpha = _bestAlpha = static_cast<betaOmegaDistribution*>(distr)->getAlpha();
	
	MDOUBLE initialGuessValueOfParamBeta; 
	if (isGamma) initialGuessValueOfParamBeta = _bestBeta = static_cast<generalGammaDistribution*>(distr)->getBeta();
	else initialGuessValueOfParamBeta = _bestBeta = static_cast<betaOmegaDistribution*>(distr)->getBeta();

	MDOUBLE initialGuessValueOfParamOmega = -1;
	MDOUBLE initialGuessValueOfParamBetaProb = -1;
	if (!isGamma) {
		initialGuessValueOfParamOmega = _bestOmega = static_cast<betaOmegaDistribution*>(distr)->getOmega();
		initialGuessValueOfParamBetaProb = _bestBetaProb = static_cast<betaOmegaDistribution*>(distr)->getBetaProb();
	}
	_bestL = likelihoodComputation2Codon::getTreeLikelihoodAllPosAlphTheSame(et,sc,spVec,distr);;
	MDOUBLE newL = _bestL;

	MDOUBLE alphaFound = 0;
	MDOUBLE kFound = 0;
	MDOUBLE betaFound = 0;
	MDOUBLE omegaFound = 0;
	MDOUBLE betaProbFound = 0;
	bool changed = false;
	int i=0;
	LOG(5,<<endl<<"Beginning optimization of parameters"<<endl<<endl);

	for (i=0; i < maxTotalIterations; ++i) {
		LOG(5,<<"Iteration Number= " << i <<endl);
		LOG(5,<<"---------------------"<<endl);		
		cout<<"Iteration number = "<< i <<endl;
		alphaFound = omegaFound = betaProbFound = kFound = betaFound=0;
		changed = false;
//ALPHA (beta or gamma distribution parameter)
		if (!isAlphaSet){
			if (isGamma) initialGuessValueOfParamAlpha = static_cast<generalGammaDistribution*>(distr)->getAlpha();
			else initialGuessValueOfParamAlpha = static_cast<betaOmegaDistribution*>(distr)->getAlpha();
			newL = -brent(lowerValueOfParamAlpha,
						initialGuessValueOfParamAlpha,
						upperBoundOnAlpha,
						evalParam(et,sc,spVec,-1,distr,isGamma),epsilonAlphaOptimization,&alphaFound); 

			LOG(5,<<"current best L= "<<_bestL<<endl<<endl);
			LOG(5,<<"new L After alpha= " << newL<<endl);
			LOG(5,<<"new alpha = " <<alphaFound<<endl<<endl);

			
			if (newL > _bestL+epsilonLikelihoodImprovment ) {// update of likelihood ,v and model.
				_bestL = newL;
				_bestAlpha = alphaFound;
				if (isGamma) static_cast<generalGammaDistribution*>(distr)->setAlpha(alphaFound);
				else static_cast<betaOmegaDistribution*>(distr)->setAlpha(alphaFound);
				for (int categor = 0; categor < spVec.size();categor++)
					static_cast<wYangModel*>(spVec[categor].getPijAccelerator()->getReplacementModel())->setW(distr->rates(categor)); 
				normalizeMatrices(spVec,distr);
				changed = true;
			} 
		}
//BETA (beta distribution parameter)
		if (!isBetaSet) {
			if (isGamma) initialGuessValueOfParamBeta = static_cast<generalGammaDistribution*>(distr)->getBeta();
			else initialGuessValueOfParamBeta = static_cast<betaOmegaDistribution*>(distr)->getBeta();
			newL = -brent(lowerValueOfParamBeta,
						initialGuessValueOfParamBeta,
						upperBoundOnBeta,
						evalParam(et,sc,spVec,-2,distr,isGamma),epsilonAlphaOptimization,&betaFound); 

			LOG(5,<<"current best L= "<<_bestL<<endl<<endl);
			LOG(5,<<"new L After beta= " << newL<<endl);
			LOG(5,<<"new beta = " <<betaFound<<endl<<endl);
		
			if (newL > _bestL+epsilonLikelihoodImprovment ) {// update of likelihood ,v and model.
				_bestL = newL;
				_bestBeta = betaFound;
				if (isGamma) static_cast<generalGammaDistribution*>(distr)->setBeta(betaFound);
				else static_cast<betaOmegaDistribution*>(distr)->setBeta(betaFound);
				for (int categor = 0; categor < spVec.size();categor++)
					static_cast<wYangModel*>(spVec[categor].getPijAccelerator()->getReplacementModel())->setW(distr->rates(categor)); 		
				normalizeMatrices(spVec,distr);
				changed = true;
			}
		}
//K parameter
		if (!isKappaSet){
			initialGuessValueOfParamTr =  static_cast<wYangModel*>(spVec[0].getPijAccelerator()->getReplacementModel())->getK();
			newL = -brent(lowerValueOfParamK,   //optimaize Tr
					initialGuessValueOfParamTr,
					upperValueOfParamK,
					evalParam(et,sc,spVec,0,distr,isGamma),epsilonKOptimization,&kFound); 
			
			LOG(5,<<"current best L= "<<_bestL<<endl<<endl);
			LOG(5,<<"new L After kappa= " << newL<<endl);
			LOG(5,<<"new kappa = " <<kFound<<endl);

			if (newL > _bestL+epsilonLikelihoodImprovment ) {// update of likelihood and model.
				_bestL = newL;
				_bestK = kFound;
				for (int categor = 0; categor < spVec.size();categor++)
					static_cast<wYangModel*>(spVec[categor].getPijAccelerator()->getReplacementModel())->setK(kFound); 
				normalizeMatrices(spVec,distr);
				changed = true;
			}
		}
//beta distribution part (betaProb and additional omega)
		if (isGamma==false && !isBetaProbSet){ //optimize  beta probs
			if (!isOmegaSet){ // optimize omega  (M8 or M8b)
				MDOUBLE omegaFound;
				newL = -brent(omegaLowerBoundary, 
						initialGuessValueOfParamOmega,
						omegaUpperBoundary,
						evalParam(et,sc,spVec,1,distr,isGamma),0.01,&omegaFound); 

				LOG(5,<<"current best L= "<<_bestL<<endl<<endl);
				LOG(5,<<"new L After additional omega caetgory = " << newL<<endl);
				LOG(5,<<"new additional omega caetgory = " <<omegaFound<<endl<<endl);
	
				if (newL > _bestL+epsilonLikelihoodImprovment ) {
					_bestL = newL;
					_bestOmega = omegaFound;
					static_cast<betaOmegaDistribution*>(distr)->setOmega(omegaFound);
					static_cast<wYangModel*>(spVec[spVec.size()-1].getPijAccelerator()->getReplacementModel())->setW(omegaFound); 	
					normalizeMatrices(spVec,distr);
					changed = true;
				}
			}
			MDOUBLE betaProbFound;	
			newL = -brent(0.0,initialGuessValueOfParamBetaProb,1.0,
					evalParam(et,sc,spVec,2,distr,isGamma),0.01,&betaProbFound); 

			LOG(5,<<"current best L= "<<_bestL<<endl<<endl);
			LOG(5,<<"new L After prob(additional omega caetgory)= " << newL<<endl);
			LOG(5,<<"new prob(additional omega caetgory)= " <<1 - betaProbFound<<endl<<endl);
			if (newL > _bestL+epsilonLikelihoodImprovment ) {// update of likelihood ,v and model.
				_bestL = newL;
				_bestBetaProb = betaProbFound;
				static_cast<betaOmegaDistribution*>(distr)->setBetaProb(betaProbFound);
				normalizeMatrices(spVec,distr);
				changed = true;
			}
		}

//BBL
		if (bblFlag==true) {
//using epsilonAlphaOptimization as the epsilon for pairwise disatnce here		
			bblEM2codon bbl(et,sc,spVec,distr,NULL,maxBBLIterations,epsilonLikelihoodImprovment,epsilonAlphaOptimization);
			newL = bbl.getTreeLikelihood();
		
			LOG(5,<<"current best L= "<<_bestL<<endl<<endl);
			LOG(5,<<"new L After BL = " << newL<<endl);
			LOG(5,<<"Tree after this BBL iteration: "<<endl);
			LOGDO(5,et.output(myLog::LogFile()));
			
			if (newL > _bestL+epsilonLikelihoodImprovment) {
				_bestL = newL;
				changed = true;
			}
		}
	
		if (changed==false)
			break;
		
	}

	LOG(5,<<endl<<"Finished optimization of parameters"<<endl<<endl);

	if (i==maxTotalIterations) {
	  LOG(5,<<"Too many iterations in function optimizeCodonModelAndBBL. The last optimized parameters are used for the calculations."<<endl<<endl);
		
	}
	
}

evalParam::~evalParam(){
  if (_distr != NULL) delete _distr;
}


evalParam::evalParam(const evalParam &other): _et(other._et),_sc(other._sc),
_spVec(other._spVec), _alphaOrKs(other._alphaOrKs),_isGamma(other._isGamma)	
{
	_distr=other._distr->clone();
}


MDOUBLE evalParam::operator()(MDOUBLE param){

	if (_alphaOrKs==-1) updateAlpha(param);
	else if (_alphaOrKs==-2) updateBeta(param);
	else if (_alphaOrKs==0) updateK(param);
	else if (_alphaOrKs==1) updateOmega(param);
	else if (_alphaOrKs==2) updateBetaProb(param);
	MDOUBLE res = likelihoodComputation2Codon::getTreeLikelihoodAllPosAlphTheSame(_et,_sc,_spVec,_distr);
	return -res;	//return -log(likelihood).
}

void evalParam::updateBeta(MDOUBLE param){
	if (_isGamma) static_cast<generalGammaDistribution*>(_distr)->setBeta(param);
	else  static_cast<betaOmegaDistribution*>(_distr)->setBeta(param);
	for (int categor = 0; categor < _spVec.size();categor++){
		static_cast<wYangModel*>(_spVec[categor].getPijAccelerator()->getReplacementModel())->setW(_distr->rates(categor)); 
		
	}
	normalizeMatrices(_spVec,_distr);
}
void evalParam::updateAlpha(MDOUBLE param){
	if (_isGamma)static_cast<generalGammaDistribution*>(_distr)->setAlpha(param);
	else static_cast<betaOmegaDistribution*>(_distr)->setAlpha(param);
	for (int categor = 0; categor < _spVec.size();categor++){
		static_cast<wYangModel*>(_spVec[categor].getPijAccelerator()->getReplacementModel())->setW(_distr->rates(categor)); 
		
	}
	normalizeMatrices(_spVec,_distr);
}

void evalParam::updateK(MDOUBLE param){
	for (int categor = 0; categor < _spVec.size();categor++){
		static_cast<wYangModel*>(_spVec[categor].getPijAccelerator()->getReplacementModel())->setK(param); 
	}
	normalizeMatrices(_spVec,_distr);
}


void evalParam::updateOmega(MDOUBLE param){
	int size = _spVec.size();
	static_cast<wYangModel*>(_spVec[size-1].getPijAccelerator()->getReplacementModel())->setW(param); 
	normalizeMatrices(_spVec,_distr);
}

void evalParam::updateBetaProb(MDOUBLE param){
	static_cast<betaOmegaDistribution*>(_distr)->setBetaProb(param);
	normalizeMatrices(_spVec,_distr);
}
