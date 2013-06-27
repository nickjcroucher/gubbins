// 	$Id: bestParamUSSRV.cpp 4951 2008-09-24 11:16:58Z osnatz $	
#include "bestParamUSSRV.h"

/* structure of this method:
(1) checks of the number of parameters to optimize, and decide how many parameters optimizations iteration,
and how many parameters+bbl iterations will be done.
(2) A loop over the parameters+bbl iterations
	(2.1) A loop over the parameters optimization iterations
		(2.1.1) Optimize alpha
		(2.1.2) Optimize nu
		(2.1.3) Optimize f
		if the likelihood wasn't changed during this loop --> parameters converged --> break
	(2.2) BBL
	if the likelihood wasn't changed during this loop --> parameters+bbl converged --> break
(3) return likelihood
*/

// ***************
// *    USSRV    *
// ***************

MDOUBLE bestParamUSSRV::operator() (tree& et,
									const sequenceContainer& sc,
									const sequenceContainer& baseSc,
									ussrvModel& model,
									const Vdouble * weights /* =NULL */,
									const MDOUBLE AlphaUpperBound /* = 15 */, 
									const MDOUBLE NuUpperBound /* = 15 */, 
									const MDOUBLE FUpperBound /* = 1 */, 
									const MDOUBLE epsilonParamOptimization /* = 0.01 */,
									const MDOUBLE epsilonLikelihoodImprovment /* = 0.01 */,
									const int maxIterations /* = 50 */,
									const int maxOfParametersAndBblIterations /* = 40 */)
{
	_bestL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;	

	bestAlphaFixedTreeUSSRV alphaOptimization;
	bestNuFixedTreeUSSRV nuOptimization;
	bestFFixedTreeUSSRV fOptimization;
	
	int it, bblIt;
	int numberOfIterations(maxIterations);
	int numberOfParametersAndBblIterations(maxOfParametersAndBblIterations);
	
	// if only one parameter is optimize (only Alpha or only Nu or only F) then we need only one iteration.
	// if we only do bbl, without any optimization of the parameters, then we don't need iterations at all.
	int countParameters2Optimize(0);
	if (_AlphaOptimizationFlag) countParameters2Optimize++;
	if (_NuOptimizationFlag) countParameters2Optimize++;
	if (_FOptimizationFlag) countParameters2Optimize++;

	if (countParameters2Optimize==0)
	{
		numberOfIterations=0;
		numberOfParametersAndBblIterations=1;
	}
	else if (countParameters2Optimize==1)
		numberOfIterations=1;
	
	if (_bblOptimizationFlag == false)
		numberOfParametersAndBblIterations = 1;
	
	_bestAlpha = model.getAlpha();
	_bestNu = model.getNu();
	_bestF = model.getF();

	bool changes(false);
	bool bblChanges(false);
	for (bblIt=0; bblIt < numberOfParametersAndBblIterations; ++bblIt)
	{
		LOG(8,<<"bestParamUSSRV, params+bbl, iteration: " << bblIt << endl);
		bblChanges = false;
		// parameters optimizations (without bbl)
		// in each iteration : optimization of Alpha and then optimization of Nu, and then of F.
		for (it=0; it < numberOfIterations; ++it)
		{
			changes = false;	
			// Alpha optimization
			if (_AlphaOptimizationFlag)
			{
				LOGDO(5,printTime(myLog::LogFile()));
				newL = alphaOptimization(et,sc,baseSc,model,weights,AlphaUpperBound,epsilonParamOptimization);

				//the improvement in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{				
					LOG(5,<<"likelihood went down in LS! (Alpha optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous alpha
					alphaOptimization.setAlpha(_bestAlpha,model);
					alphaOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else 
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					LOG(9,<<"newL = " << newL << " _bestL = " << _bestL << " epsilonLikelihoodImprovment = " << epsilonLikelihoodImprovment << endl);
					_bestL = newL;
					_bestAlpha = alphaOptimization.getBestAlpha();
					LOG(5,<<"new L = " << _bestL<<"  new Alpha = " << _bestAlpha<<endl);		
				}
			}
		
			// Nu optimization
			if (_NuOptimizationFlag)
			{
				LOGDO(5,printTime(myLog::LogFile()));
				newL = nuOptimization(et,sc,baseSc,model,weights,NuUpperBound,epsilonParamOptimization);
			
				//the improvement in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{
					LOG(5,<<"likelihood went down in LS! (Nu optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous Nu
					nuOptimization.setNu(_bestNu,model);
					nuOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					LOG(9,<<"newL = " << newL << " _bestL = " << _bestL << " epsilonLikelihoodImprovment = " << epsilonLikelihoodImprovment << endl);
					_bestL = newL;
					_bestNu = nuOptimization.getBestNu();
					LOG(5,<<"new L = " << _bestL<<"  new Nu = " << _bestNu<<endl);		
				}
			}

			// F optimization
			if (_FOptimizationFlag)
			{
				LOGDO(5,printTime(myLog::LogFile()));
				newL = fOptimization(et,sc,baseSc,model,weights,FUpperBound,epsilonParamOptimization);

				//the improvement in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{
					LOG(5,<<"likelihood went down in LS! (F optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous F
					fOptimization.setF(_bestF,model);
					fOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else 
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment ) 
					{
						changes = true;
						bblChanges = true;
					}
					LOG(9,<<"newL = " << newL << " _bestL = " << _bestL << " epsilonLikelihoodImprovment = " << epsilonLikelihoodImprovment << endl);
					_bestL = newL;
					_bestF = fOptimization.getBestF();
					LOG(5,<<"new L = " << _bestL<<"  new F = " << _bestF<<endl);						
				}
			}
			if (changes == false)
			{
				LOG(5,<<"bestParamUSSRV parameters alpha,nu,f converged!"<<endl);
				break;
			}
		}

		if (changes == true)
			LOG(5,<<"bestParamUSSRV parameters alpha, nu, f, did not converge after " << numberOfIterations << " iterations"<<endl);


		// BBL
		if (_bblOptimizationFlag == true)
		{
			LOGDO(5,printTime(myLog::LogFile()));
			bblEM2USSRV bbl(et,sc,baseSc,model,weights,maxIterations);
			newL = bbl.getTreeLikelihood();
			LOG(5,<<"current best L= "<<_bestL<<endl);
			LOG(5,<<"new L After BBL = " << newL<< " = "<< bbl.getTreeLikelihood() <<endl);
			LOG(5,<<"The new tree is: " << endl);
			if (5 <= myLog::LogLevel()) 
				et.output(myLog::LogFile());
			LOG(5,<<endl);
			if (newL > _bestL+epsilonLikelihoodImprovment)
				bblChanges = true;
			if (newL < _bestL){
				LOG(5,<<"likelihood went down in LS! (BBL)"<<endl<<"oldL = "<<_bestL);
				LOG(5,<<" newL= "<<newL<<endl) ;
			}
			else
				_bestL = newL;
		}

		if (bblChanges == false)
		{
			LOG(5,<<"bestParamUSSRV bbl and parameters converged!"<<endl);
			break;
		}
	}

	if (bblIt == numberOfParametersAndBblIterations)
		LOG(5,<<"bestParamUSSRV bbl and parameters alpha did not converge after " << numberOfParametersAndBblIterations << "iterations"<<endl);	

	LOGDO(5,printTime(myLog::LogFile()));
	return _bestL;
}



// ***************
// *    SSRV    *
// ***************

MDOUBLE bestParamSSRV::operator() (tree& et,
								   const sequenceContainer& sc,
								   stochasticProcessSSRV& ssrvSp,
								   const Vdouble * weights /* =NULL */,
								   const MDOUBLE AlphaUpperBound /* = 15 */, 
								   const MDOUBLE NuUpperBound /* = 15 */, 
								   const MDOUBLE TrTvUpperBound /* = 10 */,
								   const MDOUBLE epsilonParamOptimization /* = 0.01 */,
								   const MDOUBLE epsilonLikelihoodImprovment /* = 0.01 */,
								   const MDOUBLE epsilonBbl /*= 0.05 */,
								   const int maxIterations /* = 50 */,
								   const int maxOfParametersAndBblIterations /* = 40 */)
{
	_bestL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;	

	bestAlphaFixedTreeSSRV alphaOptimization;
	bestNuFixedTreeSSRV nuOptimization;
	bestTamura92ParamFixedTreeSSRV tamura92Optimization;
	
	int it, bblIt;
	int numberOfIterations(maxIterations);
	int numberOfParametersAndBblIterations(maxOfParametersAndBblIterations);

	// if only one parameter is optimize (only Alpha or only Nu or only tamura92) then we need only one iteration.
	// if we only do bbl, without any optimization of the parameters, then we don't need iterations at all.
	int countParameters2Optimize(0);
	if (_AlphaOptimizationFlag) countParameters2Optimize++;
	if (_NuOptimizationFlag) countParameters2Optimize++;
	if (_tamura92OptimizationFlag) countParameters2Optimize++;


	if (countParameters2Optimize==0)
	{
		numberOfIterations=0;
		numberOfParametersAndBblIterations=1;
	}
	else if (countParameters2Optimize==1)
		numberOfIterations=1;

	if (_bblOptimizationFlag == false)
		numberOfParametersAndBblIterations = 1;

	replacementModelSSRV* pMulRM = static_cast<replacementModelSSRV*>(ssrvSp.getPijAccelerator()->getReplacementModel());
	gammaDistribution* gammaDist = static_cast<gammaDistribution*>(pMulRM->getDistribution()); 
	_bestAlpha = gammaDist->getAlpha();
	_bestNu = pMulRM->getRateOfRate();


	bool changes(false);
	bool bblChanges(false);

	for (bblIt=0; bblIt < numberOfParametersAndBblIterations; ++bblIt)
	{
		bblChanges = false;

		// Set initial values of lower/upper bounds for params
		MDOUBLE AlphaLowerBoundCur = 0.0;
		MDOUBLE AlphaUpperBoundCur = AlphaUpperBound;
		MDOUBLE NuLowerBoundCur = 0.0;
		MDOUBLE NuUpperBoundCur = NuUpperBound;
		MDOUBLE TrTvLowerBoundCur = 0.0;
		MDOUBLE TrTvUpperBoundCur = TrTvUpperBound;
		MDOUBLE ThetaLowerBoundCur = 0.0;
		MDOUBLE ThetaUpperBoundCur = 1.0;
		// And for epsilon
		MDOUBLE epsilonParamOptimizationCur = epsilonParamOptimization;

		// parameters optimizations (without bbl)
		// in each iteration : optimization of Alpha and then optimization of Nu, and then of F.
		for (it=0; it < numberOfIterations; ++it)
		{
			LOG(8,<<"bestParamUSSRV, params+bbl, iteration: " << bblIt << endl);
			changes = false;	
			// Alpha optimization
			if (_AlphaOptimizationFlag)
			{
				LOGDO(5,printTime(myLog::LogFile()));
				newL = alphaOptimization(et,sc,ssrvSp,weights,AlphaLowerBoundCur,AlphaUpperBoundCur,epsilonParamOptimizationCur);

				//the improvement in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{				
					LOG(5,<<"likelihood went down in LS! (Alpha optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous alpha
					alphaOptimization.setAlpha(_bestAlpha,ssrvSp);
					alphaOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else 
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					LOG(9,<<"newL = " << newL << " _bestL = " << _bestL << " epsilonLikelihoodImprovment = " << epsilonLikelihoodImprovment << endl);
					_bestL = newL;
					_bestAlpha = alphaOptimization.getBestAlpha();
					LOG(5,<<"new L = " << _bestL<<"  new Alpha = " << _bestAlpha<<endl);		
				}

				// Narrow search range between lower/upper bounds
				AlphaLowerBoundCur = (AlphaLowerBoundCur + 2*_bestAlpha) / 3;
				AlphaUpperBoundCur = (AlphaUpperBoundCur + 2*_bestAlpha) / 3;
			}

			// Nu optimization
			if (_NuOptimizationFlag)
			{
				LOGDO(5,printTime(myLog::LogFile()));
				newL = nuOptimization(et,sc,ssrvSp,weights,NuLowerBoundCur,NuUpperBoundCur,epsilonParamOptimizationCur);

				//the improvement in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{
					LOG(5,<<"likelihood went down in LS! (Nu optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
					//go back to previous Nu
					nuOptimization.setNu(_bestNu,ssrvSp);
					nuOptimization.setBestL(_bestL); // @@@@ maybe this is unnecessary 
					//break;
				}
				else
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					LOG(9,<<"newL = " << newL << " _bestL = " << _bestL << " epsilonLikelihoodImprovment = " << epsilonLikelihoodImprovment << endl);
					_bestL = newL;
					_bestNu = nuOptimization.getBestNu();
					LOG(5,<<"new L = " << _bestL<<"  new Nu = " << _bestNu<<endl);		
				}

				// Narrow search range between lower/upper bounds
				NuLowerBoundCur = (NuLowerBoundCur + 2*_bestNu) / 3;
				NuUpperBoundCur = (NuUpperBoundCur + 2*_bestNu) / 3;
			}

			// tamura92 optimization
			if (_tamura92OptimizationFlag)
			{
				LOGDO(5,printTime(myLog::LogFile()));
				newL = tamura92Optimization(
					et,sc,ssrvSp,weights,5,epsilonLikelihoodImprovment,
					TrTvLowerBoundCur,TrTvUpperBoundCur,ThetaLowerBoundCur,ThetaUpperBoundCur,
					epsilonParamOptimizationCur,epsilonParamOptimizationCur);
				MDOUBLE bestTrTv = tamura92Optimization.getBestTrTv();
				MDOUBLE bestTheta = tamura92Optimization.getBestTheta();

				//the improvement in Likelihood is smaller than epsilon
				if (newL < _bestL)
				{
					LOG(5,<<"likelihood went down in LS! (tamura92 optimization)"<<endl<<"oldL = "<<_bestL<<" newL= "<<newL<<endl);
				}
				else
				{// update of likelihood and model.
					if (newL > _bestL+epsilonLikelihoodImprovment) 
					{
						changes = true;
						bblChanges = true;
					}
					LOG(9,<<"newL = " << newL << " _bestL = " << _bestL << " epsilonLikelihoodImprovment = " << epsilonLikelihoodImprovment << endl);
					_bestL = newL;
					LOG(5,<<"new L = " << _bestL
						<<"  new TrTv = " << bestTrTv
						<<"  new Theta = " << bestTheta <<endl);
				}

				// Narrow search range between lower/upper bounds
				TrTvLowerBoundCur = (TrTvLowerBoundCur + 2*bestTrTv) / 3;
				TrTvUpperBoundCur = (TrTvUpperBoundCur + 2*bestTrTv) / 3;

				ThetaLowerBoundCur = (ThetaLowerBoundCur + 2*bestTheta) / 3;
				ThetaUpperBoundCur = (ThetaUpperBoundCur + 2*bestTheta) / 3;
			}

			if (changes == false)
			{
				LOG(5,<<"bestParamSSRV parameters alpha,nu, and tamura92 params converged!"<<endl);
				break;
			}

			// Reduce epsilonParamOptimizationCur
			epsilonParamOptimizationCur /= 2;
		}

		if (changes == true)
			LOG(5,<<"bestParamSSRV parameters alpha, nu, and tamura92 params did not converge after " << numberOfIterations << " iterations"<<endl);


		// BBL
		if (_bblOptimizationFlag == true)
		{
			LOGDO(5,printTime(myLog::LogFile()));
			bblEM bbl(et,sc,ssrvSp,weights,maxIterations,epsilonBbl);
			newL = bbl.getTreeLikelihood();
			LOG(5,<<" current best L= "<<_bestL<<endl);
			LOG(5,<<"new L After BBL = " << newL<< " = "<< bbl.getTreeLikelihood() <<endl);
			LOG(5,<<"The new tree is: " << endl);
			if (5 <= myLog::LogLevel()) 
				et.output(myLog::LogFile());
			LOG(5,<<endl);
			if (newL > _bestL+epsilonLikelihoodImprovment)
				bblChanges = true;
			if (newL < _bestL){
				LOG(5,<<"likelihood went down in LS! (BBL)"<<endl<<"oldL = "<<_bestL);
				LOG(5,<<" newL= "<<newL<<endl) ;
			}
			else
				_bestL = newL;
		}

		if (bblChanges == false)
		{
			LOG(5,<<"bestParamSSRV bbl and parameters converged!"<<endl);
			break;
		}
	}

	if (bblIt == numberOfParametersAndBblIterations)
		LOG(5,<<"bestParamSSRV bbl and parameters alpha did not converge after " << numberOfParametersAndBblIterations << "iterations"<<endl);	

	LOGDO(5,printTime(myLog::LogFile()));
	return _bestL;
}



// Variant that can work on a const tree - only if we're not doing BBL
// WARNING: Running this with bblOptimization==true will give a fatal error
MDOUBLE bestParamSSRV::operator() (const tree& et,
								   const sequenceContainer& sc,
								   stochasticProcessSSRV& ssrvSp,
								   const Vdouble * weights /* =NULL */,
								   const MDOUBLE AlphaUpperBound /* = 15 */, 
								   const MDOUBLE NuUpperBound /* = 15 */, 
								   const MDOUBLE TrTvUpperBound /* = 10 */,
								   const MDOUBLE epsilonParamOptimization /* = 0.01 */,
								   const MDOUBLE epsilonLikelihoodImprovment /* = 0.01 */,
								   const MDOUBLE epsilonBbl /*= 0.05 */,
								   const int maxIterations /* = 50 */,
								   const int maxOfParametersAndBblIterations /* = 40 */)
{
	if (_bblOptimizationFlag == true)
		errorMsg::reportError("bestParamSSRV::operator(): Can't work on const tree if bblOptimization was requested");

	tree etNotConst(et);
	return operator()(etNotConst, sc, ssrvSp, weights,
					  AlphaUpperBound, NuUpperBound, 
					  epsilonParamOptimization, epsilonLikelihoodImprovment,
					  epsilonBbl, maxIterations,
					  maxOfParametersAndBblIterations);
}


