// $Id: bestGtrModelparams.cpp 2008-29-04 10:57:00Z nimrod $

#include "bestGtrModelParams.h"
#include <iostream>
using namespace std;

#include "bblEM.h"
#include "numRec.h"
#include "logFile.h"
#include "bestAlpha.h"

bestGtrModel::bestGtrModel(tree& et, // find best Gtr Model Params
										const sequenceContainer& sc,
										stochasticProcess& sp,
										const Vdouble * weights,
										const int maxTotalIterations,
										const MDOUBLE epsilonLikelihoodImprovment,
										const MDOUBLE epsilonLoglikelihoodForGTRParam,
										const MDOUBLE upperBoundGTRParam,
										const bool optimizeTree,
										const bool optimizeAlpha){
	LOG(5,<<"Starting bestGtrModel: find Best replacement matrix parameters"<<endl);
	MDOUBLE oldL = VERYSMALL;
	MDOUBLE newL = VERYSMALL;
	_bestL = likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(et,sc,sp,weights);

	MDOUBLE prev_a2c = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_a2c();
	MDOUBLE prev_a2g = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_a2g();
	MDOUBLE prev_a2t = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_a2t();
	MDOUBLE prev_c2g = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_c2g();
	MDOUBLE prev_c2t = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_c2t();
	MDOUBLE prev_g2t = (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->get_g2t();

	MDOUBLE prevAlpha = epsilonLoglikeForBBL;

	for (int i=0; i < maxTotalIterations; ++i) {
		//optimize a2c
		newL = -brent(0.0, prev_a2c, upperBoundGTRParam,
					  C_evalGTRParam(a2c,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_a2c);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2c(_best_a2c);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2c(prev_a2c);
			LOG(5,<<"likelihood went down in optimizing a2c"<<endl<<"oldL = "<<_bestL);
		}

		//optimize a2t
		newL = -brent(0.0, prev_a2t, upperBoundGTRParam,
					  C_evalGTRParam(a2t,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_a2t);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2t(_best_a2t);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2t(prev_a2t);
			LOG(5,<<"likelihood went down in optimizing a2t"<<endl<<"oldL = "<<_bestL);
		}

		//optimize a2g
		newL = -brent(0.0, prev_a2g, upperBoundGTRParam,
					  C_evalGTRParam(a2g,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_a2g);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2g(_best_a2g);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_a2g(prev_a2g);
			LOG(5,<<"likelihood went down in optimizing a2g"<<endl<<"oldL = "<<_bestL);
		}

		//optimize c2g
		newL = -brent(0.0, prev_c2g, upperBoundGTRParam,
					  C_evalGTRParam(c2g,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_c2g);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2g(_best_c2g);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2g(prev_c2g);
			LOG(5,<<"likelihood went down in optimizing c2g"<<endl<<"oldL = "<<_bestL);
		}

		//optimize c2t
		newL = -brent(0.0, prev_c2t, upperBoundGTRParam,
					  C_evalGTRParam(c2t,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_c2t);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2t(_best_c2t);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_c2t(prev_c2t);
			LOG(5,<<"likelihood went down in optimizing c2t"<<endl<<"oldL = "<<_bestL);
		}

		//optimize g2t
		newL = -brent(0.0, prev_g2t, upperBoundGTRParam,
					  C_evalGTRParam(g2t,et,sc,sp,weights),
					  epsilonLoglikelihoodForGTRParam,
					  &_best_g2t);
		if (newL >= _bestL) 
		{
			_bestL = newL;
			(static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_g2t(_best_g2t);//safety
		} 
		else
		{//likelihood went down!
            (static_cast<gtrModel*>(sp.getPijAccelerator()->getReplacementModel()))->set_g2t(prev_g2t);
			LOG(5,<<"likelihood went down in optimizing g2t"<<endl<<"oldL = "<<_bestL);
		}
	
		if(optimizeAlpha)
		{
			newL = -brent(0.0, prevAlpha, upperBoundForAlpha,
					  C_evalAlpha(et,sc,sp,weights),
					  epsilonLoglikeForAlphaOptimization,
					  &_bestAlpha);
			(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha); 

			if (newL >= _bestL) 
			{
				_bestL = newL;
				(static_cast<gammaDistribution*>(sp.distr()))->setAlpha(_bestAlpha); //safety
			} 
			else
			{//likelihood went down!
	            (static_cast<gammaDistribution*>(sp.distr()))->setAlpha(prevAlpha);
				LOG(5,<<"likelihood went down in optimizing alpha"<<endl<<"oldL = "<<_bestL);
			}
		}

		if(optimizeTree)
		{
			bblEM bblEM1(et,sc,sp,weights,maxBBLIt,epsilonLoglikeForBBL);
			_bestL = bblEM1.getTreeLikelihood();
		}


		// check for improvement in the likelihood
		if (_bestL > oldL+epsilonLikelihoodImprovment) {
			oldL = _bestL;
			prev_a2c = _best_a2c;
			prev_a2g = _best_a2g;
			prev_a2t = _best_a2t;
			prev_c2g = _best_c2g;
			prev_c2t = _best_c2t;
			prev_g2t = _best_g2t;
			prevAlpha = _bestAlpha;
		} else {
			break;
		}
	}
}
