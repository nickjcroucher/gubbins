// 	$Id: likelihoodComputation2USSRV.cpp 962 2006-11-07 15:13:34Z privmane $	
#include "likelihoodComputation2USSRV.h"


using namespace likelihoodComputation2USSRV;

//compute likelihood for the ssrv model and the base model.

MDOUBLE likelihoodComputation2USSRV::getTreeLikelihoodAllPosAlphTheSame(const tree& et,
							const sequenceContainer& sc, const sequenceContainer& baseSc,
							const ussrvModel& model,const Vdouble * const weights){
							
	
	computePijHom piSSRV;
	piSSRV.fillPij(et,model.getSSRVmodel());
	
	computePijGam piBase;
	piBase.fillPij(et,model.getBaseModel());
	
	MDOUBLE res =0.0;
	MDOUBLE f = model.getF();
	doubleRep LofPosSSRV(0.0),LofPosBase(0.0);
	MDOUBLE lnL(0.);
	int k;
	for (k=0; k < sc.seqLen(); ++k) {
		if (f<1.0)
			LofPosBase = likelihoodComputation::getLofPos(k,et,baseSc,piBase,model.getBaseModel());
		if (f>0.0) {
			LofPosSSRV = likelihoodComputation::getLofPos(k,et,sc,piSSRV,model.getSSRVmodel());
			if (f<1.0) 
				lnL = log(LofPosSSRV*f+(1-f)*LofPosBase);
			else // f == 1.0
				lnL = log(LofPosSSRV);
		}
		else // f == 0.0
			lnL = log(LofPosBase);

		LOG(9,<<"pos= "<<k<<" lnL= "<<lnL<<endl);
		LOG(10,<<"logLofPosBase= "<< log(LofPosBase) << " logLofPosSSRV= " << log(LofPosSSRV) << " f= " << f <<endl);
		res += lnL * (weights?(*weights)[k]:1);
	}
	return res;
}

				
				
MDOUBLE likelihoodComputation2USSRV::getTreeLikelihoodFromUp2(const tree& et,
						const sequenceContainer& sc,
						const sequenceContainer& baseSc,
						const ussrvModel & model,
						const suffStatGlobalGam& cupBase,
						const suffStatGlobalHom& cupSSRV,
						VdoubleRep& posLike, // fill this vector with each position likelihood but without the weights.
						const Vdouble * weights) {
	posLike.clear();
	MDOUBLE like = 0;
	MDOUBLE f = model.getF();
	//computing the likelihood from up:
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRep tmp=0;
		
		doubleRep tmp2 = 0; //like for the SSRV part
		// SSRV
		for (int let =0; let < model.getSSRVmodel().alphabetSize(); ++let) {
				tmp2+=cupSSRV.get(pos,et.getRoot()->id(),let) * model.getSSRVmodel().freq(let);
		}
		// Base model
		for (int categor = 0; categor < model.noOfCategor(); ++categor) {
			doubleRep veryTmp =0;
			for (int let =0; let < model.getBaseModel().alphabetSize(); ++let) {
				veryTmp+=cupBase.get(pos,categor,et.getRoot()->id(),let) * model.getBaseModel().freq(let);
			}
			tmp += veryTmp*model.getCategorProb(categor);
		}

		if(tmp<0.0) errorMsg::reportError("like< 0 in likelihoodComputation2USSRV::getTreeLikelihoodFromUp2");

		like += log((1-f)*tmp+f*tmp2) * (weights?(*weights)[pos]:1);
		posLike.push_back((1-f)*tmp+f*tmp2);
	}
	return like;
}
