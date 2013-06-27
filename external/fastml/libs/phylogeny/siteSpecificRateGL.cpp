// $Id: siteSpecificRate.cpp 3658 2008-03-05 09:25:46Z cohenofi $

#include "siteSpecificRateGL.h"
#include "numRec.h"
#include "checkcovFanctors.h"
#include "definitions.h"

using namespace siteSpecificRateGL;

MDOUBLE siteSpecificRateGL::computeML_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & likelihoodsV,
								   const sequenceContainer& sc,
								   const stochasticProcess& sp,
								   const tree& et,
								   const MDOUBLE maxRate,//20.0f
								   const MDOUBLE tol){//=0.0001f;

	ratesV.resize(sc.seqLen());
	likelihoodsV.resize(sc.seqLen());
	MDOUBLE Lsum = 0.0;

	for (int pos=0; pos < sc.seqLen(); ++pos) {
		siteSpecificRateGL::computeML_siteSpecificRate(pos,sc,sp,et,ratesV[pos],likelihoodsV[pos],maxRate,tol);
		assert(likelihoodsV[pos]>0.0);
		Lsum += log(likelihoodsV[pos]);
		LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
	return Lsum;
}

MDOUBLE siteSpecificRateGL::computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const vector<const stochasticProcess *> & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol){
	MDOUBLE Lsum = 0.0;
	ratesV.resize(sc.seqLen()); // the rates themselves
	likelihoodsV.resize(sc.seqLen()); // the log likelihood of each position
	
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		LOG(5,<<".");
		MDOUBLE bestR=-1.0; // tree1
		//		MDOUBLE LmaxR1=0;
		
		// getting the right tree for the specific position:
		const tree*  treeForThisPosition=NULL;
		if ((etVec.size() >0 ) && (treeAttributesVec[pos]>0)) {
			treeForThisPosition = & etVec[ treeAttributesVec[pos] -1];
		} else {
			errorMsg::reportError("tree vector is empty, or treeAttribute is empty, or treeAttribute[pos] is zero (it should be one)");
		}

		// getting the right stochastic process for the specific position:

		const stochasticProcess* spForThisPosition=NULL;

		if ((spVec.size() >0 ) && (spAttributesVec[pos]>0)) {
			spForThisPosition = spVec[ spAttributesVec[pos] -1];
		} else {
			errorMsg::reportError("stochastic process vector is empty, or spAttributesVec is empty, or spAttribute[pos] is zero (it should be one)");
		}

		siteSpecificRateGL::computeML_siteSpecificRate(pos,sc,*spForThisPosition,*treeForThisPosition,bestR,likelihoodsV[pos],maxRate,tol);
		ratesV[pos] = bestR;
		assert(likelihoodsV[pos]>0.0);
		Lsum += log(likelihoodsV[pos]);
		LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
	return Lsum;
}

// note that this places the likelihood, rather then the *log*likelihood into posL
void siteSpecificRateGL::computeML_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE& posL,
								 const MDOUBLE maxRate,
								 const MDOUBLE tol) {
	LOG(5,<<".");
	MDOUBLE ax=0.00001f,bx=maxRate*0.25,cx=maxRate;	// MN
	posL=-brent(ax,bx,cx,Cevaluate_L_given_r(sc,et,sp,pos),tol,&bestRate);
}

MDOUBLE siteSpecificRateGL::computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& treeAttributesVec,
						const vector<tree> & etVec,
						const stochasticProcess& sp,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol) {
	Vint spAttributesVec(sc.seqLen(),1);
	vector<const stochasticProcess* >  spVec;
	spVec.push_back(&sp);
	return computeML_siteSpecificRate(ratesV,likelihoodsV,
		spAttributesVec,treeAttributesVec,etVec,spVec,sc,maxRate,tol);
}

MDOUBLE siteSpecificRateGL::computeML_siteSpecificRate(Vdouble & ratesV,
						Vdouble & likelihoodsV,
						const Vint& spAttributesVec,
						const tree & et,
						const vector<const stochasticProcess* > & spVec,
						const sequenceContainer& sc,
						const MDOUBLE maxRate,
						const MDOUBLE tol){
	Vint treeAttributesVec(sc.seqLen(),1);				
	vector<tree>  etVec;
	etVec.push_back(et);
	return siteSpecificRateGL::computeML_siteSpecificRate(ratesV,likelihoodsV,
		spAttributesVec,treeAttributesVec,etVec,spVec,sc,maxRate,tol);
}


// THE BAYESIAN EB_EXP PART OF RATE ESTIMATION. //

void siteSpecificRateGL::computeEB_EXP_siteSpecificRate(int pos,
								 const sequenceContainer& sc,
								 const stochasticProcess& sp,
								 const computePijGam& cpg,
								 const tree &et,
								 MDOUBLE& bestRate,
								 MDOUBLE & stdRate,
								 MDOUBLE & lowerConf,
								 MDOUBLE & upperConf,
								 const MDOUBLE alphaConf, // alpha of 0.05 is considered 0.025 for each side.
								 VVdouble* LpostPerCat,
								 Vdouble* pLforMissingDataPerCat)
{
	// here we compute P(r | data)
	VdoubleRep pGivenR(sp.categories(),0.0);
	doubleRep sum=0;
	MDOUBLE LofPos_givenRateCat;
	LOG(8,<<pos+1<<"\t"); //DEBUG
	for (int cat=0; cat < sp.categories(); ++cat) {
		LofPos_givenRateCat = convert(likelihoodComputation::getLofPos(pos,et,sc,cpg[cat],sp));
		if(pLforMissingDataPerCat){
			LofPos_givenRateCat = LofPos_givenRateCat/(1- (*pLforMissingDataPerCat)[cat]);
		}
		pGivenR[cat] = LofPos_givenRateCat*sp.ratesProb(cat);
		LOG(8,<<cat<<"\t"<<LofPos_givenRateCat<<"\t"); //DEBUG
		sum+=pGivenR[cat];
	}
	LOG(8,<<"\n"); //DEBUG
	assert(sum!=0);
	
	// here we compute sigma r * P(r | data)

	doubleRep sumOfSquares(0.0);
	doubleRep bestRate_dblRep(0.0);
	
	LOG(5,<<"Pos "<<pos<<" content = "<<sc[0][pos]<<" ,total likelihood = "<<sum<<endl); //DEBUG
	
	for (int j=0; j < sp.categories(); ++j) {
		pGivenR[j]/=sum; // So that pGivenR is probability.
		                 // From here on we can convert it back
		                 // to MDOUBLE because it's not a very
		                 // small likelihood any more
		if (LpostPerCat){
			(*LpostPerCat)[j][pos]= convert(pGivenR[j]);
		}
		doubleRep tmp = pGivenR[j]*sp.rates(j);
		bestRate_dblRep += tmp;
		sumOfSquares += (tmp*sp.rates(j));

	}

	bestRate = convert(bestRate_dblRep);
	MDOUBLE varRate = convert(sumOfSquares) - convert(bestRate*bestRate);
	MDOUBLE tolerance = 0.0001; // tolerance for variance is not very exact, and also exact computation not very important
	if (varRate<-tolerance)
	    errorMsg::reportError("Error in computeEB_EXP_siteSpecificRate, varRate < 0");
	if ((varRate<0) && (varRate>=-tolerance))
	    varRate = 0;
	stdRate = sqrt(varRate);

	// detecting the confidence intervals.
	MDOUBLE oneSideConfAlpha = alphaConf/2.0; // because we are computing the two tail.
	doubleRep cdf = 0.0; // cumulative density function.
	int k=0;
	while (k < sp.categories()){
		cdf += convert(pGivenR[k]);
		if (cdf >oneSideConfAlpha) {
			lowerConf = sp.rates(k);
			break;
		} 
		k++;
	}
	while (k < sp.categories()) {
		if (cdf >(1.0-oneSideConfAlpha)) {
			upperConf = sp.rates(k);
			break;
		}
		++k;
		cdf += convert(pGivenR[k]);
	}
	if (k==sp.categories()) upperConf = sp.rates(k-1);
}

void siteSpecificRateGL::computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const sequenceContainer& sc,
								   const stochasticProcess& sp,
								   const tree& et,
								   const MDOUBLE alphaConf,
								   VVdouble* LpostPerCat,
								   Vdouble* pLforMissingDataPerCat)
{
	ratesV.resize(sc.seqLen());
	stdV.resize(sc.seqLen());
	lowerBoundV.resize(sc.seqLen());
	upperBoundV.resize(sc.seqLen());

	computePijGam cpg;
	cpg.fillPij(et,sp);
	for (int pos=0; pos < sc.seqLen(); ++pos) {
		siteSpecificRateGL::computeEB_EXP_siteSpecificRate(pos,sc,sp,cpg, et,ratesV[pos],stdV[pos],lowerBoundV[pos],upperBoundV[pos],alphaConf,LpostPerCat,pLforMissingDataPerCat);
		LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
}

void computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf){
	ratesV.resize(sc.seqLen());
	stdV.resize(sc.seqLen());
	lowerBoundV.resize(sc.seqLen());
	upperBoundV.resize(sc.seqLen());
	for (int treeNum=0; treeNum<etVec.size(); ++treeNum) {
		for (int spNum = 0; spNum<spVec.size(); ++spNum) {
            computePijGam cpg;
	    cpg.fillPij(etVec[treeNum],*(spVec[spNum]));
			for (int pos=0; pos < sc.seqLen(); ++pos) {
				if (((spAttributesVec[pos]-1)!=spNum ) || ((treeAttributesVec[pos]-1)!=treeNum )) continue;
				const tree*  treeForThisPosition=NULL;
				assert ((etVec.size() >0 ) && (treeAttributesVec[pos]>0));
				treeForThisPosition = & etVec[ treeAttributesVec[pos] -1];
				const stochasticProcess* spForThisPosition=NULL;
				assert ((spVec.size() >0 ) && (spAttributesVec[pos]>0));
				spForThisPosition = spVec[ spAttributesVec[pos] -1];
				siteSpecificRateGL::computeEB_EXP_siteSpecificRate(pos,sc,*spForThisPosition,cpg, *treeForThisPosition,ratesV[pos],stdV[pos],lowerBoundV[pos],upperBoundV[pos],alphaConf);
				LOG(5,<<" rate of pos: "<<pos<<" = "<<ratesV[pos]<<endl);
			}
		}
	}
	LOG(5,<<" number of sites: "<<sc.seqLen()<<endl);
}

// one tree many sps
void siteSpecificRateGL::computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& spAttributesVec,
							       const sequenceContainer& sc,
								   const tree & et,
								   const vector<const stochasticProcess *> & spVec,
								   const MDOUBLE alphaConf){
	Vint etAttributesVec(sc.seqLen(),1);				
	vector<tree>  etVec;
	etVec.push_back(et);
	siteSpecificRateGL::computeEB_EXP_siteSpecificRate(ratesV,stdV,lowerBoundV,upperBoundV,spAttributesVec,etAttributesVec,sc,etVec,spVec,alphaConf);
}

// one sp many trees

void siteSpecificRateGL::computeEB_EXP_siteSpecificRate(Vdouble & ratesV,
								   Vdouble & stdV,
								   Vdouble & lowerBoundV,
								   Vdouble & upperBoundV,
								   const Vint& treeAttributesVec,
							       const sequenceContainer& sc,
								   const vector<tree> & etVec,
								   const stochasticProcess & sp,
								   const MDOUBLE alphaConf){
	Vint spAttributesVec(sc.seqLen(),1);
	vector<const stochasticProcess* >  spVec;
	spVec.push_back(&sp);
	siteSpecificRateGL::computeEB_EXP_siteSpecificRate(ratesV,stdV,lowerBoundV,upperBoundV,spAttributesVec,treeAttributesVec,sc,etVec,spVec,alphaConf);
}

