#include "likelihoodComputationGL.h"

#include "definitions.h"
#include "tree.h"
#include "likelihoodComputation.h"
#include <cmath>
#include <cassert>

using namespace likelihoodComputationGL;

// account for RateCat, GainCat,LossCat 
// - For each RateCat an "external" multiplication is conducted - copy_et.multipleAllBranchesByFactor
// - the GainCat*LossCat SPs are covered by the "internal" mechanism of PijGam

/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(const tree& tr,
							const sequenceContainer& sc,
							const vector<vector<stochasticProcess*> >& spVVec,
							const distribution * distGain, const distribution * distLoss,
							unObservableData *unObservableData_p)
{	
	int numOfRateCategories = spVVec[0][0]->categories();
	vector<computePijGam> pi_vec(numOfRateCategories);
	vector<suffStatGlobalGam> ssc_vec(numOfRateCategories);
	vector<computeUpAlg> cup_vec(numOfRateCategories);
	
	likelihoodComputationGL::fillPijAndUp(tr,sc,spVVec,distGain,distLoss,pi_vec,ssc_vec,cup_vec);
	MDOUBLE res = 0.0;
	for (int k=0; k < sc.seqLen(); ++k) {
		MDOUBLE lnL = 0;
		MDOUBLE resGivenRate = 0.0;
		for(int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){
			lnL = log(likelihoodComputationGL::getProbOfPosUpIsFilledSelectionGam(k,//pos,
				tr,//const tree& 
				sc,// sequenceContainer& sc,
				spVVec,	// only needed for sp.freq(let)
				ssc_vec[rateIndex][k],//const computePijGam& ,
				distGain, distLoss)); // distributions ,
			resGivenRate += lnL * spVVec[0][0]->ratesProb(rateIndex);
		}
		LOG(20,<<"pos= "<<k+1<<" lnL= "<<lnL<<endl);
		//res += lnL;
		res += resGivenRate;
	}
	if(unObservableData_p){
		res = res - sc.seqLen()*log(1- exp(unObservableData_p->getlogLforMissingData()));
	}
	return res;
}
/********************************************************************************************
*********************************************************************************************/
void likelihoodComputationGL::fillPijAndUp(const tree& tr,
										   const sequenceContainer& sc,
										   const vector<vector<stochasticProcess*> >& spVVec,
										   const distribution * distGain, const distribution * distLoss,
										   vector<computePijGam>& pi_vec,
										   vector<suffStatGlobalGam>& ssc_vec, // info filled into suffStat
										   vector<computeUpAlg>& cup_vec)
{	
	int numOfSPs = distGain->categories()*distLoss->categories();
	int numOfRateCategories = spVVec[0][0]->categories();
	for (int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){		
		tree copy_et = tr;
		copy_et.multipleAllBranchesByFactor(spVVec[0][0]->rates(rateIndex));
		pi_vec[rateIndex]._V.resize(numOfSPs);
		//Pij
		for (int i=0; i < numOfSPs; ++i) {
			int gainIndex =fromIndex2gainIndex(i,distGain->categories(),distLoss->categories());
			int lossIndex =fromIndex2lossIndex(i,distGain->categories(),distLoss->categories());
			pi_vec[rateIndex]._V[i].fillPij(copy_et,*spVVec[gainIndex][lossIndex]);
		}
		//ComputeUp
		cup_vec[rateIndex].fillComputeUp(copy_et,sc,pi_vec[rateIndex],ssc_vec[rateIndex]);
	}
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputationGL::getProbOfPosUpIsFilledSelectionGam(const int pos,const tree& tr,
						const sequenceContainer& sc,
						const vector<vector<stochasticProcess*> >& spVVec,// only needed for sp.freq(let)
						const suffStatGlobalGamPos& cup,
						const distribution * distGain, const distribution * distLoss)
{

	doubleRep res =0;
	int numOfSPs = distGain->categories()*distLoss->categories();
	for (int categor = 0; categor < numOfSPs; ++categor) {
		doubleRep veryTmp =0.0;
		int gainCategor = fromIndex2gainIndex(categor,distGain->categories(),distLoss->categories());
		int lossCategor = fromIndex2lossIndex(categor,distGain->categories(),distLoss->categories());
		for (int let =0; let < sc.alphabetSize(); ++let) {
			veryTmp+=cup.get(categor,tr.getRoot()->id(),let) * spVVec[gainCategor][lossCategor]->freq(let);			
		}
		res += veryTmp*(distGain->ratesProb(gainCategor)*distLoss->ratesProb(lossCategor));
	}
	if ((res<-EPSILON)){
		string err = "Error in likelihoodComputationGL::getProbOfPosUpIsFilledSelectionGam, non probability value (<0) Res=";
		err+=double2string(convert(res));
		errorMsg::reportError(err);
	};
	return convert(res);
}
/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputationGL::getTreeLikelihoodFromUp2(const tree& tr,
						const sequenceContainer& sc,
						const vector<vector<stochasticProcess*> >& spVVec,// only needed for sp.freq(let)
						const suffStatGlobalGam& cup, 	//computing the likelihood from up:
						const distribution * distGain, const distribution * distLoss,
						unObservableData *unObservableData_p,
						Vdouble* posLike,
						const Vdouble * weights) 
{
	if(posLike)
		posLike->clear();
	MDOUBLE like = 0;

	int numOfSPs = distGain->categories()*distLoss->categories();
	for (int pos = 0; pos < sc.seqLen(); ++pos) {
		doubleRep tmp=0;
		for (int categor = 0; categor < numOfSPs; ++categor) {
			doubleRep veryTmp =0;
			int gainCategor = fromIndex2gainIndex(categor,distGain->categories(),distLoss->categories());
			int lossCategor = fromIndex2lossIndex(categor,distGain->categories(),distLoss->categories());
			for (int let =0; let < sc.alphabetSize(); ++let) {
				veryTmp+=cup.get(pos,categor,tr.getRoot()->id(),let) * spVVec[gainCategor][lossCategor]->freq(let);
			}
			tmp += veryTmp*(distGain->ratesProb(gainCategor)*distLoss->ratesProb(lossCategor));
		}
		if(unObservableData_p)
		    tmp = tmp/(1- exp(unObservableData_p->getlogLforMissingData()));
		if(posLike)
			posLike->push_back(log(tmp));
		like += log(tmp) * (weights?(*weights)[pos]:1);

	}
	return like;
}

/********************************************************************************************
*********************************************************************************************/
MDOUBLE likelihoodComputationGL::getTreeLikelihoodFromUp2(const tree& tr,
														  const sequenceContainer& sc,
														  const vector<vector<stochasticProcess*> >& spVVec,// only needed for sp.freq(let)
														  const vector<suffStatGlobalGam>& cup_vec, 	//computing the likelihood from up:
														  const distribution * distGain, const distribution * distLoss,
														  unObservableData *unObservableData_p,
														  Vdouble* posLike,
														  const Vdouble * weights) 
{
	if(posLike)
		posLike->resize(sc.seqLen());
	MDOUBLE like = 0;
	int numOfRateCategories = spVVec[0][0]->categories();
	for(int rateIndex=0 ; rateIndex<numOfRateCategories; ++rateIndex){
		Vdouble posLikePerCat;
		like += likelihoodComputationGL::getTreeLikelihoodFromUp2(tr,sc,spVVec,cup_vec[rateIndex], distGain,distLoss,unObservableData_p,&posLikePerCat,weights)
			* spVVec[0][0]->ratesProb(rateIndex);
		if(posLike){
			for (int k=0; k < sc.seqLen(); ++k) {			
				(*posLike)[k]+= (posLikePerCat[k]* spVVec[0][0]->ratesProb(rateIndex));
			}
		}
	}
	return like;
}

/********************************************************************************************
*********************************************************************************************/
//MDOUBLE likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSameNoComputeUp(const tree& tr,
//																			   const sequenceContainer& sc,
//																			   const vector<vector<stochasticProcess*> >& spVVec,
//																			   const distribution * distGain, const distribution * distLoss,
//																			   unObservableData *unObservableData_p)
//{	
//	MDOUBLE res = 0.0;
//	int numOfSPs = distGain->categories()*distLoss->categories();
//	for (int i=0; i < numOfSPs; ++i) {
//		int gainIndex =fromIndex2gainIndex(i,distGain->categories(),distLoss->categories());
//		int lossIndex =fromIndex2lossIndex(i,distGain->categories(),distLoss->categories());
//		res += likelihoodComputation::getTreeLikelihoodAllPosAlphTheSame(tr,sc,*spVVec[gainIndex][lossIndex])* distGain->ratesProb(gainIndex)*distLoss->ratesProb(lossIndex);
//	}
//	if(unObservableData_p){
//		res = res - sc.seqLen()*log(1- exp(unObservableData_p->getlogLforMissingData()));
//	}
//	return res;
//}




/********************************************************************************************
un-obervable data
*********************************************************************************************/

/********************************************************************************************
 used to fill the likelihood for the unobervable for each category
*********************************************************************************************/
//doubleRep likelihoodComputationGL::getLofPos(const int pos,
//										   const tree& tr,
//										   const sequenceContainer& sc,
//										   const computePijGam& pi,
//										   const stochasticProcess& sp,
//										   Vdouble& likePerCat)	// all the likdelhoodsPerCat and rateProb are filled
//{
//	//	with the pi already computed.
//	int numOfCat = sp.categories();
//	doubleRep tmp=0;
//	for (int i=0; i < numOfCat;++i) {
//		likePerCat[i] = getLofPos(pos,tr,sc,pi[i],sp)*sp.ratesProb(i);
//		likePerCat[i+numOfCat] = sp.ratesProb(i);
//		tmp += likePerCat[i];
//	}
//	return tmp;
//}
///********************************************************************************************
//likelihood computation - full data (1)
//*********************************************************************************************/
//MDOUBLE likelihoodComputationGL::getTreeLikelihoodAllPosAlphTheSame(const tree& tr,
//																	const sequenceContainer& sc,
//																	const stochasticProcess& sp,
//																	const Vdouble * const weights,
//																	Vdouble *pLforMissingDataPerCat)
//{
//	computePijGam pi;
//	pi.fillPij(tr,sp);
//	MDOUBLE res =0;
//	doubleRep LofPos;
//	int k;
//	for (k=0; k < sc.seqLen(); ++k) {
//		LofPos = likelihoodComputationGL::getLofPos(k,//pos,
//			tr,//const tree& 
//			sc,// sequenceContainer& sc,
//			pi,//const computePijGam& ,
//			sp,
//			pLforMissingDataPerCat);
//		res += log(LofPos) * (weights?(*weights)[k]:1);//const stochasticProcess& );
//	}
//	return res;
//}
//
///********************************************************************************************
//likelihood computation - per pos (1.1)
//*********************************************************************************************/
//doubleRep likelihoodComputationGL::getLofPos(const int pos,
//										   const tree& tr,
//										   const sequenceContainer& sc,
//										   const computePijGam& pi,
//										   const stochasticProcess& sp,
//										   Vdouble *pLforMissingDataPerCat)
//{
////	with the pi already computed.
//	doubleRep tmp=0;
//	int numOfCat = sp.categories();
//	Vdouble tmpPerCat;
//	tmpPerCat.resize(numOfCat);	
//	
//	for (int i=0; i < sp.categories();++i) {
//		tmpPerCat[i] = getLofPos(pos,tr,sc,pi[i],sp);
//		if(pLforMissingDataPerCat){
//			LOG(11,<<"res before MissingData correction= "<<tmpPerCat[i]);
//			tmpPerCat[i] = tmpPerCat[i]/(1- (*pLforMissingDataPerCat)[i]);
//			LOG(11,<<" after= "<<tmpPerCat[i]<<endl);
//		}
//		tmp += tmpPerCat[i]*sp.ratesProb(i);
//	}
//	return tmp;
//}
//
///********************************************************************************************
//likelihood computation - per pos, per cat (1.1.1)
//*********************************************************************************************/
//doubleRep likelihoodComputationGL::getLofPos(const int pos,
//					  const tree& tr,
//					  const sequenceContainer& sc,
//					  const computePijHom& pi,
//					  const stochasticProcess& sp)
//{
//	computeUpAlg cup;
//	suffStatGlobalHomPos ssc;
//	cup.fillComputeUp(tr,sc,pos,pi,ssc);
//
//	doubleRep tmp = 0.0;
//	for (int let = 0; let < sp.alphabetSize(); ++let) {
//		doubleRep tmpLcat=
//				ssc.get(tr.getRoot()->id(),let)*
//				sp.freq(let);
//		if (!DBIG_EQUAL(convert(tmpLcat), 0.0))
//		{
//			cerr<<"tmpLcat = "<<tmpLcat<<endl;
//			errorMsg::reportError("error in likelihoodComputation::getLofPos. likelihood is smaller than zero");
//		}
//		
//		//assert(tmpLcat>=0.0);
//		tmp+=tmpLcat;
//	}
////	cout<<"likelihoodComputation::getLofPos: tmp = "; tmp.outputn(cout);	// DEBUG EP
//	if (!(tmp>0.0)){
//		LOG(5,<<"likelihoodComputation::getLofPos: "<< tmp<<endl;);
//		LOG(5,<<"pos = "<< pos <<endl;);
//		tmp = EPSILON;
//		//errorMsg::reportError("likelihoodComputation::getLofPos: likelihood of pos was zero!",1);
//
//	}
//	return tmp;
//}
//
//Vdouble likelihoodComputationGL::getLofPosPerCat(const int pos,
//										   const tree& tr,
//										   const sequenceContainer& sc,
//										   const computePijGam& pi,
//										   const stochasticProcess& sp)
//{
////	with the pi already computed.
//    int numOfCat = sp.categories();
//	Vdouble tmp;
//	tmp.resize(numOfCat*2);
//	for (int i=0; i < numOfCat;++i) {
//		tmp[i] = getLofPos(pos,tr,sc,pi[i],sp)*sp.ratesProb(i);
//		tmp[i+numOfCat] = sp.ratesProb(i);
//	}
//	return tmp;
//}

