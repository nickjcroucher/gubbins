#include "mainbb.h"
#include "logFile.h"


int main(int argc, char* argv[]) {
	myLog::setLog("",10);
	mainbb mainbb1(argc,argv);
	return 0;
}

/*
//------------------------------------------------


#include "bbAlg.h"
#include "sequenceDataDiff.h"
sequenceContainer main1(const string& seqFile,
				   char format,
		  const string& treeFile,
		  const string& reportFileName,
		  const string& ancestralSequencesFileName,
		  const MDOUBLE alpha,
		  const int categor,
		  time_t& timeTaken,
		  clock_t& ctimeTaken,
		  const MDOUBLE recalculateExactVal); //0 never recalculate...

int veryMainLysSmallCheck() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\seq\\lys6\\junk\\seqF1.txt";
	const string treeFile1 = "C:\\tal\\seq\\lys6\\junk\\tree.txt";
	const string treeFile2 = "C:\\tal\\seq\\lys6\\junk\\tree.txt";
	const string reportFileHom = "C:\\tal\\seq\\lys6\\junk\\tmp\\reportFileHom.txt";
	const string reportFileGam = "C:\\tal\\seq\\lys6\\junk\\tmp\\reportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\seq\\lys6\\junk\\tmp\\reportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\seq\\lys6\\junk\\tmp\\ancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\seq\\lys6\\junk\\tmp\\ancstralSeqHom.txt";
	time_t time1;
	time_t time2;
clock_t ctime1;
clock_t ctime2;

	sequenceContainer sd1 = main1(seqFile,'m',treeFile1,reportFileGam,ancstralSeqGam,0.924884,4,time1,ctime1,0); // gam
	sequenceContainer sd2 = main1(seqFile,'m',treeFile2,reportFileHom,ancstralSeqHom,-3,1,time2,ctime2,0); // hom
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time2<<endl;
	out<<" time taken for gam was: "<<time1<<endl;
	out.close();
	return 0;
}

int veryMainLys() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\lys71.ngap.mase";
	const string treeFile1 = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\treehom.txt";
	const string treeFile2 = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\treegam.txt";
	const string reportFileHom = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\reportFileHom.txt";
	const string reportFileGam = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\reportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\reportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\ancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\activeProjects\\ancbb\\seq\\lys71\\ancstralSeqHom.txt";
	time_t time1;
	time_t time2;
	clock_t ctime1;
	clock_t ctime2;
	sequenceContainer sd1 = main1(seqFile,'m',treeFile1,reportFileGam,ancstralSeqGam,0.924884,4,time1,ctime1,0); // gam
	sequenceContainer sd2 = main1(seqFile,'m',treeFile2,reportFileHom,ancstralSeqHom,-3,1,time2,ctime2,0); // hom
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time2<<endl;
	out<<" time taken for gam was: "<<time1<<endl;
	out.close();
	return 0;
}

int veryMainCo1() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\co1.ngap.aln";
	const string treeFile1 = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\treehom.txt";
	const string treeFile2 = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\treegam.txt";
	const string reportFileHom = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\reportFileHom.txt";
	const string reportFileGam = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\reportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\reportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\ancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\activeProjects\\ancbb\\seq\\co1\\ancstralSeqHom.txt";
	time_t time1;
	time_t time2;
	clock_t ctime1;
	clock_t ctime2;
	sequenceContainer sd1 = main1(seqFile,'a',treeFile1,reportFileGam,ancstralSeqGam,0.257432,4,time1,ctime1,0); // gam
	sequenceContainer sd2 = main1(seqFile,'a',treeFile2,reportFileHom,ancstralSeqHom,-3,1,time2,ctime2,0); // hom
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time2<<endl;
	out<<" time taken for gam was: "<<time1<<endl;
	out.close();
	return 0;
}

int veryMainCo2() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\co2ngap.aln";
	const string treeFile1 = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\treehom.txt";
	const string treeFile2 = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\treegam.txt";
	const string reportFileHom = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\reportFileHom.txt";
	const string reportFileGam = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\reportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\reportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\ancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\activeProjects\\ancbb\\seq\\co2\\ancstralSeqHom.txt";
	time_t time1;
	time_t time2;
	clock_t ctime1;
	clock_t ctime2;
	sequenceContainer sd1 = main1(seqFile,'a',treeFile1,reportFileGam,ancstralSeqGam,0.476490,4,time1,ctime1,0); // gam
	sequenceContainer sd2 = main1(seqFile,'a',treeFile2,reportFileHom,ancstralSeqHom,-3,1,time2,ctime2,0); // hom
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time2<<endl;
	out<<" time taken for gam was: "<<time1<<endl;
	out.close();
	return 0;
}

int veryMainOpsin() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\opsin.mase";
	const string treeFile1 = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\treehom.txt";
	const string treeFile2 = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\treegam.txt";
	const string reportFileHom = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\reportFileHom.txt";
	const string reportFileGam = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\reportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\reportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\ancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\activeProjects\\ancbb\\seq\\opsin\\ancstralSeqHom.txt";
	time_t time1;
	time_t time2;
	clock_t ctime1;
	clock_t ctime2;
	sequenceContainer sd1 = main1(seqFile,'m',treeFile1,reportFileGam,ancstralSeqGam,0.331405,4,time1,ctime1,0); // gam
	sequenceContainer sd2 = main1(seqFile,'m',treeFile2,reportFileHom,ancstralSeqHom,-3,1,time2,ctime2,0); // hom
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time2<<endl;
	out<<" time taken for gam was: "<<time1<<endl;
	out.close();
	return 0;
}


int veryMainSteroid() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\noGaps.mase";
	const string treeFile1 = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\treehom.txt";
	const string treeFile2 = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\treegam.txt";
	const string reportFileHom = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\reportFileHom.txt";
	const string reportFileGam = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\reportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\reportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\ancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\ancstralSeqHom.txt";
	time_t time1;
	time_t time2;
	sequenceContainer sd1 = main1(seqFile,'m',treeFile1,reportFileGam,ancstralSeqGam,1.534586,4,time1,0); // gam
	sequenceContainer sd2 = main1(seqFile,'m',treeFile2,reportFileHom,ancstralSeqHom,-3,1,time2,0); // hom
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time2<<endl;
	out<<" time taken for gam was: "<<time1<<endl;
	out.close();
	return 0;
}


int veryMainSteroid() {// the non command line version for debugging and checking.
	const string seqFile = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\dataPreperation\\B4remGap\\ster73.snames.correct.ngap.aln";
	const string treeFile1 ="C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\buildingTree\\topologyHom.ph";
	const string treeFile2 ="C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\buildingTree\\topologyGam.ph";

	
	
	const string reportFileHom = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\NreportFileHom.txt";
	const string reportFileGam = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\NreportFileGam.txt";
	const string reportFileDiffAndTime = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\NreportFileDif.txt";
	const string ancstralSeqGam = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\NancstralSeqGam.txt";
	const string ancstralSeqHom = "C:\\tal\\activeProjects\\ancbb\\seq\\steroid\\NancstralSeqHom.txt";
	time_t time1;
	time_t time2;
	clock_t ctime1;
	clock_t ctime2;
	sequenceContainer sd1 = main1(seqFile,'a',treeFile1,reportFileHom,ancstralSeqHom,-600,1,time1,ctime1,0); // hom
	sequenceContainer sd2 = main1(seqFile,'a',treeFile2,reportFileGam,ancstralSeqGam,1.29,4,time2,ctime2,0); // gam
	sequenceDataDiff sequenceDataDiff1f(&sd1,&sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out<<" time taken for hom was: "<<time1<<endl;
	out<<" time taken for gam was: "<<time2<<endl;
	out<<" ctime taken for hom was: "<<ctime1<<endl;
	out<<" ctime taken for gam was: "<<ctime2<<endl;
	out.close();
	return 0;
}

MDOUBLE totalBranchLengh(const tree& t1) {
	MDOUBLE sum=0;
	vector<tree::nodeP> vec;
	t1.getAllNodes(vec,t1.getRoot());
	for (int i=0; i< vec.size(); ++i) {
		if (vec[i]->father != NULL) sum += vec[i]->dis2father();
		cerr<<sum<<"  "<<vec[i]->dis2father()<<endl;
	}
	return sum;
}




*/
#include "sequenceDataDiff.h"
#include "amino.h"
#include <ctime>
#include "recognizeFormat.h"
#include "uniDistribution.h"
#include "gammaDistribution.h"
#include "replacementModel.h"
#include "readDatMatrix.h"
#include "chebyshevAccelerator.h"
#include "bbAlg.h"
/*
sequenceContainer main1(const string& seqFile,
				   char format,
		  const string& treeFile,
		  const string& reportFileName,
		  const string& ancestralSequencesFileName,
		  const MDOUBLE alpha,
		  const int categor,
		  time_t& timeTaken,
		  clock_t& ctimeTaken,
		  const MDOUBLE recalculateExactVal) { // gamma distribution

	alphabet* _alph = new amino;
    ifstream f(seqFile.c_str());
	sequenceContainer original = recognizeFormat::read(f,_alph);;
	tree t1(treeFile); // with sequence data
//	t1.multipleAllBranchesByFactor(10);
	// stochastic process:

//	cerr<<" total br-len is:"<<totalBranchLengh(t1)<<endl;
//	return *sd;


	distribution *dist1 = NULL;
	if (categor ==1 ) dist1 = new uniDistribution; 
	else dist1 = new gammaDistribution(alpha,categor); 

	replacementModel *probMod=new pupAll(datMatrixHolder::jones);
	pijAccelerator *pijAcc1 = new chebyshevAccelerator(probMod); 

//	replacementModel *probMod1=new nucJC;
//	replacementModel *probMod1=new pupJTT;
//	pijAccelerator *pijAcc1= new chebyshevAccelerator(probMod1);
//	pijAccelerator *pijAcc1= new trivialAccelerator(probMod1);
	stochasticProcess* _s1 = new stochasticProcess(dist1, pijAcc1);
	bbAlg bbAlg1(t1,*_s1,original,bbAlg::both,reportFileName,recalculateExactVal);//computeAgainExactTreshold
//	bbAlg bbAlg1(&t1,_s1,bbAlg::sum,0);//computeAgainExactTreshold
//	bbAlg bbAlg1(&t1,_s1,bbAlg::max,0);//computeAgainExactTreshold
	time_t time1,time2;
	clock_t ctime1, ctime2;
	time(&time1);
	ctime1 = clock();
	cerr<<"starting time is: "<<time1<<endl;
	cerr<<"starting clock is: "<<ctime1<<endl;
	MDOUBLE res = bbAlg1.bbReconstructAllPositions(original);
	time(&time2);
	ctime2 = clock();
	cerr<<"ending time is: "<<time2<<endl;
	cerr<<"ending clock is: "<<ctime2<<endl;
	timeTaken=time2-time1;
	ctimeTaken=ctime2-ctime1;

	ofstream outi;
	outi.open(reportFileName.c_str(),ios::app);
	outi<<" the likelihood of the reconstruction is:"<<res<<endl;
	outi.close();
	sequenceContainer recS= bbAlg1.fromAncestralSequenceToSeqData();

	delete pijAcc1;
	delete dist1;
	return recS;
}
*/
/*
int mainNoCommandLine() {

//	veryMainLysSmallCheck(); // just to check that everything is working...
//	veryMainLys();
//	veryMainCo1();
//	veryMainCo2();
//	veryMainOpsin();
	veryMainSteroid();
	return 0;
}
//	const string seqFile = "C:\\tal\\seq\\lys6\\junk\\seq.txt";
//	const string treeFile = "C:\\tal\\seq\\lys6\\junk\\tree.txt";
//	const string seqFile = "C:\\tal\\seq\\lys6\\seq.txt";
//	const string treeFile = "C:\\tal\\seq\\lys6\\tree.txt";
//	main1(seqFile,treeFile,-3,1,time1);// hom

*/

//int main() {
int FindDifferencesBetween2SequenceContainerFiles() {
	const string seqFile1 = "D:\\tal\\yaep15\\fastml2.01\\originalDataForPaper\\seq_joint.txt";
	const string seqFile2 = "D:\\tal\\yaep15\\fastml2.01\\originalDataForPaper\\seq_marginal.txt";
	const string reportFileDiffAndTime = "D:\\tal\\yaep15\\fastml2.01\\originalDataForPaper\\reportFileDif.txt";

	alphabet* _alph = new amino;
    ifstream f(seqFile1.c_str());
	sequenceContainer sd1 = recognizeFormat::read(f,_alph);
	f.close();

    ifstream f2(seqFile2.c_str());
	sequenceContainer sd2 = recognizeFormat::read(f2,_alph);
	f2.close();

	sequenceDataDiff sequenceDataDiff1f(sd1,sd2);
	sequenceDataDiff1f.computeDifferences();
	ofstream outdiff(reportFileDiffAndTime.c_str(),ios::app);
	sequenceDataDiff1f.printDiff(outdiff);
	outdiff.close();
	ofstream out;
	out.open(reportFileDiffAndTime.c_str(),ios::app);
	out.close();
	return 0;
}