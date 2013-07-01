#if !defined ___BB__OPTION__T__
#define ___BB__OPTION__T__


#ifndef __STDC__
#define __STDC__ 1
#include "getopt.h"
#undef __STDC__
#else
#include "getopt.h"
#endif

#include "definitions.h"
#include <iostream>
#include <fstream>
using namespace std;

class bb_options {
public:
	MDOUBLE computeAgainExactTreshold;
	mutable bool optimizeBrLenOnStartingTree;
	bool doJoint;
	string treefile;
	string seqfile;
	enum SeqFileFormat {mase,clustal,fasta,molphy,phylip,nexus};
	SeqFileFormat seqOutputFormat;
	string treeOutFile;
	bool userProvideAlpha;
	enum distributionsNames {hom,gam};
	distributionsNames distributionName;
	enum boundMethods {max,sum,both};
	boundMethods boundMethod;
	bool verbose; // if true: print starting tree to the file: start_tree
//	tree::TREEformats outputFormat;
	enum modelNameOptions {day,jtt,lg,rev,wag,cprev,nucjc,nucgtr,aajc,nyCodon,empiriCodon};
	modelNameOptions modelName;
	int alphabet_size;
	bool removeGapsPosition;
	bool useChebyshev;
	string outTreeFileNewick;
	string outTreeFileAncestor;
	string outFile_prob_joint;
	string outFile_prob_marginal;
	string outFile_seq_joint;
	string outFile_seq_marginal;

	MDOUBLE gammaPar;
	int gammaCategies;
	string reportFile;
private:
	ostream* outPtr;
	ofstream out_f;
public:
	ostream& out() const {return *outPtr;};
	string modelNameStr() const;
	explicit bb_options(int& argc, char *argv[]);
};

#include "bb_options_list.h"
#include <string>
using namespace std;
static const string usege_splash_screen() {
	string tmp = usage();
	return tmp;
};


#endif

