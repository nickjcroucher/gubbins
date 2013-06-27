#include <cstdlib>
#include "bb_options.h"
#include "logFile.h"
#include "errorMsg.h"

bb_options::bb_options(int& argc, char *argv[]):
	computeAgainExactTreshold(0.9),
	optimizeBrLenOnStartingTree(true),
	doJoint(true),
	treefile(""),
	reportFile("log.txt"),
	outFile_seq_joint("seq.joint.txt"),
	outFile_seq_marginal("seq.marginal.txt"),
	outFile_prob_joint("prob.joint.txt"),
	outFile_prob_marginal("prob.marginal.txt"),
	seqfile(""),
    distributionName(hom),
	seqOutputFormat(clustal),
    outTreeFileNewick("tree.newick.txt"),
    outTreeFileAncestor("tree.ancestor.txt"),
	boundMethod(both),
	gammaPar(1.0),
	userProvideAlpha(false),
	gammaCategies(8),
	modelName(jtt),
	alphabet_size(20), 
	removeGapsPosition(true),
	useChebyshev(true),
	treeOutFile("TheTree.txt"),
    outPtr(&cout){ 
  	static struct option long_options[] = {{0, 0, 0, 0}}; 
	int option_index = 0;
	int c=0;
	while (c >= 0) {
		c = getopt_long(argc, argv,"a:bc:d:e:fghj:k:m:p:q:R:s:t:ux:y:z:", long_options,&option_index);
    
		switch (c) {
			case 'a': computeAgainExactTreshold=atof(optarg); break;
			case 'b': optimizeBrLenOnStartingTree=false; break;
			case 'c': gammaCategies=atoi(optarg); break;
			case 'd': outFile_prob_joint=optarg; break;
			case 'e': outFile_prob_marginal=optarg; break;
			case 'f': doJoint=false; break;
			case 'g': distributionName=gam; break;
			case 'h' :  {
						cout << "USAGE:	"<<argv[0]<<" [-options] "<<endl;
						cout << usage()<<endl;
						exit (0);
					  }	break;
			case 'j': outFile_seq_joint=optarg; break;
			case 'k': outFile_seq_marginal=optarg; break;
			case 'm': {
				switch (optarg[0]) {
					case 'd': case 'D':  modelName=day;alphabet_size=20; break;
					case 'j': case 'J':  modelName=jtt;alphabet_size=20; break;
					case 'l': case 'L':  modelName=lg;alphabet_size=20; break;
					case 'r': case 'R':  modelName=rev;alphabet_size=20; break;
					case 'w': case 'W':  modelName=wag;alphabet_size=20; break;
					case 'c': case 'C':  modelName=cprev;alphabet_size=20; break;
					case 'a': case 'A':  modelName=aajc;alphabet_size=20; break;
					case 'n': case 'N':  modelName=nucjc;alphabet_size=4; break;
					case 'g': case 'G':  modelName=nucgtr;alphabet_size=4; break;
					case 'e': case 'E':  modelName=empiriCodon;alphabet_size=61; break;
					case 'y': case 'Y':  modelName=nyCodon;alphabet_size=61; break;
					default:modelName=jtt;alphabet_size=20;
					break;
				}
			} break;
			case 'p': {
						userProvideAlpha = true;
						gammaPar=atof(optarg);
						distributionName=gam;

					  } break;
			case 'q': {
						switch (optarg[0]) {
							case 'c': seqOutputFormat=clustal; break;
							case 'f': seqOutputFormat=fasta; break;
							case 'm': seqOutputFormat=molphy; break;
							case 's': seqOutputFormat=mase; break;
							case 'p': seqOutputFormat=phylip; break;
							case 'n': seqOutputFormat=nexus; break;
							default: seqOutputFormat=clustal; break;
						}
					  } break;
			case 'R': reportFile=optarg; break;
			case 's': seqfile=optarg; break;
			case 't': treefile=optarg; break;
			case 'u': useChebyshev=false; break;
			case 'x': outTreeFileNewick=optarg; break;
			case 'y': outTreeFileAncestor=optarg; break;
			case 'z': {
				switch (optarg[0]) {
					case 's': case 'S': boundMethod=sum;  break;
					case 'm': case 'M': boundMethod=max;  break;
					case 'b': case 'B': boundMethod=both;  break;
					default:boundMethod=both;break;
				}
			} break;
			 
			//default: printf ("?? getopt returned character code 0%o ??\n", c);
		} // end of switch c
	} // end of while (c)
	if (seqfile=="") {
	  cout << "USAGE:	"<<argv[0]<<" [-options] "<<endl;
	  //cout << "cat SeqFile |"<<argv[0]<<" [-options]"<<endl <<endl;
	  cout << usage();
	  cout << endl;
	  exit (0);
	}
}


string bb_options::modelNameStr() const
{

	string res = "";
	switch (modelName)
	{
	case day:
		res = "DAY";
		break;
	case jtt:
		res = "JTT";
		break;
	case wag:
		res = "WAG";
		break;
	case lg:
		res = "LG";
		break;
	case nyCodon:
		res = "NY_CODON";
		break;
	case rev:
		res = "REV";
		break;
	case cprev:
		res = "CPREV";
		break;
	case nucjc:
		res = "NUC_JC";
		break;
	case nucgtr:
		res = "NUC_GTR";
		break;
	case aajc:
		res = "AA_JC";
		break;
	case empiriCodon:
		res = "EMPIRICAL_CODON";
		break;
	default:
		errorMsg::reportError("unknown type in  bb_options::modelNameStr");
	}
	return res;

}

