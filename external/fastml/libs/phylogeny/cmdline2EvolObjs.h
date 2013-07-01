// 	$Id: cmdline2EvolObjs.h 5928 2009-02-25 16:30:50Z privmane $	

#ifndef ___CREATESPFROMARGSINFO_H
#define ___CREATESPFROMARGSINFO_H

#include <cstdlib>
#include "amino.h"
#include "nucleotide.h"
#include "codon.h"
#include "sequenceContainer.h"
#include "tree.h"
#include "stochasticProcess.h"
#include "replacementModel.h"
#include "uniDistribution.h"
#include "trivialAccelerator.h"
#include "alphaTrivialAccelerator.h"
#include "chebyshevAccelerator.h"
#include "talRandom.h"
#include "nucJC.h"
#include "aaJC.h"
#include "hky.h"
#include "tamura92.h"
#include "gtrModel.h"
#include "logFile.h"
#include "readDatMatrix.h"
#include "gammaDistribution.h"
#include "recognizeFormat.h"
#include "replacementModelSSRV.h"
#include "stochasticProcessSSRV.h"
#include "someUtil.h"

#define DEFAULT_VALUE_FOR_ALPAH 1.0

template <class args_infoT>
class cmdline2EvolObjs {
private:
  args_infoT _args_info;
public:
  const args_infoT& getArgsInfo(void) {return(_args_info);}
  // constructors
  cmdline2EvolObjs(args_infoT &args_info) : _args_info(args_info) {
    checkParameterConsistancy();
  }
  cmdline2EvolObjs(args_infoT &args_info, bool DontChack) : _args_info(args_info) {
	//    if (!DontChack) checkParameterConsistancy();
  }
  explicit cmdline2EvolObjs(void){}; // do nothing
  void installArgsInfo(args_infoT &args_info){
	_args_info = args_info;
    checkParameterConsistancy();
  }
private:
	void checkParameterConsistancy() {
		if (!_args_info.homogeneous_flag) { // using Gamma ASRV
			if (!_args_info.alpha_given && !_args_info.optimizeAlpha_flag)
				errorMsg::reportError("Must use either 'alpha' or 'optimizeAlpha' when using Gamma ASRV");
		} else {			// using homogeneous rates
			if (_args_info.categories_given ||_args_info.alpha_given || _args_info.optimizeAlpha_given)
				errorMsg::reportError("Can't use 'categories' or 'alpha' or 'optimizeAlpha' with homogeneous rates model");
			// more tests may come here
		}

		// check compatibility of alphabet and model
		if (_args_info.alphabet_arg == 4
			&& !(_args_info.nucjc_given || _args_info.k2p_given || _args_info.tamura92_given || _args_info.gtr_given))
			errorMsg::reportError("Model type is not suitable for nucleotide alphabet");
		if (_args_info.alphabet_arg == 20
			&& (_args_info.nucjc_given || _args_info.k2p_given || _args_info.tamura92_given || _args_info.gtr_given))
			errorMsg::reportError("Model type is not suitable for amino-acid alphabet");

		if (_args_info.nu_given) {
			_args_info.ssrv_flag = true;
		}
	}

public:
	void initializeRandomSeed() {
		if	(_args_info.seed_given) {
			talRandom::setSeed(_args_info.seed_arg);
		}
	}
	void initializeLogFile() {
		myLog::setLog(_args_info.Logfile_arg, _args_info.verbose_arg);
	}

	// NOTE: Unlike other cmdline2*** classes, here a pointer to an allocated obj 
	// is returned and the user is responsible for doing delete. This is because 
	// alphabet is an abstract class, so we can't return it by value
	alphabet* cmdline2Alphabet() {
		alphabet* alphPtr = NULL;
		switch (_args_info.alphabet_arg) 
		{ // allwayes defined, with default
		case 4:	
			alphPtr = new nucleotide; 
			break;
		case 20: 
			alphPtr = new amino; 
			break;
		case 64: case 61: case 60: case 62:
			alphPtr = new codon; 
			break;
		default: errorMsg::reportError("alphabet size not supported");
		}

		// Handle mulAlphabet needed in case we use an SSRV model
		if (_args_info.ssrv_flag) {
			alphabet* mulAlphPtr = new mulAlphabet(alphPtr, _args_info.categories_arg);
			delete alphPtr;
			alphPtr = mulAlphPtr;
		}

		return alphPtr;
	}

	sequenceContainer cmdline2SequenceContainer(const alphabet * const alphPtr) {
		ifstream ins;
		istream* inPtr = &cin;	
		string sequenceFileName(_args_info.sequence_arg);
		if (sequenceFileName != "" && sequenceFileName != "-") {
			ins.open(sequenceFileName.c_str());
			if (! ins.is_open())
				errorMsg::reportError(string("Can not open sequence file ")+sequenceFileName);
			inPtr = &ins;
		}
		istream& in = *inPtr;

		sequenceContainer sc;
		if (!_args_info.ssrv_flag) {
			sc = recognizeFormat::read(in, alphPtr);
		} else {
			sequenceContainer scBase(recognizeFormat::read(in, (static_cast<const mulAlphabet*>(alphPtr))->getBaseAlphabet()));
			sc = sequenceContainer(scBase, alphPtr);
		}
		return sc;
	}
	
	void takeCareOfGaps (sequenceContainer &sc) {
		if (_args_info.gaps_flag) {
			sc.removeGapPositions();
		} else {
			sc.changeGaps2MissingData();
		}
	}

  // NOTE: Unlike other cmdline2*** classes, here a pointer to an allocated obj 
  // is returned and the user is responsible for deleting it. This is because 
  // we need to return a NULL pointer if we are not given a tree
  tree *cmdline2Tree() {
    tree *treePtr = NULL;
    if (_args_info.tree_given) { // did we get a tree
      string treeFileName(_args_info.tree_arg);
      treePtr = new tree(treeFileName);
    }
    return treePtr;
  }

  // NOTE: Unlike other cmdline2*** classes, here a pointer to an allocated obj 
  // is returned and the user is responsible for deleting it. This is because 
  // we need to return a NULL pointer if we are not given a tree
  tree *cmdline2ConstraintTree() {
    tree *constraintTreePtr = NULL;
    if (_args_info.constraint_given) { // did we get a tree
      string constraintTreeFileName(_args_info.constraint_arg);
      constraintTreePtr = new tree(constraintTreeFileName);
    }
    return constraintTreePtr;
  }

  replacementModel *cmdline2ReplacementModel() {
    replacementModel *probModPtr=NULL;
    MDOUBLE ratio =_args_info.ratio_arg;
    MDOUBLE Ap(0.25), Cp(0.25), Gp(0.25), Tp(0.25);
    sscanf(_args_info.ACGprob_arg,"%lf,%lf,%lf", &Ap, &Cp, &Gp);
    Tp=1.0-(Ap+Cp+Gp);

    if (_args_info.day_given) {
      LOG(5,<<"Using Dayhoff replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::dayhoff);
    } else if (_args_info.rev_given) {
      LOG(5,<<"Using rev replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::mtREV24);
    } else if (_args_info.wag_given) {
      LOG(5,<<"Using wag replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::wag);
    } else if (_args_info.cprev_given) {
      LOG(5,<<"Using cprev replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::cpREV45);
    } else if (_args_info.nucjc_given) {
      LOG(5,<<"Using JC for nucleotide"<<endl);
      probModPtr=new nucJC;
    } else if (_args_info.aaJC_given) {
      LOG(5,<<"Using JC for amino acids"<<endl);
      probModPtr=new aaJC;
    } else if ((_args_info.hky_given) || (_args_info.k2p_given)) {
      LOG(5,<<"Using hky replacement matrix"<<endl);
      probModPtr=new hky(Ap,Cp,Gp,Tp,ratio);
    } else if (_args_info.tamura92_given) {
      LOG(5,<<"Using the Tamura 92 replacement matrix"<<endl);
      MDOUBLE theta = Cp+Gp;
      probModPtr=new tamura92(theta, ratio);
    } else if (_args_info.gtr_given) {
      LOG(5,<<"Using the GTR replacement matrix"<<endl);
	  //Vdouble freqs = evaluateCharacterFreq(_sc);
	  Vdouble freqs;
	  freqs.push_back(0.25);
	  freqs.push_back(0.25);
	  freqs.push_back(0.25);
	  freqs.push_back(0.25);
	  probModPtr=new gtrModel(freqs);
    } else if ((_args_info.alphabet_arg == 20) && 
	       (_args_info.modelfile_given)) { // try to read the name as a file name
      LOG(5,<<"Using user supplied replacement matrix from the file "<<_args_info.modelfile_arg<<endl);
      probModPtr=new pupAll(_args_info.modelfile_arg);
    } else { /* default = if (strcmp(_args_info.model_arg,"jtt")==0) */
      probModPtr=new pupAll(datMatrixHolder::jones);
    }
		
    return probModPtr;
  }

  replacementModel *cmdline2ReplacementModelAAOnly() {
    replacementModel *probModPtr=NULL;

    if (_args_info.day_given) {
      LOG(5,<<"Using Dayhoff replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::dayhoff);
    } else if (_args_info.rev_given) {
      LOG(5,<<"Using rev replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::mtREV24);
    } else if (_args_info.wag_given) {
      LOG(5,<<"Using wag replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::wag);
    } else if (_args_info.cprev_given) {
      LOG(5,<<"Using cprev replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::cpREV45);
    } else if (_args_info.aaJC_given) {
      LOG(5,<<"Using JC for amino acids"<<endl);
	  probModPtr=new aaJC;
	} else if (_args_info.modelfile_given) { // try to read the name as a file name
      LOG(5,<<"Using user supplied replacement matrix from the file "<<_args_info.modelfile_arg<<endl);
      probModPtr=new pupAll(_args_info.modelfile_arg);
    } else { /* default = if (strcmp(_args_info.model_arg,"jtt")==0) */
      probModPtr=new pupAll(datMatrixHolder::jones);
    }
		
    return probModPtr;
  }

  bool useGamma()
  {
    return (!_args_info.homogeneous_flag);
  }

  // this function is ment for cases where a "mature" stochastic Process
  // can be produced.  If there is a chance that the user may ask for
  // alpha optimisation use the
  // "cmdline2StochasticProcessThatRequiresAlphaOptimization" version
  // instead
  inline stochasticProcess *cmdline2StochasticProcess() {
    distribution *distP = NULL;
    if (useGamma()) {
      if (_args_info.alpha_given) 
	distP =  new gammaDistribution(_args_info.alpha_arg,_args_info.categories_arg);
      else
	errorMsg::reportError("Can not create stochastic process with ASRV if no alpha is given, when working without alpha optimization");
      LOG(5,<<"Using Gamma ASRV with "<<_args_info.categories_arg<<" bins"<<endl);
    } else {
      distP =  new uniDistribution;
      LOG(5,<<"Using uniform rates"<<endl);
    }
    stochasticProcess *spPtr = cmdline2StochasticProcessInternal(*distP);
    if (distP) delete distP;
    return(spPtr);
  }
  
  // Assuming that the user asked to optimize Alpha (by bestAlphaAndBBL)
  inline stochasticProcess *cmdline2StochasticProcessThatRequiresAlphaOptimization () {
    distribution *distP = NULL;
    if (!_args_info.optimizeAlpha_given)
      errorMsg::reportError("Can't use function cmdline2StochasticProcessThatRequiresAlphaOptimization if the optimizeAlpha flag was not turned on - please inform the programmer of this error.");
    // else
    if (_args_info.alpha_given) 
      distP =  new gammaDistribution(_args_info.alpha_arg,_args_info.categories_arg);
    else
      distP =  new gammaDistribution(DEFAULT_VALUE_FOR_ALPAH,_args_info.categories_arg);
    LOG(5,<<"Using Gamma ASRV with "<<_args_info.categories_arg<<" bins"<<endl);
    stochasticProcess *spPtr = cmdline2StochasticProcessInternal(*distP);
    if (distP) delete distP;
    return(spPtr);
  }

  inline stochasticProcess *cmdline2HomogenuisStochasticProcess() {
    uniDistribution dist;  
    LOG(5,<<"Creating homogeneous rate based stochastic Process "<<endl);
    return (cmdline2StochasticProcessInternal(dist));
  }

  inline stochasticProcess cmdline2HomogenuisStochasticProcessAAOnly() {
    uniDistribution dist;  
    LOG(5,<<"Creating homogeneous rate based stochastic Process "<<endl);
    return (cmdline2StochasticProcessInternalAAOnly(dist));
  }

  inline stochasticProcess *cmdline2StochasticProcessSafe()
  {
	if (_args_info.homogeneous_flag) {
	  return cmdline2StochasticProcess();
	} else {					// we use Gamma
	  if (_args_info.optimizeAlpha_flag) {
		return cmdline2StochasticProcessThatRequiresAlphaOptimization();
	  } else if (_args_info.alpha_given) {
		return cmdline2StochasticProcess();
	} else {
	  errorMsg::reportError("Gamma ASRV requiers either --alpha or --optimizeAlpha or both.",1);
	  }
	}
	exit(1);					// should never be reached
  }
  
private:
	stochasticProcess *cmdline2StochasticProcessInternal(distribution& dist) {
		replacementModel *probModPtr=NULL;
		pijAccelerator *pijAcc=NULL;
		MDOUBLE ratio =_args_info.ratio_arg;
		MDOUBLE Ap(0.25), Cp(0.25), Gp(0.25), Tp(0.25);
		sscanf(_args_info.ACGprob_arg,"%lf,%lf,%lf", &Ap, &Cp, &Gp);
		Tp=1.0-(Ap+Cp+Gp);

		if (_args_info.day_given) {
			LOG(5,<<"Using Dayhoff replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::dayhoff);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.rev_given) {
			LOG(5,<<"Using rev replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::mtREV24);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.wag_given) {
			LOG(5,<<"Using wag replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::wag);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.cprev_given) {
			LOG(5,<<"Using cprev replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::cpREV45);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.nucjc_given) {
			LOG(5,<<"Using JC for nucleotide"<<endl);
			probModPtr=new nucJC;
			pijAcc = new trivialAccelerator(probModPtr);
		} else if (_args_info.aaJC_given) {
			LOG(5,<<"Using JC for amino acids"<<endl);
			probModPtr=new aaJC;
			pijAcc = new trivialAccelerator(probModPtr);
		} else if ((_args_info.hky_given) || (_args_info.k2p_given)) {
			LOG(5,<<"Using hky replacement matrix"<<endl);
			probModPtr=new hky(Ap,Cp,Gp,Tp,ratio);
			pijAcc = new trivialAccelerator(probModPtr);
		} else if (_args_info.tamura92_given) {
			LOG(5,<<"Using the Tamura 92 replacement matrix"<<endl);
			MDOUBLE theta = Cp+Gp;
			probModPtr=new tamura92(theta, ratio);
			pijAcc = new trivialAccelerator(probModPtr);
		} else if (_args_info.gtr_given) {
			LOG(5,<<"Using the GTR replacement matrix"<<endl);
			//Vdouble freqs = evaluateCharacterFreq(_sc);
			Vdouble freqs;
			freqs.push_back(0.25);
			freqs.push_back(0.25);
			freqs.push_back(0.25);
			freqs.push_back(0.25);
			probModPtr=new gtrModel(freqs);
			pijAcc = new trivialAccelerator(probModPtr);
		} else if ((_args_info.alphabet_arg == 20) && 
				   (_args_info.modelfile_given)) { // try to read the name as a file name
			LOG(5,<<"Using user supplied replacement matrix from the file "<<_args_info.modelfile_arg<<endl);
			probModPtr=new pupAll(_args_info.modelfile_arg);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else { /* default = if (strcmp(_args_info.model_arg,"jtt")==0) */
			probModPtr=new pupAll(datMatrixHolder::jones);
			pijAcc = new chebyshevAccelerator(probModPtr);
		}		

		stochasticProcess *spPtr = NULL;
		if (!_args_info.ssrv_flag) {
			spPtr = new stochasticProcess(&dist, pijAcc);

		} else {
			// Using a Site-Specific Rate Variation model
			replacementModelSSRV probModSsrv(&dist,probModPtr,_args_info.nu_arg);
			if (pijAcc) delete pijAcc;
			pijAcc = new trivialAccelerator(&probModSsrv);
			spPtr = new stochasticProcessSSRV(pijAcc);
			LOG(5,<<"cmdline2StochasticProcessInternal: Created stochasticProcessSSRV"<<endl);
		}

		// if rate is given in input, set it.
		if (_args_info.inputRate_given)
			spPtr->setGlobalRate(_args_info.inputRate_arg);

		if (probModPtr) delete probModPtr;
		if (pijAcc) delete pijAcc;
		return spPtr;
	}

	stochasticProcess cmdline2StochasticProcessInternalAAOnly(distribution& dist) {
		replacementModel *probModPtr=NULL;
		pijAccelerator *pijAcc=NULL;

		if (_args_info.day_given) {
			LOG(5,<<"Using Dayhoff replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::dayhoff);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.rev_given) {
			LOG(5,<<"Using rev replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::mtREV24);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.wag_given) {
			LOG(5,<<"Using wag replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::wag);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.cprev_given) {
			LOG(5,<<"Using cprev replacement matrix"<<endl);
			probModPtr=new pupAll(datMatrixHolder::cpREV45);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else if (_args_info.aaJC_given) {
			LOG(5,<<"Using JC for amino acids"<<endl);
			probModPtr=new aaJC;
			pijAcc = new trivialAccelerator(probModPtr);
		} else if (_args_info.modelfile_given) { // try to read the name as a file name
			LOG(5,<<"Using user supplied replacement matrix from the file "<<_args_info.modelfile_arg<<endl);
			probModPtr=new pupAll(_args_info.modelfile_arg);
			pijAcc = new chebyshevAccelerator(probModPtr);
		} else { /* default = if (strcmp(_args_info.model_arg,"jtt")==0) */
			probModPtr=new pupAll(datMatrixHolder::jones);
			pijAcc = new chebyshevAccelerator(probModPtr);
		}		
		stochasticProcess sp(&dist, pijAcc);

		// if rate is given in input, set it.
		// if (_args_info.inputRate_given)
		// sp.setGlobalRate(_args_info.inputRate_arg);

		if (probModPtr) delete probModPtr;
		if (pijAcc) delete pijAcc;
		return sp;
	}

public:
  stochasticProcess cmdline2ExactGammaStochasticProcess() {
    uniDistribution dist;  
    LOG(5,<<"Creating exact Gamma based stochastic Process "<<endl);
    if(!_args_info.alpha_given)
      errorMsg::reportError("Using exact Gamma requires alpha to be set");
    pupAll *probModPtr=NULL;
    //		pijAccelerator *pijAcc=NULL;
    alphaTrivialAccelerator *pijAcc=NULL;
 
    if (_args_info.day_given) {
      LOG(5,<<"Using Dayhoff replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::dayhoff);
      pijAcc = new alphaTrivialAccelerator(probModPtr,_args_info.alpha_arg);
    } else if (_args_info.rev_given) {
      LOG(5,<<"Using rev replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::mtREV24);
      pijAcc = new alphaTrivialAccelerator(probModPtr,_args_info.alpha_arg);
    } else if (_args_info.wag_given) {
      LOG(5,<<"Using wag replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::wag);
      pijAcc = new alphaTrivialAccelerator(probModPtr,_args_info.alpha_arg);
    } else if (_args_info.cprev_given) {
      LOG(5,<<"Using cprev replacement matrix"<<endl);
      probModPtr=new pupAll(datMatrixHolder::cpREV45);
      pijAcc = new alphaTrivialAccelerator(probModPtr,_args_info.alpha_arg);
    } else if ((_args_info.alphabet_arg == 20) && 
	       (_args_info.modelfile_given)) { // try to read the name as a file name
      LOG(5,<<"Using user supplied replacement matrix from the file "<<_args_info.modelfile_arg<<endl);
      probModPtr=new pupAll(_args_info.modelfile_arg);
      pijAcc = new alphaTrivialAccelerator(probModPtr,_args_info.alpha_arg);
    } else if (_args_info.nucjc_given || 
	       _args_info.aaJC_given || 
	       _args_info.hky_given || 
	       _args_info.k2p_given || 
	       _args_info.tamura92_given || 
	       _args_info.gtr_given) {
      errorMsg::reportError("Exact Gamma stochastic process only works with pupAll model");
    } else { /* default = if (strcmp(_args_info.model_arg,"jtt")==0) */
      probModPtr=new pupAll(datMatrixHolder::jones);
      pijAcc = new alphaTrivialAccelerator(probModPtr,_args_info.alpha_arg);
    }		
    stochasticProcess sp(&dist, pijAcc);
 		
    // if rate is given in input, set it.
    if (_args_info.inputRate_given)
      sp.setGlobalRate(_args_info.inputRate_arg);
 
    if (probModPtr) delete probModPtr;
    if (pijAcc) delete pijAcc;
    return sp;
  }
 
public:
  // NOTE: the user must check:
  // if the returned stream is an ofstream object (an actual file) it should be deleted
  // if the returned stream is an ostream object (cout) do nothing
  ostream *cmdline2OutputStream() {
    ostream *outPtr;
    string outFileName(_args_info.outputfile_arg);
    if (outFileName == "") outFileName="-";
    if (outFileName == "-") {
      outPtr = &cout;
    } else {
      outPtr = new ofstream(outFileName.c_str());
      if (!outPtr->good()) errorMsg::reportError(string("Can't open for writing the file ")+outFileName);
    }
    return outPtr;
  }

  // NOTE: the user must check:
  // if the returned stream is an ofstream object (an actual file) it should be deleted
  // if the returned stream is an ostream object (cout) do nothing
  ostream *cmdline2TreeOutputStream() {
    ostream *outPtr;
    string outFileName(_args_info.treeoutputfile_arg);
    if (outFileName == "") outFileName="-";
    if (outFileName == "-") {
      outPtr = &cout;
    } else {
      outPtr = new ofstream(outFileName.c_str());
      if (!outPtr->good()) errorMsg::reportError(string("Can't open for writing the file ")+outFileName);
    }
    return outPtr;
  }

  void consistencyCheck (tree *treePtr, tree *constraintTreePtr) {
    if (treePtr!=NULL) {
      if (constraintTreePtr !=NULL) {
	/*				constraints c1(*constraintTreePtr);
					c1.setTree(*treePtr);
					if (!c1.fitsConstraints()){
					LOG(1,<<"Input tree does not fit constraints!"<<endl);
					LOGDO(1,c1.outputMissingClads(myLog::LogFile()));
					errorMsg::reportError("Please enter a starting tree that fits the constraints");
					}
	*/			}
    }
  }

public:
	// Read from file the posterior distribution of rates for each sequence site
	VVdoubleRep cmdline2PosteriorRates() {
		if (!_args_info.posteriorRates_given)
			errorMsg::reportError("cmdline2EvolObjs::cmdline2PosteriorRates: This method shouldn't be used if --posteriorRates was not given");
		ifstream in(_args_info.posteriorRates_arg);
		if (!in.is_open())
			errorMsg::reportError(string("Can not open sequence file ")+string(_args_info.posteriorRates_arg));

		string line, number, rest; // For splitting the line into separate numbers
		VdoubleRep posterior(_args_info.categories_arg);  // Current line
		VVdoubleRep posteriorRates; // Accumulate all lines
		getline(in, line);

		// Each loop reads one line of numbers
		while (in) {
			// split line into numbers
			for(int cat=0; cat<_args_info.categories_arg; ++cat) {
				splitString2(line, " ", number, rest);
				if (number.size() == 0)
					errorMsg::reportError(string("cmdline2EvolObjs::cmdline2PosteriorRates: Bad line with too few numbers in file ")
										  +_args_info.posteriorRates_arg+": "+line);
				posterior[cat] = atof(number.c_str());
			}
			posteriorRates.push_back(posterior);
			getline(in, line);
		}
		return posteriorRates;
	}
};

#endif
