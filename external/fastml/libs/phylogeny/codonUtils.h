#ifndef CODON_UTILS_H
#define CODON_UTILS_H

#include <iostream>
#include "nucleotide.h"
#include "codon.h"
#include "amino.h"
#include "logFile.h"
#include "fastaFormat.h"
#include "clustalFormat.h"
#include "recognizeFormat.h"
#include "someUtil.h"
#include "definitions.h"
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "wYangModel.h"
#include "evaluateCharacterFreq.h"
#include "geneticCodeHolder.h"
#include "codon.h"
using namespace std;


void printHelp();
void checkInputSeqLength(string codonFile);
sequenceContainer convertCodonToAmino(sequenceContainer &codonSc,codon *codonAlph);
vector<vector<int> > create7ColorValues();
void outToRasmolFile(string fileName,vector<int>& color4Site);

void normalizeMatrices(vector<stochasticProcess> & spVec,const distribution * forceDistr);

Vdouble freqCodonF3x4(const sequenceContainer &nucSc,codon *coAlph);

void kaks2Color(const Vdouble & kaksVec,const Vdouble &lowerBoundV, 
		const sequence & refSeq, string fileName,codon *co);

#endif
