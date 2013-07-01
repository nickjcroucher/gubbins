// $Id: evaluateCharacterFreq.h 3895 2008-04-21 07:38:32Z itaymay $

#ifndef __Evaluate_Character_Freq_h
#define __Evaluate_Character_Freq_h

#include <iostream>
using namespace std;

#include "sequenceContainer.h"
#include "definitions.h"

vector<MDOUBLE> sumAlphabetCounts(const sequenceContainer & sc);
vector<MDOUBLE> evaluateCharacterFreq(const sequenceContainer & sc);
VVdouble evaluateCharacterFreqOneForEachGene(const vector<sequenceContainer> & scVec);
vector<MDOUBLE> evaluateCharacterFreqBasedOnManyGenes(const vector<sequenceContainer> & scVec);

void changeCountsToFreqs(vector<MDOUBLE>& charFreq);
void makeSureNoZeroFreqs(vector<MDOUBLE> & charFreq);

//returns the number of each character in each position
void getCharacterCounts(const sequenceContainer & sc, VVint& counts4pos);
//returns the number of different character types in each position
void getCharacterType4pos(const sequenceContainer & sc, Vint& charactersType4pos);
//returns the distribution of the different character types in each position along the whole alignment
void getCharacterTypeDistribution(const sequenceContainer & sc, Vint& charactersTypeDist);
#endif
