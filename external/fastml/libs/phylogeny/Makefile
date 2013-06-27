#! /usr/local/bin/gmake
# $Id: Makefile 6111 2009-04-26 14:08:30Z cohenofi $
# makfile for yaep5



# use LOGREP=t or DOUBLEREP=t to activate logRep or doubleRep respectively (or setenv DOUBLEREP in the shell)
#DOUBLEREP=t
#LOGREP=t


Libsources= AddLog.cpp NNiProp.cpp NNiSep.cpp Nni.cpp aaJC.cpp				\
   allTrees.cpp allTreesSeparateModel.cpp alphabet.cpp amino.cpp			\
   bestAlpha.cpp bestAlphaManyTrees.cpp bestHKYparam.cpp bootstrap.cpp			\
   bblEM.cpp bblEMfixRoot.cpp bblEMProprtional.cpp bblEMSeperate.cpp	\
   chebyshevAccelerator.cpp clustalFormat.cpp codon.cpp codonJC.cpp			\
   computeCounts.cpp computeDownAlg.cpp computeMarginalAlg.cpp				\
   computePijComponent.cpp computeUpAlg.cpp computeUpAlgFactors.cpp			\
   computeSubstitutionCounts.cpp											\
   computePosteriorExpectationOfSubstitutions.cpp							\
   computePosteriorExpectationOfSubstitutions_nonReversibleSp.cpp			\
   ConversionUtils.cpp countTableComponent.cpp datMatrixHolder.cpp distanceTable.cpp	\
   distribution.cpp errorMsg.cpp evaluateCharacterFreq.cpp				\
   fastStartTree.cpp fastaFormat.cpp findRateOfGene.cpp					\
   fromCountTableComponentToDistance.cpp fromCountTableComponentToDistancefixRoot.cpp 	\
   fromCountTableComponentToDistanceProp.cpp fromQtoPt.cpp				\
   generalGammaDistributionFixedCategories.cpp gammaDistribution.cpp gammaUtilities.cpp	\
   generalGammaDistribution.cpp getRandomWeights.cpp goldmanYangModel.cpp		\
   granthamChemicalDistances.cpp hky.cpp khTest.cpp likeDist.cpp likeDistfixRoot.cpp \
   likeDistProp.cpp likelihoodComputation.cpp						\
   likelihoodComputationFactors.cpp logFile.cpp maseFormat.cpp				\
   molphyFormat.cpp nexusFormat.cpp nj.cpp njConstrain.cpp				\
   nucJC.cpp nucleotide.cpp numRec.cpp Parameters.cpp phylipFormat.cpp			\
   pijAccelerator.cpp readDatMatrix.cpp readTree.cpp recognizeFormat.cpp		\
   replacementModel.cpp searchStatus.cpp seqContainerTreeMap.cpp			\
   sequence.cpp sequenceContainer.cpp simulateTree.cpp					\
   siteSpecificRate.cpp someUtil.cpp split.cpp splitMap.cpp				\
   splitTreeUtil.cpp stochasticProcess.cpp suffStatComponent.cpp			\
   talRandom.cpp tree.cpp treeIt.cpp treeUtil.cpp uniDistribution.cpp			\
   uniformDistribution.cpp cmdline2EvolObjs.cpp						\
   generalGammaDistributionLaguerre.cpp gammaDistributionLaguerre.cpp GLaguer.cpp					\
   givenRatesMLDistance.cpp distanceBasedSeqs2Tree.cpp					\
   posteriorDistance.cpp pairwiseGammaDistance.cpp doubleRep.cpp			\
   logRep.cpp indel.cpp indelModel.cpp mulAlphabet.cpp					\
   replacementModelSSRV.cpp stochasticProcessSSRV.cpp bestAlphaAndNu.cpp		\
   C_evalParamUSSRV.cpp matrixUtils.cpp betaOmegaDistribution.cpp			\
   betaUtilities.cpp betaDistribution.cpp geneticCodeHolder.cpp				\
   samplingSequences.cpp bblEM2USSRV.cpp bestParamUSSRV.cpp				\
   likeDist2USSRV.cpp ussrvModel.cpp likelihoodComputation2USSRV.cpp			\
   fromCountTableComponentToDistance2USSRV.cpp normalDist.cpp				\
   tamura92.cpp bestTamura92param.cpp phylipSequentialFormat.cpp			\
   simulateCodonsJumps.cpp													\
   simulateJumpsAbstract.cpp												\
   ssrvDistanceSeqs2Tree.cpp multipleStochasticProcess.cpp distributionPlusInvariant.cpp\
   gammaDistributionFixedCategories.cpp  generalGammaDistributionPlusInvariant.cpp  gammaDistributionPlusInvariant.cpp   \
   distributionPlusCategory.cpp simulateJumps.cpp computeJumps.cpp seqeuncesFilter.cpp 					\
   optGammaMixtureLS.cpp mixtureDistribution.cpp  suffStatGammaMixture.cpp GamMixtureOptimizer.cpp optGammaMixtureEM.cpp gainLossAlphabet.cpp \
   wYangModel.cpp codonUtils.cpp likelihoodComputation2Codon.cpp likeDist2Codon.cpp unObservableData.cpp likelihoodComputationGL.cpp	\
	threeStateModel.cpp threeStateAlphabet.cpp betaDistributionFixedCategories.cpp betaDistributionFixedCategoriesWithOmegaUniform.cpp	\
	bblEM2codon.cpp bestAlphaAndK.cpp fromCountTableComponentToDistance2Codon.cpp\
	gtrModel.cpp bestGtrModelParams.cpp

# do not use: fromInstructionFile.cpp, simulateSequnce.cpp  split.save.cpp


# LibCsources= cmdline.c
# LibCsources += getopt.c getopt1.c 

EXEC = 
#TEST_EXEC_SUB = split_test splitMap_test bootstrap_test
TEST_EXEC = $(addprefix tests/,$(TEST_EXEC_SUB))
LIB = libEvolTree.a
DEBUGLIB = $(LIB:.a=Debug.a)
DOUBLEREPLIB = $(LIB:.a=DoubleRep.a)


#CC=g++
CXX=g++
CC=$(CXX)

#requres 2.13, but may work with 2.11
GENGETOPT = gengetopt
# osX/tiger
#GENGETOPT = /opt/local/bin/gengetopt

.SECONDARY: semphy_cmdline.c semphy_cmdline.h

#LDFLAGS=  

CPPFLAGS=  -O3 -Wall -Wno-sign-compare -I. -DLOG -ftemplate-depth-32
CPPFLAGSDEBUG= -g -Wall -Wno-sign-compare -I. -DLOG -ftemplate-depth-32 -DVERBOS
#CPPFLAGSDOU=  $(CPPFLAGS) 
#-pg


#CPPFLAGS+= -I/usr/include/g++-v3
#CPPFLAGS+= -DLOG -DLOGCLS -DMEMCHK

# sources
sources= $(Libsources) $(LibCsources) $(addsuffix .cpp,$(EXEC) $(TEST_EXEC))

.PHONY: tests lib test debug %.debug
.PHONY:  dat DOUBLEREP doubleRep

all: lib $(EXEC)

test: all tests
	+cd tests; make -k test

ifdef DOUBLEREP
CPPFLAGS+= -DLOGREP
CPPFLAGSDEBUG += -DLOGREP
LDFLAGSDEBUG +=  -DLOGREP
endif


debug: CPPFLAGS = -g -Wall -Wno-sign-compare -I. -DLOG -ftemplate-depth-32
debug:  $(DEBUGLIB)
pl:
	@echo "lib  ="$(LIB)
	@echo "debug="$(DEBUGLIB)
#debug: all 
#	cp libEvolTree.a libEvolTreeDebug.a

# <<<<<<< Makefile
# %.debug: CPPFLAGS = -g -Wall -Wno-sign-compare -I. -DLOG -ftemplate-depth-25
# % debug: LIB = libEvolTreeDebug.a
# %.debug: %
# 	@echo "made \""$(*)"\" in debug mode"

# =======
#>>>>>>> 2.34

lib: $(LIB)

$(LIB): $(Libsources:.cpp=.o) $(LibCsources:.c=.o)
	ar rv $@ $?
	ranlib $@

tags: *.cpp *.h
	etags --members --language=c++ $^ 

$(EXEC) $(TEST_EXEC): $(LIB)
tests: $(TEST_EXEC)

-include make.dep

install:
	cd ../fast; make -f Makefile.lib install_do


clean:
	-rm -f $(LIB) $(DEBUGLIB) $(DOUBLEREPLIB) $(EXEC) $(TEST_EXEC) *.o


ifneq ($(wildcard make.dep), make.dep)
 make.dep: depend
endif

depend makedep: _make.dep
	@mv -f _make.dep make.dep 

_make.dep: $(sources)
	@echo making depend
#	$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $^ | sed '\''s/\($*\)\.o[ :]*/\1.o $@ : /g'\'' > $@ ; [ -s $@ ] || rm -f $@' 
	@$(SHELL) -ec '$(CC) -MM $(CPPFLAGS) $^ | sed "s/\(^[^.]*\)\.o/\1.o \1.debug.o/g" > $@' 
_fast:
	cd ../fast;make -f Makefile.lib -k all

fast.% _fast.%:
	cd ../fast;make -f Makefile.lib -k $(*)


simulateSequnce: simulateSequnce_cmdline.o


evolObjsTest.ggo: evolObjs.header evolObjs.args
	cat $^ > $@


# commandline (gengetopts)
%_cmdline.h %_cmdline.c: %.ggo
	$(GENGETOPT)  -i$< -F$(*)_cmdline

%.dat.q: %.dat
	awk 'BEGIN{RS="[\n\r]+";};{print "\" "$$0" \"\r"}' $< >  $@
#	cat $@

DAT = cpREV45.dat.q dayhoff.dat.q jones.dat.q mtREV24.dat.q wag.dat.q HIVb.dat.q HIVw.dat.q

dat: $(DAT)

cleandat:
	rm $(DAT)

datMatrixHolder.o: $(DAT)

.PRECIOUS: $(DAT)

debug: LIB = $(DEBUGLIB)

%.debug: CPPFLAGS = $(CPPFLAGSDEBUG)
%.debug: %
	@echo "made \""$(*)"\" in debug mode"


%.debug.o: %.c
	$(CC) -c $(CPPFLAGSDEBUG) $(CFLAGS) $< -o $@

%.debug.o: %.cpp
	$(CXX) -c $(CPPFLAGSDEBUG) $(CXXFLAGS) $< -o $@

$(DEBUGLIB): $(Libsources:.cpp=.debug.o) $(LibCsources:.c=.debug.o)
	ar rv $@ $?
	ranlib $@

doubleRep: LOGREP=t
doubleRep: CPPFLAGS+= -DLOGREP
doubleRep: $(DOUBLEREPLIB)

%.doubleRep.o: %.c
	$(CC) -c $(CPPFLAGS) $(CFLAGS) $< -o $@

%.doubleRep.o: %.cpp
	$(CXX) -c $(CPPFLAGS) $(CXXFLAGS) $< -o $@

$(DOUBLEREPLIB): $(Libsources:.cpp=.doubleRep.o) $(LibCsources:.c=.doubleRep.o)
	ar rv $@ $?
	ranlib $@

# DO NOT DELETE
