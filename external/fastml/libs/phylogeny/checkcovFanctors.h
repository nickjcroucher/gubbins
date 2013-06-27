// $Id: checkcovFanctors.h 1732 2007-02-26 13:45:41Z itaymay $

#ifndef ____CHECKCOV__FANCTORS
#define ____CHECKCOV__FANCTORS
#include "definitions.h"
#include "tree.h"

#include "likelihoodComputation.h"
using namespace likelihoodComputation;
#include "sequenceContainer.h"
#include "stochasticProcess.h"
#include "logFile.h"

#include <cmath>

//#define VERBOS

#ifdef VERBOS
#include <iostream>
using namespace std;
#endif

class Cevaluate_L_given_r{
public:
  explicit Cevaluate_L_given_r(	const sequenceContainer& sd,
								const tree& t1,
								const stochasticProcess& sp,
								const int pos)
			:_sd(sd),_t1(t1),_pos(pos), _sp(sp) {}
private:
	const sequenceContainer& _sd;
	const tree& _t1;
	const int _pos;
	const stochasticProcess& _sp;
public:
	MDOUBLE operator() (const MDOUBLE r) {
		
		MDOUBLE tmp1= convert(getLofPos(_pos,_t1,_sd,_sp,r));
#ifdef VERBOS
			LOG(5,<<" r = "<<r<<" l = "<<tmp1<<endl);
#else
			LOG(12,<<" r = "<<r<<" l = "<<tmp1<<endl);
#endif
		return -tmp1;
	}
};

// THIS FUNCTION IS USED ONLY BY ITAY MAYROSE AND ONLY HE KNOWS WHAT IS INSIDE...
// ONE DAY HE WILL WRITE .DOC FILES...
class Cevaluate_Posterior_given_r {
public:
  explicit Cevaluate_Posterior_given_r(	const sequenceContainer& seqContainer,
								const tree& t1,
								const stochasticProcess& sp,
								const MDOUBLE alpha,
								const int pos)
			:m_seqContainer(seqContainer), m_alpha(alpha),m_tree(t1), m_pos(pos), m_sp(sp) {}
public:
	MDOUBLE operator() (const MDOUBLE r) 
	{
		
		MDOUBLE l= convert(getLofPos(m_pos, m_tree, m_seqContainer, m_sp, r));
		#ifdef VERBOS
			LOG(5,<<" r = "<<r<<" l = "<<l<<endl);
		#endif
		MDOUBLE prior = exp((-m_alpha) * r) * pow(r, m_alpha - 1);
		return -(l * prior);
	}

private:
	const sequenceContainer& m_seqContainer;
	const MDOUBLE m_alpha;
	const tree& m_tree;
	const int m_pos;
	const stochasticProcess& m_sp;

};

// WHEN YOU WANT TWO TREE TO HAVE THE SAME RATE AT A SPECIFIC POSITION.
class Cevaluate_L_sum_given_r{
public:
	explicit Cevaluate_L_sum_given_r(const stochasticProcess& sp,
									 const sequenceContainer& sd,
									 const tree &inLTree1,
									 const tree &inLTree2,
									 const int pos)
			:_sp(sp), _sd(sd), _tree1(inLTree1),_tree2(inLTree2), _pos(pos){};

private:
	const stochasticProcess _sp;
	const sequenceContainer _sd;
	const tree& _tree1;
	const tree& _tree2;
	const int _pos;
public:
	MDOUBLE operator() (const MDOUBLE r) {
		MDOUBLE tmp1= convert(getLofPos(_pos,_tree1,_sd,_sp,r));
		MDOUBLE tmp2= convert(getLofPos(_pos,_tree2,_sd,_sp,r));
		MDOUBLE tmp= tmp1*tmp2;
		return -tmp;
	}
};

#endif
