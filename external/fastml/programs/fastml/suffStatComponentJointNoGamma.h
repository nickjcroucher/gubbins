#ifndef SUFF_STAT_COMPONENT_JOINT_NO_GAMMA_H___
#define SUFF_STAT_COMPONENT_JOINT_NO_GAMMA_H___

#include "definitions.h"
#include <vector>
#include <cassert>
using namespace std;

class suffStatSpecHomPosJointNoGamma{ // this is for a specific node.
	public:
		void set(const int letterInFather,const int val) {
			_V[letterInFather]=val;
		}
		
		int get(const int letterInFather) const	{
			return _V[letterInFather];
		}

		void allocatePlace(const int alphabetSize) {
			_V.resize(alphabetSize);
		}
		bool isEmpty (){return (_V.empty());};
		size_t size() {return _V.size();}
	private:
		Vint _V;//size = alphabet size
};

class suffStatGlobalHomPosJointNoGamma{ // this is for all nodes
	public:
		void set(const int nodeId,const int letterInFather,const int val) {
			_V[nodeId].set(letterInFather,val);
		}
		
		int get(const int nodeId,const int letterInFather) const	{
			return _V[nodeId].get(letterInFather);
		}

		void allocatePlace(const int numOnNodes,const int alphabetSize) {
			_V.resize(numOnNodes);
			for (int i=0;i<_V.size();++i) {_V[i].allocatePlace(alphabetSize);}
		}
		bool isEmpty (){return (_V.empty());}
		size_t size() {return _V.size();}

	private:
		vector<suffStatSpecHomPosJointNoGamma> _V;//size = letter
};


#endif
