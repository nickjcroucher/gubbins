// $Id: countTableComponent.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___COUNT_TABLE_COMPONENT
#define ___COUNT_TABLE_COMPONENT

#include "definitions.h"
#include <cassert>

class countTableComponentHom{
public:  

	void setCount( const int letter1,
					const int letter2,
					const MDOUBLE val) {
		_countValues[letter1][letter2]=val;
	}
	int alphabetSize() const {return _countValues.size();}
	void zero();
	MDOUBLE getCounts(	const int letter1,
						const int letter2) const	{
		return _countValues[letter1][letter2];
	}
	void addToCounts(const int let1,const int let2,const MDOUBLE val) {
		_countValues[let1][let2]+=val;
	}
	bool isEmpty (){return (_countValues.empty());};
	void countTableComponentAllocatePlace(const int alphabetSize);
	void printTable(ostream & out) const;
private:					
	VVdouble _countValues;//letter1,letter2

};

class countTableComponentGam{
public:  

	void setCount( const int letter1,
					const int letter2,
					const int rateCategor,
					const MDOUBLE val) {
		_countValues[rateCategor].setCount(letter1,letter2,val);
	}
	
	int alphabetSize() const {return _countValues.empty()?0:_countValues[0].alphabetSize();}
	void zero(){
		for (int i=0; i < _countValues.size(); ++i ) _countValues[i].zero();
	}


	MDOUBLE getCounts(	const int letter1,
						const int letter2,
						const int rateCategor) const {
		assert(_countValues[rateCategor].getCounts(letter1,letter2)>=0);
		return _countValues[rateCategor].getCounts(letter1,letter2);
	}

	void addToCounts(const int let1,const int let2,
		const int rate,const MDOUBLE val) {
		_countValues[rate].addToCounts(let1,let2,val);
	}

	bool isEmpty (){return (_countValues.empty());};

	void countTableComponentAllocatePlace(const int alphabetSize,
		const int numberOfrateCategories) {
		_countValues.resize(numberOfrateCategories);
		for (int i=0; i < _countValues.size(); ++i ){
			_countValues[i].countTableComponentAllocatePlace(alphabetSize);
		}
	}
		void printTable(ostream & out) const {
			for (int i=0; i < _countValues.size(); ++i) {
				_countValues[i].printTable(out);
			}
		}
	countTableComponentHom& operator[] (int i) {return _countValues[i];}
	const countTableComponentHom& operator[] (int i) const {return _countValues[i];}
private:
	vector<countTableComponentHom> _countValues;//letter1,letter2,rateCategor 

};

#endif

