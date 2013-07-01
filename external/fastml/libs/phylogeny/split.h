// $Id: split.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___SPLIT
#define ___SPLIT

#include "definitions.h"
#include <set>
#include <vector>
#include <iostream>
#include <cassert>
using namespace std;


// this split always has the member "1" in it.  
// if not, it will take the reverse of the split, so that it dose have the "1" member.

class split {
public:
  explicit split (const int max=0); // empty split

// get an itarator of members and the max member. 

template<class Iterator>
split (Iterator& i,
		 Iterator& end,
		 int max):_max(max),  _reverse(true){ 
	for(int j=0;j<max;++j)
		_set[1].insert(j);
  
	for (;i!=end;++i){
		assert((*i)<_max && (*i) >= 0);
		_set[0].insert(*i);
		_set[1].erase(*i);
		if (*i==0)		// if we add "0", we may need to reverse the split
		reverse();
	}
}

  bool isMember(const int key) const;
  int size() const ;
  void print(ostream& sout = cout) const;
  bool lessThen(const split& other) const;
  bool compatible(const split& other) const ;

  // remove the key from the active set to the non-active set or vice versa.
  // for example if the split is {0,1 | 2}
  // reverseMembership(1) will change the split to this one: {0 | 1,2 }
  void reverseMembership(const int key); 

  void getId(vector<int> & id) const  {
    id.clear();
    bool small(_set[0].size()>_set[1].size());
    for (set<int>::const_iterator i=_set[small].begin();i!=_set[small].end();++i)
      id.push_back(*i);
  }

private:
  void reverse();


  int _max;			// max element.  all elements are asumed to be in the range [1..max]
  set<int> _set[2];
  bool _reverse;
};

bool operator<(const split& a, 
	       const split& b) ;

  

ostream& operator<< (ostream &sout,  const split& split) ;



#endif // ___SPLIT
