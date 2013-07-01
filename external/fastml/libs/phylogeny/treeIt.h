// $Id: treeIt.h 962 2006-11-07 15:13:34Z privmane $

#ifndef ___TREE_IT
#define ___TREE_IT
#include "definitions.h"
#include "errorMsg.h"
#include "tree.h"


class treeIterTopDown{
public:
	treeIterTopDown(tree& t) : _t(t) , _current(_t.getRoot()) {
		_childCheck.push_back(0);
	}
	tree::nodeP first() {
		_childCheck.clear();
		_childCheck.push_back(0);
		_current = _t.getRoot();
		return _t.getRoot();
	}
	tree::nodeP next() {
		if (_childCheck.empty()) return NULL;
		if (_childCheck[_childCheck.size()-1]<_current->getNumberOfSons()) {
			_current = _current->getSon(_childCheck[_childCheck.size()-1]);
			_childCheck[_childCheck.size()-1]++;
			_childCheck.push_back(0);
		}
		else {
			_current = _current->father();
			_childCheck.pop_back();
			return next();
		}
		return _current;
	}
	tree::nodeP operator++(int) {return next();}
	tree::nodeP operator++() {return next();}
	tree::nodeP end(){ return NULL;}
	tree::nodeP operator-> (){ return _current;}
	tree::TreeNode& operator* (){return *_current;}
	bool operator!= (tree::nodeP t) {return (t != this->_current);}
private:
	vector<int> _childCheck;
	tree& _t;
	tree::nodeP _current;
};

class treeIterTopDownConst{
public:
	treeIterTopDownConst(const tree& t) : _t(t) , _current(_t.getRoot()) {
		_childCheck.push_back(0);
	}
	tree::nodeP first() {
		_childCheck.clear();
		_childCheck.push_back(0);
		_current = _t.getRoot();
		return _t.getRoot();
	}
	tree::nodeP next() {
		if (_childCheck.empty()) return NULL;
		if (_childCheck[_childCheck.size()-1]<_current->getNumberOfSons()) {
			_current = _current->getSon(_childCheck[_childCheck.size()-1]);
			_childCheck[_childCheck.size()-1]++;
			_childCheck.push_back(0);
		}
		else {
			_current = _current->father();
			_childCheck.pop_back();
			return next();
		}
		return _current;
	}
	tree::nodeP operator++(int) {return next();}
	tree::nodeP operator++() {return next();}
	tree::nodeP end(){ return NULL;}
	tree::nodeP operator-> (){ return _current;}
	tree::TreeNode& operator* (){return *_current;}
	bool operator!= (tree::nodeP t) {return (t != this->_current);}
private:
	vector<int> _childCheck;
	const tree& _t;
	tree::nodeP _current;
};

class treeIterDownTopConst{
public:
	treeIterDownTopConst(const tree& t) : _t(t) , _current(_t.getRoot()) {
		_childCheck.push_back(0);
	}
	const tree::nodeP first() {
		_childCheck.clear();
		_childCheck.push_back(0);
		_current = _t.getRoot();
		return next();
	}
	const tree::nodeP next() {
		if (_childCheck[_childCheck.size()-1]>_current->getNumberOfSons()) {//checked
			_current = _current->father();
			if (!_current) return NULL;
			_childCheck.pop_back();
			_childCheck[_childCheck.size()-1]++;
			return next();
		}
		else if (_childCheck[_childCheck.size()-1]<_current->getNumberOfSons()) {
			_current = _current->getSon(_childCheck[_childCheck.size()-1]);
			_childCheck.push_back(0);
			return next();
		}
//		else //if (_childCheck[_childCheck.size()-1]==_current->getNumberOfSons()) 
//		{
				_childCheck[_childCheck.size()-1]++;
				return _current;
//		}
		
//		return next();
	}
	const tree::nodeP operator++(int) {return next();}
	const tree::nodeP operator++() {return next();}
	const tree::nodeP end(){ return NULL;}
	const tree::nodeP operator-> (){ return _current;}
	const tree::TreeNode& operator* (){return *_current;}
	bool operator!= (tree::nodeP t) {return (t != this->_current);}
private:
	vector<int> _childCheck;
	const tree& _t;
	tree::nodeP _current;
};

#endif
