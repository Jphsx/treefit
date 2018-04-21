#ifndef _TREE_
#define _TREE_
#include "Node.h"
#include <vector>
#include <iostream>
#include <string.h>
using namespace std;



class Tree{
	public:
//tree construction methods
	Node* newNode(int id);

	vector<double> castVector_double(vector<string> v);
	vector<int> castVector_int(vector<string> v);
	vector<string> splitString(string str, string delimiter);
	Node* constructTree(vector<string> serial, int* serialIndex);
	void preOrderTraverse(Node* tree);
	void postOrderTraverse(Node* tree);
	void setTreeMasses(Node* root, vector<double> masses, int* massptr);
	void setTreePdgCodes(Node* root, vector<int> pdgs, int* pdgptr );
	void markTreeLeaves(Node* root);
	void findLeaves(Node* root, Node* originalParent);
	void populateNLeaves(Node* root);
	void setParents(Node* root, Node* parentptr);
	void printTree(Node* root);
	Node* treeInit(string pdg, string serial, string mass, string delimiter, int TESTNUM);

	Node* getParentNextNonLeafChild();

	//TODO add deconstructor for NODE

};
//printing is generic template, so define it externally as a templated too, therefore tree is not templated
//these methods are defined in the cpp, with no specific class attachment
//void printvector(vector<type> v);
//void print2dvec(vector<vector<type> > v);

#endif





