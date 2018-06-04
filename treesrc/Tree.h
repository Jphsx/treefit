#ifndef _TREE_
#define _TREE_
#include "Node.h"
#include <vector>
#include <iostream>
#include <string.h>
using namespace std;



class Tree{
	public:
	~Tree();
	void deletePostTraverse(Node* root);

	int lastNonLeafNodeId;
	Node* Root;
	//tree construction methods
	Node* newNode(int id);
	

	vector<double> castVector_double(vector<string> v);
	vector<int> castVector_int(vector<string> v);
	vector<string> splitString(string str, string delimiter);
	Node* constructTree(vector<string> serial, int* serialIndex);
	void preOrderTraverse(Node* tree);
	void postOrderTraverse(Node* tree);
	void setTreeMasses(Node* root, vector<float> masses, int* massptr);
	void setTreePdgCodes(Node* root, vector<int> pdgs, int* pdgptr );
	void markTreeLeaves(Node* root);
	void findLeaves(Node* root, Node* originalParent);
	void populateNLeaves(Node* root);
	void setParents(Node* root, Node* parentptr);
	void getLastNonLeafNodeId(Node* root, int* id);
	void printTree(Node* root);
	//void treeInit(string pdg, string serial, string mass, string delimiter, int TESTNUM);
	void treeInit(vector<int> pdg, string serial, vector<float>  mass, string delimiter, int TESTNUM);


	static Node* locateAncestorNearestNonLeafChild(Node* root);
	//returns the pointer to the closest non leaf sibling
	static Node* getParentNextNonLeafChild(Node* parent, int callingNodeID);
	//gets the first child that is not a leaf
	static Node* getFirstNonLeafChild(Node* root);

	//TODO add deconstructor for NODE
	template <typename type>
	static void printvector(vector<type> v);
	template <typename type>
	static void print2dvec(vector<vector<type> > v);
	
	//iterate and print all non leaf current combos
	static void printfit(Node* root);
	

};


#endif





