#ifndef _NODE_
#define _NODE_
#include <vector>
#include <cstddef> //for NULL definition
using namespace std;
//Tree node class definition
class Node{

	public:
	~Node(){
	parent = NULL;
	for(int i=0; i<children.size(); i++){
		if(children.at(i) != NULL){
			children.at(i)=NULL;
		}
	}
	};

	Node* parent = NULL;
	vector<Node*> children{};
	int pdg = -1;
	double mass = -1.0; // -1 for leaf or no 			constraint at this node
	bool isLeaf = 0;

	//unique ID from key
	int nodeId = -1;
	
	//store all combinations of final state particles
	//indices of initial particle array
	vector<vector<int> > combinations{};
	vector<vector<int> > combinations_pdgs{};
	
	//place holder for current set of particles for this node
	vector<int> currentunusedparts{};
	vector<int> currentcombination{};
	vector<int> currentcombination_pdgs{};

	//the number of final state particles from this particular decay node
	int nLeaves=0;
	//the pdgcodes of all the final state particles for this
	vector<int> leafpdgs{};
	
};

/*
Node::~Node(){
	for(int i=0; i<children.size(); i++){
		if(children.at(i) != NULL){
			children.at(i)=NULL;
		}
	}
}*/

#endif

