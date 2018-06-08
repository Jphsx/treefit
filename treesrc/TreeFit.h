#ifndef _TREEFIT_
#define _TREEFIT_


#include "Tree.h"
#include "Particle.h"
#include "Combinatorics.h"


using namespace std;


class TreeFit{

	public:
	TreeFit();	
	//set of reconstructed particles containers
	vector<Particle*> recoparts{};
	//ID is index of particle on recoparts
	vector<int> recoIDs{};

	//when moving to a new event clear the old reco particles
	//and old fit combinations
	void clearEvent();

	int* LASTNONLEAFID;
	Tree* ParticleTree;

	//Fit table is a 3d vector
	//first index is ith Node
	//second index is the jth fit combination for that node
	//third index is the index of the particle on recoparts in the jth combination for ith node
	std::vector< std::vector< std::vector<int>>> fitTable{};
	std::vector<  
	void initTable();
	void printTable();
	
	void addrecopart(Particle* pc);

	void printParticles(vector<Particle*> parts);
	vector<vector<int> > makepdgcombinations(vector<vector<int> > combinations);
	vector<int> getpdgcombo(vector<int> combo);
	void generatefitcombinations(Node* root, vector<int> parentcombo);

	void addFitToTable(Node* root);
		

};




#endif
