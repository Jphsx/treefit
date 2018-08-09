#ifndef _TREEFIT_
#define _TREEFIT_


#include "Tree.h"
#include "Particle.h"
#include "Combinatorics.h"


//using namespace std;


class TreeFit{

	public:
	TreeFit();	
	//set of reconstructed particles containers
	std::vector<Particle*> recoparts{};

	//set of fit particles containers
	std::vector<Particle*> fitparts{};

	//ID is index of particle on recoparts
	std::vector<int> recoIDs{};

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
	//duplicate vector, but instead of indices it holds pdgs
	//good for debugging purposes (printing out)
	std::vector< std::vector< std::vector<int>>> fitPdgs{};
	 
	void initTable();
	void printTable();
	
	void addrecopart(Particle* pc);

	void printParticles(std::vector<Particle*> parts);
	std::vector<std::vector<int> > makepdgcombinations(std::vector<std::vector<int> > combinations);
	std::vector<int> getpdgcombo(std::vector<int> combo);
	void generatefitcombinations(Node* root, std::vector<int> parentcombo);

	void addFitToTable(Node* root);

	//determine the subset of particles in set to share a common vertex
	std::vector<int> getVertexSet(std::vector<int> combo, int nodeId, std<vector<std::vector<int> > fit); 
		

};




#endif
