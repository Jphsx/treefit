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
	void clearEvent();

	int LASTNONLEAFID;
	Tree* ParticleTree;

	//TODO add fitTable

	
	void addrecopart(Particle* pc);

	void printParticles(vector<Particle*> parts);
	vector<vector<int> > makepdgcombinations(vector<vector<int> > combinations);
	vector<int> getpdgcombo(vector<int> combo);
	void generatefitcombinations(Node* root, vector<int> parentcombo);
		

};




#endif
