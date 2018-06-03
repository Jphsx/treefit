#ifndef _TREEFIT_
#define _TREEFIT_


#include "Tree.h"
#include "Particle.h"
#include "Combinatorics.h"


using namespace std;


class TreeFit{

	public:
		
	//global set of reconstructed particlesX
	vector<Particle*> recoparts{};
	//ID is index of particle on recoparts
	vector<int> recoIDs{};
	int LASTNONLEAFID;
	Tree* globalTree;

	void initializerecoparts(vector<int> recopdgs);
	void printParticles(vector<Particle*> parts);
	vector<vector<int> > makepdgcombinations(vector<vector<int> > combinations);
	vector<int> getpdgcombo(vector<int> combo);
	void generatefitcombinations(Node* root, vector<int> parentcombo)
		

};




#endif
