


#include "Tree.h"
#include "Particle.h"
#include "Combinatorics.h"

using namespace std;



//global set of reconstructed particles
vector<Particle*> recoparts{};
vector<int> recoIDs{};
//this function should probably take in a TLV or ReconstructedParticle in addition to pdg array to fully populate the object
void initializerecoparts(vector<int> recopdgs){
	for(int i=0; i<recopdgs.size(); i++){
		Particle* p = new Particle();
		p->recopdg=recopdgs.at(i);
		recoparts.push_back(p);
		//also creat a vector of unique 
		//recoparticle id's
		recoIDs.push_back(i);
		//the unique ID is just the index on the 
		//reco parts array
	}
}
void printParticles(vector<Particle*> parts){
	for(int i=0; i<parts.size(); i++){
		cout<<"Particle Index "<<i<<endl;
		cout<<parts.at(i)->recopdg<<endl;
		cout<<"Status "<<parts.at(i)->used<<endl;
		cout<<endl;
	}
}
/*void generatefitcombinations(Node* root, vector<int> parentcombo, const int last_non_leaf_id){

	//if(root->isLeaf) return; stop returning on leaf, generate combos at leaves (Mark used there)
	//first generate all combinations for this node
	//arguments(n,k,vector) n choose k and put index combos onto vector

	
	generateIndicesCombinations(parentcombo.size(), root->nLeaves, root->combinations, parentcombo);
	//filtercombinations(root->leafpdgs, root->combinations, root->combinations_pdgs, parentcombo);
	filtercombinations(root->leafpdgs, root->combinations);//, parentcombo);
	
	//iterate through combinations and through all children recursively
	//for(int i=0; i
	//++*CURRENT_NON_LEAF_NODES;
	for(int i=0; i<root->combinations.size(); i++){
		//first set the current combination
		root->currentcombination = root->combinations.at(i);
		//a current combination has been set, tally the global number of combinations
		//this is contributing to the number of potenial consraints for a fit
		//++*CURRENT_NON_LEAF_NODES;

		//if(*CURRENT_NON_LEAF_NODES == TOTAL_NON_LEAF_NODES){
		//check the node id, if its the last non leaf node id need to fit
		if(root->nodeId == last_non_leaf_id){
			//do fit
			//TODO: validate the fit (ancestry)
			cout<<"FIT"<<endl;
			printfit(theTree);
			cout<<endl;
			//--*CURRENT_NON_LEAF_NODES;
			
		}

		
	
		for(int k=0; k<root->children.size(); k++){
			//recurse through the children
			//generatefitcombinations(root->children.at(k));
			//generatefitcombinations(root->children.at(k), root->currentcombination, CURRENT_NON_LEAF_NODES, TOTAL_NON_LEAF_NODES);
			generatefitcombinations(root->children.at(k), root->currentcombination, last_non_leaf_id);
		}
		//before we proceed to the next combination clear the previously used
		root->currentcombination.clear();
		//--*CURRENT_NON_LEAF_NODES;
		
	}
	//--*CURRENT_NON_LEAF_NODES;
	root->combinations.clear();
	
	return;
}*/

int main(){ 

	Tree* tree;
	string delimiter = " ,";
	vector<int> recopdgs{};
	
/*	string preorder_pdg = " 443 331 321 -321 221 211 -211 111 22 22 ";
	string preorder_key = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial = "0 1 2 ) 3 ) ) 4 5 ) 6 ) 7 8 ) 9 ) ) ) )";
	string preorder_mass = " 3.096 1.0195 -1 -1 0.547 -1 -1 0.135 -1 -1 ";
	
	
	tree = testTree(preorder_pdg, preorder_serial, preorder_mass, delimiter, 1);
	delete tree;
	tree=NULL;

/*
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//test tree #2
	string preorder_pdg2 = " 443 221 111 22 22 211 -211 331 321 -321";
	string preorder_key2 = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial2 = "0 1 2 3 ) 4 ) ) 5 ) 6 ) ) 7 8 ) 9 ) ) )";
	string preorder_mass2 = " 3.096 0.547 0.135 -1 -1 -1 -1 1.0195 -1 -1";
	
	tree = testTree(preorder_pdg2, preorder_serial2, preorder_mass2, delimiter, 2);
	delete tree;
	tree=NULL;
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg3 = " 223 221 111 22 22 111 22 22 111 22 22 22";
	string preorder_key3 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial3 ="0 1 2 3 ) 4 ) ) 5 6 ) 7 ) ) 8 9 ) 10 ) ) ) 11 ) )";
	string preorder_mass3 = "0.782 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1 -1";
	
	tree = testTree(preorder_pdg3, preorder_serial3, preorder_mass3, delimiter, 3);
	delete tree;
	tree=NULL;

/*
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"CASE 3pi0 + photon "<<std::endl;
	string preorder_pdg4 = " 223 22 221 111 22 22 111 22 22 111 22 22";
	string preorder_key4 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial4 = "0 1 ) 2 3 4 ) 5 ) ) 6 7 ) 8 ) ) 9 10 ) 11 ) ) ) )";
	string preorder_mass4 = "0.782 -1 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1";
	//try best case
	recopdgs = { 22, 22, 22, 22, 22, 22, 22};
	initializerecoparts(recopdgs);

	tree = testTree(preorder_pdg4, preorder_serial4, preorder_mass4, delimiter, 4);
	delete tree;
	tree=NULL;
	theTree = NULL;
	recopdgs.clear();
	recoparts.clear();
	cout<<std::endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
*/
	string preorder_pdg5 = " 221 211 -211 111 22 22";
	string preorder_key5 = "0 1 2 3 4 5";
	string preorder_serial5 = "0 1 ) 2 ) 3 4 ) 5 ) ) )";
	string preorder_mass5 = "0.547 -1 -1 0.135 -1 -1";

	recopdgs = { 211,-211, 13, -13,-211, 22, 22};
	//initializeusedvector(recopdgs.size());
	initializerecoparts(recopdgs);
	printParticles(recoparts);

	tree->treeInit(preorder_pdg5, preorder_serial5, preorder_mass5, delimiter, 5);
	delete tree;
	tree = NULL;
	//theTree = NULL; 
	//recopdgs.clear();
	//recoparts.clear();

/*
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg6 = " 221 111 22 22 211 -211";
	string preorder_key6 = "0 1 2 3 4 5";
	string preorder_serial6 = "0 1 2 ) 3 ) ) 4 ) 5 ) )";
	string preorder_mass6 = "0.537 0.135 -1 -1 -1 -1";

	//preorder combination generation
	//define the reconstucted particles (by pdg)
	// using tree 6, eta -> pi+ pi- pi0
	//arbitrary reconstructed particles: 4 tracks 3 photons 7 particles total
	recopdgs = { 211, -211, 13, -13, 22, 22, 22};
	
	//start by setting global array for all particles unused
	//initializeusedvector(recopdgs.size());
	initializerecoparts(recopdgs);	

	//tree fitting combination starting with simplest cases
	tree = testTree(preorder_pdg6, preorder_serial6, preorder_mass6, delimiter, 6);
	delete tree;
	tree = NULL;
	theTree = NULL; //null global tree after deleting the contents that it points to
	recopdgs.clear(); 
	recoparts.clear();
//////////////////////////////////////////////////////////////////////////////////////////////////////////////
//test 7 same as 6 but only eta constraint
/////////NOTE://///////////
// the mass = -1 nodes should be dealt with on the fit object generating level
// eliminating and reorganizing combinations based on -1 flag in mass is too tricky at tree level
// all combinations will still be generatated at this -1 node regardless of its constraint condition
///////////////////////////
	string preorder_pdg7 = " 221 111 22 22 211 -211";
	string preorder_key7 = "0 1 2 3 4 5";
	string preorder_serial7 = "0 1 2 ) 3 ) ) 4 ) 5 ) )";
	string preorder_mass7 = "0.537 -1 -1 -1 -1 -1";
	
	recopdgs = { 211, -211, 13, -13, 22, 22, 22};
	
	//start by setting global array for all particles unused
	//initializeusedvector(recopdgs.size());
	initializerecoparts(recopdgs);
	//tree fitting combination starting with simplest cases
	tree = testTree(preorder_pdg7, preorder_serial7, preorder_mass7, delimiter, 7);
	delete tree;
	tree = NULL;
	theTree = NULL; 
	recopdgs.clear();
	recoparts.clear();
*/
}

