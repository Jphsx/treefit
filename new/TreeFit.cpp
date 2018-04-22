


#include "Tree.h"
#include "Particle.h"
#include "Combinatorics.h"


using namespace std;



//global set of reconstructed particles
vector<Particle*> recoparts{};
//ID is index of particle on recoparts
vector<int> recoIDs{};
int LASTNONLEAFID;

Tree* globalTree;
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
		//cout<<"Status "<<parts.at(i)->used<<endl;
		cout<<endl;
	}
}
vector<vector<int> > makepdgcombinations(vector<vector<int> > combinations){
	vector<vector<int> > pdgcombinations;
	vector<int> pdgcombo;
	for(int i=0; i<combinations.size();i++){
		for(int j=0; j<combinations.at(i).size(); j++){
			pdgcombo.push_back(recoparts.at(combinations.at(i).at(j))->recopdg);
		}
		pdgcombinations.push_back(pdgcombo);
		pdgcombo.clear();
	}
	return pdgcombinations;
}
vector<int> getpdgcombo(vector<int> combo){
	vector<int> pdgcombo;
	for(int i=0; i<combo.size(); i++){
		pdgcombo.push_back(recoparts.at(combo.at(i))->recopdg);
	}
	return pdgcombo;	
}
void generatefitcombinations(Node* root, vector<int> parentcombo){

	
	//no reason to ever visit a leaf
	//if this is a non leaf do combination stuff
	if(root->isLeaf) return;

	//first generate all combinations for this node
	//arguments(n,k,vector) n choose k and put index combos onto vector
	Combinatorics::generateParticlesCombinations(parentcombo.size(), root->nLeaves, root->combinations, parentcombo);
	//map root combinations to a combination vector containing the pdgs instead of recoids
	
	Combinatorics::filtercombinations(root->leafpdgs, root->combinations, makepdgcombinations(root->combinations));
	
	//iterate through combinations and through all children recursively
	
	for(int i=0; i<root->combinations.size(); i++){
		//first set the current combination
		root->currentcombination = root->combinations.at(i);
		
		root->currentunusedparts = root->combinations.at(i);
		Combinatorics::quickSort(root->currentcombination, 0, root->currentcombination.size()-1);
		Combinatorics::quickSort(root->currentunusedparts,0, root->currentunusedparts.size()-1);

		//do after sorting now pdgs follow sort indices
		root->currentcombination_pdgs= getpdgcombo(root->combinations.at(i));

		//modify parent unused parts
		//these sets need to be sorted
		if(root->parent != NULL){
			root->parent->currentunusedparts = Combinatorics::subtractSets(root->parent->currentunusedparts, root->currentcombination);
		}

		

		//move down to first child if possible
		if( Tree::getFirstNonLeafChild(root) != NULL){
			generatefitcombinations(Tree::getFirstNonLeafChild(root),  root->currentcombination);
		
		}//if we cant move down
		// we travel upwards toward the next
		//ancestor non leaf //what if this is only 2 particle?
		else if(Tree::locateAncestorNearestNonLeafChild(root) != NULL){
			//the nearest ancestor will inherit set from its parent
			generatefitcombinations(Tree::locateAncestorNearestNonLeafChild(root), Tree::locateAncestorNearestNonLeafChild(root)->parent->currentunusedparts);
		}

		if(root->nodeId == LASTNONLEAFID){
			//do fit
			cout<<"FIT"<<endl;
			Tree::printfit(globalTree->Root);
			cout<<endl;
			
		}
		
		//before we move to the next combination
		//add current combo into parent set
		if(root->parent != NULL){
			root->parent->currentunusedparts = Combinatorics::addSets(root->parent->currentunusedparts, root->currentcombination);
		}

		
		
	}//end combinations loop
		
	////we are done back up to the previous call
	root->currentcombination.clear();
	root->currentcombination_pdgs.clear();
	root->currentunusedparts.clear();
	root->combinations.clear();
	
	return;
}

int main(){ 

	//Tree* tree;
	//Node* root;
	string delimiter = " ,";
	vector<int> recopdgs{};
	
	cout<<"Testing Set Operations with sets A, B "<<endl;
	vector<int> A = { 1, 2, 3, 40, 50, 60 };
	vector<int> B = { 3, 40, 45, 50 };
	cout<<" A: ";
	Tree::printvector(A);
	cout<<endl;
	cout<<" B: ";
	Tree::printvector(B);
	cout<<endl;

	vector<int> solution;
	solution = Combinatorics::addSets(A,B);
	cout<<" A+B = ";
	Tree::printvector(solution);
	cout<<endl;
	
	solution = Combinatorics::subtractSets(A,B);
	cout<<" A-B = ";
	Tree::printvector(solution);
	cout<<endl;
	
	
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

*/
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"CASE 3pi0 + photon "<<std::endl;
	Tree* tree4;

	string preorder_pdg4 = " 223 22 221 111 22 22 111 22 22 111 22 22";
	string preorder_key4 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial4 = "0 1 ) 2 3 4 ) 5 ) ) 6 7 ) 8 ) ) 9 10 ) 11 ) ) ) )";
	string preorder_mass4 = "0.782 -1 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1";
	//try best case
	recopdgs = { 22, 22, 22, 22, 22, 22, 22};
	initializerecoparts(recopdgs);
	printParticles(recoparts);

	tree4->treeInit(preorder_pdg4, preorder_serial4, preorder_mass4, delimiter, 4);
	cout<<"tree constructed"<<endl;
	globalTree=tree4;
	generatefitcombinations(tree4->Root, recoIDs);
	//delete tree;
	//tree=NULL;
	//theTree = NULL;
	recopdgs.clear();
	recoparts.clear();
	cout<<std::endl;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
	Tree* tree5;
	string preorder_pdg5 = " 221 211 -211 111 22 22";
	string preorder_key5 = "0 1 2 3 4 5";
	string preorder_serial5 = "0 1 ) 2 ) 3 4 ) 5 ) ) )";
	string preorder_mass5 = "0.547 -1 -1 0.135 -1 -1";

	recopdgs = { 211,-211, 13, -13,-211, 22, 22};
	initializerecoparts(recopdgs);
	//printParticles(recoparts);

	 tree5->treeInit(preorder_pdg5, preorder_serial5, preorder_mass5, delimiter, 5);
	LASTNONLEAFID = tree5->lastNonLeafNodeId;
	globalTree=tree5;
	generatefitcombinations(tree5->Root, recoIDs);

	//delete tree;
	//tree = NULL;
	//globalTree = NULL; 
	recopdgs.clear();
	recoparts.clear();

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Tree* tree6;
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
	initializerecoparts(recopdgs);	
	printParticles(recoparts);
	
	//tree fitting combination starting with simplest cases
	tree6->treeInit(preorder_pdg6, preorder_serial6, preorder_mass6, delimiter, 6);
	LASTNONLEAFID = tree6->lastNonLeafNodeId;
	globalTree=tree6;
	generatefitcombinations(tree6->Root, recoIDs);
	
	//delete tree;
	//tree = NULL;
	//globalTree = NULL; 
	recopdgs.clear();
	recoparts.clear();
*/


}

