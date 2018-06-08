
#include "TreeFit.h"

//this function should probably take in a TLV or ReconstructedParticle in addition to pdg array to fully populate the object
TreeFit::TreeFit(){
	ParticleTree = new Tree();
	LASTNONLEAFID = &(ParticleTree->lastNonLeafNodeId);
}

void TreeFit::addrecopart(Particle* pc){
	recoparts.push_back(pc);
	recoIDs.push_back(recoparts.size()-1);
}
void TreeFit::printParticles(vector<Particle*> parts){
	for(int i=0; i<parts.size(); i++){
		cout<<"Particle Index/recoID "<<i<<endl;
		Particle::printParticle(parts.at(i));
	}
}
vector<vector<int> > TreeFit::makepdgcombinations(vector<vector<int> > combinations){
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
vector<int> TreeFit::getpdgcombo(vector<int> combo){
	vector<int> pdgcombo;
	for(int i=0; i<combo.size(); i++){
		pdgcombo.push_back(recoparts.at(combo.at(i))->recopdg);
	}
	return pdgcombo;	
}
void TreeFit::generatefitcombinations(Node* root, vector<int> parentcombo){

	
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

		
		if(root->nodeId == *LASTNONLEAFID){
			//do fit
			//cout<<"FIT"<<endl;
			Tree::printfit(ParticleTree->Root);
			cout<<endl;
			
		}
		//move down to first child if possible
		else if( Tree::getFirstNonLeafChild(root) != NULL){
			generatefitcombinations(Tree::getFirstNonLeafChild(root),  root->currentcombination);
		
		}//if we cant move down
		// we travel upwards toward the next
		//ancestor non leaf //what if this is only 2 particle?
		else if(Tree::locateAncestorNearestNonLeafChild(root) != NULL){
			//the nearest ancestor will inherit set from its parent
			generatefitcombinations(Tree::locateAncestorNearestNonLeafChild(root), Tree::locateAncestorNearestNonLeafChild(root)->parent->currentunusedparts);
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
void TreeFit::clearEvent(){
	recoparts.clear();
	recoIDs.clear();
}
void TreeFit::addFitToTable(Node* root){
/*	if(root->isLeaf) return;

	//cout<<"Node "<<root->nodeId;
	//cout<<" FitRecoIDs: ";
	printvector(root->currentcombination);
	cout<<" FitRecoPDGs: ";
	printvector(root->currentcombination_pdgs);
	cout<<endl;
	for(int i=0; i< root->children.size(); i++){
		printfit(root->children.at(i));
	}
*/
}
void TreeFit::initTable(){
	//fitTable.reserve(*LASTNONLEAFID +1);
	std::vector< std::vector< std::vector<int>>> table(*LastNONLEAFID + 1); 
	fitTable = table;
	std::cout<<fitTable.size()<<" size "<<std::endl;
}
//Testing framework////////////////
/*int main(){ 

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
*/	
	
	
/*
	Tree* tree1;
	string preorder_pdg1 = " 443 331 321 -321 221 211 -211 111 22 22 ";
	string preorder_key1 = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial1 = "0 1 2 ) 3 ) ) 4 5 ) 6 ) 7 8 ) 9 ) ) ) )";
	string preorder_mass1 = " 3.096 1.0195 -1 -1 0.547 -1 -1 0.135 -1 -1 ";
	
	
	recopdgs = {211, -211, 22, 22, 321, -321, 22 , 13, -321};

	
	initializerecoparts(recopdgs);
	printParticles(recoparts);

	tree1->treeInit(preorder_pdg1, preorder_serial1, preorder_mass1, delimiter, 1);
	tree1->printTree(tree1->Root);
	cout<<"tree constructed"<<endl;
	LASTNONLEAFID = tree1->lastNonLeafNodeId;
	ParticleTree=tree1;
	generatefitcombinations(tree1->Root, recoIDs);
	//delete tree;
	//tree=NULL;
	//theTree = NULL;
	recopdgs.clear();
	recoparts.clear();
	cout<<"DONE";
	cout<<std::endl;
*/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//test tree #2
/*
	Tree* tree2;
	string preorder_pdg2 = " 443 221 111 22 22 211 -211 331 321 -321";
	string preorder_key2 = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial2 = "0 1 2 3 ) 4 ) ) 5 ) 6 ) ) 7 8 ) 9 ) ) )";
	string preorder_mass2 = " 3.096 0.547 0.135 -1 -1 -1 -1 1.0195 -1 -1";
	
	recopdgs = {211, -211, 22, 22, 321, -321, 22 , 13, -321};

	
	initializerecoparts(recopdgs);
	printParticles(recoparts);

	tree2->treeInit(preorder_pdg2, preorder_serial2, preorder_mass2, delimiter, 2);
	tree2->printTree(tree2->Root);
	cout<<"tree constructed"<<endl;
	LASTNONLEAFID = tree2->lastNonLeafNodeId;
	ParticleTree=tree2;
	generatefitcombinations(tree2->Root, recoIDs);
	//delete tree;
	//tree=NULL;
	//theTree = NULL;
	recopdgs.clear();
	recoparts.clear();
	cout<<"DONE";
	cout<<std::endl;
	
*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*
	Tree* tree3;
	string preorder_pdg3 = " 223 221 111 22 22 111 22 22 111 22 22 22";
	string preorder_key3 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial3 ="0 1 2 3 ) 4 ) ) 5 6 ) 7 ) ) 8 9 ) 10 ) ) ) 11 ) )";
	string preorder_mass3 = "0.782 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1 -1";
	
	recopdgs = { 22, 22, 22, 22, 22, 22, 22};
	initializerecoparts(recopdgs);
	printParticles(recoparts);

	tree3->treeInit(preorder_pdg3, preorder_serial3, preorder_mass3, delimiter, 3);
	tree3->printTree(tree3->Root);
	cout<<"tree constructed"<<endl;
	LASTNONLEAFID = tree3->lastNonLeafNodeId;
	ParticleTree=tree3;
	generatefitcombinations(tree3->Root, recoIDs);
	//delete tree;
	//tree=NULL;
	//theTree = NULL;
	recopdgs.clear();
	recoparts.clear();
	cout<<"DONE";
	cout<<std::endl;
*/
//	

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*	cout<<"CASE 3pi0 + photon "<<std::endl;
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
	tree4->printTree(tree4->Root);
	cout<<"tree constructed"<<endl;
	LASTNONLEAFID = tree4->lastNonLeafNodeId;
	ParticleTree=tree4;
	generatefitcombinations(tree4->Root, recoIDs);
	//delete tree;
	//tree=NULL;
	//theTree = NULL;
	recopdgs.clear();
	recoparts.clear();
	cout<<std::endl;
*/
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
	tree5->printTree(tree5->Root);
	LASTNONLEAFID = tree5->lastNonLeafNodeId;
	ParticleTree=tree5;
	generatefitcombinations(tree5->Root, recoIDs);

	//delete tree;
	//tree = NULL;
	//ParticleTree = NULL; 
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
	//printParticles(recoparts);
	
	//tree fitting combination starting with simplest cases
	tree6->treeInit(preorder_pdg6, preorder_serial6, preorder_mass6, delimiter, 6);
	std::cout<<std::endl;
	tree6->printTree(tree6->Root);
	LASTNONLEAFID = tree6->lastNonLeafNodeId;
	ParticleTree=tree6;
	generatefitcombinations(tree6->Root, recoIDs);
	
	//delete tree;
	//tree = NULL;
	//ParticleTree = NULL; 
	recopdgs.clear();
	recoparts.clear();

if(ParticleTree != NULL){
//delete ParticleTree;
}

return 0;
}
*/


