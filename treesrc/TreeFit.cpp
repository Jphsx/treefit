
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
void TreeFit::printParticles(std::vector<Particle*> parts){
	for(int i=0; i<parts.size(); i++){
		if(parts.at(i) != NULL){
			std::cout<<"Particle Index/recoID "<<i<<std::endl;
			Particle::printParticle(parts.at(i));
		}
	}
}
std::vector<std::vector<int> > TreeFit::makepdgcombinations(std::vector<std::vector<int> > combinations){
	std::vector<std::vector<int> > pdgcombinations;
	std::vector<int> pdgcombo;
	for(int i=0; i<combinations.size();i++){
		for(int j=0; j<combinations.at(i).size(); j++){
			pdgcombo.push_back(recoparts.at(combinations.at(i).at(j))->recopdg);
		}
		pdgcombinations.push_back(pdgcombo);
		pdgcombo.clear();
	}
	return pdgcombinations;
}
std::vector<int> TreeFit::getpdgcombo(std::vector<int> combo){
	std::vector<int> pdgcombo;
	for(int i=0; i<combo.size(); i++){
		pdgcombo.push_back(recoparts.at(combo.at(i))->recopdg);
	}
	return pdgcombo;	
}
void TreeFit::generatefitcombinations(Node* root, std::vector<int> parentcombo){

	
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
			addFitToTable(ParticleTree->Root);
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
	fitTable.clear();
	fitPdgs.clear();
	fitparts.clear();
}
void TreeFit::addFitToTable(Node* root){
	if(root->isLeaf) return;// testing new imp, if is leaf  put on a empty table

	//if(root->isLeaf){
		
	//}
	//else{
		fitTable.at(root->nodeId).push_back(root->currentcombination);
		fitPdgs.at(root->nodeId).push_back(root->currentcombination_pdgs);
	//}
	for(int i=0; i< root->children.size(); i++){
		addFitToTable(root->children.at(i));
	}

}
void TreeFit::initTable(){
	std::vector< std::vector< std::vector<int>>> idtable(*LASTNONLEAFID + 1); 
	std::vector< std::vector< std::vector<int>>> pdgtable(*LASTNONLEAFID + 1);

	fitTable = idtable;
	fitPdgs = pdgtable;
}
void TreeFit::printTable(){
	int tpdg=0;
	for(int j=0; j< fitTable.at(0).size(); j++){
		std::cout<<" FIT: "<<j <<std::endl;
		for(int i=0; i< fitTable.size(); i++){
		//we require at least 1 node to fit (node 0)
			
			if(fitTable.at(i).size() > 0){
				//this is not a leaf node
				Tree::getNodePdg(ParticleTree->Root, i, &tpdg);
				std::cout<<"Node: "<< i<<" Pdg: "<< tpdg <<" ";
				
				std::cout<<"Final State RecoIDs: { ";
				for(int k=0; k<fitTable.at(i).at(j).size(); k++){
					std::cout<<fitTable.at(i).at(j).at(k)<<" ";
				}	
				std::cout<<"} Pdgs: {";
				
				for(int k=0; k<fitPdgs.at(i).at(j).size(); k++){
					std::cout<<fitPdgs.at(i).at(j).at(k)<<" ";
				}
				std::cout<<"}"<<std::endl;
			}
			
		}
		std::cout<<std::endl;				
	}
	

}
std::vector<int> TreeFit::getVertexSet(std::vector<int> combo, int nodeId, std::vector<std::vector<int> > fit){
	
	
	std::vector<int> subset = combo;
	//loop over the fit for elements after the given node
	for(int i=nodeId+1; i<fit.size(); i++){
		//if this node is non leaf do set math
		if( fit.at(i).size() > 0 ){
			subset = Combinatorics::subtractSets(subset, fit.at(i));
		}
	}

	return subset;

}
void TreeFit::printComparison(std::vector<Particle*> reco, std::vector<Particle*> fit, std::vector<int> combo){
	//TODO write this
	
	//loop over combo
	for(int i=0; i<combo.size(); i++){
		//print particles at combo indices
		std::cout<<"Particle: "<<reco.at(combo.at(i))->part->getType() <<std::endl;
		std::cout<<"Reco: ";
		Particle::printReconstructedParticle(reco.at(combo.at(i))->part );
		std::cout<<"Fit : ";
		Particle::printReconstructedParticle(fit.at(combo.at(i))->part);
		if(reco.at(combo.at(i))->isTrack){
			std::cout<<"Reco: ";
			Particle::printTrack(reco.at(combo.at(i))->track);
			Particle::printLocalParameters(reco.at(combo.at(i))->localParams);
			Particle::printLocalErrors(reco.at(combo.at(i))->localErrors);
			std::cout<<"Fit : ";
			Particle::printTrack(fit.at(combo.at(i))->track);
			Particle::printLocalParameters(fit.at(combo.at(i))->localParams);
			Particle::printLocalErrors(fit.at(combo.at(i))->localErrors);
			std::cout<<"Reco Error Matrix "<<std::endl;
			Particle::printCovarianceMatrix(reco.at(combo.at(i))->track->getCovMatrix(), reco.at(combo.at(i))->localParams.size());				
			std::cout<<"Fit Error Matrix " <<std::endl;
			Particle::printCovarianceMatrix(fit.at(combo.at(i))->track->getCovMatrix(), fit.at(combo.at(i))->localParams.size());
		}
		else{
			std::cout<<"Reco Error Matrix "<<std::endl;
			Particle::printCovarianceMatrix(reco.at(combo.at(i))->part->getCovMatrix(), reco.at(combo.at(i))->localParams.size());				
			std::cout<<"Fit Error Matrix " <<std::endl;
			Particle::printCovarianceMatrix(fit.at(combo.at(i))->part->getCovMatrix(), fit.at(combo.at(i))->localParams.size());
		}
	
	}
	
}
void TreeFit::printParentComparison(TLorentzVector recoparent, TLorentzVector fitparent, int nodeId, Node* root){
	Node* parentNode = Tree::getNode(root, nodeId);
	std::cout<<"Parent Particle: "<<parentNode->pdg<<std::endl;
	std::cout<<"Reco: ";
	Particle::printTLorentzVector(recoparent);	
	std::cout<<"Fit : ";
	Particle::printTLorentzVector(fitparent);
	std::cout<<std::endl;

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


