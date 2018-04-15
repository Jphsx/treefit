#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h> //atof atoi

using namespace std;
//this number should be set when the tree is built
//int TOTAL_NON_LEAF_NODES = -1;

//this number is incremented when generating fit combinations
//int CURRENT_NON_LEAF_NODES=0;

//When TOTAL == CURRENT we should perform a fit

//////globals///////////////////////////////////////////////////////////////////////////////////
//a general particle container that aggregates all possible information or a reconstructed particle
class Particle{
	public:
	int recopdg;
	//some particle object should be in here as well
	//ReconstructedParticle/Track*
	//TLorentzVector RECO
	//MCParticle
	//TLorentzVector MC
	//LocalParamterization 
	//LocalParameterization Errors
	//used flag
	bool used = 0;

};
//the global object that contains all of the reconstructed particles

vector<Particle*> recoparts{};
//this function should probably takin in a TLV or ReconstructedParticle
void initializerecoparts(vector<int> recopdgs){
	for(int i=0; i<recopdgs.size(); i++){
		Particle* p = new Particle();
		p->recopdg=recopdgs.at(i);
		recoparts.push_back(p);
	}
}
///////////////////////////////////////////////////////////////////////////////////////////////

class Node{
	public:
	Node* parent = NULL;
	vector<Node*> children{};
	int pdg = -1;
	double mass = -1.0; // -1 for leaf or no constraint at this node
	bool isLeaf = 0;

	//unique ID from key
	int nodeId = -1;
	
	//store all combinations of final state particles
	//indices of initial particle array
	vector<vector<int> > combinations{};
	vector<vector<int> > combinations_pdgs{};
	
	//place holder for current set of particles for this node
	vector<int> currentcombination{};
	vector<int> currentcombination_pdgs{};

	//the number of final state particles from this particular decay node
	int nLeaves=0;
	//the pdgcodes of all the final state particles for this
	vector<int> leafpdgs{};

	//helper data structures
	//quick access to see which children are leaves
	//the values are 0 for no a leaf, otherwise the pdgcode is the element
	//TODO change the values from 1's
	vector<int> childrenleaves{};
	//for each fit combination, need to store which indices that are flagged used
	//so they can be unflagged in each iteration
	//the incices used is the INDEX OF THE FLAG 
	//vector<int> usedpartsindicesflagged{};
	
};

/////the global tree
Node* theTree = NULL;
////////////
Node* newNode (int id){
    Node* temp = new Node();
    temp->nodeId=id;
    return temp;
}

template <class type>
void printvector(vector<type> v){
	for(int i=0; i<v.size(); i++){
			cout<<v.at(i)<<" ";
	}
	//cout<<endl;
}
template <class type>
void print2dvec(vector<vector<type> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}
}
vector<double> castVector_double(vector<string> v){
	vector<double> v_double;
	for(int i=0; i<v.size(); i++){
		v_double.push_back(atof( v.at(i).c_str() ));
	}
	return v_double;
}
vector<int> castVector_int(vector<string> v){
	vector<int> v_int;	
	for(int i=0; i<v.size(); i++){
		v_int.push_back(atoi( v.at(i).c_str() ));
	}
	return v_int;
}
vector<string> splitString(string str, string delimiter){
	vector<string> vectokens{};
	char* dup = strdup(str.c_str());
	char* dup_delim = strdup(delimiter.c_str());
	char* tokens = strtok(dup,dup_delim);
	while( tokens != NULL){
		vectokens.push_back(string(tokens));
		tokens = strtok(NULL,dup_delim);
	}
	return vectokens;
}
Node* constructTree(vector<string> serial, int* serialIndex){
		
	//cout<<serial.at(*serialIndex)<<" "<<*serialIndex<<endl;
	//need a base case for return
	if(*serialIndex > serial.size() ) return NULL;
	
//	cout<<atoi(serial.at(*serialIndex).c_str())<<endl;
	Node* n = newNode ( atoi(serial.at(*serialIndex).c_str()) );
    	++*serialIndex;
	//cout<<"made a node "<< n->nodeId<<endl;
	
	//if this is not the end of a branch then make a child and recurse
	while( serial.at(*serialIndex).find(")") == std::string::npos ){
	//	cout<<"made a child for "<< n->nodeId<<endl;
		
		n->children.push_back( constructTree(serial, serialIndex) );
	//	cout<<"child pushed onto "<<n->nodeId<<endl;
	}
	
	++*serialIndex;
	return n;
	
	
}
void preOrderTraverse(Node* tree){
	cout<<tree->nodeId<<" ";
	for(int i=0; i<tree->children.size(); i++){
		preOrderTraverse(tree->children.at(i));
	}
	return;
}
//goes right to left
void postOrderTraverse(Node* tree){
	for(int i=0; i<tree->children.size(); i++){
		postOrderTraverse(tree->children.at(i));
	}
	cout<<tree->nodeId<<" ";
	return;
}
//a preorder traversal that sets the mass constraint mass for each node based on the input string
void setTreeMasses(Node* root, vector<double> masses, int* massptr ){
	root->mass = masses.at(*massptr);
	++*massptr;
	for(int i=0; i< root->children.size(); i++){
		setTreeMasses(root->children.at(i), masses, massptr);
	}
	return;
}
//a preorder traversal that sets the pdg codes for each node based on the input string
void setTreePdgCodes(Node* root, vector<int> pdgs, int* pdgptr ){
	root->pdg = pdgs.at(*pdgptr);
	++*pdgptr;
	for(int i=0; i< root->children.size(); i++){
		setTreePdgCodes(root->children.at(i), pdgs, pdgptr);
	}
	return;
}
void markTreeLeaves(Node* root){
	if(root->children.size() == 0){
		root->isLeaf=true;
		return;
	}
	else{
		root->isLeaf=false;
	}
	for(int i=0; i< root->children.size(); i++){
		markTreeLeaves(root->children.at(i));
	}
	return;
}
//find how many leaves, and the pdg of the leaves of a particular node
void findLeaves(Node* root, Node* originalParent){
	if(root->isLeaf){
		 originalParent->nLeaves++;
		 originalParent->leafpdgs.push_back( root->pdg );
	}
	for(int i=0; i<root->children.size(); i++){
		findLeaves(root->children.at(i), originalParent);
	}
	return;
		
}
void populateNLeaves(Node* root){
	//int counter= 0;
	findLeaves(root, root);
	//cout<<"found the leaves counter = "<<counter<<endl;
	//root->nLeaves=counter;
	for(int i=0; i<root->children.size(); i++){
		populateNLeaves(root->children.at(i));
	}
	return;
}
void countnonleafnodes(Node* root, int* numnonleafnodes){
		if(root->isLeaf) return;
		++*numnonleafnodes;
		for(int i=0; i<root->children.size(); i++){
			countnonleafnodes(root->children.at(i), numnonleafnodes);
		}
		return;
}
//TODO this function
void setChildrenLeaves(Node* root){

	//do a preorder traversal, check all the children at a node to see if hey are leaves
	for(int i=0; i<root->children.size(); i++){
		if(root->children.at(i)->isLeaf){
			root->childrenleaves.push_back(root->children.at(i)->pdg);
		}
		else{
			root->childrenleaves.push_back(0);
		}
		setChildrenLeaves(root->children.at(i));
	}
}
void getLastNonLeafNodeId(Node* root, int* id){
	if(!root->isLeaf){
		cout<<"this node "<<root->nodeId<<endl;
		*id = root->nodeId;
	}
	for(int i=0; i<root->children.size(); i++){
		getLastNonLeafNodeId(root->children.at(i), id);
	}
}
//make sure the root is predefined by passing in NULL
void setParents(Node* root, Node* parentptr){
	root->parent = parentptr;
	for(int i=0; i<root->children.size(); i++){
		setParents(root->children.at(i), root);
	}
}
void printNodeAncestors(Node* root){
	
}
//preorder print all tree information
void printTree(Node* root){
	cout<< "NodeId: "<< root->nodeId <<" Pdg: "<< root->pdg <<" Mass: "<< root->mass << " isLeaf= "<<root->isLeaf <<" Children { ";
	for(int i=0; i<root->children.size(); i++){
		cout<< root->children.at(i)->nodeId << " ";
	}
	cout<< "}"<<" nLeaves= "<< root->nLeaves << " Leaf Pdgs: { " ;
	printvector(root->leafpdgs);
	cout<< "}"<<" childrenleaves: ";
	printvector(root->childrenleaves);
	cout<<endl;
	for(int i=0; i<root->children.size(); i++){
		printTree(root->children.at(i));
	}
}

////////////////////////////////////////////////////END TREE STUFF///////////////////////////////////////////////
////////////////////////////////////////////////////Begin Combinatorics//////////////////////////////////////////


template <typename type>
void generateSubsets(std::vector<type> v, int k, int start, int currLen, std::vector<bool> used, std::vector<std::vector<type> >& combinations){
	if(currLen == k) {
		std::vector<type> items;
		for(unsigned int i=0; i<v.size(); i++){
			if(used.at(i) == true) {
				items.push_back(v.at(i));
			}
		}
		combinations.push_back(items);
		return;
	}
	if (start == v.size()){
		return;
	}
	used[start] = true;
	generateSubsets(v,k,start+1,currLen+1,used, combinations);

	used[start] = false;
	generateSubsets(v,k,start+1,currLen,used, combinations);
}
//n choose k
//this will not work with complex objects like string
template<typename type>
void generateIndicesCombinations(int vectorsize, int nparticles, std::vector<std::vector<type> >& combinations, vector<type> recoset){
	int n = vectorsize;
	int k = nparticles;

		std::vector<bool> used(n);	

		for(int i=0; i<n; i++){
			used.at(i)=false;
		}

	generateSubsets(recoset,k,0,0,used,combinations);
}
/* C implementation QuickSort */
//implementation modified and copied from 
//http://www.geeksforgeeks.org/quick-sort/
 
// A utility function to swap two elements
template <class type>
void swap(type* a, type* b)
{
    type t = *a;
    *a = *b;
    *b = t;
}
 
/* This function takes last element as pivot, places
   the pivot element at its correct position in sorted
    array, and places all smaller (smaller than pivot)
   to left of pivot and all greater elements to right
   of pivot */
template <class type>
int partition (vector<type>& arr, int low, int high)
{
    int pivot = arr.at(high);    // pivot
    int i = (low - 1);  // Index of smaller element
 
    for (int j = low; j <= high- 1; j++)
    {
        // If current element is smaller than or
        // equal to pivot
        if (arr.at(j) <= pivot)
        {
            i++;    // increment index of smaller element
            swap(&arr.at(i), &arr.at(j));
        }
    }
    swap(&arr.at(i + 1), &arr.at(high));
    return (i + 1);
}
 
/* The main function that implements QuickSort
 arr --> Array to be sorted,
  low  --> Starting index,
  high  --> Ending index */
template <class type>
void quickSort(vector<type>& arr, int low, int high)
{
    if (low < high)
    {
        /* pi is partitioning index, arr[p] is now
           at right place */
        int pi = partition(arr, low, high);
 
        // Separately sort elements before
        // partition and after partition
        quickSort(arr, low, pi - 1);
        quickSort(arr, pi + 1, high);
    }
}
//does the combo have particles in it that are not in the needed final state? input is pdgs mapped from indices
bool containsfinalstateparticle(vector<int> fsp_pdgs, vector<int> combo){
	
	for(int i=0; i<combo.size(); i++){
		if(std::find(fsp_pdgs.begin(), fsp_pdgs.end(), recoparts.at(combo.at(i))->recopdg) != fsp_pdgs.end()) {
    			/* fsp contains particle */
		} else {
    			/* fsp does not contain particle */
			return false;
		}
	}
	//if we havent returned false yet then combo is good, return true
	return true;
}
//same as containsfinal state particle, but diffrenent names for code clarity
/*bool subsetofparentcombo(vector<int>& parentcurrentcombo,vector<int>& combination){
	for(int i=0; i<combination.size(); i++){
		if(std::find(parentcurrentcombo.begin(), parentcurrentcombo.end(), combination.at(i)) != parentcurrentcombo.end()) {
    			// fsp contains particle 
		} else {
    			// fsp does not contain particle 
			return false;
		}
	}
	return true;
}*/
//combo argument is  pdgs, is the combination of pdgs the same as the required combination?
bool finalstatepdgmatch(vector<int> fsp_pdgs, vector<int> combo){

	//create a copy of both arrays, we don't want the originals sorted because we need to preserve indices
	vector<int> fsp_copy{};
	fsp_copy = fsp_pdgs;	//this is a deep copy	

	//TODO Make an array of pdgs by extracting the pdg values from the combo array
	vector<int> combo_copy{};// combo copy is the pdgs from each of the combo indices
	//the vectors had better be the same size
	for(int i=0; i<combo.size(); i++){ 
		//fsp_copy.push_back(fsp_pdgs.at(i));
		//combo_copy.push_back(combo.at(i));
		combo_copy.push_back(recoparts.at(combo.at(i))->recopdg);
		
	}
	
	//create a matched vector of boolean values
	//get fsp pdgs and combo pdgs
	//sort both arrays using quicksort
	//if both arrays have the same pdg content then the two arrays should be identical after sorting
	quickSort(fsp_copy, 0, fsp_copy.size()-1);
	quickSort(combo_copy, 0, combo_copy.size()-1);

	for(int i=0; i< fsp_copy.size(); i++){
		if( fsp_copy.at(i) != combo_copy.at(i) ) return false;
	}
	
	return true;
		
}
void filtercombinations(vector<int> fsp, vector<vector<int> >& combinations){//, vector<int> parentcombo ){

	vector<vector<int> > filteredcombos;
	//vector<vector<int> > filteredpdgcombos;
//the three main conditions for selecting a valid combination of particles is

	for(int i=0; i<combinations.size(); i++){

		// 1- the set contains the same pdgs as the decay is supposed to have 
		// this may be redundant (can only use 1-b to save time)
		if(!containsfinalstateparticle(fsp,combinations.at(i))) continue;

		// 1-b the set contains exactly the same composition of pdgs
		if(!finalstatepdgmatch(fsp,combinations.at(i))) continue;

		// 1-c the combinations are a subset of the parent set of indices
		//his is obsolete because combinations are generated from he parent combo
	//	if(!subsetofparentcombo(parentcombo,combinations.at(i))) continue;
		
		// 2- none of the particles in the set are marked as used from previous sets
		//if(containsusedparticle(combinations.at(i))) continue;//used particle found reject combo

		// 3- none of the particles are duplicates
			//require E/p uniqueness
		//if all the tests are passed push the combination onto the new object
		filteredcombos.push_back(combinations.at(i));
		//filteredpdgcombos.push_back(pdgcombinations.at(i));

	}

	combinations.clear();
	//pdgcombinations.clear();//
	//maybe combinations deleted when out of scope,so, try a hard copy
	//TODO change this back i think it was okay
	for(int i=0; i<filteredcombos.size(); i++){
		combinations.push_back(filteredcombos.at(i));
		//pdgcombinations.push_back(filteredpdgcombos.at(i));
	}
}

///////////////////////////////////////////////////End Combinatorics///////////////////////////////////////////////
//////////////////////////////////////////////////Begin generating fit combinations from reco particles///////////


//TODO
//function that counts node combinations chosen
//preorder traverse and print current combinations and node information
void printfit(Node* root){
	
	if(root->isLeaf) return;
	cout<< "NodeId: "<< root->nodeId <<" Pdg: "<< root->pdg << " FSPs Indices: ";
	//these are indices of the global paricle vector
	printvector(root->currentcombination);
	cout<< " FSPs Pdgs: ";
	for(int i=0; i< root->currentcombination.size(); i++){
		cout<< recoparts.at(root->currentcombination.at(i))->recopdg <<" ";
	}
	cout<<endl;

	for(int i=0; i<root->children.size(); i++){
		printfit(root->children.at(i));
	}
	return;
	
}
//preorder combination generation recursively
//first call pass in reco indices
//void generatefitcombinations(Node* root, vector<int> parentcombo, int* CURRENT_NON_LEAF_NODES, const int TOTAL_NON_LEAF_NODES ){
void generatefitcombinations(Node* root, vector<int> parentcombo, const int last_non_leaf_id){

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
}

vector<int> makerecoindices(){
	vector<int> recoindices{};
	for(int i=0; i<recoparts.size(); i++){
		recoindices.push_back(i);
	}
	return recoindices;
}
Node* testTree(string pdg, string serial, string mass, string delimiter, int TESTNUM){
	int index = 0;
	Node* root = new Node();
	root = constructTree( splitString(serial,delimiter), &index);
	theTree = root;
	cout<<"TREE "<<TESTNUM<<endl;
	preOrderTraverse(root);
	cout<<endl;
	postOrderTraverse(root);
	cout<<endl;

	index=0;
	setTreeMasses(root, castVector_double(splitString(mass,delimiter)), &index);
	index=0;
	setTreePdgCodes(root, castVector_int(splitString(pdg,delimiter)), &index);
	markTreeLeaves(root);
	populateNLeaves(root);
	setChildrenLeaves(root);
	printTree(root);
	//setmaxusedparts(root);
	//int currentnodes = 0;
	//int totalnodes = 0;//TODO make a count of non leaf nodes
	//countnonleafnodes(root,&totalnodes);

	
	int last_node_id=0;
	cout<<last_node_id<<" last "<<std::endl;
	getLastNonLeafNodeId(root , &last_node_id);
	cout<<last_node_id<<" last "<<std::endl;
	vector<int> recoindices = makerecoindices();
        //generatefitcombinations(root, recoindices, &currentnodes, totalnodes);
	generatefitcombinations(root, recoindices, last_node_id);
	return root;
	
}
///////////////////////////////////////////////////End fitting combinations//////////////////////////////////////
int main(){ 

	Node* tree;
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

	tree = testTree(preorder_pdg5, preorder_serial5, preorder_mass5, delimiter, 5);
	delete tree;
	tree = NULL;
	theTree = NULL; 
	recopdgs.clear();
	recoparts.clear();

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
*?
}
