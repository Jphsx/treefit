#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h> //atof atoi

using namespace std;

//////globals///////////////////////////////////////////////////////////////////////////////////
//vector marking particles as used or free
vector<int> usedparts{};

//the set of pdgs from the reconstructed set of particles
vector<int> recopdgs{};
int maxusedparts=-1;



///////////////////////////////////////////////////////////////////////////////////////////////

class Node{
	public:
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
	vector<int> usedpartsindicesflagged{};
	
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
//input the root of the fully constructed tree to determine the maximum number of particles going to be used in a fit
void setmaxusedparts(Node* root){
	maxusedparts = root->nLeaves;
	return;
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


void generateSubsets(std::vector<int> v, int k, int start, int currLen, std::vector<bool> used, std::vector<std::vector<int> >& combinations){
	if(currLen == k) {
		std::vector<int> indices;
		for(unsigned int i=0; i<v.size(); i++){
			if(used.at(i) == true) {
				indices.push_back(i);
			}
		}
		combinations.push_back(indices);
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
void generateIndicesCombinations(int vectorsize, int nparticles, std::vector<std::vector<int> >& combinations){
	int n = vectorsize;
	int k = nparticles;

		std::vector<int> master;
		std::vector<bool> used(n);
	
		for(int i=0; i<n; i++){
			master.push_back(i);
			used.at(i)=false;
		}


	generateSubsets(master,k,0,0,used,combinations);

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
//combo is index vector
bool containsusedparticle(vector<int> combo){
	for(int i=0; i<combo.size(); i++){
		//if a particle in the combo is used, return true
		if(usedparts.at(combo.at(i))) return true; 
	}
	return false;
}
//does the combo have particles in it that are not in the needed final state? input is pdgs mapped from indices
bool containsfinalstateparticle(vector<int> fsp_pdgs, vector<int> combo_pdgs){
	
	for(int i=0; i<combo_pdgs.size(); i++){
		if(std::find(fsp_pdgs.begin(), fsp_pdgs.end(), combo_pdgs.at(i)) != fsp_pdgs.end()) {
    			/* fsp contains particle */
		} else {
    			/* fsp does not contain particle */
			return false;
		}
	}
	//if we havent returned false yet then combo is good, return true
	return true;
}
bool subsetofparentcombo(parentcombo,combinations.at(i)){

}
//combo argument is  pdgs, is the combination of pdgs the same as the required combination?
bool finalstatepdgmatch(vector<int>& fsp_pdgs, vector<int>& combo_pdgs){

	//create a copy of both arrays, we don't want the originals sorted because we need to preserve indices
	vector<int> fsp_copy{};
	vector<int> combo_copy{};
	//the vectors had better be the same size
	for(int i=0; i<fsp_pdgs.size(); i++){ 
		fsp_copy.push_back(fsp_pdgs.at(i));
		combo_copy.push_back(combo_pdgs.at(i));
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
void filtercombinations(vector<int> fsp, vector<vector<int> >& combinations, vector<vector<int> >& pdgcombinations, vector<int> parentcombo ){

	vector<vector<int> > filteredcombos;
	vector<vector<int> > filteredpdgcombos;
//the three main conditions for selecting a valid combination of particles is

	for(int i=0; i<combinations.size(); i++){

		// 1- the set contains the same pdgs as the decay is supposed to have 
		// this may be redundant (can only use 1-b to save time)
		if(!containsfinalstateparticle(fsp,pdgcombinations.at(i))) continue;

		// 1-b the set contains exactly the same composition of pdgs
		if(!finalstatepdgmatch(fsp,pdgcombinations.at(i))) continue;

		//the combinations are a subset of the parent set of indices
		if(!subsetofparentcombo(parentcombo,combinations.at(i))) continue;
		
		// 2- none of the particles in the set are marked as used from previous sets
		if(containsusedparticle(combinations.at(i))) continue;//used particle found reject combo

		// 3- none of the particles are duplicates
			//require E/p uniqueness
		//if all the tests are passed push the combination onto the new object
		filteredcombos.push_back(combinations.at(i));
		filteredpdgcombos.push_back(pdgcombinations.at(i));

	}

	combinations.clear();
	pdgcombinations.clear();//
	//maybe combinations deleted when out of scope,so, try a hard copy
	//TODO change this back i think it was okay
	for(int i=0; i<filteredcombos.size(); i++){
		combinations.push_back(filteredcombos.at(i));
		pdgcombinations.push_back(filteredpdgcombos.at(i));
	}
}
//map indices to pdgs
//vector<vector<int> >
void mapindextopdg( vector<vector<int> >& combinations, vector<vector<int> >& pdgcombinations){
	vector<vector<int> > pdgcombos;

		for(int i=0; i<combinations.size(); i++){
			vector<int> subarray;
			for(int j=0; j<combinations.at(i).size(); j++){
				subarray.push_back( recopdgs.at( combinations.at(i).at(j) ) );
			}
			pdgcombos.push_back(subarray);
			subarray.clear();
		}
		//return pdgcombos;
		pdgcombinations = pdgcombos;	
	
}
///////////////////////////////////////////////////End Combinatorics///////////////////////////////////////////////
//////////////////////////////////////////////////Begin generating fit combinations from reco particles///////////
void initializeusedvector(int size){
	for(int i=0; i<size; i++){
		usedparts.push_back(0);
	}
	return;
}
bool allpartsused(){
	//count the number of used parts, if it is equal to maxused parts return true
	int counter=0;
	for(int i=0; i<usedparts.size(); i++){
		if(usedparts.at(i) == 1) counter++;
	}
	if(counter == maxusedparts){
		return true;
	}
	else{
		return false;
	}
}
//preorder traverse and print current combinations and node information
void printfit(Node* root){
	
	if(root->isLeaf) return;

	cout<< "NodeId: "<< root->nodeId <<" Pdg: "<< root->pdg << " FSPs Indices: ";
	printvector(root->currentcombination);
	cout<< " FSPs Pdgs: ";
	printvector(root->currentcombination_pdgs);
	cout<<endl;

	for(int i=0; i<root->children.size(); i++){
		printfit(root->children.at(i));
	}
	return;
	
}
//preorder combination generation recursively
void generatefitcombinations(Node* root, vector<int> parentcombo){

	if(root->isLeaf) return;
	//first generate all combinations for this node
	//arguments(n,k,vector) n choose k and put index combos onto vector
	//TODO what if mass=-1?? we dont want to make combinations

	/*generateIndicesCombinations(recopdgs.size(),root->nLeaves,root->combinations);
	mapindextopdg( root->combinations, root->combinations_pdgs);
	filtercombinations(root->leafpdgs, root->combinations, root->combinations_pdgs );
	*/
	generateIndicesCombinations(recopdgs.size(), root->nLeaves, root->combinations);
	mapindextopdg( root->combinations, root->combinations_pdgs);
	filtercombinations(root->leafpdgs, root->combinations, root->combinations_pdgs, parentcombo);
	
	//iterate through combinations and through all children recursively
	//for(int i=0; i
	for(int i=0; i<root->combinations.size(); i++){
		//first set the current combination
		root->currentcombination = root->combinations.at(i);
		root->currentcombination_pdgs = root->combinations_pdgs.at(i);

		//for the children that are leaves, we will mark them as used for this combo
		for(int j=0; j<root->childrenleaves.size(); j++){
			//if element==false no flag
			//if the element is nonzero it is a leaf, so we need to mark
			if(root->childrenleaves.at(j) != 0){
				
				//find the leaf that corresponds to index in combinations and usedparts
				for(int x=0; x<root->currentcombination.size(); x++){
					if(!usedparts.at(root->currentcombination.at(x)) && recopdgs.at(root->currentcombination.at(x)) == root->childrenleaves.at(j) ){				
						usedparts.at(root->currentcombination.at(x))=1;
						root->usedpartsindicesflagged.push_back(root->currentcombination.at(x));

						
					}
				}
				
				//when a particle is marked check and see if we need to fit
				//all parts used needs changed, we need to check and see if n flags = n leaves of the root
				if(allpartsused()){
					cout<<"performing fit"<<endl;
					//TODO printfit function;
					printfit(theTree);

					for(int l=0; l<root->usedpartsindicesflagged.size(); l++){
						usedparts.at(root->usedpartsindicesflagged.at(l))=0;
					}
					root->usedpartsindicesflagged.clear();
					break;
				}
			}
		}
	
		for(int k=0; k<root->children.size(); k++){
			//recurse through the children
			//generatefitcombinations(root->children.at(k));
			generatefitcombinations(root->children.at(k), root->currentcombination);
		}
		//before we proceed to the next combination unmark the leaves previously used
		root->currentcombination.clear();
		root->currentcombination_pdgs.clear();
		
	//	for(int l=0; l<root->usedpartsindicesflagged.size(); l++){
	//		usedparts.at(root->usedpartsindicesflagged.at(l))=0;
	//	}
	//	root->usedpartsindicesflagged.clear();
		
	}
	root->combinations.clear();
	root->combinations_pdgs.clear();
	return;
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
	setmaxusedparts(root);
        generatefitcombinations(root, recopdgs);

	return root;
	
}
///////////////////////////////////////////////////End fitting combinations//////////////////////////////////////
int main(){ 

	Node* tree;
	string delimiter = " ,";
/*	string preorder_pdg = " 443 331 321 -321 221 211 -211 111 22 22 ";
	string preorder_key = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial = "0 1 2 ) 3 ) ) 4 5 ) 6 ) 7 8 ) 9 ) ) ) )";
	string preorder_mass = " 3.096 1.0195 -1 -1 0.547 -1 -1 0.135 -1 -1 ";
	
	
	tree = testTree(preorder_pdg, preorder_serial, preorder_mass, delimiter, 1);
	delete tree;
	tree=NULL;
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
///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	string preorder_pdg4 = " 223 22 221 111 22 22 111 22 22 111 22 22";
	string preorder_key4 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial4 = "0 1 ) 2 3 4 ) 5 ) ) 6 7 ) 8 ) ) 9 10 ) 11 ) ) ) )";
	string preorder_mass4 = "0.782 -1 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1";

	tree = testTree(preorder_pdg4, preorder_serial4, preorder_mass4, delimiter, 4);
	delete tree;
	tree=NULL;
///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg5 = " 221 211 -211 111 22 22";
	string preorder_key5 = "0 1 2 3 4 5";
	string preorder_serial5 = "0 1 ) 2 ) 3 4 ) 5 ) ) )";
	string preorder_mass5 = "0.547 -1 -1 0.135 -1 -1";

	tree = testTree(preorder_pdg5, preorder_serial5, preorder_mass5, delimiter, 5);
	delete tree;
	tree = NULL;

*/
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg6 = " 221 111 22 22 211 -211";
	string preorder_key6 = "0 1 2 3 4 5";
	string preorder_serial6 = "0 1 2 ) 3 ) ) 4 ) 5 ) )";
	string preorder_mass6 = "0.537 0.135 -1 -1 -1 -1";

	
	
	//delete tree;
	//tree = NULL;

	//preorder combination generation
	//define the reconstucted particles (by pdg)
	// using tree 6, eta -> pi+ pi- pi0
	//arbitrary reconstructed particles: 4 tracks 3 photons 7 particles total
	recopdgs = { 211, -211, 13, -13, 22, 22, 22};
	//start by setting global array for all particles unused
	initializeusedvector(recopdgs.size());
		

	//tree fitting combination starting with simplest cases
	tree = testTree(preorder_pdg6, preorder_serial6, preorder_mass6, delimiter, 6);
	//take the tree and fit combinations
		
}
