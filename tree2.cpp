#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> //atof atoi

using namespace std;

//////globals///////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////////////////////////////////////////////////////////////////

class Node{
	public:
	vector<Node*> children;
	int pdg;
	double mass; // -1 for leaf or no constraint at this node
	bool isLeaf;

	//unique ID from key
	int nodeId;
	
	//store all combinations of final state particles
	//indices of initial particle array
	vector<int> particles;

	//the number of final state particles from this particular decay node
	int nLeaves=0;
	//the pdgcodes of all the final state particles for this
	vector<int> leafpdgs;
	
};
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
	vector<string> vectokens;
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
//preorder print all tree information
void printTree(Node* root){
	cout<< "NodeId: "<< root->nodeId <<" Pdg: "<< root->pdg <<" Mass: "<< root->mass << " isLeaf= "<<root->isLeaf <<" Children { ";
	for(int i=0; i<root->children.size(); i++){
		cout<< root->children.at(i)->nodeId << " ";
	}
	cout<< "}"<<" nLeaves= "<< root->nLeaves << " Leaf Pdgs: { " ;
	printvector(root->leafpdgs);
	cout<< "}"<<endl;
	for(int i=0; i<root->children.size(); i++){
		printTree(root->children.at(i));
	}
}
//should test return a tree?
void testTree(string pdg, string serial, string mass, string delimiter, int TESTNUM){
	int index = 0;
	Node* root = constructTree( splitString(serial,delimiter), &index);
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
	printTree(root);
	
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
//here remove combinations that contain used, particles, and wrong pdgs
vector<vector<int> > filterCombinations(){

	//the three main conditions for selecting a valid combination of particles is

		// 1- the set contains the same pdgs as the decay is supposed to have 

		// 2- none of the particles in the set are marked as used from previous sets

		// 3- none of the particles are duplicates
			//require E/p uniqueness

	


}
///////////////////////////////////////////////////End Combinatorics///////////////////////////////////////////////
int main(){ 
	string preorder_pdg = " 443 331 321 -321 221 211 -211 111 22 22 ";
	string preorder_key = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial = "0 1 2 ) 3 ) ) 4 5 ) 6 ) 7 8 ) 9 ) ) ) )";
	string preorder_mass = " 3.096 1.0195 -1 -1 0.547 -1 -1 0.135 -1 -1 ";
	string delimiter = " ,";
	
	testTree(preorder_pdg, preorder_serial, preorder_mass, delimiter, 1);
	
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//test tree #2
	string preorder_pdg2 = " 443 221 111 22 22 211 -211 331 321 -321";
	string preorder_key2 = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial2 = "0 1 2 3 ) 4 ) ) 5 ) 6 ) ) 7 8 ) 9 ) ) )";
	string preorder_mass2 = " 3.096 0.547 0.135 -1 -1 -1 -1 1.0195 -1 -1";
	
	testTree(preorder_pdg2, preorder_serial2, preorder_mass2, delimiter, 2);
	
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg3 = " 223 221 111 22 22 111 22 22 111 22 22 22";
	string preorder_key3 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial3 ="0 1 2 3 ) 4 ) ) 5 6 ) 7 ) ) 8 9 ) 10 ) ) ) 11 ) )";
	string preorder_mass3 = "0.782 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1 -1";
	
	testTree(preorder_pdg3, preorder_serial3, preorder_mass3, delimiter, 3);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	
	string preorder_pdg4 = " 223 22 221 111 22 22 111 22 22 111 22 22";
	string preorder_key4 = "0 1 2 3 4 5 6 7 8 9 10 11";
	string preorder_serial4 = "0 1 ) 2 3 4 ) 5 ) ) 6 7 ) 8 ) ) 9 10 ) 11 ) ) ) )";
	string preorder_mass4 = "0.782 -1 0.547 0.135 -1 -1 0.135 -1 -1 0.135 -1 -1";

	testTree(preorder_pdg4, preorder_serial4, preorder_mass4, delimiter, 4);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg5 = " 221 211 -211 111 22 22";
	string preorder_key5 = "0 1 2 3 4 5";
	string preorder_serial5 = "0 1 ) 2 ) 3 4 ) 5 ) ) )";
	string preorder_mass5 = "0.547 -1 -1 0.135 -1 -1";

	testTree(preorder_pdg5, preorder_serial5, preorder_mass5, delimiter, 5);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	string preorder_pdg6 = " 221 111 22 22 211 -211";
	string preorder_key6 = "0 1 2 3 4 5";
	string preorder_serial6 = "0 1 2 ) 3 ) ) 4 ) 5 ) )";
	string preorder_mass6 = "0.537 0.135 -1 -1 -1 -1";

	testTree(preorder_pdg6, preorder_serial6, preorder_mass6, delimiter, 6);

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
	//preorder combination generation

}
