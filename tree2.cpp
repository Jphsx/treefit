#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> //atof atoi

using namespace std;
class Node{
	//Node* parent;
	public:
	vector<Node*> children;
	int pdg;
	double mass;
	bool isLeaf;
	//unique ID from key
	int nodeId;
};
Node* newNode (int id){
    Node* temp = new Node();//(struct Node*)malloc(sizeof(Node));
//    cout<<"making node"<<endl;
    temp->nodeId=id;
 //  cout<<"assigned id "<<id<<endl;
    return temp;
}

template <class type>
void printvector(vector<type> v){
	for(int i=0; i<v.size(); i++){
			cout<<v.at(i)<<" ";
	}
	cout<<endl;
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
		
	cout<<serial.at(*serialIndex)<<" "<<*serialIndex<<endl;
	//need a base case for return
	if(*serialIndex > serial.size() ) return NULL;
	
//	cout<<atoi(serial.at(*serialIndex).c_str())<<endl;
	Node* n = newNode ( atoi(serial.at(*serialIndex).c_str()) );
    	++*serialIndex;
	cout<<"made a node "<< n->nodeId<<endl;
	
	//if this is not the end of a branch then make a child and recurse
	while( serial.at(*serialIndex).find(")") == std::string::npos ){
		cout<<"made a child for "<< n->nodeId<<endl;
		
		n->children.push_back( constructTree(serial, serialIndex) );
		cout<<"child pushed onto "<<n->nodeId<<endl;
	}
	
	++*serialIndex;
	return n;
	
	
}
void preOrderTraverse(Node* tree){
	cout<<tree->nodeId<<" ";
	for(int i=0; i<tree->children.size(); i++){
		preOrderTraverse(tree->children.at(i));
	}
}
//goes right to left
void postOrderTraverse(Node* tree){
	for(int i=0; i<tree->children.size(); i++){
		postOrderTraverse(tree->children.at(i));
	}
	cout<<tree->nodeId<<" ";
}
int main(){ 
	string preorder = " 443 331 321 -321 221 211 -211 111 22 22 ";
	string preorder_key = "0 1 2 3 4 5 6 7 8 9";
	string preorder_serial = "0 1 2 ) 3 ) ) 4 5 ) 6 ) 7 8 ) 9 ) ) ) )";
	string preorder_mass = " 3.096 1.0195 -1 -1 0.547 -1 -1 0.135 -1 -1 ";
	string delimiter = " ,";
	
	int index=0;
	Node* root = constructTree( splitString(preorder_serial, delimiter), &index ); 
	
	preOrderTraverse(root);
	cout<<endl;
	postOrderTraverse(root);
	
}
