

#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> //atof atoi

using namespace std;

//global stuff from xml
//use for mapping pdgs with masses in tree
vector<int> m_pdgs ={ 443, 333, 221, 111 };
vector<double> masses = { 3.096, 1.020, 0.547, 0.135}; 

struct decayTree{
	decayTree* parent;
	vector<decayTree*> children;
	int pdg;
	double mass;
};

void print2Dvector(vector<vector<string> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}

	
}
void print2Dvector(vector<vector<int> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}

	
}
decayTree* buildRootNode(vector<int> firstSubTree ){
	//assign pdg of root
	decayTree* root = new decayTree();
	root->pdg = firstSubTree.at(0);

	//loop over the first children and make decay trees
	for(int i=1; i<firstSubTree.size(); i++){
		decayTree* gen1child = new decayTree();
		gen1child = firstSubTree.at(i);
		root->children.push_back(gen1child);
	}

	return root

}
decayTree* buildChildrenNode(decayTree* root, vector<vector<int>> tree){
	//look at the children if a childs child vecsize = 0 then look for children and make a new node
	
	

}
decayTree* buildTree(vector<string> dstring){
	decayTree* root = new decayTree();
	//delimit array further at [,],][
	vector<vector<string> > topology;
	 std::string delim ("[");
	int i=0;
	while( i< dstring.size()-1){
		std::size_t found = dstring.at(i).find(delim);
  		if (found!=std::string::npos){
			vector<string> subtree;
			i++;
			found = dstring.at(i).find(delim);
		
			while(found==std::string::npos && i < dstring.size()-1){
				subtree.push_back(dstring.at(i));
				i++;
				found = dstring.at(i).find(delim);
				
			}
		topology.push_back(subtree);
		subtree.clear();
		//i++;
		}
	}
	//cout<<"strings"<<endl;
	//print2Dvector(topology);
	//copy onto int vector
	vector<vector<int> > topology_int;
	for(int i=0; i<topology.size(); i++){
		vector<int> subvec;
		for(int j=0; j< topology.at(i).size(); j++){
			
			subvec.push_back(atoi(topology.at(i).at(j).c_str()));
		}
		topology_int.push_back(subvec);
		subvec.clear();
	}
	cout<<"ints"<<endl;
	print2Dvector(topology_int);
	
	return NULL;
}
int main(){

	vector<string> neutrals = {"g1","g2","g3"};
	vector<string> charged = {"pi+1", "pi-2", "pi+3","pi-4","k+1","k-2","k+3","k-4"};
	string decay = "[ 443 331 221 ][ 333,321 , -321 ][ 221 211 -211 111 ][ 111 22 22 ]";

	
	
	vector<string> decaytokens;
	char * dup = strdup(decay.c_str());
	 char * tokens = strtok(dup, " ,");
	
	while(tokens != NULL){
		 decaytokens.push_back(string(tokens));
        // the call is treated as a subsequent calls to strtok:
        // the function continues from where it left in previous invocation
        tokens = strtok(NULL, " ,");
	}
	//for(int i=0; i<decaytokens.size(); i++){
	//	cout<<decaytokens.at(i)<<endl;
	//}
	decayTree* d = buildTree(decaytokens);
	

}
