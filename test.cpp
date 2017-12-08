

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

template <class type>
void print2Dvector(vector<vector<type> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}

}
template <class type>
void printvector(vector<type> v){
	for(int i=0; i<v.size(); i++){
		cout<<v.at(i)<<" ";
	}
	cout<<endl;
}
/*decayTree* buildRootNode(vector<int> firstSubTree ){
	//assign pdg of root
	decayTree* root = new decayTree();
	root->pdg = firstSubTree.at(0);

	//loop over the first children and make decay trees
	for(int i=1; i<firstSubTree.size(); i++){
		decayTree* gen1child = new decayTree();
		gen1child = firstSubTree.at(i);
		root->children.push_back(gen1child);
	}

	return root;

}*/
vector<decayTree*> createNodes(vector<vector<int> > subtrees, vector<double> masses, decayTree* node, int subindex){
	//look at the children if a childs child vecsize = 0 then look for children and make a new node
	//vector<decayTree*> nodes;
	//build each node based on the subtree arrays
	//for(int i=0; i<subtrees.size(); i++){
	//	decayTree* node = new decayTree();
		
	//}
	
	//the children should start as null
	//assign children to node
	if(node->children.size()==0){
		//make the children
		for(int i=0; i<subtrees.at(subindex)
	}
	

}
decayTree* buildTree(vector<string> dstring){
	decayTree* root = new decayTree();
	//delimit array further at [,],][
	vector<vector<string> > topology;
	cout<<"input dstring ";
	printvector(dstring);


	//the first character should always be [ so start i=1
	vector<string> subtree;	
	for(int i=1; i<dstring.size(); i++){
		if (dstring.at(i).find("[") != std::string::npos || dstring.at(i).find("]") != std::string::npos ){
			topology.push_back(subtree);
			subtree.clear();	
			continue;
			
		}
		subtree.push_back(dstring.at(i));
	}

	//copy onto int vector
	vector<vector<int> > topology_int;
	vector<double> topology_mass;
	for(int i=0; i<topology.size(); i++){
		vector<int> subvec;
		//size-1 because last element is reserved for mass
		for(int j=0; j< topology.at(i).size()-1; j++){
			
			subvec.push_back(atoi(topology.at(i).at(j).c_str()));
		}	
		//get masses from vector
		topology_mass.push_back(atof( topology.at(i).at( topology.at(i).size()-1 ).c_str() ) );

		topology_int.push_back(subvec);
		subvec.clear();
	}


	cout<<"pdg ints"<<endl;
	print2Dvector(topology_int);
	cout<<"masses"<<endl;
	printvector(topology_mass);

	//make all the nodes of the tree
	createNodes(topology_int,topology_mass);
	
	return NULL;
}
int main(){

	vector<string> neutrals = {"g1","g2","g3"};
	vector<string> charged = {"pi+1", "pi-2", "pi+3","pi-4","k+1","k-2","k+3","k-4"};
	string decay = "[ 443 331 221 3.096 ][ 331,321 , -321 1.0195 ][ 221 211 -211 111 0.547 ][ 111 22 22 0.135 ]";

	
	
	vector<string> decaytokens;
	char * dup = strdup(decay.c_str());
	 char * tokens = strtok(dup, " ,");
	
	while(tokens != NULL){
		 decaytokens.push_back(string(tokens));
        // the call is treated as a subsequent calls to strtok:
        // the function continues from where it left in previous invocation
        tokens = strtok(NULL, " ,");
	}
	
	decayTree* d = buildTree(decaytokens);
	

}
