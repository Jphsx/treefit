#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h> //atof atoi

using namespace std;
struct Node{
	//Node* parent;
	vector<Node*> children;
	int pdg;
	double mass;

	bool isLeaf;

	//unique ID from key
	int nodeId;
};
struct Node* newNode (int id, int pdg, double mass, bool leaf)
{
     Node* temp = new Node();//malloc(sizeof(Node));
    temp->nodeId=id;
    temp->pdg=pdg;
    temp->mass=mass;
    temp->isLeaf=leaf;
    return temp;
}

template <class type>
void printvector(vector<type> v){
	for(int i=0; i<v.size(); i++){
			cout<<v.at(i)<<" ";
	}
cout<<endl;
}
vector<int> vectorizeString_int(string str, string delimiter){
	vector<int> vectokens;
	char* dup = strdup(str.c_str());
	char* dup_delim = strdup(delimiter.c_str());
	char* tokens = strtok(dup,dup_delim);
	while( tokens != NULL){
		vectokens.push_back(atoi(tokens));
		tokens = strtok(NULL,dup_delim);
	}
	return vectokens;
}
vector<double> vectorizeString_double(string str, string delimiter){
	vector<double> vectokens;
	char* dup = strdup(str.c_str());
	char* dup_delim = strdup(delimiter.c_str());
	char* tokens = strtok(dup,dup_delim);
	while( tokens != NULL){
		vectokens.push_back(atof(tokens));
		tokens = strtok(NULL,dup_delim);
	}
	return vectokens;
}
// A recursive function to construct Full from pre[] and post[]. 
// preIndex is used to keep track of index in pre[].
// l is low index and h is high index for the current subarray in post[]
struct Node* constructTreeUtil (vector<int> pre, vector<int> post, int* preIndex, int l, int h){
    // Base case
    if (*preIndex >= pre.size() || l > h)
        return NULL;
 
    // The first node in preorder traversal is root. So take the node at
    // preIndex from preorder and make it root, and increment preIndex
    struct node* root = newNode ( pre[*preIndex] );
    ++*preIndex;
 
    // If the current subarry has only one element, no need to recur
    if (l == h)
        return root;
 
    // Search the next element of pre[] in post[]
    int i;
    for (i = l; i <= h; ++i)
        if (pre[*preIndex] == post[i])
            break;
 
    // Use the index of element found in postorder to divide postorder array in
    // two parts. Left subtree and right subtree
    if (i <= h)
    {
        root->left = constructTreeUtil (pre, post, preIndex, l, i, size);
        root->right = constructTreeUtil (pre, post, preIndex, i + 1, h, size);
    }
 
    return root;
}
struct Node *constructTree (vector<int> pre,vector<int> post[])
{
    int preIndex = 0;
    return constructTreeUtil (pre, post, &preIndex, 0, size - 1);
}
 
int main(){ 
	
	//traversals in string form
	//
	string preorder = " 443 331 321 -321 221 211 -211 111 22 22 ";
	//since pdgs are not unique we must define unique nodes to map to pdg values and masses
	string preorder_key = "0 1 2 3 4 5 6 7 8 9";

	string postorder = " 321 -321 331 211 -211 22 22 111 221 443 ";
	string postorder_key ="2 3 1 5 6 8 9 7 4 0";

	string preorder_mass = " 3.096 1.0195 -1 -1 0.547 -1 -1 0.135 -1 -1 ";
	string preorder_leaves = " 0 0 1 1 0 1 1 0 1 1 ";
	
	string delimiter = " ,";
	printvector(vectorizeString_int(preorder,delimiter));
	printvector(vectorizeString_int(postorder,delimiter));
	printvector(vectorizeString_double(preorder_mass,delimiter));
	printvector(vectorizeString_int(preorder_leaves,delimiter));

	//create unique key value pairs for each particle using indices
	vector<int> pre = vectorizeString_int(preorder, delimiter);
	vector<int> post = vectorizeString_int(postorder, delimiter);
	
	
	vector<int> pre_key = vectorizeString_int(preorder_key, delimiter);
	vector<int> post_key = vectorizeString_int(postorder_key, delimiter);

	vector<double> pre_mass = vectorizeString_double(preorder_mass, delimiter);
	vector<int> pre_leaf = vectorizeString_int(preorder_leaves, delimiter);
	






	

}
