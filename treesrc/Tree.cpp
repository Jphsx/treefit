#include "Tree.h"
Tree::Tree(){
	Root = NULL;
	lastNonLeafNodeId=-1;
}
Tree::~Tree(){
	
	deletePostTraverse(Root);
	Root=NULL;
}
void Tree::deletePostTraverse(Node* root){
	
	for(int i=0; i<root->children.size(); i++){
		deletePostTraverse(root->children.at(i));
	}	
	delete root;

}
Node* Tree::newNode (int id){
    Node* temp = new Node();
    temp->nodeId=id;
    return temp;
}

template <typename type>
void Tree::printvector(vector<type> v){
	for(int i=0; i<v.size(); i++){
			cout<<v.at(i)<<" ";
	}
	//cout<<endl;
}
template <typename type>
void Tree::print2dvec(vector<vector<type> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}
}
vector<double> Tree::castVector_double(vector<string> v){
	vector<double> v_double;
	for(int i=0; i<v.size(); i++){
		v_double.push_back(atof( v.at(i).c_str() ));
	}
	return v_double;
}
vector<int> Tree::castVector_int(vector<string> v){
	vector<int> v_int;	
	for(int i=0; i<v.size(); i++){
		v_int.push_back(atoi( v.at(i).c_str() ));
	}
	return v_int;
}
vector<string> Tree::splitString(string str, string delimiter){
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
Node* Tree::constructTree(vector<string> serial, int* serialIndex){
		
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
void Tree::preOrderTraverse(Node* tree){
	cout<<tree->nodeId<<" ";
	for(int i=0; i<tree->children.size(); i++){
		preOrderTraverse(tree->children.at(i));
	}
	return;
}
//goes right to left
void Tree::postOrderTraverse(Node* tree){
	for(int i=0; i<tree->children.size(); i++){
		postOrderTraverse(tree->children.at(i));
	}
	cout<<tree->nodeId<<" ";
	return;
}
//a preorder traversal that sets the mass constraint mass for each node based on the input string
void Tree::setTreeMasses(Node* root, vector<float> masses, int* massptr ){
	root->mass = masses.at(*massptr);
	++*massptr;
	for(int i=0; i< root->children.size(); i++){
		setTreeMasses(root->children.at(i), masses, massptr);
	}
	return;
}
//a preorder traversal that sets the vertex constraint for each node based on the input vector
void Tree::setTreeVertexConstraints(Node* root, vector<int> VCs, int* vcptr){
	root->VC = VCs.at(*vcptr);
	++*vcptr;
	for(int i=0; i< root->children.size(); i++){
		setTreeVertexConstraints(root->children.at(i), VCs, vcptr);
	}
	return;
}
//a preorder traversal that sets the pdg codes for each node based on the input string
void Tree::setTreePdgCodes(Node* root, vector<int> pdgs, int* pdgptr ){
	root->pdg = pdgs.at(*pdgptr);
	++*pdgptr;
	for(int i=0; i< root->children.size(); i++){
		setTreePdgCodes(root->children.at(i), pdgs, pdgptr);
	}
	return;
}
void Tree::markTreeLeaves(Node* root){
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
void Tree::findLeaves(Node* root, Node* originalParent){
	if(root->isLeaf){
		 originalParent->nLeaves++;
		 originalParent->leafpdgs.push_back( root->pdg );
	}
	for(int i=0; i<root->children.size(); i++){
		findLeaves(root->children.at(i), originalParent);
	}
	return;
		
}
void Tree::populateNLeaves(Node* root){
	//int counter= 0;
	findLeaves(root, root);
	//cout<<"found the leaves counter = "<<counter<<endl;
	//root->nLeaves=counter;
	for(int i=0; i<root->children.size(); i++){
		populateNLeaves(root->children.at(i));
	}
	return;
}


//make sure the root is predefined by passing in NULL
void Tree::setParents(Node* root, Node* parentptr){
	root->parent = parentptr;
	for(int i=0; i<root->children.size(); i++){
		setParents(root->children.at(i), root);
	}
}
Node* Tree::getParentNextNonLeafChild(Node* parent, int callingNodeID){
	if(parent != NULL){
		for(int i=0; i<parent->children.size(); i++){
			if( !parent->children.at(i)->isLeaf && parent->children.at(i)->nodeId > callingNodeID){
				return parent->children.at(i);
			}	
		}
	}
	return NULL;
}
Node* Tree::locateAncestorNearestNonLeafChild(Node* root){
	if(root->parent == NULL) return NULL;
	if(getParentNextNonLeafChild(root->parent, root->nodeId) == NULL){
		return locateAncestorNearestNonLeafChild(root->parent);
	}
	else{
		return getParentNextNonLeafChild(root->parent, root->nodeId);
	}
	
}
Node* Tree::getFirstNonLeafChild(Node* root){

	for(int i=0; i<root->children.size(); i++){
		if(!root->children.at(i)->isLeaf){
			return root->children.at(i);
		}
	}
	//this node is all leaves
	return NULL;

}
void Tree::printfit(Node* root){
	if(root->isLeaf) return;

	cout<<"Node "<<root->nodeId;
	cout<<" FitRecoIDs: ";
	printvector(root->currentcombination);
	cout<<" FitRecoPDGs: ";
	printvector(root->currentcombination_pdgs);
	cout<<endl;
	for(int i=0; i< root->children.size(); i++){
		printfit(root->children.at(i));
	}
	
}
//this method should be removed later
void Tree::getNodePdg(Node* root, int id, int* pdg){
	
	if(root->nodeId == id){
		*pdg = root->pdg;
		return;
	}
	for(int i=0; i<root->children.size(); i++){
		getNodePdg(root->children.at(i), id, pdg);
	}
	return ;
}
Node* Tree::getNode(Node* root, int id){
	//bubble up the node ptr through returns
	Node* node = NULL;
	if(root->nodeId == id){
		node = root;
		return node;
	}
	else{
		//this isnt the node so keep searching
		for(int i=0; i<root->children.size();i++){
			node = getNode(root->children.at(i), id);
			//if we end up finding it in a child
			// break out, dont iterate through children
			if(node != NULL) return node;
		}
	}
	//worst case here, we cant find it
	//null will get returned
	return node;
}
void Tree::getLastNonLeafNodeId(Node* root, int* id){
	if(!root->isLeaf){
		//cout<<"this node "<<root->nodeId<<endl;
		*id = root->nodeId;
	}
	for(int i=0; i<root->children.size(); i++){
		getLastNonLeafNodeId(root->children.at(i), id);
	}
}
//preorder print all tree information
void Tree::printTree(Node* root){
	cout<<endl;
	cout<< "NodeId: "<< root->nodeId <<" Pdg: "<< root->pdg <<" Mass: "<< root->mass << " VertexConstraint = "<< root->VC << " isLeaf= "<<root->isLeaf <<" Children { ";
	for(int i=0; i<root->children.size(); i++){
		cout<< root->children.at(i)->nodeId << " ";
	}
	cout<< "}"<<" nLeaves= "<< root->nLeaves << " Leaf Pdgs: { " ;
	printvector(root->leafpdgs);
	cout<< "}";
	
	cout<< " parentID ";
	if(root->parent != NULL) cout<< root->parent->nodeId;
	
	cout<<endl;
	cout<<" next non leaf sibling: ";
	Node* temp = getParentNextNonLeafChild(root->parent, root->nodeId);
	if(temp != NULL){
	cout<< temp->nodeId;
	}
	cout<<endl;

	cout<<" Nearest ancestor non-leaf child: ";
	temp = locateAncestorNearestNonLeafChild(root);
	if(temp != NULL){
	cout<< temp->nodeId;
	}
	cout<<endl;

	cout<<" First non-leaf child: ";
	temp = getFirstNonLeafChild(root);
	
	if(temp != NULL){
	cout<< temp->nodeId;
	}
	
	cout<<endl;
	for(int i=0; i<root->children.size(); i++){
		printTree(root->children.at(i));
	}
}
void Tree::treeInit(vector<int> pdg, string serial, vector<float> mass, vector<int> vertex, string delimiter){
	int index = 0;
	Node* root = new Node();
	root = constructTree( splitString(serial,delimiter), &index);
	//theTree = root;
	cout<<"PreOrder Traversal: ";
	preOrderTraverse(root);
	cout<<endl;
	cout<<"PostOrder Traversal: ";
	postOrderTraverse(root);
	cout<<endl;

	index=0;
	setTreeMasses(root, mass, &index);
	index=0;
	setTreePdgCodes(root,pdg, &index);
	index=0;
	setTreeVertexConstraints(root,vertex,&index);
	
	markTreeLeaves(root);
	
	populateNLeaves(root);
	
	setParents(root,NULL);
	
	int lastNonLeaf=0;
	getLastNonLeafNodeId(root,&lastNonLeaf);
	lastNonLeafNodeId = lastNonLeaf;
	
	cout<<"Last non-leaf NodeID: "<<lastNonLeafNodeId<<endl;
	
	Root = root;
}

