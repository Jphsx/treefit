//this is the templated prototype for combination generation
//this produces combinations of objects rather than indices based on a single array
#include <vector>
#include <iostream>
#include <string.h>

using namespace std;

template <typename type>
void print2dvec(vector<vector<type> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}
}
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
int main(){
	vector<vector<int> > combos{};
	vector<int> recoparts = {1, 10, 13, 15, 2, 7};
	//(n, k, combos) n choose k
	generateIndicesCombinations(6, 3, combos, recoparts);
	print2dvec(combos);
	cout<<"number combos "<< combos.size()<<endl;
	cout<<endl;
	
/*
	vector<vector<int> > combos2{};
	vector<char> recoobjs = {'a','b','c','d'};
	generateIndicesCombinations(4,2,combos2, recoobjs);
	print2dvec(combos2);
	cout<<"number combos "<< combos2.size()<<endl;
*/
	
}
