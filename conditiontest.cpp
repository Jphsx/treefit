#include <vector>
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include <algorithm>

using namespace std;
//global vector of reconstructed particles
vector<int> recopdgs ={ 211, -211, 13, -13, 22, 22, 22, 22 };
vector<int> usedparts = {0, 0, 0, 0, 0, 0, 0, 0};






template <class type>
void print2dvec(vector<vector<type> > v){
	for(int i=0; i<v.size(); i++){
		for(int j=0; j<v.at(i).size(); j++){
			cout<<v.at(i).at(j)<<" ";
		}
		cout<<endl;
	}
}
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
//combo argument is unique particle index not pdgs, is the combination of pdgs the same as the required combination?
bool finalstatepdgmatch(vector<int> fsp, vector<int> combo){

	//create a matched vector of boolean values
	vector<bool> used;
	for(int i=0;fsp.size();i++  ){
		used.push_back(false);
	}
	
		
	return false;
		
}
vector<vector<int> > filtercombinations(vector<int> fsp, vector<vector<int> >& combinations, vector<vector<int> > pdgcombinations ){

	vector<vector<int> > filteredcombos;
//the three main conditions for selecting a valid combination of particles is

	for(int i=0; i<combinations.size(); i++){

		// 1- the set contains the same pdgs as the decay is supposed to have 
		if(!containsfinalstateparticle(fsp,pdgcombinations.at(i))) continue;
		
		// 2- none of the particles in the set are marked as used from previous sets
		if(containsusedparticle(combinations.at(i))) continue;//used particle found reject combo

		// 3- none of the particles are duplicates
			//require E/p uniqueness
		//if all the tests are passed push the combination onto the new object
		filteredcombos.push_back(combinations.at(i));

	}

	return filteredcombos;
}
//map indices to pdgs
vector<vector<int> > mapindextopdg( vector<vector<int> >& combinations){
	vector<vector<int> > pdgcombos;

		for(int i=0; i<combinations.size(); i++){
			vector<int> subarray;
			for(int j=0; j<combinations.at(i).size(); j++){
				subarray.push_back( recopdgs.at( combinations.at(i).at(j) ) );
			}
			pdgcombos.push_back(subarray);
			subarray.clear();
		}
		return pdgcombos;
	
}


int main(){
	//test vector vector
	vector<vector<int> > combos;
	//assess 4particle final state from the reconstructed set
	//8 choose 4
	generateIndicesCombinations(8,4,combos);
	print2dvec(combos);
	cout<<"n combinations ="<< combos.size()<<endl;

	vector<vector<int> > pdgcombos;
	pdgcombos = mapindextopdg(combos);
	print2dvec(pdgcombos);
	
	vector<int> finalstatepdgs ={ 211, -211, 22, 22 };

	vector<vector<int> > filteredcombos;
	filteredcombos= filtercombinations(finalstatepdgs, combos, pdgcombos);

	cout<<"filtered n combinations = "<< filteredcombos.size()<<endl;
	print2dvec(filteredcombos);
	
	vector<vector<int> > filteredpdgcombos;
	filteredpdgcombos = mapindextopdg(filteredcombos);
	print2dvec(filteredpdgcombos);
	

}
