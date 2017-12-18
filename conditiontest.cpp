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
template <class type>
void printvector(vector<type> v){
	for(int i=0; i<v.size(); i++){
			cout<<v.at(i)<<" ";
	}
	//cout<<endl;
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
//combo argument is  pdgs, is the combination of pdgs the same as the required combination?
bool finalstatepdgmatch(vector<int>& fsp_pdgs, vector<int>& combo_pdgs){

	//create a matched vector of boolean values
	//get fsp pdgs and combo pdgs
	//sort both arrays using quicksort
	//if both arrays have the same pdg content then the two arrays should be identical after sorting
	quickSort(fsp_pdgs, 0, fsp_pdgs.size()-1);
	quickSort(combo_pdgs, 0, combo_pdgs.size()-1);

	for(int i=0; i< fsp_pdgs.size(); i++){
		if( fsp_pdgs.at(i) != combo_pdgs.at(i) ) return false;
	}
	
	return true;
		
}
vector<vector<int> > filtercombinations(vector<int> fsp, vector<vector<int> >& combinations, vector<vector<int> > pdgcombinations ){

	vector<vector<int> > filteredcombos;
//the three main conditions for selecting a valid combination of particles is

	for(int i=0; i<combinations.size(); i++){

		// 1- the set contains the same pdgs as the decay is supposed to have 
		// this may be redundant (can only use 1-b to save time)
		if(!containsfinalstateparticle(fsp,pdgcombinations.at(i))) continue;

		// 1-b the set contains exactly the same composition of pdgs
		if(!finalstatepdgmatch(fsp,pdgcombinations.at(i))) continue;
		
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


	//sorting test: sort this random array
	vector<int> numbers = {10 , -90, 400, 2, -7, 13, 999, -400, 10, 10 ,10 ,2, 14, 399, -399};
	
	quickSort(numbers, 0, numbers.size()-1);
    	cout<<"Sorted array: ";
   	printvector(numbers);
   	cout<<endl;

	//now test marking a particle used
	//last photon is now used
	cout<<"same combination generating but with last photon marked used"<<endl;
	usedparts.at(usedparts.size()-1)=1;
	filteredcombos = filtercombinations(finalstatepdgs, combos, pdgcombos);
	filteredpdgcombos = mapindextopdg(filteredcombos);
	print2dvec(filteredcombos);
	print2dvec(filteredpdgcombos);
	cout<<endl;
	

}
