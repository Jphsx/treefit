#include "Combinatorics.h"

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
void generateParticlesCombinations(int vectorsize, int nparticles, std::vector<std::vector<type> >& combinations, vector<type> recoset){
	int n = vectorsize;
	int k = nparticles;

		std::vector<bool> used(n);	

		for(int i=0; i<n; i++){
			used.at(i)=false;
		}

	generateSubsets(recoset,k,0,0,used,combinations);
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
//does the combo have particles in it that are not in the needed final state? input is pdgs mapped from indices
bool containsfinalstateparticle(vector<int> fsp_pdgs, vector<int> combo){
	
	for(int i=0; i<combo.size(); i++){
		if(std::find(fsp_pdgs.begin(), fsp_pdgs.end(), recoparts.at(combo.at(i))->recopdg) != fsp_pdgs.end()) {
    			/* fsp contains particle */
		} else {
    			/* fsp does not contain particle */
			return false;
		}
	}
	//if we havent returned false yet then combo is good, return true
	return true;
}
//same as containsfinal state particle, but diffrenent names for code clarity
/*bool subsetofparentcombo(vector<int>& parentcurrentcombo,vector<int>& combination){
	for(int i=0; i<combination.size(); i++){
		if(std::find(parentcurrentcombo.begin(), parentcurrentcombo.end(), combination.at(i)) != parentcurrentcombo.end()) {
    			// fsp contains particle 
		} else {
    			// fsp does not contain particle 
			return false;
		}
	}
	return true;
}*/
//combo argument is  pdgs, is the combination of pdgs the same as the required combination?
bool finalstatepdgmatch(vector<int> fsp_pdgs, vector<int> combo){

	//create a copy of both arrays, we don't want the originals sorted because we need to preserve indices
	vector<int> fsp_copy{};
	fsp_copy = fsp_pdgs;	//this is a deep copy	

	//TODO Make an array of pdgs by extracting the pdg values from the combo array
	vector<int> combo_copy{};// combo copy is the pdgs from each of the combo indices
	//the vectors had better be the same size
	for(int i=0; i<combo.size(); i++){ 
		//fsp_copy.push_back(fsp_pdgs.at(i));
		//combo_copy.push_back(combo.at(i));
		combo_copy.push_back(recoparts.at(combo.at(i))->recopdg);
		
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
void filtercombinations(vector<int> fsp, vector<vector<int> >& combinations){//, vector<int> parentcombo ){

	vector<vector<int> > filteredcombos;
	//vector<vector<int> > filteredpdgcombos;
//the three main conditions for selecting a valid combination of particles is

	for(int i=0; i<combinations.size(); i++){

		// 1- the set contains the same pdgs as the decay is supposed to have 
		// this may be redundant (can only use 1-b to save time)
		if(!containsfinalstateparticle(fsp,combinations.at(i))) continue;

		// 1-b the set contains exactly the same composition of pdgs
		if(!finalstatepdgmatch(fsp,combinations.at(i))) continue;

		// 1-c the combinations are a subset of the parent set of indices
		//his is obsolete because combinations are generated from he parent combo
	//	if(!subsetofparentcombo(parentcombo,combinations.at(i))) continue;
		
		// 2- none of the particles in the set are marked as used from previous sets
		//if(containsusedparticle(combinations.at(i))) continue;//used particle found reject combo

		// 3- none of the particles are duplicates
			//require E/p uniqueness
		//if all the tests are passed push the combination onto the new object
		filteredcombos.push_back(combinations.at(i));
		//filteredpdgcombos.push_back(pdgcombinations.at(i));

	}

	combinations.clear();
	//pdgcombinations.clear();//
	//maybe combinations deleted when out of scope,so, try a hard copy
	//TODO change this back i think it was okay
	for(int i=0; i<filteredcombos.size(); i++){
		combinations.push_back(filteredcombos.at(i));
		//pdgcombinations.push_back(filteredpdgcombos.at(i));
	}
}
