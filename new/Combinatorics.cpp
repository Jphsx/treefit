#include "Combinatorics.h"

template <typename type>
void Combinatorics::generateSubsets(std::vector<type> v, int k, int start, int currLen, std::vector<bool> used, std::vector<std::vector<type> >& combinations){
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
void Combinatorics::generateParticlesCombinations(int vectorsize, int nparticles, std::vector<std::vector<type> >& combinations, vector<type> recoset){
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
void Combinatorics::swap(type* a, type* b)
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
int Combinatorics::partition (vector<type>& arr, int low, int high)
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
void Combinatorics::quickSort(vector<type>& arr, int low, int high)
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

//combo argument is  pdgs, is the combination of pdgs the same as the required combination?
bool Combinatorics::finalstatepdgmatch(vector<int> fsp_pdgs, vector<int> pdgcombo){


	//create a matched vector of boolean values
	//get fsp pdgs and combo pdgs
	//sort both arrays using quicksort
	//if both arrays have the same pdg content then the two arrays should be identical after sorting
	quickSort(fsp_pdgs, 0, fsp_pdgs.size()-1);
	quickSort(pdgcombo, 0, pdgcombo.size()-1);

	for(int i=0; i< fsp_pdgs.size(); i++){
		if( fsp_pdgs.at(i) != pdgcombo.at(i) ) return false;
	}
	
	return true;
		
}
//create a combinations vector to work alongside combinations
void Combinatorics::filtercombinations(vector<int> fsp, vector<vector<int> >& combinations, vector<vector<int> >& pdgcombinations){//, vector<int> parentcombo ){

	vector<vector<int> > filteredcombos;
	//vector<vector<int> > filteredpdgcombos;
//the three main conditions for selecting a valid combination of particles is

	for(int i=0; i<combinations.size(); i++){

		
		// 1 the set contains exactly the same composition of pdgs
		if(!finalstatepdgmatch(fsp,pdgcombinations.at(i))) 			continue;


		// 2- none of the particles are duplicates
			//require E/p uniqueness
		//if all the tests are passed push the combination onto the new object
		filteredcombos.push_back(combinations.at(i));
		

	}

	combinations.clear();
	combinations = filteredcombos;

}
/////////////set operations///////////////////////

//the sets need to be sorted before being combined

vector<int> Combinatorics::addSets(vector<int> a, vector<int> b){
	
  	vector<int> setSum(a.size()+b.size());                     
  	std::vector<int>::iterator it;

 	it=std::set_union (a.begin(), a.end(), b.begin(), b.end(), setSum.begin());
                                               
	setSum.resize(it-setSum.begin());
	return setSum;
}

vector<int> Combinatorics::subtractSets(vector<int> a, vector<int> b){
	
	vector<int> setDiff(a.size()+b.size());
	std::vector<int>::iterator it;
	
	it=std::set_difference (a.begin(), a.end(), b.begin(), b.end(), setDiff.begin());
                                              
  	setDiff.resize(it-setDiff.begin()); 
	return setDiff;

}
