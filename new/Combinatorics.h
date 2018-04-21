#ifndef _COMBIN_
#define _COMBIN_
#include <vector>
#include <algorithm>
#include <iostream>
using namespace std;

class Combinatorics{

	public:
	template<typename type>
	static void generateSubsets(std::vector<type> v, int k, int start, int currLen, std::vector<bool> used, std::vector<std::vector<type> >& combinations);
	
	template<typename type>
	static void generateParticlesCombinations(int vectorsize, int nparticles, std::vector<std::vector<type> >& combinations, vector<type> recoset);

	//sorting (quicksort)
	template<typename type>
	static void swap(type* a, type* b);
	template<typename type>
	static int partition (vector<type>& arr, int low, int high);
	template<typename type>
	static void quickSort(vector<type>& arr, int low, int high);
	
	static bool finalstatepdgmatch(vector<int> fsp_pdgs, vector<int> pdgcombo);
	static void filtercombinations(vector<int> fsp, vector<vector<int> >& combinations, vector<vector<int> >& pdgcombinations);

	//add elements of set b into a
	static vector<int> addSets(vector<int> a, vector<int> b);
	
	//remove elements of b from a
	static vector<int> subtractSets(vector<int> a, vector<int> b);
	
};

#endif
