#ifndef _COMBIN_
#define _COMBIN_
#include <vector>
template <typename type>
class Combinatorics{

	public:
	static void generateSubsets(std::vector<type> v, int k, int start, int currLen, std::vector<bool> used, std::vector<std::vector<type> >& combinations);
	
	static void generateParticlesCombinations(int vectorsize, int nparticles, std::vector<std::vector<type> >& combinations, vector<type> recoset);

	//sorting (quicksort)
	static void swap(type* a, type* b);
	static int partition (vector<type>& arr, int low, int high);
	static void quickSort(vector<type>& arr, int low, int high);
	static bool containsfinalstateparticle(vector<int> fsp_pdgs, vector<int> combo);
	static bool finalstatepdgmatch(vector<int> fsp_pdgs, vector<int> combo);
	static void filtercombinations(vector<int> fsp, vector<vector<int> >& combinations);
	
	
};

#endif
