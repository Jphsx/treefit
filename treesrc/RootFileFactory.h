#ifndef _FACTORY_
#define _FACTORY_
#include "Particle.h"
#include <vector>
#include <string>

class RootFileFactory{


	public:
	//initialize the TTree with constructor
	RootFileFactory(int nodeId, int pdg);
	
	TFile* rootFile{};	
	TTree* tree{};

	std::vector<std::vector<double> > recoLocalParams{};
  	std::vector<std::vector<double> > recoLocalErrors{};
	std::vector<std::vector<double> > fitLocalParams{};
	std::vector<std::vector<double> > fitLocalErrors{};

	//this is the indexing array
	std::vector<double> pdgs{};
	
	//this nodes resonance mass
	double RecoMass{};
	double FitProbability{};
	double Chisq{};
	
	//this nodes total energy
	double RecoEnergy{};
	double FitEnergy{};
	
	void addFittedParticles(Particle* fitcontainer);
	void addReconstructedParticles(Particle* recocontainer);
	void addParticleSets(std::vector<Particle*> fitcontainer, std::vector<Particle*> recocontainer);
	void addpdg( int pdg );
	void addFitDetails(double fitprob, double chisq); 
	//TODO MonteCarlo stuff

	//tree management functions
	void TreeFillAndClear();
	


};
#endif
