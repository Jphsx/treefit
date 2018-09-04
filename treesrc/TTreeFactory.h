#ifndef _FACTORY_
#define _FACTORY_
#include "Particle.h"
#include <vector>
#include <string>

#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include "ThreeVector.h"


class TTreeFactory{


	public:
	//initialize the TTree with constructor
	TTreeFactory(int nodeId, int pdg, TFile* f);

	//vertex variable management
	void initVertexVars();
	bool usingVertex{};
	void clearVertex();
		
	TTree* tree{};
	

	std::vector<std::vector<double> > recoLocalParams{};
  	std::vector<std::vector<double> > recoLocalErrors{};
	std::vector<std::vector<double> > fitLocalParams{};
	std::vector<std::vector<double> > fitLocalErrors{};
	
	//vertex information
	std::vector<double> vertex{};
	std::vector<double> vertexErrors{};

	//this is the indexing array
	std::vector<double> pdgs{};

	//the resonance for the particles contained in this class
	//Params : { Px Py Pz E }
	std::vector<double> recoParentParams{};
	std::vector<double> fitParentParams{};
	//Errors : { dPx dPy dPz dE }
	std::vector<double> recoParentErrors{}; //TODO these
	std::vector<double> fitParentErrors{};
	
	
	//this nodes resonance mass
	double RecoMass{};
	double FitMass{};
	double FitProbability{};
	double Chisq{};
	
	//this nodes total energy
	double RecoEnergy{};
	double FitEnergy{};
	
	void addFittedParticle(Particle* fitcontainer);
	void addReconstructedParticle(Particle* recocontainer);
	void addParticleSets(std::vector<Particle*> fitcontainer, std::vector<Particle*> recocontainer);
	void addpdg( int pdg );
	void addFitDetails(double fitprob, double chisq); 
	void addVertexDetails(VertexFitObject* vfo);
	//TODO MonteCarlo stuff?

	//parent error population
	void addFitParentErrors(float* cov);
	void addRecoParentErrors(float* cov);

	//tree management functions
	void TreeFillAndClear();
	
	//debugging
	void printParams(std::vector<double> params);
	


};
#endif
