#include "RootFileFactory.h"


RootFileFactory::RootFileFactory(int nodeId, int pdg){
	
	//to_string is c++11
	std::string filename = std::to_string(pdg) + "resonance_node" + std::tostring(nodeId);
	
	rootFile = new TFile(filename.c_str(),"RECREATE");
	tree = new TTree("tree", "tree");
  	tree->SetDirectory(rootFile);
	

	tree->Branch("RecoEnergy", &RecoEnergy);
  	tree->Branch("FitEnergy", &FitEnergy);
        tree->Branch("RecoMass", &RecoMass);
  	tree->Branch("FitProbability", &FitProbability );
	tree->Branch("Chisq", &Chisq);
        tree->Branch("recoLocalParams.", &recoLocalParams);
	tree->Branch("recoLocalErrors.", &recoLocalErrors);
	tree->Branch("fitLocalParams.",&fitLocalParams);
        tree->Branch("fitLocalErrors.", &fitLocalErrors);

}
void RootFileFactory::addFittedParticle(Particle* fitcontainer){
	fitLocalParams.push_back(fitcontainer->localParams);
	fitLocalErrors.push_back(fitcontainer->localErrors);
	//TODO fit pulls
	
}
void RootFileFactory::addReconstructedParticle(Particle* recocontainer){
	recoLocalParams.push_back(recocontainer->localParams);
	recoLocalErrors.push_back(recocontainer->localErrors);
	
}
void RootFileFactory::addpdg( int pdg){
	pdgs.pushback(pdg);
}

void RootFileFactory::addFitDetails(double fitprob, double chisq){
	FitProbability = fitprob ;
	Chisq = chisq ;
}
void RootFileFactory::addParticleSets(std::vector<Particle*> fitcontainer, std::vector<Particle*> recocontainer){
	//containers better be the same size, so just loop once
	TLorentzVector fitsum, recosum;
	for(unsigned int i=0; i< fitcontainer.size(); i++){
		addFittedParticle(fitcontainer.at(i));
		addReconstructedParticle(recocontainer.at(i));
		addpdg(recocontainer.at(i)->part->getType());
		RecoMass = recocontainer.at(i)->part->getMass();
		fitsum += fitcontainer->v;
		recosum += recocontainer->v;
	}
	RecoEnergy = recosum.E();
	FitEnergy = fitsum.E();

}
void RootFileFactory::TreeFillAndClear(){
	//fill the tree
	tree->Fill();
	//clear out all the vectors for the next event
	for(unsigned int i=0; i<recoLocalParams.size(); i++){
		recoLocalParams.at(i).clear();
		recoLocalErrors.at(i).clear();
		fitLocalParams.at(i).clear();
		fitLocalErrors.at(i).clear();
	}
	recoLocalParams.clear();
	recoLocalErrors.clear();
	fitLocalParams.clear();
	fitLocalErrors.clear();
	pdgs.clear();
	
}

