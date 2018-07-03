#include "TTreeFactory.h"


TTreeFactory::TTreeFactory(int nodeId, int pdg, TFile* f){
	
	//to_string is c++11
	std::string treeid = "resonance" + std::to_string(pdg) + "_node" + std::to_string(nodeId);
	
	tree = new TTree(treeid.c_str(), treeid.c_str());
  	tree->SetDirectory(f);
	

	tree->Branch("RecoEnergy", &RecoEnergy);
  	tree->Branch("FitEnergy", &FitEnergy);
        tree->Branch("RecoMass", &RecoMass);
  	tree->Branch("FitProbability", &FitProbability );
	tree->Branch("Chisq", &Chisq);
        tree->Branch("recoLocalParams.", &recoLocalParams);
	tree->Branch("recoLocalErrors.", &recoLocalErrors);
	tree->Branch("fitLocalParams.",&fitLocalParams);
        tree->Branch("fitLocalErrors.", &fitLocalErrors);
	tree->Branch("pdgs.",&pdgs);

}
void TTreeFactory::addFittedParticle(Particle* fitcontainer){
	fitLocalParams.push_back(fitcontainer->localParams);
	fitLocalErrors.push_back(fitcontainer->localErrors);
	//TODO fit pulls
	
}
void TTreeFactory::addReconstructedParticle(Particle* recocontainer){
	recoLocalParams.push_back(recocontainer->localParams);
	recoLocalErrors.push_back(recocontainer->localErrors);
	printParams(recocontainer->localParams);
	printParams(recocontainer->localErrors);
	
}
void TTreeFactory::addpdg( int pdg){
	pdgs.push_back(pdg);
}

void TTreeFactory::addFitDetails(double fitprob, double chisq){
	FitProbability = fitprob ;
	Chisq = chisq ;
}
void TTreeFactory::addParticleSets(std::vector<Particle*> fitcontainer, std::vector<Particle*> recocontainer){
	//containers better be the same size, so just loop once
	TLorentzVector fitsum, recosum;
	for(unsigned int i=0; i< fitcontainer.size(); i++){
		addFittedParticle(fitcontainer.at(i));
		addReconstructedParticle(recocontainer.at(i));
		addpdg(recocontainer.at(i)->part->getType());

		fitsum += *fitcontainer.at(i)->v;
		recosum += *recocontainer.at(i)->v;
		std::cout<<"FIT ";
     std::cout<<"TLV: (Px,Py,Pz,P,E,M) "<<
		fitcontainer.at(i)->v->Px()<< " " <<
		fitcontainer.at(i)->v->Py()<< " " <<
		fitcontainer.at(i)->v->Pz()<< " " <<
		fitcontainer.at(i)->v->P() << " " <<
		fitcontainer.at(i)->v->E() << " " <<
		fitcontainer.at(i)->v->M() << " " <<std::endl;
	std::cout<<"RECO ";
	std::cout<<"TLV: (Px,Py,Pz,P,E,M) "<<
		recocontainer.at(i)->v->Px()<< " " <<
		recocontainer.at(i)->v->Py()<< " " <<
		recocontainer.at(i)->v->Pz()<< " " <<
		recocontainer.at(i)->v->P() << " " <<
		recocontainer.at(i)->v->E() << " " <<
		recocontainer.at(i)->v->M() << " " <<std::endl;
	}
	RecoEnergy = recosum.E();
	std::cout<<"RECO ";
     std::cout<<"TLV: (Px,Py,Pz,P,E,M) "<<
		recosum.Px()<< " " <<
		recosum.Py()<< " " <<
		recosum.Pz()<< " " <<
		recosum.P() << " " <<
		recosum.E() << " " <<
		recosum.M() << " " <<std::endl;
	FitEnergy = fitsum.E();
	std::cout<<"FIT ";
	std::cout<<"TLV: (Px,Py,Pz,P,E,M) "<<
		fitsum.Px()<< " " <<
		fitsum.Py()<< " " <<
		fitsum.Pz()<< " " <<
		fitsum.P() << " " <<
		fitsum.E() << " " <<
		fitsum.M() << " " <<std::endl;
	RecoMass = recosum.M();

}
void TTreeFactory::TreeFillAndClear(){
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
void TTreeFactory::printParams(std::vector<double> params){
	for(int i=0; i<params.size(); i++){
		std::cout<<params.at(i)<<" ";
	}
	std::cout<<std::endl;
}

