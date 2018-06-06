#ifndef _PARTICLE_
#define _PARTICLE_

#include "lcio.h"
#include <iostream>
#include <vector>
#include "TLorentzVector.h"
#include "
////////////////////
//a general particle container that aggregates all possible information or a reconstructed particle
class Particle{
	public:
	int recopdg;

	//distinguish recoparts from tracks
	bool isTrack;

	//some particle object should be in here as well
	//ReconstructedParticle/Track*
	//TLorentzVector RECO
	//MCParticle
	//TLorentzVector MC
	//LocalParamterization 
	//LocalParameterization Errors


	//printing utilities
	//print 5 track parameters
	static void printTrack(Track* t);
	static void printTrackPxPyPz(Track* t, double B);
	//print pdg/charge/px py pz E M
	static void printReconstructedParticle(ReconstructedParticle* p);
	//print px py pz E M
	static void printTLorentzVector(TLorentzVector v);

	//from the 5 track parameters return a vector of 
	//momentum components px,py,pz
	static std::vector<double> getTrackPxPyPz(Track* t, double BField);
	
	

};


#endif
