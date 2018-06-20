#ifndef _PARTICLE_
#define _PARTICLE_

#include "lcio.h"
#include <iostream>
#include <vector>
#include "TLorentzVector.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
//#include "TrackParticleFitObject.h"
//#include "JetFitObject.h"

typedef lcio::Track Track ;
typedef lcio::ReconstructedParticle ReconstructedParticle ;
////////////////////
//a general particle container that aggregates all possible information or a reconstructed particle
class Particle{
	public:
	Particle();
	//build from recos/tracks
	Particle(ReconstructedParticle* p, Track* t, double B);

	//build from fitobjects
	Particle(JetFitObject* jfo, TrackParticleFitObject tpfo, int pdg, float mass);

	//need to write a good destructor here?
	//TODO ~Particle()
	
	int recopdg;	
	double Bfield;
	//distinguish recoparts from tracks
	bool isTrack;

	//ReconstructedParticle/Track*
	Track* track;
	ReconstructedParticle* part;

	//TLorentzVector RECO
	TLorentzVector* v;

	//MCParticle
	//TLorentzVector MC

	//LocalParamterization 
	std::vector<double> localParams{};
	//LocalParameterization Errors sigma _not_ sigma^2
	std::vector<double> localErrors{};


	//printing utilities
	//print 5 track parameters
	static void printTrack(Track* t);
	//print track momentum only
	static void printTrackPxPyPz(Track* t, double B);
	//print pdg/charge/px py pz E M
	static void printReconstructedParticle(ReconstructedParticle* p);
	//print px py pz E M
	static void printTLorentzVector(TLorentzVector* v);
	//print whatever the local parameters are
	static void printLocalParameters(std::vector<double> params);
	//print whatever the local errors are (sigmas)
	static void printLocalErrors(std::vector<double> errors);

	//from the 5 track parameters return a vector of 
	//momentum components px,py,pz
	static std::vector<double> getTrackPxPyPz(Track* t, double BField);
	
	//make TLorentzVector from ReconstructedParticle/Track
	static TLorentzVector* getTLorentzVector(ReconstructedParticle* p);
	static TLorentzVector* getTLorentzVector(Track* t, double Mass, double B);

	static void printParticle(Particle* pc);
	
	

};


#endif
