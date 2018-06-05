#ifndef
#define _MATCHING_


//class that matches particles between collections

#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
#include "vector"
#include "iostream"

class Matching{
	

	//from the 5 track parameters return a vector of 
	//momentum components px,py,pz
	static std::vector<double> getTrackPxPyPz(Track* t, double BField);
	
	//simple matching make sure each momentum component is at most 1 sigma ~0.001 GeV apart
	static Track* MatchParticleToTrack(ReconstructedParticle* p, std::vector<Track*> tvec );
	
	

};
#endif
