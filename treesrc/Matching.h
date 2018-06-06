#ifndef _MATCHING_
#define _MATCHING_


//class that matches particles between collections
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
#include <math.h>
#include <vector>
#include <iostream>
#include "Particle.h"
using namespace lcio;
class Matching{
	

	public:
	//simple matching make sure each momentum component is at most 1 sigma ~0.001 GeV apart
	static Track* MatchParticleToTrack(ReconstructedParticle* p, std::vector<Track*> tvec , double BField);
	
	

};
#endif
