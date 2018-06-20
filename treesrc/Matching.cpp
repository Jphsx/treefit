#include "Matching.h"
using namespace lcio;


Track* Matching::MatchParticleToTrack(ReconstructedParticle* p, std::vector<Track*> tvec, double BField){

	if(p->getCharge() == 0) return NULL;
	const double* pxpypz = p->getMomentum();
	std::vector<double> txtytz;
	double dpx,dpy,dpz;
	for(unsigned int i=0; i< tvec.size(); i++){
		txtytz = Particle::getTrackPxPyPz(tvec.at(i), BField);
		dpx = abs(pxpypz[0] - txtytz.at(0));
		dpy = abs(pxpypz[1] - txtytz.at(1));
		dpz = abs(pxpypz[2] - txtytz.at(2));
		if( (dpx<0.015) && (dpy<0.015) && (dpz<0.015)){
			//track is matched
			return tvec.at(i);
		}
		txtytz.clear();
	}
	//if no match is found
	std::cout<<"Particle "<<p->getType()<< " No Matching Track"<<std::endl;
	return NULL;	
};
