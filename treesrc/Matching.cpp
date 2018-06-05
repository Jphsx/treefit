#include "Matching.h"
using namespace lcio;


std::vector<double> getTrackPxPyPz(Track* t, double BField){
	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = BField*c*mm2m*eV2GeV;
 
 	
	double cosLambda = 1 / sqrt(1 + t->getTanLambda()*t->getTanLambda() );
	double P = (eB/fabs(t->getOmega()))/cosLambda;
	double sinLambda = t->getTanLambda()*cosLambda;
	double cosPhi = cos(t->getPhi());
	double sinPhi = sin(t->getPhi());
	double px = P*cosLambda*cosPhi;
	double py = P*cosLambda*sinPhi;
	double pz = P*sinLambda;
	std::vector<double> txtytz;
	txtytz.push_back(px);
	txtytz.push_back(py);
	txtytz.push_back(pz);
	return txtytz;
}
Track* MatchParticleToTrack(ReconstructedParticle* p, std::vector<Track*> tvec, double BField){

	const double* pxpypz = p->getMomentum();
	std::vector<double> txtytz;
	double dpx,dpy,dpz;
	for(unsigned int i=0; i< tvec.size(); i++){
		txtytz = getTrackPxPyPz(tvec.at(i), BField);
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
	std::cout<<"Particle "<<p->getType()<< "No Matching Track"<<std::endl;
	return NULL;	
}
