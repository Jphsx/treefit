#include "Particle.h"


Particle::Particle( ){}
Particle::Particle(ReconstructedParticle* p, Track* t, double B ){
	//if no track this is just a neutral particle	
	if(t==NULL){
		isTrack = false;
	}
	else{
		isTrack = true;
	}
	//set pdgcode
	recopdg = p->getType();
	part = p;
	track = t;// t can be null
	Bfield = B;
	//if its a track populate tlv from track not reco
	if(isTrack){
		v = getTLorentzVector(t,p->getMass(),B);
	}
	else{
		v = getTLorentzVector(p);
	}

}
void Particle::printTrack(Track* t){
	std::cout<<"Track: (D0,Z0,ome,tanL,phi) "<< 
		t->getD0()<<" "<<
		t->getZ0()<<" "<<
		t->getOmega()<<" "<<
		t->getTanLambda()<<" "<<
		t->getPhi()<<std::endl;
}
void Particle::printReconstructedParticle(ReconstructedParticle* p){
	if(p !=NULL) std::cout<<p<<std::endl;
	if(p == NULL) std::cout<<"p is null"<<std::endl;
	const double* mom = p->getMomentum();	
	std::cout<<"Particle "<< p->getType() <<": "<<
	"(Px,Py,Pz,E,M,q) "<<
	mom[0]<< " "<<mom[1]<< " "<<mom[2]<< " "
	<<p->getEnergy()<<" "<<p->getMass()<<" "
	<<p->getCharge()<<std::endl;
		
	
}
void Particle::printTLorentzVector(TLorentzVector* v){

	std::cout<<"TLV: (Px,Py,Pz,P,E,M) "<<
		v->Px()<< " " <<
		v->Py()<< " " <<
		v->Pz()<< " " <<
		v->P() << " " <<
		v->E() << " " <<
		v->M() << " " <<std::endl;
}
std::vector<double> Particle::getTrackPxPyPz(Track* t, double BField){
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
void Particle::printTrackPxPyPz(Track* t, double B){
	std::vector<double> txtytz = getTrackPxPyPz(t,B);
	std::cout<<"Track: (Px,Py,Pz) "
		<<txtytz.at(0)<<" "
		<<txtytz.at(1)<<" "
		<<txtytz.at(2)<<std::endl;
}
TLorentzVector* Particle::getTLorentzVector(ReconstructedParticle* p){
	TLorentzVector* tlv = new TLorentzVector();
	const double* mom = p->getMomentum();
	tlv->SetXYZM(mom[0],mom[1],mom[2],p->getMass());
	return tlv;
}
TLorentzVector* Particle::getTLorentzVector(Track* t, double Mass, double B){
	TLorentzVector* tlv = new TLorentzVector();
	std::vector<double> txtytz = getTrackPxPyPz(t, B);
	tlv->SetXYZM(txtytz.at(0),txtytz.at(1),txtytz.at(2), Mass);
	return tlv;
}
void Particle::printParticle(Particle* pc){
	std::cout<<std::endl;
	//printReconstructedParticle(pc->part);
	printTLorentzVector(pc->v);
	if(pc->isTrack){
		std::cout<<"this is a track "<<std::endl;
		printTrackPxPyPz(pc->track, pc->Bfield);
		printTrack(pc->track);	
	}
	
}
