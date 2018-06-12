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
	//also populate the local parameterization are error matrices
	if(isTrack){
		v = getTLorentzVector(t,p->getMass(),B);
		localParams.push_back(t->getD0());
		localParams.push_back(t->getPhi());
		localParams.push_back(t->getOmega());
		localParams.push_back(t->getZ0());
		localParams.push_back(t->getTanLambda());
		localErrors.push_back(std::sqrt(t->getCovMatrix()[0]));//d0 
            localErrors.push_back(std::sqrt(t->getCovMatrix()[2]));//phi
            localErrors.push_back(std::sqrt(t->getCovMatrix()[5]));//ome
             localErrors.push_back(std::sqrt(t->getCovMatrix()[9]));//z0
            localErrors.push_back(std::sqrt(t->getCovMatrix()[14]));//tanL
	}
	else{
		v = getTLorentzVector(p);
		localParams.push_back(p->getEnergy());//E
		localParams.push_back(v.getTheta());//theta
		localParams.push_back(v.getPhi());//phi
		//reconstructed particle covmatrix must be modified 
		//to use the E,theta,phi error model
		localErrors.push_back(std::sqrt(p->getCovMatrix()[0]));//dE
		localErrors.push_back(std::sqrt(p->getCovMatrix()[2]));//dtheta
		localErrors.push_back(std::sqrt(p->getCovMatrix()[5]));//dphi
	}
	

}
void Particle::printTrack(Track* t){
	std::cout<<"Track: (D0,phi,ome,tanL,phi) "<< 
		t->getD0()<<" "<<
		t->getPhi()<<" "<<
		t->getOmega()<<" "<<
		t->getZ0()<<" "<<
		t->getTanLambda()<<std::endl;
}
void Particle::printReconstructedParticle(ReconstructedParticle* p){
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
void Particle::printLocalParameters(std::vector<double> params){
	std::cout<<"Local Params: "<<
	for(int i=0; i<localParams.size(); i++){
		std::cout<< params.at(i) <<" ";
	}
	std::cout<<std::endl;
}
void Particle::printLocalErrors(std::vector<double> errors){
	std::cout<<"Local Errors: "<<
	for(int i=0; i<localErrors.size(); i++){
		std::cout<< errors.at(i) <<" ";
	}
	std::cout<<std::endl;
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

	printReconstructedParticle(pc->part);
	printTLorentzVector(pc->v);

	if(pc->isTrack){
		printTrackPxPyPz(pc->track, pc->Bfield);
		printTrack(pc->track);	
	}

	printLocalParameters(pc->localParams);
	printLocalErrors(pc->localErrors);
	
	std::cout<<std::endl;
	
}
