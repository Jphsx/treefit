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
		localParams.push_back(v->Theta());//theta
		localParams.push_back(v->Phi());//phi
		//reconstructed particle covmatrix must be modified 
		//to use the E,theta,phi error model
		localErrors.push_back(std::sqrt(p->getCovMatrix()[0]));//dE
		localErrors.push_back(std::sqrt(p->getCovMatrix()[2]));//dtheta
		localErrors.push_back(std::sqrt(p->getCovMatrix()[5]));//dphi
	}
	

}

Particle::Particle(JetFitObject* jfo, TrackParticleFitObject* tpfo, int pdg, float mass ){
	//can either be jfo or tfo only
	if(tpfo==NULL){
		isTrack = false;
		//this is not a track, make ReconstructedParticle
		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		ParticleIDImpl* newPDG = new ParticleIDImpl();
		newPDG->setPDG(pdg);
		newPDG->setLikelihood(1.0);
		//for readability add local param variables
		float E = jfo->getParam(0);
		float Theta = jfo->getParam(1);
		float Phi = jfo->getParam(2);
		//calculate px,py,pz
		float* mom = new float[3];
		mom[0] = sqrt( E*E - mass*mass)*cos(Theta)*sin(Phi);
		mom[1] = sqrt( E*E - mass*mass)*sin(Theta)*sin(Phi);
 		mom[2] = sqrt( E*E - mass*mass)*cos(Theta);	
		
		p->setMomentum(mom);
		p->setEnergy(E);
		//give the reco part an E,theta,phi cov matrix
		//we need to construct the lower diagonal manually
		float* cov = new float[6];
		int index = 0;
		for(int i=0; i<=2; i++){
			for(int j=0; j<=i; j++){
				cov[index]=jfo->getCov(i,j);
				index++;	
			}
		}
		p->setCovMatrix(cov);
		p->setMass(mass);
		p->setCharge(0.0);
		p->addParticleID(newPDG);
		p->setParticleIDUsed(newPDG);
		p->setType(pdg);

		part = p;
	}
	else{
		isTrack = true;
		//this is a track, make a Track*
		TrackImpl* t = new TrackImpl();
		//we need to rescale all the parameters
		//in TrackParticleFitObject*
		std::vector<double> scaleFactor{1.e-2, 1., 1.e-3, 1.e-2, 1., 1., 1.};

		t->setD0(tpfo->getParam(0)*scaleFactor.at(0)); //Impact parameter in r-phi
		t->setPhi(tpfo->getParam(1)*scaleFactor.at(1)); //phi of track at reference point (primary vertex)
		t->setOmega(tpfo->getParam(2)*scaleFactor.at(2));// signed curvature in 1/mm 
		t->setZ0(tpfo->getParam(3)*scaleFactor.at(3)); //Impact parameter in r-z
		t->setTanLambda(tpfo->getParam(4)*scaleFactor.at(4));// dip of the track in r-z at primary vertex
		//manually make the lower diagonal covariance matrix 
		float* cov = new float[15];	
		int index = 0;
		for(int i=0; i<=4; i++){
			for(int j=0; j<=i; j++){
				cov[index]=tpfo->getCov(i,j)*scaleFactor.at(i)*scaleFactor.at(j);
				index++;	
			}
		}
	
		t->setCovMatrix(cov);
		track = t;
		//also set bfield
		Bfield = tpfo->bfield;

		//now make a reconstructed particle to go with
		//with the track and store additional details
		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		ParticleIDImpl* newPDG = new ParticleIDImpl();
		newPDG->setPDG(pdg);
		newPDG->setLikelihood(1.0);
		
		float* mom = new float[3];
		std::vector<double> mom_vec = getTrackPxPyPz( t, tpfo->bfield);
		mom[0] = mom_vec.at(0);
		mom[1] = mom_vec.at(1);
 		mom[2] = mom_vec.at(2);	
		
		float P = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
		p->setMomentum(mom);
		p->setEnergy( sqrt(P*P + mass*mass ) );

		p->setMass(mass);
		p->setCharge(tpfo->getCharge());
		p->addParticleID(newPDG);
		p->setParticleIDUsed(newPDG);
		p->setType(pdg);
		//dont worry about setting the cov in the 
		//reconstructedparticle, just only use the
		//track covariance matrix
		part = p;	
	}

	//do tlv and error arrays
	if(isTrack){
		v = getTLorentzVector(track,part->getMass(),Bfield);
		localParams.push_back(track->getD0());
		localParams.push_back(track->getPhi());
		localParams.push_back(track->getOmega());
		localParams.push_back(track->getZ0());
		localParams.push_back(track->getTanLambda());
		localErrors.push_back(std::sqrt(track->getCovMatrix()[0]));//d0 
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[2]));//phi
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[5]));//ome
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[9]));//z0
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[14]));//tanL
	}
	else{
		v = getTLorentzVector(part);
		localParams.push_back(part->getEnergy());//E
		localParams.push_back(v->Theta());//theta
		localParams.push_back(v->Phi());//phi
		//reconstructed particle covmatrix must be modified 
		//to use the E,theta,phi error model
		localErrors.push_back(std::sqrt(part->getCovMatrix()[0]));//dE
		localErrors.push_back(std::sqrt(part->getCovMatrix()[2]));//dtheta
		localErrors.push_back(std::sqrt(part->getCovMatrix()[5]));//dphi
	}

}
//need to supply the old impact parameters
Particle::Particle(JetFitObject* jfo, LeptonFitObject* lfo, int pdg, float mass , float d0, float z0, double B){
	//can either be jfo or lfo only
	if(lfo==NULL){
		isTrack = false;
		//this is not a track, make ReconstructedParticle
		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		ParticleIDImpl* newPDG = new ParticleIDImpl();
		newPDG->setPDG(pdg);
		newPDG->setLikelihood(1.0);
		//for readability add local param variables
		float E = jfo->getParam(0);
		float Theta = jfo->getParam(1);
		float Phi = jfo->getParam(2);
		//calculate px,py,pz
		float* mom = new float[3];
		mom[0] = sqrt( E*E - mass*mass)*cos(Theta)*sin(Phi);
		mom[1] = sqrt( E*E - mass*mass)*sin(Theta)*sin(Phi);
 		mom[2] = sqrt( E*E - mass*mass)*cos(Theta);	
		
		p->setMomentum(mom);
		p->setEnergy(E);
		//give the reco part an E,theta,phi cov matrix
		//we need to construct the lower diagonal manually
		float* cov = new float[6];
		int index = 0;
		for(int i=0; i<=2; i++){
			for(int j=0; j<=i; j++){
				cov[index]=jfo->getCov(i,j);
				index++;	
			}
		}
		p->setCovMatrix(cov);
		p->setMass(mass);
		p->setCharge(0.0);
		p->addParticleID(newPDG);
		p->setParticleIDUsed(newPDG);
		p->setType(pdg);

		part = p;
	}
	else{
		isTrack = true;
		//this is a track, make a Track*
		TrackImpl* t = new TrackImpl();
		std::vector<double> trackparams{};
		trackparams = getTrackHelix(lfo, d0, z0, B);

		t->setD0(d0); //Impact parameter in r-phi
		t->setPhi(lfo->getParam(2)); //phi of track at reference point (primary vertex)
		t->setOmega(trackparams.at(2));// signed curvature in 1/mm 
		t->setZ0(z0); //Impact parameter in r-z
		t->setTanLambda(trackparams.at(4));// dip of the track in r-z at primary vertex
		//manually make the lower diagonal covariance matrix 
		float* cov = new float[15];	
		int index = 0;
		for(int i=0; i<=2; i++){
			for(int j=0; j<=i; j++){
				//cov[index]=tpfo->getCov(i,j)*scaleFactor.at(i)*scaleFactor.at(j);
				cov[index] = lfo->getCov(i,j);
				index++;	
			}
		}
	
		t->setCovMatrix(cov);
		track = t;
		//also set bfield
		Bfield = B;

		//now make a reconstructed particle to go with
		//with the track and store additional details
		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		ParticleIDImpl* newPDG = new ParticleIDImpl();
		newPDG->setPDG(pdg);
		newPDG->setLikelihood(1.0);
		
		float* mom = new float[3];
		std::vector<double> mom_vec = getTrackPxPyPz( t, Bfield);
		mom[0] = mom_vec.at(0);
		mom[1] = mom_vec.at(1);
 		mom[2] = mom_vec.at(2);	
		
		float P = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
		p->setMomentum(mom);
		p->setEnergy( sqrt(P*P + mass*mass ) );

		p->setMass(mass);
		p->setCharge(lfo->getParam(0)*sqrt(mom[0]*mom[0] + mom[1]*mom[1]));
		p->addParticleID(newPDG);
		p->setParticleIDUsed(newPDG);
		p->setType(pdg);
		//dont worry about setting the cov in the 
		//reconstructedparticle, just only use the
		//track covariance matrix
		part = p;
	}
	if(isTrack){
		v = getTLorentzVector(track,part->getMass(),Bfield);
		localParams.push_back(lfo->getParam(0));
		localParams.push_back(lfo->getParam(1));
		localParams.push_back(lfo->getParam(2));
		
		localErrors.push_back(std::sqrt(track->getCovMatrix()[0]));//k 
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[2]));//theta
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[5]));//phi
            	
	}
	else{
		v = getTLorentzVector(part);
		localParams.push_back(part->getEnergy());//E
		localParams.push_back(v->Theta());//theta
		localParams.push_back(v->Phi());//phi
		//reconstructed particle covmatrix must be modified 
		//to use the E,theta,phi error model
		localErrors.push_back(std::sqrt(part->getCovMatrix()[0]));//dE
		localErrors.push_back(std::sqrt(part->getCovMatrix()[2]));//dtheta
		localErrors.push_back(std::sqrt(part->getCovMatrix()[5]));//dphi
	}

}
void Particle::printTrack(Track* t){
	std::cout<<"Track: (d0,phi,ome,z0,tanL) "<< 
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
std::vector<double> Particle::getTrackHelix(LeptonFitObject* lfo, double d0, double z0, double B){
	//const double c = 2.99792458e8; // m*s^-1        
  	//const double mm2m = 1e-3;
  	//const double eV2GeV = 1e-9;
  	//const double eB = BField*c*mm2m*eV2GeV;
	std::vector<double> helixparams{};	
	
	double tanlambda = tan(lfo->getParam(1)); //does this angle need adjusted?
	//double omega = eBField/(fitp.P()*coslambda);
	double omega = lfo->getParam(0)*B;
	/*if(meast->getOmega() < 0){
		omega = -omega;
	}*/
	helixparams.push_back(d0);
	helixparams.push_back(lfo->getParam(2));
	helixparams.push_back(omega);
	helixparams.push_back(z0);
	helixparams.push_back(tanlambda);

	
	return helixparams;
}
void Particle::printTrackPxPyPz(Track* t, double B){
	std::vector<double> txtytz = getTrackPxPyPz(t,B);
	std::cout<<"Track: (Px,Py,Pz) "
		<<txtytz.at(0)<<" "
		<<txtytz.at(1)<<" "
		<<txtytz.at(2)<<std::endl;
}
void Particle::printLocalParameters(std::vector<double> params){
	std::cout<<"Local Params: ";
	for(unsigned int i=0; i<params.size(); i++){
		std::cout<< params.at(i) <<" ";
	}
	std::cout<<std::endl;
}
void Particle::printLocalErrors(std::vector<double> errors){
	std::cout<<"Local Errors: ";
	for(unsigned int i=0; i<errors.size(); i++){
		std::cout<< errors.at(i) <<" ";
	}
	std::cout<<std::endl;
}
void Particle::printCovarianceMatrix(std::vector<float> cov, int npar){
	//the number of elements in the lower diagonal is
	//sum( i ) where i= 0->npar
	int nelem = 0;
	int k=0;
	while(k<=npar){
		nelem = nelem + k;
		k++;
	}
	//reuse k for linebreaking
	k=0;
	int kinc=2;
	for(int i=0; i<nelem; i++){
		std::cout<<cov[i]<<" ";
		if(i==k){
			std::cout<<std::endl;
			k = k + kinc;
			kinc++;
		}
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
	std::cout<<"Covariance Matrix: "<<std::endl;
	if(pc->isTrack){
		printCovarianceMatrix(pc->track->getCovMatrix(),pc->localParams.size());
	}
	else{
		printCovarianceMatrix(pc->part->getCovMatrix(),pc->localParams.size());
	}
	std::cout<<std::endl;
	
}
