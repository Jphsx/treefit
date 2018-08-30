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
		localParams.push_back(v.Theta());//theta
		localParams.push_back(v.Phi());//phi
		//reconstructed particle covmatrix must be modified 
		//to use the E,theta,phi error model
		localErrors.push_back(std::sqrt(p->getCovMatrix()[0]));//dE
		localErrors.push_back(std::sqrt(p->getCovMatrix()[2]));//dtheta
		localErrors.push_back(std::sqrt(p->getCovMatrix()[5]));//dphi
	}
	

}

Particle::Particle(JetFitObject* jfo, TrackParticleFitObject* tpfo, VertexFitObject* vfo, int pdg, float mass ){
	//can either be jfo or tfo only
	if(tpfo==NULL){
		isTrack = false;
		//this is not a track, make ReconstructedParticle
		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		//ParticleIDImpl* newPDG = new ParticleIDImpl();
		//newPDG->setPDG(pdg);
		//newPDG->setLikelihood(1.0);
		
		
		//calculate px,py,pz
		float* mom = new float[3];
		mom[0] = jfo->getPx();
		mom[1] = jfo->getPy();
 		mom[2] = jfo->getPz();	
		
		p->setMomentum(mom);
		p->setEnergy(jfo->getE());
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
		//p->addParticleID(newPDG);
		//p->setParticleIDUsed(newPDG);
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

		//set the reference point to the vertex from vfo
		ThreeVector vtx = vfo->getVertex();
		float* ref = new float[3];
		ref[0] = vtx.getX();
		ref[1] = vtx.getY();
		ref[2] = vtx.getZ();

		t->setReferencePoint(ref);		

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
		//ParticleIDImpl* newPDG = new ParticleIDImpl();
		//newPDG->setPDG(pdg);
		//newPDG->setLikelihood(1.0);
		
		float* mom = new float[3];
		//std::vector<double> mom_vec = getTrackPxPyPz( t, tpfo->bfield);
		mom[0] = tpfo->getPx();
		mom[1] = tpfo->getPy();
 		mom[2] = tpfo->getPz();	
		
		//float P = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
		p->setMomentum(mom);
		p->setEnergy( tpfo->getE() );

		p->setMass(mass);
		p->setCharge(tpfo->getCharge());
		//p->addParticleID(newPDG);
		//p->setParticleIDUsed(newPDG);
		p->setType(pdg);
		//dont worry about setting the cov in the 
		//reconstructedparticle, just only use the
		//track covariance matrix
		part = p;	
	}

	//do tlv and error arrays
	if(isTrack){
		//v = getTLorentzVector(track,part->getMass(),Bfield);
		//TODO change this to direct creation from getPx getPy etc..
		TLorentzVector tlv;
		//tlv.SetPxPyPzE(tpfo->getPx(), tpfo->getPy(), tpfo->getPz(), tpfo->getE());
		tlv.SetXYZM(tpfo->getPx(), tpfo->getPy(), tpfo->getPz(), mass);
		v = tlv;
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
		//v = getTLorentzVector(part);
		TLorentzVector tlv;
		//tlv.SetPxPyPzE(jfo->getPx(), jfo->getPy(), jfo->getPz(), jfo->getE());
		tlv.SetXYZM(jfo->getPx(), jfo->getPy(), jfo->getPz(), mass);
		v = tlv;
		localParams.push_back(part->getEnergy());//E
		localParams.push_back(v.Theta());//theta
		localParams.push_back(v.Phi());//phi
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
		//ParticleIDImpl* newPDG = new ParticleIDImpl();
		//newPDG->setPDG(pdg);
		//newPDG->setLikelihood(1.0);
		
		
		//calculate px,py,pz
		float* mom = new float[3];
		mom[0] = jfo->getPx();
		mom[1] = jfo->getPy();
 		mom[2] = jfo->getPz();	
		
		p->setMomentum(mom);
		p->setEnergy(jfo->getE());
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
		//p->addParticleID(newPDG);
		//p->setParticleIDUsed(newPDG);
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
		//ParticleIDImpl* newPDG = new ParticleIDImpl();
		//newPDG->setPDG(pdg);
		//newPDG->setLikelihood(1.0);
		
		float* mom = new float[3];
		//std::vector<double> mom_vec = getTrackPxPyPz( t, Bfield);
		mom[0] = lfo->getPx();
		mom[1] = lfo->getPy();
 		mom[2] = lfo->getPz();	
		
		//float P = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
		p->setMomentum(mom);
		p->setEnergy( lfo->getE() );

		p->setMass(mass);
		p->setCharge(lfo->getParam(0)*sqrt(mom[0]*mom[0] + mom[1]*mom[1]));
		//p->addParticleID(newPDG);
		//p->setParticleIDUsed(newPDG);
		p->setType(pdg);
		//dont worry about setting the cov in the 
		//reconstructedparticle, just only use the
		//track covariance matrix
		part = p;
	}
	if(isTrack){
		//v = getTLorentzVector(track,part->getMass(),Bfield);
		TLorentzVector tlv;
		//tlv.SetPxPyPzE(lfo->getPx(), lfo->getPy(), lfo->getPz(),lfo->getE() );
		tlv.SetXYZM(lfo->getPx(), lfo->getPy(), lfo->getPz(), mass);
		v = tlv;
		localParams.push_back(lfo->getParam(0));
		localParams.push_back(lfo->getParam(1));
		localParams.push_back(lfo->getParam(2));
		
		localErrors.push_back(std::sqrt(track->getCovMatrix()[0]));//k 
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[2]));//theta
            	localErrors.push_back(std::sqrt(track->getCovMatrix()[5]));//phi
            	
	}
	else{ 
		//v = getTLorentzVector(part); 
		TLorentzVector tlv;
		//tlv.SetPxPyPzE(jfo->getPx(), jfo->getPy(), jfo->getPz(), jfo->getE() );
		tlv.SetXYZM(jfo->getPx(), jfo->getPy(), jfo->getPz(), mass);
		v = tlv;
		localParams.push_back(part->getEnergy());//E   
		localParams.push_back(v.Theta());//theta
		localParams.push_back(v.Phi());//phi
		//reconstructed particle covmatrix must be modified 
		//to use the E,theta,phi error model
		localErrors.push_back(std::sqrt(part->getCovMatrix()[0]));//dE
		localErrors.push_back(std::sqrt(part->getCovMatrix()[2]));//dtheta
		localErrors.push_back(std::sqrt(part->getCovMatrix()[5]));//dphi
	}

}
//TODO change oldtrk to new oldPart
std::vector<double> Particle::constructSameTrackJacobian(Track* t1, Track* t2 ){
 
	//p1 is original p2 is p'

	double d01 = t1->getD0();
 	double omega1 = t1->getOmega();
	double q1 = omega1/fabs(omega1);
	double phi1 = t1->getPhi();
	double z01 = t1->getZ0();
	double tanLambda1 = t1->getTanLambda();	

	
	double d02 = t2->getD0();
	double phi2 = t2->getPhi();
	double s = -(phi2-phi1)/omega1;

	std::vector<double> jacobian{};
	
	jacobian.push_back( cos(phi2-phi1) ); //dd0'/dd0
	jacobian.push_back( -(q1/omega1 - d01)*sin(phi2-phi1) ); //dd0'/dphi
	jacobian.push_back( ((q1*q1)/(omega1*omega1))*(cos(phi2-phi1)-1) ); //dd0'/domega
	jacobian.push_back( 0 ); //dd0'/dz0
	jacobian.push_back( 0 ); //dd0'/dtl
	
	jacobian.push_back( sin(phi2-phi1)/(q1/omega1 - d02 )); //dphi'/dd0
	jacobian.push_back( (q1/omega1 - d01)*cos(phi2-phi1)/ (q1/omega1 - d02) ); //dphi'/dphi
	jacobian.push_back( ((q1*q1)/(omega1*omega1))*sin(phi2-phi1)/ (q1/omega1 - d02) ); //dphi'/domega
	jacobian.push_back( 0 ); //dphi'/dz0
	jacobian.push_back( 0 ); //dphi'/dtl

	jacobian.push_back( 0 ); //domega'/dd0
	jacobian.push_back( 0 ); //domega'/dphi
	jacobian.push_back( 1 ); //domega'/domega
	jacobian.push_back( 0 ); //domega'/dz0
	jacobian.push_back( 0 ); //domega'/dtl
	
	jacobian.push_back( 0 ); //dz0'/dd0
	jacobian.push_back( 0 ); //dz0'/dphi
	jacobian.push_back( 0 ); //dz0'/domega
	jacobian.push_back( 1 ); //dz0'/dz0
	jacobian.push_back( s ); //dz0'/dtl

	jacobian.push_back( 0 ); //dtl'/dd0
	jacobian.push_back( 0 ); //dtl'/dphi
	jacobian.push_back( 0 ); //dtl'/domega
	jacobian.push_back( 0 ); //dtl'/dz0
	jacobian.push_back( 1 ); //dtl'/dtl

	return jacobian;
}
//p1 is unprimed p2 is primed
//return LowerDiagonal
float* Particle::transformSameTrackCov(double* oldcov, Track* t1, Track* t2){
	
	std::vector<double> jacobian = constructSameTrackJacobian(t1, t2 );
	//convert this to a 1d vec
	double* jac = new double[25];
	//double* oldcov_ = new double[25];
	for(int i=0; i<jacobian.size(); i++){
		jac[i] = jacobian.at(i);
		//oldcov_[i] = oldcov.at(i);
	}
	
		
	TMatrixD Dmatrix(5,5,jac,"F");
	TMatrixD Vmatrix(5,5, oldcov, "F");
	TMatrixD Covmatrix(5,5); 
	Covmatrix.Mult( TMatrixD( Dmatrix, TMatrixD::kTransposeMult, Vmatrix) ,Dmatrix);

	float* newcov = new float[15];
	int index =0;
	for(int i=0; i<=4; i++){
		for(int j=0; j<=i; j++){
			newcov[index] = (float) Covmatrix(i,j);
			index++;
		}
	}

	return newcov;
}
void Particle::printTrack(Track* t){
	std::cout<<"Track: (d0,phi,ome,z0,tanL) "<< 
		t->getD0()<<" "<<
		t->getPhi()<<" "<<
		t->getOmega()<<" "<<
		t->getZ0()<<" "<<
		t->getTanLambda()<<std::endl;
	std::cout<<"Reference Point (x,y,z): ";
	const float* ref = t->getReferencePoint();
	std::cout<<ref[0]<<" "<<ref[1]<<" "<<ref[2]<<std::endl; 
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
Particle::Particle(Particle* oldPart, std::vector<double> vtx){
	//this constructor creates a new particle by updating a track to new reference point
	//and creating an updated reconstructed particle for this track
double factor = 1e-2;

	Track* oldtrk = oldPart->track;
	TrackImpl* t = new TrackImpl();\
	//with the new reference vertex update the parameters

	//get the original reference point for easy manipulation
	const float* oldref = oldtrk->getReferencePoint();
	std::vector<double> ref(3);
	for(int i=0; i<ref.size(); i++){
		ref.at(i) = oldref[i];
	}


	//readability variables
	double dx = vtx.at(0) - ref.at(0);
	double dy = vtx.at(1) - ref.at(1);
	double q = oldtrk->getOmega()/fabs(oldtrk->getOmega());
	double R = q/oldtrk->getOmega();
	double d0 = oldtrk->getD0();
	double phi = oldtrk->getPhi();
	double z0 = oldtrk->getZ0();
	double omega = oldtrk->getOmega();
	double tanLambda = oldtrk->getTanLambda();
	

	//phi2
	double phiNew = atan2(sin(phi)- dx/(R-d0) , cos(phi) + dy/(R-d0)  );
	
	//double d0New = d0 + dx*sin(phi) - dy*cos(phi) + (dx*cos(phi) + dy*sin(phi))*tan( (phiNew-phi)/2 );
	double a = (-2*dx*sin(phi)+2*dy*cos(phi))/(R-d0);
	double b = (dx*dx + dy*dy)/( (R-d0)*(R-d0) );
	//double d0New = R - (R-d0)*sqrt( 1+a+b  );
	double d0New = 0;

	double s = -(phiNew-phi)/omega;
	//double z0New = z0 + s*tanLambda;
	double z0New = 0;

	t->setD0(d0New); //Impact parameter in r-phi
	t->setPhi(phiNew); //phi of track at reference point (primary vertex)
	t->setOmega(omega);// signed curvature in 1/mm 
	t->setZ0(z0New); //Impact parameter in r-z
	t->setTanLambda(tanLambda);// dip of the track in r-z at primary vertex
	//set the reference point as the new vertex
	float* newref = new float[3];
	newref[0] = vtx.at(0);
	newref[1] = vtx.at(1);
	newref[2] = vtx.at(2);

	t->setReferencePoint(newref);

	//take the cov matrix, make it square
	std::vector<float> oldcov = oldtrk->getCovMatrix();

	//easiest quick fix, justmake a 2d array
	//do memory allocation
	std::vector<std::vector<double> > cov(5);
	std::vector<double> covcol(5);
	for(int i=0; i<cov.size(); i++){
		cov.at(i) = covcol;
	}
	int index=0;
	//make the square matrix
	for(int i=0; i<cov.size(); i++){
		for(int j=0; j<=i; j++){
			cov.at(i).at(j) = (double) oldcov[index];
			//we can also make it symmetric here
			cov.at(j).at(i) = (double) oldcov[index];
			index++;
		}
	}
	//convert the cov matrix back to 1d for transformation
	index = 0;
	double* cov1d = new double[25];
	for(int i=0; i<cov.size(); i++){
		for(int j=0; j<cov.at(i).size(); j++){
			cov1d[index] = cov.at(i).at(j);
			index++;
		}
	}

	//with the full square cov we can now apply the jacobian transformation
	t->setCovMatrix(transformSameTrackCov(cov1d, oldtrk, t) );

	track = t;
	//also set bfield
	Bfield = oldPart->Bfield;

	//print testing for track comparison
	std::cout<<"OldTrack ";
	printTrack(oldtrk);
	std::cout<<"NewTrack ";
	printTrack(t);
	std::cout<<"old cov"<<std::endl;
	printCovarianceMatrix(oldcov, 5);
	std::cout<<"new cov"<<std::endl;
	printCovarianceMatrix(t->getCovMatrix(),5);

	//now make a reconstructed particle to go with
	//with the track and store additional details
	ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
	//ParticleIDImpl* newPDG = new ParticleIDImpl();
	//newPDG->setPDG(pdg);
	//newPDG->setLikelihood(1.0);

	//recalculate E/P
		
	float* mom = new float[3];
	std::vector<double> mom_vec = getTrackPxPyPz( t, Bfield);
	mom[0] = mom_vec[0];
	mom[1] = mom_vec[1];
 	mom[2] = mom_vec[2];	

	//define mass
	double mass = oldPart->part->getMass();
		
	double P = sqrt(mom[0]*mom[0] + mom[1]*mom[1] + mom[2]*mom[2]);
	p->setMomentum(mom);
	p->setEnergy( sqrt(P*P + mass*mass ));

	p->setMass(mass);
	p->setCharge(q);
	//p->addParticleID(newPDG);
	//p->setParticleIDUsed(newPDG);
	p->setType(oldPart->part->getType());
	//dont worry about setting the cov in the 
	//reconstructedparticle, just only use the
	//track covariance matrix
	part = p;
	
	
	TLorentzVector tlv;

	tlv.SetXYZM(mom_vec.at(0), mom_vec.at(1), mom_vec.at(2), mass);
	v = tlv;
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

void Particle::printReconstructedParticle(ReconstructedParticle* p){
	const double* mom = p->getMomentum();	
	std::cout<<"Particle "<< p->getType() <<": "<<
	"(Px,Py,Pz,E,M,q) "<<
	mom[0]<< " "<<mom[1]<< " "<<mom[2]<< " "
	<<p->getEnergy()<<" "<<p->getMass()<<" "
	<<p->getCharge()<<std::endl;
		
	
}
void Particle::printTLorentzVector(TLorentzVector v){

	std::cout<<"TLV: (Px,Py,Pz,P,E,M) "<<
		v.Px()<< " " <<
		v.Py()<< " " <<
		v.Pz()<< " " <<
		v.P() << " " <<
		v.E() << " " <<
		v.M() << " " <<std::endl;
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
	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = B*c*mm2m*eV2GeV;
	std::vector<double> helixparams{};	
	
	double tanlambda = 1.0/tan(lfo->getParam(1)); //does this angle need adjusted?
	//double omega = eBField/(fitp.P()*coslambda);
	double omega = lfo->getParam(0)*eB;
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

void Particle::printCovarianceMatrix(float* cov, int npar){
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
TLorentzVector Particle::getTLorentzVector(ReconstructedParticle* p){
	TLorentzVector tlv; 
	const double* mom = p->getMomentum();
	tlv.SetXYZM(mom[0],mom[1],mom[2],p->getMass());
	//tlv.SetPxPyPzE(mom[0],mom[1],mom[2],p->getEnergy());
	
	return tlv;
}
TLorentzVector Particle::getTLorentzVector(Track* t, double Mass, double B){
	TLorentzVector tlv;
	std::vector<double> txtytz = getTrackPxPyPz(t, B);
	double P = sqrt(txtytz.at(0)*txtytz.at(0) + txtytz.at(1)*txtytz.at(1) +txtytz.at(2)*txtytz.at(2) );
	double E = sqrt(P*P + Mass*Mass);
	tlv.SetXYZM(txtytz.at(0),txtytz.at(1),txtytz.at(2), Mass);
	//tlv.SetPxPyPzE(txtytz.at(0),txtytz.at(1),txtytz.at(2),E);
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
