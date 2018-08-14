#include "Covariance.h"

int Covariance::getNparams(std::vector<Particle*> parts, std::vector<int> combo){
	std::vector<int> nparams{};
	for(int i=0; i<combo.size(); i++){
		nparams.push_back( parts.at( combo.at(i) )->localParams.size() );
	}
	int Nparams = 0;
	for(int i =0; i<nparams.size(); i++){
		Nparams += nparams.at(i);
	}
	return Nparams;
	
}
std::vector<int> Covariance::getnparamsvec(std::vector<Particle*> parts, std::vector<int> combo){
	std::vector<int> nparams{};
	
	for(int i=0; i<combo.size(); i++){
		nparams.push_back( parts.at( combo.at(i) )->localParams.size() );
	}
	return nparams;
}
int Covariance::getNparts(std::vector<int> combo){
	int Nparts = combo.size();
	return Nparts;
}
double* Covariance::get4VecLD(double* cov){

	/*xx xy xz xE
	  yx yy yz yE
	  zx zy zz zE
	  Ex Ey Ez EE */
	double* lowerDiagonal = new double[10];
	lowerDiagonal[0] = cov[0];
	lowerDiagonal[1] = cov[4];
	lowerDiagonal[2] = cov[5];
	lowerDiagonal[3] = cov[8];
	lowerDiagonal[4] = cov[9];
	lowerDiagonal[5] = cov[10];
	lowerDiagonal[6] = cov[12];
	lowerDiagonal[7] = cov[13];
	lowerDiagonal[8] = cov[14];
	lowerDiagonal[9] = cov[15];

	return lowerDiagonal;	
}
void Covariance::printCovarianceMatrix(double* cov, int rows, int columns){
	//std::cout<<"the cov size "<< cov.size() << std::endl;
	for(int i=0; i<rows*columns; i++){
		if(i%columns == 0 ){ std::cout<<std::endl;}
		std::cout<<cov[i]<<" ";
		
	}
	std::cout<<std::endl;

}

void Covariance::printSectoredCovarianceMatrix(std::vector<std::vector<std::vector<double> > > cov ){
	for(int i=0; i<cov.size(); i++){
		for(int j=0; j<cov.at(i).size(); j++){
			std::cout<<"SECTOR "<<i<<" "<<j<<std::endl;
			for(int k=0; k<cov.at(i).at(j).size(); k++){
				std::cout<<cov.at(i).at(j).at(k)<<" ";
			}	
			std::cout<<std::endl;
		}
	}

}
double* Covariance::matrix2DTo1D( std::vector<std::vector<double>  > mat, std::vector<int> nparams ){
	//make vector and copy it? no
	//add all sizes together
	
	//this is the total number of elements
	int size = 0;
	for(int i=0; i<mat.size(); i++){
		size += mat.at(i).size();
		
	}	
	//assuming there is 1 element get N parts
	int Nparts = mat.size();

	double* cov1d = new double[size];
	//we need to know the jth particles number of params

	//do iterator trick
	std::vector<std::vector<double>::iterator> its;
	for(unsigned int i=0; i<mat.size(); i++){
		std::vector<double>::iterator it = mat.at(i).begin();
		its.push_back(it);
	}

	std::vector<double> mat_1d{};
	while(its.at(Nparts-1) < mat.at(Nparts-1).end()){
		
			for(int i=0; i<mat.size(); i++){
				
				if(its.at(i) < mat.at(i).end()){
					mat_1d.insert(mat_1d.end(), its.at(i), its.at(i)+nparams.at(i));
					its.at(i) = its.at(i) + nparams.at(i);
				}
				
			}
	}
	//copy vector onto double*
	for(int i=0; i<mat_1d.size(); i++){
		cov1d[i] = mat_1d.at(i);
	}
	
	return cov1d;

}
double* Covariance::matrix3DTo1D( std::vector<std::vector<std::vector<double> > > mat , std::vector<int> nparams){

	int size = 0;
	for(int i=0; i<mat.size(); i++){
		for(int j=0; j<mat.at(i).size(); j++){
			size += mat.at(i).at(j).size();
		}
	}	
	//assuming there is 1 element get N parts
	int Nparts = mat.size();

	double* cov1d = new double[size];
	//we need to know the jth particles number of params

	//do iterator trick
	std::vector<std::vector<std::vector<double>::iterator > > its(Nparts);
	//need to do memory allocation for 3d vector
	std::vector<std::vector<double>::iterator> itcol(Nparts);
	for(int i=0; i<Nparts; i++){
		its.at(i) = itcol;
		for( int j=0; j<Nparts; j++){
			std::vector<double>::iterator it = mat.at(i).at(j).begin();
			its.at(i).at(j) = it;
		}
	}

	

	std::vector<double> mat_1d{};
	while(its.at(Nparts-1).at(Nparts-1) < mat.at(Nparts-1).at(Nparts-1).end()){
		
			std::cout<<"parsed "<<std::endl;
		
		for(int i=0; i<mat.size(); i++){

			while(its.at(i).at(Nparts-1) < mat.at(i).at(Nparts-1).end()){
				for( int j=0; j<mat.at(i).size(); j++){
					if(its.at(i).at(j) < mat.at(i).at(j).end()){
						mat_1d.insert(mat_1d.end(), its.at(i).at(j), its.at(i).at(j)+nparams.at(j));
						its.at(i).at(j) = its.at(i).at(j) + nparams.at(j);
					}
				}
			}//while
				
		}
	}
	//copy vector onto double*
	for(int i=0; i<mat_1d.size(); i++){
		cov1d[i] = mat_1d.at(i);
	}
	
	return cov1d;

}
std::vector<double> Covariance::constructJFOJacobian(Particle* p){
	std::cout<<"in the JFO jac "<<std::endl;
	//std::vector< std::vector<string> > jacobian{};
	/*
	  dPx/de dPx/dtheta dPx/dphi
	  dPy/de dPy/dtheta dPy/dphi
	  dPz/de dPz/dtheta dPz/dphi
	  dE/de  dE/dtheta  dE/dphi
	*/ 
	
	std::vector<double> jacobian{};
	jacobian.push_back( p->v.Px()/p->v.P());// dPx/de
	jacobian.push_back( p->v.Pz()*p->v.Px()/p->v.Perp());// dPx/dtheta
	jacobian.push_back( -p->v.Py() );// dPx/dphi
	
	jacobian.push_back( p->v.Py()/p->v.P());// dPy/de
	jacobian.push_back( p->v.Pz()*p->v.Py()/p->v.Perp());// dPy/dtheta
	jacobian.push_back( p->v.Px());// dPy/dphi
	
	jacobian.push_back( p->v.Pz()/p->v.P());// dPz/de
	jacobian.push_back( -p->v.Perp());// dPz/dtheta
	jacobian.push_back( 0.0 );// dPz/dphi
	
	jacobian.push_back( 1.0);// dE/de
	jacobian.push_back( 0.0);// dE/dtheta
	jacobian.push_back( 0.0);// dE/dphi
	
	std::cout<<"made jfo "<<std::endl;
	return jacobian;
	


}
std::vector<double> Covariance::constructLFOJacobian(Particle* p){
	std::cout<<"in the lfo jac "<<std::endl;
	
	/* dPx/dk dPx/dtheta dPx/dphi
	  dPy/dk dPy/dtheta dPy/dphi
	  dPz/dk dPz/dtheta dPz/dphi
	  dE/dk dE/dtheta dE/dphi */
	
	std::vector<double> jacobian{};
	
	jacobian.push_back( -p->v.Perp() * p->v.Px() );//"dPx/dk" );
	jacobian.push_back( p->v.Pz()*p->v.Px()/ p->v.Perp());//"dPx/dtheta" );
	jacobian.push_back( -p->v.Py());//"dPx/dphi" );
	
	jacobian.push_back( -p->v.Perp() * p->v.Py() );//"dPy/dk" );
	jacobian.push_back( p->v.Pz()*p->v.Py()/ p->v.Perp() );//"dPy/dtheta" );
	jacobian.push_back( p->v.Px() );//"dPy/dphi" );
	
	jacobian.push_back( -p->v.Perp()*p->v.P() );//"dPz/dk" );
	jacobian.push_back( -p->v.Perp());//"dPz/dtheta" );
	jacobian.push_back( 0.0);//"dPz/dphi" );
	
	jacobian.push_back( 1.0);//"dE/dk" );
	jacobian.push_back( 0.0);//"dE/dtheta" );
	jacobian.push_back( 0.0);//"dE/dphi" );
	
	std::cout<<"made LFO "<<std::endl;
	return jacobian;

}
std::vector<double> Covariance::constructTPFOJacobian(Particle* p){
	
	const double c = 2.99792458e8; // m*s^-1        
  	const double mm2m = 1e-3;
  	const double eV2GeV = 1e-9;
  	const double eB = p->Bfield*c*mm2m*eV2GeV;
 
 	
	double cosLambda = 1 / sqrt(1 + p->track->getTanLambda()*p->track->getTanLambda() );
	double P = (eB/fabs(p->track->getOmega()))/cosLambda;
	double sinLambda = p->track->getTanLambda()*cosLambda;
//test	double cosPhi = cos(p->track->getPhi());
//test	double sinPhi = sin(p->track->getPhi());
	
	double adjustedphi = p->track->getPhi();
	if(adjustedphi<0){
		adjustedphi = 2*3.14159 + adjustedphi;
	}
	double cosPhi = cos(adjustedphi);
	double sinPhi = sin(adjustedphi);

	double px = P*cosLambda*cosPhi;
	double py = P*cosLambda*sinPhi;
	double pz = P*sinLambda;
	double E = sqrt( P*P + p->part->getMass() * p->part->getMass() );

	double omega = p->track->getOmega();
	double q = p->part->getCharge();
	//omega = fabs(omega);
	

	std::vector<double> jacobian{};
	
	jacobian.push_back(0); //dpx/dd0
	jacobian.push_back(-py); //dpx/dphi ...
 	jacobian.push_back(-px/omega  * q); //dpx/dome
	jacobian.push_back(0); //dpx/dz0
	jacobian.push_back(0); //dpx/dtanl

	jacobian.push_back(0); //dpy/dd0
	jacobian.push_back(px); //dpy/dphi
	jacobian.push_back(-py/omega  * q); //dpy/dome
	jacobian.push_back(0);//dpy/dz0
	jacobian.push_back(0); //dpy/dtanl

	jacobian.push_back(0);//dpz/dd0'
	jacobian.push_back(0);//dpz/dphi
	jacobian.push_back(-pz/omega  * q); //dpz/dome
	jacobian.push_back(0);//dpz/dz0
	jacobian.push_back(eB/omega   * q);//dpz/dtanl

	jacobian.push_back(0);//de/dd0
	jacobian.push_back(0);//de/dphi
	jacobian.push_back( -(P*P)/(omega*E)  * q );//de/dome 
	jacobian.push_back(0);//de/dz0
	jacobian.push_back( (pz*P*cosLambda)/E );//de/dtanl
	

	return jacobian;
}
double* Covariance::constructJacobian(std::vector<Particle*> parts, std::vector<int> combo,  int FO_Option){

	//local param size in particle will give the number of params
	std::vector<int> nparams = getnparamsvec(parts,combo);
	int Nparts = getNparts(combo);
	

	int Nparams = getNparams(parts,combo);

	
	//this will preserve its position in the matrix
	std::vector<Particle*> p{};
	for(int i=0; i<combo.size(); i++){
		p.push_back(parts.at(combo.at(i)));
	}
	
	

	std::vector<std::vector<double> > jacobian(Nparts);
	

	for(int i = 0; i < Nparts; i++){
		if(p.at(i)->isTrack){
			if( FO_Option == 1){
				jacobian.at(i) =  (constructLFOJacobian(p.at(i) ));
			}
			if( FO_Option == 2){
				jacobian.at(i) = (constructTPFOJacobian(p.at(i) ));
			}
		}
		if(!p.at(i)->isTrack){
			//JFO
			jacobian.at(i) =(constructJFOJacobian(p.at(i) ));
		}
	}


	double* jac1d = matrix2DTo1D(jacobian, nparams);




	

	
	return jac1d;


}
std::vector<std::vector<std::vector<double> > > Covariance::matrix1DTo3D(double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> combo){

	//put the global cov onto an vector so we can use iterators on it
	std::vector<double> _globalcov{};
	for(int i=0; i<dim*dim; i++){
		_globalcov.push_back(globalcov[i]);
	}
	
	std::vector<int> nparams = getnparamsvec(parts,combo);
	int Nparts = getNparts(combo);
	int Nparams = getNparams(parts,combo);

	std::vector<std::vector<std::vector<double> > > cov(Nparts);
	std::vector<std::vector<double> > covcol(Nparts);

	//memory management //covsubmatix is a sector in the 2d mat
	std::vector<double> covsubmatrix{};
	for(int i =0; i< cov.size(); i++){
		cov.at(i) = covcol;
		for(int j=0; j<cov.at(i).size(); j++){
			cov.at(i).at(j) = covsubmatrix;
		}
	}

	//extrapolate the gencov into the 3d matrix
	//keep looping over parameters over entire cov
	std::vector<double>::iterator it=_globalcov.begin();
	
	int R = 0;
	//param pointer
	std::vector<int>::iterator param_it = nparams.begin();
	
	int paramThreshold = *(param_it);

	for(int i=0; i<Nparams; i++){
		
		if( param_it < nparams.end() && i == paramThreshold  ){
				param_it++;
				paramThreshold += *(param_it);
				R++;
		}


		for(int j=0; j<nparams.size(); j++){
			//ith row jth column
			//extract params from nparams and put in ij sector
			covsubmatrix.clear();
			if(it < _globalcov.end()){
				covsubmatrix.insert(covsubmatrix.end(),it,it+nparams.at(j)); 
				cov.at(R).at(j).insert(cov.at(R).at(j).end(), covsubmatrix.begin(), covsubmatrix.end()); 
				it+= nparams.at(j);
			}
		}
	}
	

	//test print of the sectored matrix
	//printSectoredCovarianceMatrix( cov);


	return cov;
	

}

double* Covariance::getSubGlobalCov( double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> globalCombo, std::vector<int> subCombo){
	
	if(globalCombo.size() == subCombo.size()){
		//this is the same thing just return the global cov
		return globalcov;
	}


	//make the manageable sectored version of the global covariance matrix
	std::vector<std::vector<std::vector<double> > > sectoredGlobalCov = matrix1DTo3D(globalcov, dim, parts, globalCombo);

	//locate the indices on globalcombo of the particles matched between combo sets
	//these indices are the locations of the particles on the sectored matrix
	std::vector<int> indices{};
	//do a simple match


	for(int i=0; i<subCombo.size(); i++){
	
		for(int j=0; j<globalCombo.size(); j++){
			if( subCombo.at(i) == globalCombo.at(j)){
				indices.push_back(j);
				break;
			}
		}		
	}

	
	//now create memory space for the new submatrix
	std::vector<std::vector<std::vector<double> > > subcov(subCombo.size());
	std::vector<std::vector<double> > subcovcol(subCombo.size());

	//memory management
	std::vector<double> covsubmatrix{};
	for(int i =0; i< subcov.size(); i++){
		subcov.at(i) = subcovcol;
		for(int j=0; j<subcov.at(i).size(); j++){
			subcov.at(i).at(j) = covsubmatrix;
		}
	}

	//now go through and extract the proper matrix sectors
	//use a selection sort type method
	//do it in two loops for more readability
	
	std::vector<std::vector<double> > tempcov{};
	for(int i=0; i<indices.size(); i++){
		for(int j=0; j<indices.size(); j++){
			tempcov.push_back( sectoredGlobalCov.at(indices.at(i)).at(indices.at(j)) );
		}
	}
	//translate the 2d to 3d for printing and testing
	int R=0; 
	int J=0;
	int threshold = subCombo.size();
	for(int i=0; i< tempcov.size(); i++){
		if( i == threshold ){
			R++;
			threshold += threshold;
			J=0;
		}
		subcov.at(R).at(J) = tempcov.at(i);
		J++;
	}

	//printSectoredCovarianceMatrix(subcov );
	
	//print the submatrix
	//transform to 1d and return
	
	double * vec = matrix3DTo1D(subcov, getnparamsvec(parts,subCombo));
	//TODO return 1d submatrix
	return vec;

}
void Covariance::rescaleSector(std::vector<double>& covSector){
	//scale factors for the parameters d0,phi,omega,z0,tanLambda	
	std::vector<double> scaleFactor{1.e-2, 1., 1.e-3, 1.e-2, 1.};
	int threshold = 5;
	int sFindex = 0;
	for(int i=0; i<covSector.size(); i++){
		covSector.at(i) = covSector.at(i) * scaleFactor.at(sFindex);
		sFindex++;
		if(sFindex == threshold) sFindex =0;
	}
}
double* Covariance::rescaleGlobalCov(double* globalCov, int dim, std::vector<Particle*> parts, std::vector<int> combo){
	
		
	//go to 3d to find the track sectors 	
	std::vector<std::vector<std::vector<double> > > cov3d = matrix1DTo3D(globalCov, dim, parts, combo);

	//scan through  parts/combo, the order the tracks appear is the order they appear on the cov matrix
	//make a new parts vector to make the indices correspond with the 3d matrix
	std::vector<Particle*> p{};
	for(int i=0; i<combo.size(); i++){
		p.push_back( parts.at(combo.at(i)) );
	}
	
	//locate the tracks on the vector p
	for(int k=0; k<p.size(); k++){
		//if it is a track rescale at the corresponding sector/index
		if(p.at(k)->isTrack){
			
			//loop through the all sectors containing k and rescale kk, kj, jk
			for(int i=0; i<cov3d.size(); i++){
				for( int j=0; j<cov3d.at(i).size(); j++){
					if( i == k){//rescale
						rescaleSector( cov3d.at(i).at(j) );
					}
					if( j == k){//rescale
						rescaleSector( cov3d.at(i).at(j) );
					}
				}
			}
					

		}
	}
	//transform 3d back to 1d and return it
	double* rescaledCov = matrix3DTo1D( cov3d , getnparamsvec(parts, combo));
	return rescaledCov;
	
	
}
double* Covariance::removeVFOSectors(double* globalCov, int dim, std::vector<Particle*> parts, std::vector<int> combo){
	//for 1 vertex constraint the vfo sector is 5 parameters and it will be the last diagonal chunk in the global cov
	//this part is a little hacky:
	//add an artificial particle to parts so matrix 1d to 3d sees vfo as another particle
	//also add this artificial particle to the combo
	//this will allow the use of 1d->3d without restructuring
	Particle* fake = new Particle();
	//the functions only look at localParams.size() so we need to make this 5
	//then fake particle parameters will be consistent with the vfo added parameters
	//i think the vfo params are start1, start2, x, y, z
	std::vector<double> fakeParams{ 1,2,3,4,5 };
	fake->localParams = fakeParams;
	std::vector<Particle*> fakeParts = parts;
	fakeParts.push_back(fake);
	
	std::vector<int> fakeCombo = combo;
	fakeCombo.push_back( fakeParts.size() -1 );

	std::vector<std::vector<std::vector<double> > > global3d = matrix1DTo3D(globalCov, dim, fakeParts, fakeCombo);
	
	//with the 3d vector remove the last row and column
	std::vector<std::vector<std::vector<double> > > trimmed3d(global3d.size()-1);
	//memory allocation
	std::vector<std::vector<double> > trimmed3d_col(global3d.size()-1);
	for(int i=0; i<trimmed3d.size(); i++){
		trimmed3d.at(i) = trimmed3d_col;
	}
	//iterate over the new 3d matrix and transfer sectors from the global3d onto it
	for(int i=0; i<trimmed3d.size(); i++){ 
		for(int j=0; j<trimmed3d.size(); j++){
			trimmed3d.at(i).at(j) = global3d.at(i).at(j);
		}
	}
	//std::cout<<"PRINTING TRIMMED 3d MATRIX"<<std::endl;
//	printSectoredCovarianceMatrix( trimmed3d );
	
	//put the matrix back to 1d
	double* trimmed1d = matrix3DTo1D( trimmed3d, getnparamsvec(parts, combo) );

	return trimmed1d;
		


}
float* Covariance::get4VecCovariance(double* globalCov, int dim, std::vector<Particle*> parts, std::vector<int> globalCombo, std::vector<int> subCombo,  int FO_Option){
	

	//get Nparams
	int Nparams = getNparams(parts, subCombo);

	//our sub matrix of the global matrix
	double* subcov;
	
	//if we are doing vertex fitting with ONLY 1 VERTEX CONSTRAINT
	//we need to trim off the vfo elements of the covariance matrix
	if(FO_Option == 2){
		double* trimmedcov = removeVFOSectors(globalCov, dim, parts, subCombo);
		subcov = getSubGlobalCov(trimmedcov,dim-5, parts, globalCombo, subCombo);
			//if we use tpfo we need to rescale globalCov
		subcov = rescaleGlobalCov(subcov, getNparams( parts, subCombo),  parts, subCombo);
	}
	else{		
		//get the sub covariance matrix
		subcov = getSubGlobalCov(globalCov,dim, parts, globalCombo, subCombo);	
	}
	//get the jacobian for this submatrix
	//the jacobian retrieved is actually the transpose
	double* jacobian = constructJacobian(parts,subCombo, FO_Option);
	
	//do the matrix calculation
	//figure out all matrix dimensions
	//TMatrixD Dmatrix(4,Nparams, jacobian, "F");
	TMatrixD Dmatrix(Nparams,4,jacobian,"F");
	TMatrixD Vmatrix(Nparams,Nparams, subcov, "F");
 
        TMatrixD Covmatrix(4,4); 
	Covmatrix.Mult( TMatrixD( Dmatrix, TMatrixD::kTransposeMult, Vmatrix) ,Dmatrix);
	//Covmatrix.Mult( Dmatrix, TMatrixD( Vmatrix, TMatrixD::kMultTranspose, Dmatrix)); 

	//turn matrix into storable double*
	double* newcov = new double[Nparams*Nparams];
	int index =0;
	for(int i=0; i<4; i++){
		for(int j=0; j<4; j++){
			newcov[index] = Covmatrix(i,j);
			index++;
		}
	}
	printCovarianceMatrix(newcov,4,4);
	//convert the matrix to lower diagonal
	double* newLDcov = get4VecLD(newcov);
	float* newLDcovf = new float[10];
	//convert to float //TODO properly change everything to floats
	for(int i=0; i<10; i++){
		newLDcovf[i] = (float) newLDcov[i];
	}

	return newLDcovf;

}



