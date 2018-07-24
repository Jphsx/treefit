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
	std::cout<<"finished print"<<std::endl;
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

	double* cov1d = new double(size);
	//we need to know the jth particles number of params

	//do iterator trick
	std::vector<std::vector<string>::iterator> its;
	for(unsigned int i=0; i<mat.size(); i++){
		std::vector<string>::iterator it = mat.at(i).begin();
		its.push_back(it);
	}

	std::vector<string> mat_1d{};
	while(its.at(Nparts-1) < mat.at(Nparts-1).end()){
		
			std::cout<<"parsed "<<std::endl;
		
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
double* Covariance::matrix3DTo1D( std::vector<std::vector<std::vector<double> > > mat ){

	int size = 0;
	for(int i=0; i<mat.size(); i++){
		for(int j=0; j<mat.at(i).size(); j++){
			size += mat.at(i).at(j).size();
		}
	}	
	//assuming there is 1 element get N parts
	int Nparts = mat.size();

	double* cov1d = new double(size);
	//we need to know the jth particles number of params

	//do iterator trick
	std::vector<std::vector<std::vector<double> > ::iterator> its(Nparts);
	//need to do memory allocation for 3d vector
	std::vector<std::vector<double> :: iterator> itcol(Nparts);
	for(int i=0; i<Nparts; i++){
		its.at(i) = itcol;
		for( j=0; j<Nparts; j++){
			std::vector<double>::iterator it = mat.at(i).at(j).begin();
			its.at(i).at(j) = it;
		}
	}

	

	std::vector<string> mat_1d{};
	while(its.at(Nparts-1).at(Nparts-1) < mat.at(Nparts-1).at(Nparts-1).end()){
		
			std::cout<<"parsed "<<std::endl;
		
			for(int i=0; i<mat.size(); i++){
				for( int j=0; j<mat.at(i).size(); j++){
					if(its.at(i) < mat.at(i).end()){
						mat_1d.insert(mat_1d.end(), its.at(i).at(j), its.at(i).at(j)+nparams.at(j));
						its.at(i).at(j) = its.at(i).at(j) + nparams.at(j);
					}
				}
				
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
	//go row by row
	//this jacobian is by columns
	std::vector<string> jacobian{};
	jacobian.push_back( 1);// "dPx/de" );
	jacobian.push_back( 1);//"dPx/dtheta" );
	jacobian.push_back( 1);//"dPx/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( 2);//"dPy/de" );
	jacobian.push_back( 2);//"dPy/dtheta" );
	jacobian.push_back( 2);//"dPy/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( 3);//"dPz/de" );
	jacobian.push_back( 3);//"dPz/dtheta" );
	jacobian.push_back( 3);//"dPz/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( 4);//"dE/de" );
	jacobian.push_back( 4);//"dE/dtheta" );
	jacobian.push_back( 4);//"dE/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	
	std::cout<<"made jfo "<<std::endl;
	return jacobian;

}
std::vector<double> Covariance::constructLFOJacobian(Particle* p){
	std::cout<<"in the lfo jac "<<std::endl;
	//std::vector< std::vector<string> > jacobian{};
	/* dPx/dk dPx/dtheta dPx/dphi
	  dPy/dk dPy/dtheta dPy/dphi
	  dPz/dk dPz/dtheta dPz/dphi
	  dE/dk dE/dtheta dE/dphi */
	//this jacobian is by columns
	std::vector<string> jacobian{};
	
	jacobian.push_back( 5);//"dPx/dk" );
	jacobian.push_back( 5);//"dPx/dtheta" );
	jacobian.push_back( 5);//"dPx/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( 6);//"dPy/dk" );
	jacobian.push_back( 6);//"dPy/dtheta" );
	jacobian.push_back( 6);//"dPy/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( 7);//"dPz/dk" );
	jacobian.push_back( 7);//"dPz/dtheta" );
	jacobian.push_back( 7);//"dPz/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( 8);//"dE/dk" );
	jacobian.push_back( 8);//"dE/dtheta" );
	jacobian.push_back( 8);//"dE/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	std::cout<<"made LFO "<<std::endl;
	return jacobian;


}
double* Covariance::constructJacobian(std::vector<Particle*> parts, std::vector<int> combo){

	//local param size in particle will give the number of params
	std::vector<int> nparams = getnparamsvec(parts,combo);
	int Nparts = getNparts(combo);
	
	std::cout<<"NPARTS "<<Nparts<<std::endl;
	int Nparams = getNparams(parts,combo);
	std::cout<<"NPARAMS "<< Nparams<<std::endl;
	//start with 4d matrix
	// then unwrap the 4d into a 2d matrix
	//need to store the fit parts onto a local array with no gaps
	//this will preserve its position in the matrix
	std::vector<Particle*> p{};
	for(int i=0; i<combo.size(); i++){
		p.push_back(parts.at(combo.at(i)));
	}
	
	

	std::vector<std::vector<string> > jacobian(Nparts);
	
	std::cout<<"going into the big matrix "<<std::endl;
	for(int i = 0; i < Nparts; i++){
		if(p.at(i)->isTrack){
			//LFO
			jacobian.at(i) =  (constructLFOJacobian(p.at(i) ));
		}
		if(!p.at(i)->isTrack){
			//JFO
			jacobian.at(i) =(constructJFOJacobian(p.at(i) ));
		}
	}
	std::cout<<"finished big matrix"<<std::endl;

	/*std::vector<std::vector<string>::iterator> its;
	for(unsigned int i=0; i<jacobian.size(); i++){
		std::vector<string>::iterator it = jacobian.at(i).begin();
		its.push_back(it);
	}

	std::vector<string> jac_1d{};
	while(its.at(Nparts-1) < jacobian.at(Nparts-1).end()){
		
			std::cout<<"parsed "<<std::endl;
		
			for(int i=0; i<jacobian.size(); i++){
				
				if(its.at(i) < jacobian.at(i).end()){
					jac_1d.insert(jac_1d.end(), its.at(i), its.at(i)+nparams.at(i));
					its.at(i) = its.at(i) + nparams.at(i);
				}
				
			}
		
	}*/
	double* jac1d = matrix2Dto1D(jacobian, nparams);



	std::cout<<"finished parsing"<<std::endl;
	

	
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
	std::cout<<"allocated mem"<<std::endl;
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
		std::cout << " R , i "<< R << " "<< i <<std::endl;

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
	printSectoredCovarianceMatrix( cov);


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

	std::cout<<"SUB SIZE, GLOBLSIZE "<< subCombo.size() << " "<< globalCombo.size() <<std::endl;
	for(int i=0; i<subCombo.size(); i++){
	
		for(int j=0; j<globalCombo.size(); j++){
			if( subCombo.at(i) == globalCombo.at(j)){
				indices.push_back(j);
				break;
			}
		}		
	}
	std::cout<<"just did index matching"<<std::endl;
	
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

	printSectoredCovarianceMatrix(subcov );
	
	//print the submatrix
	//transform to 1d and return
	
	double * vec = matrix3DTo1D(subcov, getnparamsvec(parts,combo));
	//TODO return 1d submatrix
	return vec;

}



