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
int Covariance::getNparts(std::vector<Particle*> parts, std::vector<int> combo){
	int Nparts = combo.size();
	return Nparts;
}
void Covariance::printCovarianceMatrix(std::vector<string> cov, int dim){
	std::cout<<"the cov size "<< cov.size() << std::endl;
	for(int i=0; i<cov.size(); i++){
		if(i%dim == 0 ){ std::cout<<std::endl;}
		std::cout<<cov.at(i)<<" ";
		
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
std::vector<string> Covariance::constructJFOJacobian(Particle* p){
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
	jacobian.push_back( "dPx/de" );
	jacobian.push_back( "dPx/dtheta" );
	jacobian.push_back( "dPx/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( "dPy/de" );
	jacobian.push_back( "dPy/dtheta" );
	jacobian.push_back( "dPy/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( "dPz/de" );
	jacobian.push_back( "dPz/dtheta" );
	jacobian.push_back( "dPz/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( "dE/de" );
	jacobian.push_back( "dE/dtheta" );
	jacobian.push_back( "dE/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	
	std::cout<<"made jfo "<<std::endl;
	return jacobian;

}
std::vector<string> Covariance::constructLFOJacobian(Particle* p){
	std::cout<<"in the lfo jac "<<std::endl;
	//std::vector< std::vector<string> > jacobian{};
	/* dPx/dk dPx/dtheta dPx/dphi
	  dPy/dk dPy/dtheta dPy/dphi
	  dPz/dk dPz/dtheta dPz/dphi
	  dE/dk dE/dtheta dE/dphi */
	//this jacobian is by columns
	std::vector<string> jacobian{};
	
	jacobian.push_back( "dPx/dk" );
	jacobian.push_back( "dPx/dtheta" );
	jacobian.push_back( "dPx/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( "dPy/dk" );
	jacobian.push_back( "dPy/dtheta" );
	jacobian.push_back( "dPy/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( "dPz/dk" );
	jacobian.push_back( "dPz/dtheta" );
	jacobian.push_back( "dPz/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	jacobian.push_back( "dE/dk" );
	jacobian.push_back( "dE/dtheta" );
	jacobian.push_back( "dE/dphi" );
	//jacobian.push_back(row);
	//row.clear();
	std::cout<<"made LFO "<<std::endl;
	return jacobian;


}
std::vector<string> Covariance::constructEmptyJacobian(int nrow, int ncol){
	std::cout<<"in the empty jac "<<std::endl;
	//std::vector< std::vector<string> > jacobian{};
	
	std::vector<string> jacobian{};
	for(int i =0; i<ncol*nrow; i++){
		jacobian.push_back("0");
	}
	//for(int i=0; i<nrow; i++){
	//	jacobian.push_back(row);
//	}
//	row.clear();
	return jacobian;
}
std::vector<string> Covariance::constructJacobian(std::vector<Particle*> fitparts, std::vector<int> fitCombo){

	//local param size in particle will give the number of params
	std::vector<int> nparams{};
	int Nparts = fitCombo.size();
	for(int i=0; i<fitCombo.size(); i++){
		nparams.push_back( fitparts.at( fitCombo.at(i) )->localParams.size() );
	}
	std::cout<<"NPARTS "<<Nparts<<std::endl;
	int Nparams = 0;
	for(int i =0; i<nparams.size(); i++){
		Nparams += nparams.at(i);
	}
	std::cout<<"NPARAMS "<< Nparams<<std::endl;
	//start with 4d matrix
	// then unwrap the 4d into a 2d matrix
	//need to store the fit parts onto a local array with no gaps
	//this will preserve its position in the matrix
	std::vector<Particle*> fitp{};
	for(int i=0; i<fitCombo.size(); i++){
		fitp.push_back(fitparts.at(fitCombo.at(i)));
	}
	
	//with fitp make the jacobian
	//4 d matrix
	//std::vector< std::vector< std::vector< std::vector<string> > > > jacobian{};
	//allocate memory for 4d array
	//std::vector< std::vector< std::vector< std::vector<string> > > > jacobian(Nparts);
	//std::vector< std::vector< std::vector<string> > > jacobianrow(Nparts);

	std::vector<std::vector<string> > jacobian(Nparts);
	
	std::cout<<"going into the big matrix "<<std::endl;
	for(int i = 0; i < Nparts; i++){
		if(fitp.at(i)->isTrack){
			//LFO
			jacobian.at(i) =  (constructLFOJacobian(fitp.at(i) ));
		}
		if(!fitp.at(i)->isTrack){
			//JFO
			jacobian.at(i) =(constructJFOJacobian(fitp.at(i) ));
		}
	}
	std::cout<<"finished big matrix"<<std::endl;

	std::vector<std::vector<string>::iterator> its;
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
		
	}

/* save all this code for global cov decomposition ////////////////////////////////////////////////////////////////
	std::vector< std::vector< std::vector<string> > > jacobian(Nparts);
	std::vector<std::vector<string> > jacobiancol(Nparts);

	std::cout<<"trying new allocation method"<<std::endl;
	for(int i =0; i<Nparts; i++){
		jacobian.at(i) = jacobiancol;
	}


	//ij are the ith and jth particle
	//kl are the ith and jth kth and lth parameters
	std::cout<<"going into the big matrix "<<std::endl;
	for(int i = 0; i < Nparts; i++){
		for(int j=0; j< Nparts; j++){
			if (i==j){
				if(fitp.at(i)->isTrack){
					//LFO
					jacobian.at(i).at(j)=  (constructLFOJacobian(fitp.at(i) ));
				}
				if(!fitp.at(i)->isTrack){
					//JFO
					jacobian.at(i).at(j) =(constructJFOJacobian(fitp.at(i) ));
				}
			}else{
			//if i!= j then make empty guy with
			// i rows and j cols
			jacobian.at(i).at(j)=(constructEmptyJacobian(nparams.at(i), nparams.at(j)));
			}
		}
	}
	std::cout<<"finished big matrix"<<std::endl;
	//turn the 4d matrix into a 1 d matrix


	std::cout<<"trying new allocation method"<<std::endl;
	for(int i =0; i<Nparts; i++){
		its.at(i) = itscol;
	}
	for(unsigned int i=0; i< jacobian.size(); i++){
		for( unsigned int j=0; j< jacobian.at(i).size(); j++){
			std::vector<string>::iterator it = jacobian.at(i).at(j).begin();
			its.at(i).at(j) = it;
		}
	}
	//parse this guy

	std::vector<string> jac_1d{};
	while(its.at(Nparts-1).at(Nparts-1) < jacobian.at(Nparts-1).at(Nparts-1).end()){
		
			std::cout<<"parsed "<<std::endl;
		
			for(int i=0; i<jacobian.size(); i++){
				for(int j=0; j<jacobian.at(i).size(); j++){
					if(its.at(i).at(j) < jacobian.at(i).at(j).end()){
						jac_1d.insert(jac_1d.end(), its.at(i).at(j), its.at(i).at(j)+nparams.at(j));
						its.at(i).at(j) = its.at(i).at(j) + nparams.at(j);
					}
				}
			}
		
	}
*/				

	std::cout<<"finished parsing"<<std::endl;
	

	
	return jac_1d;


}
std::vector<std::vector<std::vector<double> > > Covariance::sectorGlobalCov(double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> combo){

	//put the global cov onto an vector so we can use iterators on it
	std::vector<double> _globalcov{};
	for(int i=0; i<dim*dim; i++){
		_globalcov.push_back(globalcov[i]);
	}
	
	std::vector<int> nparams{};
	int Nparts = combo.size();
	for(int i=0; i<combo.size(); i++){
		nparams.push_back( parts.at( combo.at(i) )->localParams.size() );
	}
	std::cout<<"NPARTS "<<Nparts<<std::endl;
	int Nparams = 0;
	for(int i =0; i<nparams.size(); i++){
		Nparams += nparams.at(i);
	}
	std::cout<<"NPARAMS "<< Nparams<<std::endl;

	std::vector<std::vector<std::vector<double> > > cov(Nparts);
	std::vector<std::vector<double> > covcol(Nparts);

	//memory management
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

			//if(i%dim == 0){ std::cout<<std::endl;}
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
	for(int i=0; i<cov.size(); i++){
		for(int j=0; j<cov.at(i).size(); j++){
			std::cout<<"SECTOR "<<i<<" "<<j<<std::endl;
			for(int k=0; k<cov.at(i).at(j).size(); k++){
				std::cout<<cov.at(i).at(j).at(k)<<" ";
			}	
			std::cout<<std::endl;
		}
	}


	return cov;
	

}
std::vector<double> Covariance::getSubGlobalCov( double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> globalCombo, std::vector<int> subCombo){
	
/*	std::vector<int> nparams{};
	int Nparts = combo.size();
	for(int i=0; i<combo.size(); i++){
		nparams.push_back( parts.at( combo.at(i) )->localParams.size() );
	}
	std::cout<<"NPARTS "<<Nparts<<std::endl;
	int Nparams = 0;
	for(int i =0; i<nparams.size(); i++){
		Nparams += nparams.at(i);
	}
	std::cout<<"NPARAMS "<< Nparams<<std::endl; */

	//make the manageable sectored version of the global covariance matrix
	std::vector<std::vector<std::vector<double> > > sectoredGlobalCov = sectorGlobalCov(globalcov, dim, parts, globalCombo);

	//locate the indices on globalcombo of the particles matched between combo sets
	//these indices are the locations of the particles on the sectored matrix
	std::vector<int> indices{};
	//do a simple match
	for(int i=0; i<subCombo.size(); i++){
		for(int j=0; j<globalCombo.size(); i++){
			if( subCombo.at(i) == globalCombo.at(i)){
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
	std::vector<double> emptyvec{};
	return emptyvec;

}



