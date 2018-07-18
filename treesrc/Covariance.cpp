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

	std::vector<std::vector<std::vector<string>::iterator> >its(Nparts);
	std::vector<std::vector<string>::iterator> itscol(Nparts);
	/*for(unsigned int i=0; i<matrices.size(); i++){
		std::vector<double>::iterator it = matrices.at(i).begin();
		its.push_back(it);
	}*/
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
					jac_1d.insert(jac_1d.end(), its.at(i).at(j), its.at(i).at(j)+nparams.at(j));
					its.at(i).at(j) = its.at(i).at(j) + nparams.at(j);
				}
			}
	}
				

	std::cout<<"finished parsing"<<std::endl;
	
	//SEG IS HERE START HERE
/*	std::vector<string> jac{};
	for(int i=0; i<Nparts; i++){
		for(int j=0; j<Nparts; j++){
			for(int k=0; k< jacobian.at(i).at(j).size(); k++){
				for(int l=0; l< jacobian.at(i).at(j).at(k).size(); j++){
					jac.push_back(jacobian.at(i).at(j).at(k).at(l));
				}

			}
		}
	}

	return jac;*/
	
	return jac_1d;


}
void Covariance::printCovarianceMatrix(std::vector<string> cov, int dim){
	for(int i=0; i<dim; i++){
		if(i%dim == 0 ){ std::cout<<std::endl;}
		std::cout<<cov.at(i)<<" ";
		
	}
	std::cout<<std::endl;
}

