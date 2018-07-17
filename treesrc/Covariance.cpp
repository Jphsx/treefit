#include "Covariance.h"

int Covariance::getNparams(std::vector<Particle*> parts, std::vector<int> combo){
	std::vector<int> nparams{};
	for(int i=0; i<fitCombo.size(); i++){
		nparams.push_back( fitparts.at( fitCombo.at(i) )->localParams.size() );
	}
	int Nparams = 0;
	for(int i =0; i<nparams.size(); i++){
		Nparams += nparams.at(i);
	}
	return Nparams;
	
}
int Covariance::getNparts(std::vector<Particle*> parts, std::vector<int> combo){
	int Nparts = fitCombo.size();
	return Nparts;
}
std::vector<std::vector<string> > Covariance::constuctJFOJacobian(Particle* p ){

	std::vector< std::vector<string> > jacobian{};
	/*
	  dPx/de dPx/dtheta dPx/dphi
	  dPy/de dPy/dtheta dPy/dphi
	  dPz/de dPz/dtheta dPz/dphi
	  dE/de  dE/dtheta  dE/dphi
	*/ 
	//go row by row
	std::vector<string> row{};
	row.push_back( "dPx/de" );
	row.push_back( "dPx/dtheta" );
	row.push_back( "dPx/dphi" );
	jacobian.push_back(row);
	row.clear();
	row.push_back( "dPy/de" );
	row.push_back( "dPy/dtheta" );
	row.push_back( "dPy/dphi" );
	jacobian.push_back(row);
	row.clear();
	row.push_back( "dPz/de" );
	row.push_back( "dPz/dtheta" );
	row.push_back( "dPz/dphi" );
	jacobian.push_back(row);
	row.clear();
	row.push_back( "dE/de" );
	row.push_back( "dE/dtheta" );
	row.push_back( "dE/dphi" );
	jacobian.push_back(row);
	row.clear();
	
	return jacobian;

}
std::vector<std::vector<string> > Covariance::constuctLFOJacobian(Particle* p ){
	
	std::vector< std::vector<string> > jacobian{};
	/* dPx/dk dPx/dtheta dPx/dphi
	  dPy/dk dPy/dtheta dPy/dphi
	  dPz/dk dPz/dtheta dPz/dphi
	  dE/dk dE/dtheta dE/dphi */
	std::vector<string> row{};
	
	row.push_back( "dPx/dk" );
	row.push_back( "dPx/dtheta" );
	row.push_back( "dPx/dphi" );
	jacobian.push_back(row);
	row.clear();
	row.push_back( "dPy/dk" );
	row.push_back( "dPy/dtheta" );
	row.push_back( "dPy/dphi" );
	jacobian.push_back(row);
	row.clear();
	row.push_back( "dPz/dk" );
	row.push_back( "dPz/dtheta" );
	row.push_back( "dPz/dphi" );
	jacobian.push_back(row);
	row.clear();
	row.push_back( "dE/dk" );
	row.push_back( "dE/dtheta" );
	row.push_back( "dE/dphi" );
	jacobian.push_back(row);
	row.clear();
	
	return jacobian;


}
std::vector<std::vector<string> > Covariance::constructEmptyJacobian(int nrow, int ncol){
	std::vector< std::vector<string> > jacobian{};
	std::vector<string> row{};
	for(int i =0; i<ncol; i++){
		row.push_back("0");
	}
	for(int i=0; i<nrow; i++){
		jacobian.push_back(row);
	}
	row.clear();
	return jacobian;
}
std::vector<string> Covariance::constructJacobian(std::vector<Particle*> fitparts, std::vector<int> fitCombo){

	//local param size in particle will give the number of params
	std::vector<int> nparams{};
	int Nparts = fitCombo.size();
	for(int i=0; i<fitCombo.size(); i++){
		nparams.push_back( fitparts.at( fitCombo.at(i) )->localParams.size() );
	}

	int Nparams = 0;
	for(int i =0; i<nparams.size(); i++){
		Nparams += nparams.at(i);
	}
	//start with 4d matrix
	// then unwrap the 4d into a 2d matrix
	//need to store the fit parts onto a local array with no gaps
	//this will preserve its position in the matrix
	std::vector<Particle*> fitp{};
	for(int i=0; i<fitCombo.size(); i++){
		fitp.push_back(fitCombo.at(i));
	}
	
	//with fitp make the jacobian
	//4 d matrix
	std::vector< std::vector< std::vector< std::vector<string> > > > jacobian{};
	//ij are the ith and jth particle
	//kl are the ith and jth kth and lth parameters
	
	for(int i = 0; i < Nparts; i++){
		for(int j=0; j< Nparts; j++){
			if (i==j){
				if(fitp.at(i)->isTrack){
					//LFO
					jacobian.at(i).push_back(constuctLFOJacobian(fitp.at(i) ));
				}
				if(!fitp.at(i)->isTrack){
					//JFO
					jacobian.at(i).push_back(constuctJFOJacobian(fitp.at(i) ));
				}
			}
			//if i!= j then make empty guy with
			// i rows and j cols
			jacobian.at(i).push_back(constructEmptyJacobian(npars.at(i), npars.at(j)));
			
		}
	}

	//turn the 4d matrix into a 1 d matrix
	std::vector<string> jac{};
	for(int i=0; i<Nparts; i++){
		for(int j=0; j<Nparts; j++){
			for(int k=0; k< jacobian.at(i).at(j).size(); k++){
				for(int l=0; l< jacobian.at(i).at(j).at(k).size(); j++){
					jac.push_back(jacobian.at(i).at(j).at(k).at(l);
				}

			}
		}
	}

	return jac;


}
void Covariance::printCovarianceMatrix(std::vector<string> cov, int dim){
	for(int i=0; i<dim; i++){
		if(i%dim == 0){ std::cout<<std::endl;
		std::cout<<cov.at(i)<<" ";
	}
	std::cout<<std::endl;
}

