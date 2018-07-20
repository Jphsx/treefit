#ifndef _COVARIANCE_
#define _COVARIANCE_


//class that calculates the 4 vector covariance matrix
//based on the k-parameters and n-particle fit global covariance
//matrix
#include "lcio.h"
#include "EVENT/ReconstructedParticle.h"
#include "EVENT/Track.h"
#include <iostream>
#include "Particle.h"
#include <vector>
#include <string>
typedef lcio::Track Track ;
typedef lcio::ReconstructedParticle ReconstructedParticle ;
typedef std::string string;
class Covariance{

	public:
	//the jacobian derivatives are completely dependent on
	//the type of fit objects used
	static int getNparams(std::vector<Particle*> parts, std::vector<int> combo);
	static int getNparts(std::vector<Particle*> parts, std::vector<int> combo);
	

	//first testing construct a jacobian with strings
	//take in all fit parts and the combo for this resonance
	static std::vector<string> constructJacobian(std::vector<Particle*> fitparts, std::vector<int> fitCombo);

	static std::vector<string> constructJFOJacobian(Particle* p);
	static std::vector<string> constructLFOJacobian(Particle* p);
	//static void printCovarianceMatrix(std::vector<std::vector<string> > cov);
	static void printCovarianceMatrix(std::vector<string> cov, int dim);//dim = Nparam
	static void printSectoredCovarianceMatrix(std::vector<std::vector<std::vector<double> > > cov );
	//TODO 
	//static std::Vector<string> constuct TPFOJacobian();
	
	//TODO delete this stupid method	
	static std::vector<string> constructEmptyJacobian(int nrow, int ncol);
	//turn 1d covariance matrix into a more manageable 3d matrix
	static std::vector<std::vector<std::vector<double> > > sectorGlobalCov(double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> combo);
	
	static std::vector<double> getSubGlobalCov( double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> globalCombo, std::vector<int> subCombo);

};
#endif
