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
#include "TMath.h"
#include "TMatrixD.h"
typedef lcio::Track Track ;
typedef lcio::ReconstructedParticle ReconstructedParticle ;
typedef std::string string;
class Covariance{

	public:
	//main method to call
	//Options 1=LFO 2=TPFO
	static float* get4VecCovariance(double* globalCov, int dim, std::vector<Particle*> parts, std::vector<int> globalCombo, std::vector<int> subCombo,  int FO_Option);
	static double* get4VecLD(double* cov);
	

	//the jacobian derivatives are completely dependent on
	//the type of fit objects used
	static int getNparams(std::vector<Particle*> parts, std::vector<int> combo);
	static int getNparts( std::vector<int> combo);
	static std::vector<int> getnparamsvec(std::vector<Particle*> parts, std::vector<int> combo);
	

	//first testing construct a jacobian with strings
	//take in all fit parts and the combo for this resonance
	static double* constructJacobian(std::vector<Particle*> parts, std::vector<int> combo, int FO_Option);
	static std::vector<double> constructJFOJacobian(Particle* p);
	static std::vector<double> constructLFOJacobian(Particle* p);
	static std::vector<double> constructTPFOJacobian(Particle* p);
	
	static void printCovarianceMatrix(double* cov, int rows, int columns);//dim = Nparam
	static void printSectoredCovarianceMatrix(std::vector<std::vector<std::vector<double> > > cov );
	//TODO 
	
	//manipulation for tpfo cov matrix rescaling and calculation
	static double* rescaleGlobalCov(double* globalCov, int dim, std::vector<Particle*> parts, std::vector<int> combo);
	static void rescaleSector(std::vector<double>& covSector);
	//currently this will only work for a single vertex fit
	//since multiple vertex fits wont work, i dont know how to manipulate the global cov which is created from multiple vertex constraints (unknown dimensions/ parameter order)
	static double* removeVFOSectors(double* globalCov, int dim, std::vector<Particle*> parts, std::vector<int> combo);	

	//produce covariance sub matrix
	static double* getSubGlobalCov( double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> globalCombo, std::vector<int> subCombo);

	//matrix management functions
	static double* matrix2DTo1D( std::vector<std::vector<double>  > mat, std::vector<int> nparams );
	static double* matrix3DTo1D( std::vector<std::vector<std::vector<double> > > mat, std::vector<int> nparams );
	//turn 1d covariance matrix into a more manageable 3d matrix
	static std::vector<std::vector<std::vector<double> > > matrix1DTo3D(double* globalcov, int dim, std::vector<Particle*> parts, std::vector<int> combo);

	//testing function force diagonalized matrix
	static double* forcediagonalmatrix(std::vector<std::vector<std::vector<double> > > 3dmat);//TODO write this

};
#endif
