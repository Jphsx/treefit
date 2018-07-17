#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackImpl.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"
#include "lcio.h"
#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "TFile.h"
#include "TTree.h"
#include "CLHEP/Vector/LorentzVector.h"
#include "TLorentzVector.h"
//typedef CLHEP::HepLorentzVector LorentzVector ;
#include "LeptonFitObject.h"
#include "TrackParticleFitObject.h"
#include "JetFitObject.h"
#include "VertexFitObject.h"
#include "OPALFitterGSL.h"
//#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"
#include "TH1D.h"


typedef lcio::Track Track ;
typedef lcio::ReconstructedParticle ReconstructedParticle ;


#include "TreeFit.h"
#include "Matching.h"
#include "Particle.h"
#include "TTreeFactory.h"
#include "Covariance.h"


class treeFitter : public marlin::Processor {
  
	public:
  
	virtual marlin::Processor*  newProcessor() { return new treeFitter ; }

	treeFitter(const treeFitter&) = delete ;
	treeFitter& operator=(const treeFitter&) = delete ;
  
	treeFitter() ;

	/** Called at the beginning of the job before anything is read.
   	*  Use to initialize the proscessor, e.g. book histograms.
   	*/
  	virtual void init() ;
  	/** Called for every run.
   	*/
  	virtual void processRunHeader( LCRunHeader* run ) ;

  	/** Called for every event - the working horse.
   	*/
  	virtual void processEvent( LCEvent * evt ) ; 


 	/** Called after data processing for clean up.
   	*/
  	virtual void end() ;

	private:
	/*************
	track the event number for printing
	*************/
	int evtNo;
	/*************
	store the Bfield
	*************/
	double B;
	
	/**************
	LCIO Collection gathering methods
	Locates the specified track/reconstrucedParticle/MCParticle 
	collections an extracts the particles into global arrays
	**************/ 
  	bool FindPFOs( LCEvent* evt );
  	bool FindTracks( LCEvent* evt );
  	bool FindMCParticles( LCEvent* evt); 

	/**************
	The global data structures that house the particles for
	a given event
	**************/
	std::vector<ReconstructedParticle*>_pfovec{};
  	std::vector<Track*> _trackvec{};
  	std::vector<MCParticle*> _mcpartvec{};

	/*************
	The strings that hold the names of the input/output LCIO collections
	*************/
	std::string _inputTrackCollectionName;
	std::string _inputParticleCollectionName;
	std::string _inputMcParticleCollectionName;
	std::string _outputParticleCollectionName;

	/*************
	Creates the fit objects and performs the fit with all
	specified constraints
	*************/
	OPALFitterGSL* fitParticles(std::vector< std::vector<int>> fit);
	/*************
	Global vector to store all fit objects for post fit
	easy access
	*************/
	std::vector<ParticleFitObject*> FitObjects{};
	
	/************
	Function to created fit particles from the fit objects
	*************/
	//TODO

	/*************
	Functions to generate the proper output collections
	**************/
	//void createLCOutputParticles(LCCollectionVec* recparcol, std::vector<std::vector<int> > fit, double fitProb);
	ReconstructedParticleImpl* createLCOutputParticleTree(LCCollectionVec* recparcol,Node* root, std::vector<std::vector<int> > fit, double fitProb);

	void FindMassConstraintCandidates( LCCollectionVec* recparcol);

	/** Processor Parameters **/
	/** tree construction parameters **/
	std::vector<int> _preorderPdgs{};
	std::vector<float> _preorderMass{};
	std::string _preorderSerial;

	/*************
	 object that stores/manages
	 the tree and generates combinations
	*************/
	TreeFit* TFit;

	/*************
	File to store all the TTrees from the treefit
	And vector to store all the tree containers
	*************/
	TFile* file;
	//all trees created in init
	std::vector<TTreeFactory*> ttrees{};

	//cut parameters
	double _fitProbabilityCut;

};


