

#include "marlin/Processor.h"
#include "EVENT/ReconstructedParticle.h"
#include "IMPL/TrackImpl.h"
#include "EVENT/Track.h"
#include "EVENT/MCParticle.h"

#include <vector>
#include "IMPL/LCCollectionVec.h"
#include "TFile.h"
#include "TTree.h"

#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "gear/BField.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "TLorentzVector.h"

#include "LeptonFitObject.h"
#include "TrackParticleFitObject.h"
#include "JetFitObject.h"
#include "VertexFitObject.h"
#include "OPALFitterGSL.h"
//#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"

#include "TH1D.h"

// Marlin stuff
#include <marlin/Global.h>
// ROOT stuff
#include "TMath.h"
#include "TMatrixD.h"

#include <cstdlib>
#include <cmath>

#include "../src/TreeFit.cpp"

using namespace lcio ;

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

	void FindMassConstraintCandidates( LCCollectionVec* recparcol);

	/** Processor Parameters **/
	/** tree construction parameters **/
	std::vector<int> _preorderPdgs{};
	std::vector<double> _preorderMass{};
	std::string _preorderSerial;

};


