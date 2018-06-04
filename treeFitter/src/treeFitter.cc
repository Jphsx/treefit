#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "gear/BField.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
typedef CLHEP::HepLorentzVector LorentzVector ;
typedef CLHEP::Hep3Vector Vector3D ;

// MarlinKinfit stuff
#include "LeptonFitObject.h"
#include "JetFitObject.h"
#include "OPALFitterGSL.h"
#include "NewFitterGSL.h"
#include "NewtonFitterGSL.h"
#include "MassConstraint.h"

// Marlin stuff
#include <marlin/Global.h>
// the event display

// ROOT stuff
#include "TMath.h"
#include "TMatrixD.h"

#include <cstdlib>
#include <cmath>
#include "treeFitter.h"

using namespace lcio;

treeFitter atreeFitter;


/********************
Register processor parameters and input/output collections
*********************/
treeFitter::treeFitter() : marlin::Processor("treeFitter") {

	//register parameters

	//tree construction parameters
	std::vector<int> preorderPdgs;
	preorderPdgs.push_back(0);
	registerProcessorParameter("preorderPdgs",
				   "Preorder traversal of pdg codes in the particle tree",
				   _preorderPdgs,
				   preorderPdgs);
/*
	std::vector<float> preorderMass;
	preorderMass.push_back(0.0);
	registerProcessorParameter("preorderMass",
				   "Preorder traversal of masses GEV in the particle tree",
				   _preorderMass,
				   preorderMass);
*/
	std::string preorderSerial = " ) ";
	registerProcessorParameter("preorderSerial",
				   "Preorder Serialization of the tree using unique node IDs",
				   _preorderSerial,
				   preorderSerial);

 
	//input collection parameters
	std::string inputParticleCollectionName = "x";
  	registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     	"InputParticleCollectionName" , 
			     	"Input Particle Collection Name "  ,
			     	_inputParticleCollectionName,
			      	inputParticleCollectionName);

 	std::string inputTrackCollectionName = "x";
  	registerInputCollection( LCIO::TRACK,
				"InputTrackCollectionName" ,
				"Input Track Collection Name " ,
				_inputTrackCollectionName,
				inputTrackCollectionName);

	std::string inputMcParticleCollectionName = "x";
	registerInputCollection( LCIO::MCPARTICLE,
				"MCParticleCollection" ,
				"Name of the MCParticle input collection" ,
				_inputMcParticleCollectionName,
				inputMcParticleCollectionName);

/*
	//output collection parameters this should be modified
	std::string outputParticleCollectionName = "x";
	registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			     	"OutputParticleCollectionName" , 
			     	"Output Particle Collection Name "  ,
			     	_outputParticleCollectionName,
			     	outputParticleCollectionName);

  	std::string outputTrackCollectionName = "x";
  	registerOutputCollection( LCIO::TRACK,
			    	"OutputTrackCollectionName",
			    	"Output Particle Collection Name" ,
			    	_outputTrackCollectionName,
			    	outputTrackCollectionName);
*/
return;
}

/********************
Print processor paramters, initalize global variables like 
event number, initialize the output TTree
*********************/
void treeFitter::init() {
//evtno
//treeinit TTree
//build particle tree
	TFit = new TreeFit();
	string delimiter = " ,";
	TFit->ParticleTree->treeInit(_preorderPdgs, _preorderSerial, _preorderMass, delimiter, 1);
  return;
}
/*******************
output header information
********************/
void treeFitter::processRunHeader( LCRunHeader* run) {
  streamlog_out(MESSAGE) << " processRunHeader "  << run->getRunNumber() << std::endl ; 
}
/*******************
call supporting functions to build reconstructed particle collection
call the fitter
********************/
void treeFitter::processEvent( LCEvent * evt ) { 
	// Make a new vector of particles
  	LCCollectionVec * recparcol = new LCCollectionVec(LCIO::RECONSTRUCTEDPARTICLE);

	//find pfos
	//find tracks
	//if(FindPFOs(evt) && FindTracks(evt)){
	//temp
	if(true){
		//TODO: if(using mcparticles){
		//find mcparticles
		//FindMCParticles(evt);
	
		//call fitter
		FindMassConstraintCandidates(recparcol);
	}


	//add collection to event
	evt->addCollection( recparcol,  _outputParticleCollectionName.c_str() );

	return;
}
/*****************
finish and write out to rootfile
******************/
void treeFitter::end(){
	//if(_fitAnalysis){
	//	rootFile->Write();
	//}
  	return;
}
/*****************
locate the pfo collection with specified name
populated the global pfo vectors with particles from that collection for this event
******************/
bool treeFitter::FindPFOs( LCEvent* evt ) {

	bool collectionFound = false;

  	// clear old global pfovector
	_pfovec.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;
	
	//iterate over collections, find the matching name
  	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){
     
		//if found print name and number of elements
    		if(*itname==_inputParticleCollectionName){ 
			LCCollection* collection = evt->getCollection(*itname);
			std::cout<< "Located Pfo Collection "<< *itname<< " with "<< collection->getNumberOfElements() << " elements " <<std::endl;
			collectionFound = true;

 			//add the collection elements to the global vector
      			for(int i=0; i<collection->getNumberOfElements(); i++){
				ReconstructedParticle* recoPart = dynamic_cast<ReconstructedParticle*>(collection->getElementAt(i));
				_pfovec.push_back(recoPart);
      			}
    		}
  	}
	
	if(!collectionFound){
		std::cout<<"Pfo Collection "<< _inputParticleCollectionName << "not found"<<std::endl;
	}

   
	return collectionFound;
}
/*****************
locate the track collection with specified name
populate the global track vectors with tracks from the collection for this event
******************/
bool treeFitter::FindTracks( LCEvent* evt ) {

	bool collectionFound = false;

	// clear old global track vector
 	_trackvec.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;

	//iterate over collections, find the matching name
 	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){   

		//if found print name and number of elements 
    		if(*itname==_inputTrackCollectionName){
      			LCCollection* collection = evt->getCollection(*itname);
      			std::cout<< "Located Track Collection "<< *itname<< " with "<< collection->getNumberOfElements() << " elements " <<std::endl;
     		
			collectionFound = true;
		
			//add the collection elements to the global vector
      			for(int i=0; i<collection->getNumberOfElements(); i++){
				Track* track = dynamic_cast<Track*>(collection->getElementAt(i));
				_trackvec.push_back(track);
      			}
    		}
  	}

  	if(!collectionFound){
		std::cout<<"Track Collection "<< _inputTrackCollectionName << "not found"<<std::endl;
	}

  	return collectionFound;
}
/*****************
locate the mcparticle collection with specified name
populate the global mcp vector with mc particles from the collection for this event
******************/
bool treeFitter::FindMCParticles( LCEvent* evt ){
   
	bool collectionFound = false;

  	// clear old global MCParticle vector
  	_mcpartvec.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;

	//iterate over collections, find the matching name
  	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){    
    
		//if found print name and number of elements 
		if(*itname==_inputMcParticleCollectionName){
      			LCCollection* collection = evt->getCollection(*itname);
     			std::cout<< "Located MC Collection "<< *itname<< " with "<< collection->getNumberOfElements() << " elements " <<std::endl;
      			collectionFound = true;
      
			//add the collection elements to the global vector
			for(int i=0;i<collection->getNumberOfElements();i++){
				MCParticle* mcPart = dynamic_cast<MCParticle*>(collection->getElementAt(i));
				_mcpartvec.push_back(mcPart);

       
      			}
    		}
  	}

  	if(!collectionFound){
		std::cout<<"MC Collection "<< _inputMcParticleCollectionName << "not found"<<std::endl;
	}
  
  	return collectionFound;
}

void treeFitter::FindMassConstraintCandidates(LCCollectionVec * recparcol) {
	//print global tree

}


