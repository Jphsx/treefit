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

	std::vector<float> preorderMass;
	preorderMass.push_back(0.0);
	registerProcessorParameter("preorderMass",
				   "Preorder traversal of masses GEV in the particle tree",
				   _preorderMass,
				   preorderMass);

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
				"McParticleCollectionName" ,
				"Name of the MCParticle input collection" ,
				_inputMcParticleCollectionName,
				inputMcParticleCollectionName);


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

return;
}

/********************
Print processor paramters, initalize global variables like 
event number, initialize the output TTree
*********************/
void treeFitter::init() {
	//evtno set to 0
	evtNo=0;
	//set the BField
	B = marlin::Global::GEAR->getBField().at(gear::Vector3D(0.,0.,0.)).z();
	//print all processor input parameters
	printParameters();

	//build particle tree
	TFit = new TreeFit();
	string SerialDelimiter = " ,";
	TFit->ParticleTree->treeInit(_preorderPdgs, _preorderSerial, _preorderMass, SerialDelimiter);
	TFit->ParticleTree->printTree(TFit->ParticleTree->Root);

	//build TTREE
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
	if(FindPFOs(evt) && FindTracks(evt)){
	//temp
	
		//TODO: if(using mcparticles){
		//find mcparticles
		FindMCParticles(evt);
	
		//call fitter
		this->FindMassConstraintCandidates(recparcol);
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
//input a fit from the fit table 
OPALFitterGSL* treeFitter::fitParticles(std::vector< std::vector<int>> fit){
		//general procedure
		OPALFitterGSL *  fitter = new OPALFitterGSL();

		//make a FO vector to contain both neutral and charged FOs, the index of the FO should match the index of the recopart in TFit
 		std::vector<ParticleFitObject*> FO_vec(TFit->recoparts.size());
		//use the first node it populate all the FOs
		//do this assignment at the end
		// FitObjects = FO_vec;

		//adding this index var for readability
		int recoindex=-1;
		//iterate over the tree roots combination (node 0)
		for(int i=0; i<fit.at(0).size(); j++){
			//make an object at index
			//fit.at(0).at(j) this is index on reco part
			recoindex = fit.at(0).at(j)
			if(TFit->recoparts.at(recoindex)->isTrack){
				//this is a track make LFO
				FO_vec.at(recoindex) = new LeptonFitObject(
				TFit->recoparts.at(recoindex)->track, 
				TFit->recoparts.at(recoindex)->Bfield, 
				TFit->recoparts.at(recoindex)->part->getMass()); 
			}
			else{
				//this is not a track make JFO
				FO->vec.at(recoindex) = new JetFitObject(
				TFit->recoparts.at(recoindex)->localParams.at(0), 
 				TFit->recoparts.at(recoindex)->localParams.at(1), 
				TFit->recoparts.at(recoindex)->localParams.at(2),  
				TFit->recoparts.at(recoindex)->localErrors.at(0), 
				TFit->recoparts.at(recoindex)->localErrors.at(1), 
				TFit->recoparts.at(recoindex)->localErrors.at(2), 
				TFit->recoparts.at(recoindex)->part->getMass() );
				
			}
			
			
		}
		for(unsigned int i=0; i< recoparts.size(); i++){
			//iterate over the fit table to 

		}
		
		//make mass constraint objects
		
		
		//have to make local parameterizations for errors on jfo
		//make arrays of JFOs and LFO

		//add corresponding FOs to MC objects

		//add all FOs to fitter

		//add all mc constraints to fitter
		return fitter;
}

void treeFitter::FindMassConstraintCandidates(LCCollectionVec * recparcol) {
	std::cout<<"EVENT "<<evtNo<<std::endl;
	//print each particle directly 
	//cout.precision(10);	
	//start populating the Particle* structure in TFit
	for(unsigned int i=0; i<_pfovec.size(); i++){
		Particle* pc = new Particle(_pfovec.at(i), Matching::MatchParticleToTrack(_pfovec.at(i), _trackvec,B),B);
		
		TFit->addrecopart(pc);
	}	
	
	TFit->printParticles(TFit->recoparts);
	std::cout<<std::endl;

	//Do the fits (print them)
	TFit->initTable();
	TFit->generatefitcombinations(TFit->ParticleTree->Root, TFit->recoIDs);
	//print the fit table	
	TFit->printTable();
			

	//advance to next event
	evtNo++;
	_pfovec.clear();
	_trackvec.clear();
	//might need to run a destructor here first
	TFit->clearEvent();
	return;
}


