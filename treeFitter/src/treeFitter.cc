
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
	registerProcessorParameter("preorderPdgs",
				   "Preorder traversal of pdg codes in the particle tree",
				   _preorderPdgs,
				   preorderPdgs);

	std::vector<int> preorderMass;
	registerProcessorParameter("preorderMass",
				   "Preorder traversal of masses [GeV] in the particle tree",
				   _preorderMass,
				   preorderMass);

	std::string preorderSerial = " ) ";
	registerProcessorParameter("preorderSerial",
				   "Preorder Serialization of the tree using unique node IDs",
				   _preorderSerial,
				   preorderSerial); 

return;
}

/********************
Print processor paramters, initalize global variables like 
event number, initialize the output TTree
*********************/
void treeFitter::init() {
  if(_printing>1)printParameters(); 
//evtno
//treeinit TTree
//build particle tree
	globalTree->treeInit(_preorderPdgs, _preorderSerial, _preorderMass, " ", 1);
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
void MassConstraintFitter::processEvent( LCEvent * evt ) { 
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
	evt->addCollection( recparcol, _outputParticleCollectionName.c_str() );

	return;
}
/*****************
finish and write out to rootfile
******************/
void MassConstraintFitter::end(){
	if(_fitAnalysis){
		rootFile->Write();
	}
  	return;
}
/*****************
locate the pfo collection with specified name
populated the global pfo vectors with particles from that collection for this event
******************/
bool MassConstraintFitter::FindPFOs( LCEvent* evt ) {

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
			std::cout<< "Located Pfo Collection "<< *itname<< " with "<< collection->getNumberofElements() << " elements " <<std::endl;
			collectionFound = true;

 			//add the collection elements to the global vector
      			for(unsigned int i=0; i<collection->getNumberofElements(); i++){
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
bool MassConstraintFitter::FindTracks( LCEvent* evt ) {

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
      			std::cout<< "Located Track Collection "<< *itname<< " with "<< collection->getNumberofElements() << " elements " <<std::endl;
     		
			collectionFound = true;
		
			//add the collection elements to the global vector
      			for(unsigned int i=0; i<collection->getNumberofElements(); i++){
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
bool MassConstraintFitter::FindMCParticles( LCEvent* evt ){
   
	bool collectionFound = false;

  	// clear old global MCParticle vector
  	_mcpartvec.clear();
  	typedef const std::vector<std::string> StringVec ;
  	StringVec* strVec = evt->getCollectionNames() ;

	//iterate over collections, find the matching name
  	for(StringVec::const_iterator itname=strVec->begin(); itname!=strVec->end(); itname++){    
    
		//if found print name and number of elements 
		if(*itname==_mcParticleCollectionName){
      			LCCollection* collection = evt->getCollection(*itname);
     			std::cout<< "Located MC Collection "<< *itname<< " with "<< collection->getNumberofElements() << " elements " <<std::endl;
      			collectionFound = true;
      
			//add the collection elements to the global vector
			for(unsigned int i=0;i<nelem;i++){
				MCParticle* mcPart = dynamic_cast<MCParticle*>(collection->getElementAt(i));
				_mcpartvec.push_back(mcPart);

       
      			}
    		}
  	}

  	if(!collectionFound){
		std::cout<<"MC Collection "<< _mcParticleCollectionName << "not found"<<std::endl;
	}
  
  	return collectionFound;
}

void MassConstraintFitter::FindMassConstraintCandidates(LCCollectionVec * recparcol) {
	//print global tree

}


