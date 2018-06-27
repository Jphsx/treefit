#include "EVENT/LCIO.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCParameters.h"
#include "EVENT/ReconstructedParticle.h"
#include "gear/BField.h"
#include "IMPL/ReconstructedParticleImpl.h"

#include "CLHEP/Vector/LorentzVector.h"
#include "CLHEP/Vector/ThreeVector.h"
//typedef CLHEP::HepLorentzVector LorentzVector ;
//typedef CLHEP::Hep3Vector Vector3D ;

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

	//build TTREE from non leaf nodes
	//iterate through the tree by retrieving each node, then create a treefactory at each nonleaf
	file = new TFile("rootFile.root","RECREATE");
	Node* nonleafnode;
	for(int i=0; i<= *(TFit->LASTNONLEAFID); i++){
		nonleafnode = Tree::getNode(TFit->ParticleTree->Root, i);
		if(!nonleafnode->isLeaf){
			//get naming details for this node
			ttrees.push_back(new TTreeFactory(i, nonleafnode->pdg, file));
		}
	}		
	

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
	file->Write();
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
		
		
		//this array should match reco parts in size, but only the final particles
		//used will be populated into the array
 		std::vector<ParticleFitObject*> FO_vec(TFit->recoparts.size());
		//use the first node it populate all the FOs
		//adding this index var for readability
		int recoindex=-1;
		//iterate over the tree roots combination (node 0)
		for(int j=0; j<fit.at(0).size(); j++){
			//make an object at index
			//fit.at(0).at(j) this is index on reco part
			recoindex = fit.at(0).at(j);
			if(TFit->recoparts.at(recoindex)->isTrack){
				//this is a track make LFO
				/*FO_vec.at(recoindex) = new LeptonFitObject(
				TFit->recoparts.at(recoindex)->track, 
				TFit->recoparts.at(recoindex)->Bfield, 
				TFit->recoparts.at(recoindex)->part->getMass()); */
				//testing track particle fit object
				FO_vec.at(recoindex) = new TrackParticleFitObject(
				TFit->recoparts.at(recoindex)->track,
				TFit->recoparts.at(recoindex)->part->getMass());
				TrackParticleFitObject* tpfo = (TrackParticleFitObject*) FO_vec.at(recoindex);
				tpfo->setBfield(TFit->recoparts.at(recoindex)->Bfield);
			}
			else{
				//this is not a track make JFO
				FO_vec.at(recoindex) = new JetFitObject(
				TFit->recoparts.at(recoindex)->localParams.at(0), 
 				TFit->recoparts.at(recoindex)->localParams.at(1), 
				TFit->recoparts.at(recoindex)->localParams.at(2),  
				TFit->recoparts.at(recoindex)->localErrors.at(0), 
				TFit->recoparts.at(recoindex)->localErrors.at(1), 
				TFit->recoparts.at(recoindex)->localErrors.at(2), 
				TFit->recoparts.at(recoindex)->part->getMass() );
				
			}
			//add all FOs to the fitter
			fitter->addFitObject( FO_vec.at(recoindex) );
			
		}
			
		
		
		//make mass constraint objects for each node in the tree
		//that has a specified mass constraint
		
		//iterate through the fit, get mass and combination
		//and add them to the corresponding constraint
		for(int i=0; i<fit.size(); i++){
			//find the node for the current fit
			Node* node = TFit->ParticleTree->getNode(TFit->ParticleTree->Root, i);
			
			
			if(node->mass != -1){
				//make a new constraint
				MassConstraint* mc = new MassConstraint(node->mass);

				//get the FOs by iterating over j
				std::vector<ParticleFitObject*>* mcFitObjects = new vector<ParticleFitObject*>();
				
				//iterating over the combo in fit i
				for(int j=0; j<fit.at(i).size(); j++){
					//add to the array of FOs
					//we have to use an array because ParticleConstraint  is weird
					mcFitObjects->push_back(FO_vec.at( fit.at(i).at(j) ));
				}//end j
				//add FOs to constraint
				mc->setFOList( mcFitObjects );
				//instead of using a mcvector try just immediately pushing onto the fitter
				fitter->addConstraint(mc);
			}//end if
		}//end i
		
		//save the FOs globally so we can easily
		//access/print the fitted particles
		 FitObjects = FO_vec;
		//do the fit
		fitter->fit();
		
		std::cout<<"Fit Probability: "<<fitter->getProbability()<<std::endl;
		
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
	std::cout<<"Reconstructed Particles "<<std::endl;
	TFit->printParticles(TFit->recoparts);
	std::cout<<std::endl;

	//Do the fits (print them)
	TFit->initTable();
	TFit->generatefitcombinations(TFit->ParticleTree->Root, TFit->recoIDs);
	//print the fit table	
	TFit->printTable();

	//use these variables to save the jth fit with the
	//highest fit probability
	double bestfitprob=0.0;
	std::vector<std::vector<int> > bestfit{};

	//do the fits
	OPALFitterGSL*  fitter; 
	//extract each fit onto a 2d fit vector
	//this is a single fit from the fit table
	for(int j = 0; j<TFit->fitTable.at(0).size(); j++){
		std::vector<std::vector<int> > fit(TFit->fitTable.size());
		for(int i=0; i<TFit->fitTable.size(); i++){
			//if it has particles to fit then proceed
			if(TFit->fitTable.at(i).size() != 0){
				fit.at(i) = TFit->fitTable.at(i).at(j);
			}
		}
		
		fitter = fitParticles(fit);
		//check and see if this is the best fit
		if(fitter->getProbability() > bestfitprob){
			bestfit = fit;
			bestfitprob = fitter->getProbability();
 		}

		//make the fit particles from the FOs
		for(int k=0; k<FitObjects.size(); k++){
			if(FitObjects.at(k)==NULL){
				continue;
			} 
			
			if(TFit->recoparts.at(k)->isTrack){
				TFit->addfitpart( new Particle(NULL, (TrackParticleFitObject*) FitObjects.at(k), TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass()) );
				
			}
			if(!TFit->recoparts.at(k)->isTrack){
				TFit->addfitpart( new Particle( (JetFitObject*) FitObjects.at(k), NULL, TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass()) );
			}
			
						
		}
		//print every fit
		std::cout<<"Fitted Particles in Fit "<< j <<std::endl;
		TFit->printParticles(TFit->fitparts);
		//TODO save each fit passing a certain probability
		//cut to output collection
		//For local plots save the best fit probability fit
		//this is the event we will plot
		
		//before we move on to the next set of combinations
		//clear the fit (fit is the fit combo from the fit table)
		fit.clear();
		//also clear fitparts after each fit and FOs
		TFit->fitparts.clear();
		FitObjects.clear();//this vector will not change capacity when cleared
	}//fitTable iteration
	
	//redo the best fit, and send the particles to the TTrees in the Rootfiles
	fitter = fitParticles(bestfit);
	//iterate through the fit, create the needed fit particles
	//put the particles on new vectors, now indexed by a pdg vector
	//re-vectoring will get rid of gaps in used reco/fit particle vector
	std::vector<Particle*> recop,fitp;
	//populate reco and fit at the same time
	//iterate over each nodeId in fit
	int index;
	//use index for iterating over treefactory vector, since we cant have gaps in its array
	for(unsigned int i=0; i<bestfit.size(); i++){
		//if this node is non leaf, we can populate its tree
		if(bestfit.at(i) == NULL) continue;
			//iterate over the fit particles
			for(unsigned int k=0; k<bestfit.at(i).size(); k++){
				reco.push_back(TFit->recoparts.at( bestfit.at(i).at(j) ));
				fit.push_back(TFit->fitparts.at( bestfit.at(i).at(j) ));
			}
			ttrees.at(index)->addParticleSets(fitp,recop);
			ttrees.at(index)->addFitDetails(fitter->getProbability(), fitter->getChi2());
			ttrees.at(index)->TreeFillAndClear();
			index++;	
		
	}
	



	//advance to next event
	evtNo++;
	_pfovec.clear();
	_trackvec.clear();
	//might need to run a destructor here first
	TFit->clearEvent();
	FitObjects.clear();//clear the bestfit
	return;
}


