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


	//cut input parameters
	registerProcessorParameter( "FitProbabilityCut" , 
			      	  "Minimum fit probability"  ,
			          _fitProbabilityCut,
			          (double)0.001);

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
		//this->FindMassConstraintCandidates(recparcol);
		this->FindMassConstraintCandidates(NULL);
	}

	std::cout<<"seg at the collection??"<<std::endl;
	//add collection to event
	evt->addCollection( recparcol,  _outputParticleCollectionName.c_str() );
	std::cout<<"seg after the collection??"<<std::endl;	
	
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
		//printing fit for testing
		for(int i=0; i<fit.size(); i++){
		std::cout<< i <<"     ";
		for(int k=0; k<fit.at(i).size(); k++){
			std::cout<<fit.at(i).at(k)<<" ";
		}
		std::cout<<std::endl;	
	}
		
		
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
				FO_vec.at(recoindex) = new LeptonFitObject(
				TFit->recoparts.at(recoindex)->track, 
				TFit->recoparts.at(recoindex)->Bfield, 
				TFit->recoparts.at(recoindex)->part->getMass()); 
				//testing track particle fit object
			/*	FO_vec.at(recoindex) = new TrackParticleFitObject(  //tpfo does not preserve the mass constraint
				TFit->recoparts.at(recoindex)->track,
				TFit->recoparts.at(recoindex)->part->getMass());
				TrackParticleFitObject* tpfo = (TrackParticleFitObject*) FO_vec.at(recoindex);
				tpfo->setBfield(TFit->recoparts.at(recoindex)->Bfield);
			*/
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
				std::cout<<" THE NODE MASS "<<node->mass<<std::endl;
				MassConstraint* mc = new MassConstraint(node->mass);
				//MassConstraint mc(double(node->mass));

				//get the FOs by iterating over j
				std::vector<ParticleFitObject*>* mcFitObjects = new vector<ParticleFitObject*>();
				
				//iterating over the combo in fit i
				for(int j=0; j<fit.at(i).size(); j++){
					//add to the array of FOs
					//we have to use an array because ParticleConstraint  is weird
					mcFitObjects->push_back(FO_vec.at( fit.at(i).at(j) ));
//trying to add each guy individually and
//make sure to do a cast for each type
					/*if(TFit->recoparts.at( fit.at(i).at(j) )->isTrack){
						//track push back casted tpfo
					//	mc->addToFOList(*(TrackParticleFitObject*)FO_vec.at( fit.at(i).at(j) ));
						mc->addToFOList(*(LeptonFitObject*)FO_vec.at(fit.at(i).at(j) ));
					}
					else{
						//not a track add jfo
					}	mc->addToFOList(*(JetFitObject*)FO_vec.at( fit.at(i).at(j) ));
					*/
				}//end j
				//add FOs to constraint
				mc->setFOList( mcFitObjects );
				//instead of using a mcvector try just immediately pushing onto the fitter
				fitter->addConstraint(mc);


			}//end if
		}//end i
		
		
		
		//do the fit
		fitter->fit();
		//save the FOs globally so we can easily
		//access/print the fitted particles
		//check fit here directly
		std::vector<TLorentzVector> fittlv{};
		TLorentzVector temp;
		for(int i=0; i<FO_vec.size(); i++){
			if(FO_vec.at(i) != NULL){
			temp.SetPxPyPzE(FO_vec.at(i)->getPx(),FO_vec.at(i)->getPy(),FO_vec.at(i)->getPz(),FO_vec.at(i)->getE());
			fittlv.push_back(temp);
			}
		}
		TLorentzVector parent;
		for(int i=0; i<fittlv.size(); i++){
			parent += fittlv.at(i);
		}
		std::cout<<"THE PARENT "<<parent.E()<<" "<<parent.M()<<std::endl;;
		//Particle::printTLorentzVector(parent);
		

		 FitObjects = FO_vec;
		std::cout<<"Fit Probability: "<<fitter->getProbability()<<std::endl;
		
		return fitter;
}
//TODO this function
//void treeFitter::createFitPartsfromFitObjects(FOVEC){}

//recursive function called by createLCOutput, creates the tree of recoparts/tracks to be put onto output LCIO
ReconstructedParticleImpl* treeFitter::createOutputParticle(Node* root, double fitProb, std::vector<std::vector<int> > fit){
	
		if(root->isLeaf) return NULL;// this should never happen
		//guarante debug print
		std::cout<<std::endl;
		std::cout<<std::endl;

		//create a reconstructed particle for non leaf node
		ReconstructedParticleImpl* p = new ReconstructedParticleImpl();
		ParticleIDImpl* newPDG = new ParticleIDImpl();
			
		newPDG->setPDG(root->pdg);
		newPDG->setLikelihood(1.0);

		//look up constituent particles, add together to get this particle
		//while we are at it, add up the total charge
		//this nodeID is index of fit combination
		TLorentzVector parentParticle;
		float charge=0.0;
		std::cout<<"SEG1"<<std::endl;
		std::cout<<"root node id "<<root->nodeId<<std::endl;
		//try printing fit
		for(int i=0; i<fit.size(); i++){
			for(int j=0; j<fit.at(i).size(); j++){
				std::cout<<"node "<<i<<" size "<< fit.at(i).size() << "element j" << fit.at(i).at(j)<<std::endl;
			}
		}
		std::cout<<"FIT PRINTED"<<std::endl;
		for(int i=0; i<fit.at(root->nodeId).size(); i++){
			//each element in the array at this fit location is an index of reco/FO/fit particle
			parentParticle += TFit->fitparts.at(fit.at(root->nodeId).at(i))->v;
			charge += TFit->recoparts.at(fit.at(root->nodeId).at(i))->part->getCharge();
		}
				std::cout<<"SEG1F"<<std::endl;
		//set px,py,pz
		float* mom = new float[3];
		mom[0] = parentParticle.Px();
		mom[1] = parentParticle.Py();
 		mom[2] = parentParticle.Pz();	
		
		p->setMomentum(mom);
		p->setEnergy(parentParticle.E());
		//give the reco part an E,theta,phi cov matrix
		//we need to construct the lower diagonal manually
		//TODO constuct matrix for every single particle resonance
		float* cov = new float[10];
		int index = 0;
	/*	for(int i=0; i<=2; i++){
			for(int j=0; j<=i; j++){
				cov[index]=jfo->getCov(i,j);
				index++;	
			}
		}*/
		for(int i=0; i<10; i++){
			cov[i] = 0.0;
		}
		p->setCovMatrix(cov);
		p->setMass(root->mass);
		p->setCharge(charge);
		p->addParticleID(newPDG);
		p->setParticleIDUsed(newPDG);
		p->setType(root->pdg);
		p->setGoodnessOfPID(fitProb);

		//go through children, identify the leaves and add them
		//if we have a nonleaf child make a reconstructed particle for that resonance
		std::vector<int> parentSet = fit.at(root->nodeId);
		std::vector<int> childSet{};
				std::cout<<"SEG2"<<std::endl;
		for(int i=0; i<root->children.size(); i++){
				//subtract all non leaf  children sets from parent, the remaining is the leaves to add 
			if(!root->children.at(i)->isLeaf){
				childSet = fit.at(root->children.at(i)->nodeId);
				parentSet = Combinatorics::subtractSets(parentSet,childSet);
			}
		}
				std::cout<<"SEG2F"<<std::endl;
		//the remaining (if any) particles on parentSet are leaves that can be added immediately
				std::cout<<"SEG3"<<std::endl;
		for(int i=0; i<parentSet.size(); i++){
			if(TFit->recoparts.at(parentSet.at(i))->isTrack){
				//this is a track dont add reconstructedParticle*
				p->addTrack( TFit->fitparts.at(parentSet.at(i))->track);
			}
			else{
				//this is a neutral add the correct object
				p->addParticle( TFit->fitparts.at(parentSet.at(i))->part);
			}
		}
				std::cout<<"SEG3F"<<std::endl;
		//now deal with the non leaves, iterate through children again and create the other particles
				std::cout<<"SEG4"<<std::endl;
		for(int i=0; i<root->children.size(); i++){
			if(!root->children.at(i)->isLeaf){
				p->addParticle( createOutputParticle(root->children.at(i),fitProb,fit) );
			}
		}
				std::cout<<"SEG4F"<<std::endl;
		std::cout<<"created this parent particle :: "<<std::endl;
		Particle::printReconstructedParticle(p);
		return p;

}
void treeFitter::createLCOutputParticles(LCCollectionVec* recparcol, std::vector<std::vector<int> > fit, double fitProb){
		//create the LCIO reconstructed particle tree
		//save the particles to an output collection		
		//calreccol->setSubset(true);	   is this needed?
		recparcol->addElement( createOutputParticle(TFit->ParticleTree->Root, fitProb, fit  ));
		  
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
	//prep fitparticles to match the size of recoparts
	//std::vector<Particle*>  fitparts(TFit->recoparts.size());
	//TFit->fitparts = fitparts;

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
	double bestfitprob=-1.0;
	std::vector<std::vector<int> > bestfit{};

	//do the fits
	OPALFitterGSL*  fitter; 
	
	for(int j = 0; j<TFit->fitTable.at(0).size(); j++){
		//extract each fit onto a 2d fit vector
		std::vector<std::vector<int> > fit(TFit->fitTable.size());
		//this is a single fit from the fit table
		for(int i=0; i<TFit->fitTable.size(); i++){
			//if it has particles to fit then proceed
			if(TFit->fitTable.at(i).size() != 0){
				fit.at(i) = TFit->fitTable.at(i).at(j);
			}
		}
		
		fitter = fitParticles(fit);
		
		//check and see if this is the best fit and exceeds the minimal probability cut
		if(fitter->getProbability() > bestfitprob && fitter->getProbability() > _fitProbabilityCut){
			bestfit = fit;
			bestfitprob = fitter->getProbability();
 		}

		//make the fit particles from the FOs
		//make a temp vec to give to tfit, to make sure fitparts in synchronized in size
		std::vector<Particle*> fit_vec(TFit->recoparts.size());
		TFit->fitparts = fit_vec;
		for(int k=0; k<FitObjects.size(); k++){
			if(FitObjects.at(k)==NULL){
				continue;
			} 
			
			if(TFit->recoparts.at(k)->isTrack){
			//	TFit->fitparts.at(k) = new Particle(NULL, (TrackParticleFitObject*) FitObjects.at(k), TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass()) ;
				TFit->fitparts.at(k) = new Particle(NULL, (LeptonFitObject*) FitObjects.at(k), TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass(), TFit->recoparts.at(k)->track->getD0() ,TFit->recoparts.at(k)->track->getZ0(), TFit->recoparts.at(k)->Bfield);
				
			}
			if(!TFit->recoparts.at(k)->isTrack){
				TFit->fitparts.at(k) = new Particle( (JetFitObject*) FitObjects.at(k), NULL, TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass()) ;
			}
			
						
		}
		if(fitter->getProbability() > _fitProbabilityCut){
			//if we pass, save this particle hypothesis and fit to the outputcollection
			createLCOutputParticles(recparcol, fit, fitter->getProbability());
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
		FitObjects.clear();//this vector will not change capacity when cleared size will be 0 though
	}//fitTable iteration
	
	//redo the best fit, and send the particles to the TTrees in the Rootfiles
	std::cout<<"is fault here"<<std::endl;
	//make sure there was at least 1 fit

	if(bestfitprob != -1.0){
		fitter = fitParticles(bestfit);
	std::cout<<"nothere "<<std::endl;
		std::cout<<FitObjects.size()<<" FO size "<<std::endl;
		std::cout<<TFit->fitparts.size()<<" fitparts size"<<std::endl;
	std::cout<<"edit"<<std::endl;
	//remake fitparticles
		std::vector<Particle*> fit_vec(TFit->recoparts.size());
		TFit->fitparts = fit_vec;
		for(int k=0; k<FitObjects.size(); k++){
			if(FitObjects.at(k)==NULL){
				continue;
			} 
			
			if(TFit->recoparts.at(k)->isTrack){
				//TFit->fitparts.at(k) = new Particle(NULL, (TrackParticleFitObject*) FitObjects.at(k), TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass()) ;
				TFit->fitparts.at(k) = new Particle(NULL, (LeptonFitObject*) FitObjects.at(k), TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass(), TFit->recoparts.at(k)->track->getD0() ,TFit->recoparts.at(k)->track->getZ0(), TFit->recoparts.at(k)->Bfield);
				
			}
			if(!TFit->recoparts.at(k)->isTrack){
				TFit->fitparts.at(k) = new Particle( (JetFitObject*) FitObjects.at(k), NULL, TFit->recoparts.at(k)->recopdg, TFit->recoparts.at(k)->part->getMass(), -1, -1, TFit->recoparts.at(k)->Bfield) ;
			}
			
						
		}

	std::cout<<"print bestfit"<<std::endl;
		for(int i=0; i<bestfit.size(); i++){
			std::cout<< i <<"     ";
		for(int k=0; k<bestfit.at(i).size(); k++){
			std::cout<<bestfit.at(i).at(k)<<" ";
		}
		std::cout<<std::endl;	
	}
	//iterate through the fit, create the needed fit particles
	//put the particles on new vectors, now indexed by a pdg vector
	//re-vectoring will get rid of gaps in used reco/fit particle vector
	std::vector<Particle*> recop{};
	std::vector<Particle*> fitp{};
	//populate reco and fit at the same time
	//iterate over each nodeId in fit
	int index=0;
	//use index for iterating over treefactory vector, since we cant have gaps in its array
	for(unsigned int i=0; i<bestfit.size(); i++){
			std::cout<<"where is fault????"<<std::endl;
			//nodeId should by construction match fit index with ttrees index
			//iterate over the fit particles
			std::cout<<"best fit size "<< bestfit.size() << std::endl;
			for(unsigned int k=0; k<bestfit.at(i).size(); k++){
				std::cout<<"is it in here??"<<std::endl;
				std::cout<<"size at i "<<i<<" "<<bestfit.at(i).size() <<std::endl;
				std::cout<<"getting recop"<<std::endl;
				recop.push_back(TFit->recoparts.at( bestfit.at(i).at(k) ));
				std::cout<<"getting fitp"<<std::endl;
				fitp.push_back(TFit->fitparts.at( bestfit.at(i).at(k) ));
				std::cout<<"got p's "<<std::endl;
				
			}
			std::cout<<"is the fault at trees"<<std::endl;
			ttrees.at(index)->addParticleSets(fitp,recop);
			std::cout<<"is the fault at trees"<<std::endl;
			ttrees.at(index)->addFitDetails(fitter->getProbability(), fitter->getChi2());
			std::cout<<"is the fault at trees"<<std::endl;
			ttrees.at(index)->TreeFillAndClear();
			std::cout<<"weve went past trees"<<std::endl;
			index++;	
		
		
	}
	recop.clear();
	fitp.clear();
	}//end bestfit   

	std::cout<<"seg at the end??"<<std::endl;

	//advance to next event
	evtNo++;
std::cout<<"seg at the end??"<<std::endl;
	_pfovec.clear();
	_trackvec.clear();
std::cout<<"seg at the end??"<<std::endl;
	//might need to run a destructor here first
	TFit->clearEvent();
	FitObjects.clear();//clear the bestfit
std::cout<<"seg at the end??"<<std::endl;
	return;
}


