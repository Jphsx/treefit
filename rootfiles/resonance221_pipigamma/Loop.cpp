#include "OutputClass.cpp"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TFile.h"

void Loop(){
	TTree* tree=0;
	TFile *f = TFile::Open("MassConstraint.root");
 	TFile *f2 = new TFile("test.root","RECREATE");
	f->GetObject("resonance221_node0",tree);
  OutputClass* oc = new OutputClass();


	TTreeReader mydoubleReader0(tree);
	TTreeReaderValue<Double_t> mydoublevar0(mydoubleReader0,"RecoEnergy");
	while(mydoubleReader0.Next()){	
		oc->hRecoEnergy->Fill(*mydoublevar0);
	}
	TTreeReader mydoubleReader1(tree);
	TTreeReaderValue<Double_t> mydoublevar1(mydoubleReader1,"FitEnergy");
	while(mydoubleReader1.Next()){	
		oc->hFitEnergy->Fill(*mydoublevar1);
	}
	TTreeReader mydoubleReader2(tree);
	TTreeReaderValue<Double_t> mydoublevar2(mydoubleReader2,"RecoMass");
	while(mydoubleReader2.Next()){	
		oc->hRecoMass->Fill(*mydoublevar2);
	}
	TTreeReader mydoubleReader3(tree);
	TTreeReaderValue<Double_t> mydoublevar3(mydoubleReader3,"FitProbability");
	while(mydoubleReader3.Next()){	
		oc->hFitProbability->Fill(*mydoublevar3);
	}
	TTreeReader mydoubleReader4(tree);
	TTreeReaderValue<Double_t> mydoublevar4(mydoubleReader4,"Chisq");
	while(mydoubleReader4.Next()){	
		oc->hChisq->Fill(*mydoublevar4);
	}
  	TTreeReader myvecReader5_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar5_0(myvecReader5_0,"pdgs");
  	while(myvecReader5_0.Next()){
		oc->hpdgs_param_0->Fill(myvecvar5_0->at(0));
	}
  	TTreeReader myvecReader5_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar5_1(myvecReader5_1,"pdgs");
  	while(myvecReader5_1.Next()){
		oc->hpdgs_param_1->Fill(myvecvar5_1->at(1));
	}
  	TTreeReader myvecReader5_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar5_2(myvecReader5_2,"pdgs");
  	while(myvecReader5_2.Next()){
		oc->hpdgs_param_2->Fill(myvecvar5_2->at(2));
	}
	f2->Write();
}