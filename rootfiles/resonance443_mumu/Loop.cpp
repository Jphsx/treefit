#include "OutputClass.cpp"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TFile.h"

void Loop(){
	TTree* tree=0;
	TFile *f = TFile::Open("MassConstraint.root");
 	TFile *f2 = new TFile("test.root","RECREATE");
	f->GetObject("resonance443_node0",tree);
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
  	TTreeReader myvecReader6_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar6_0(myvecReader6_0,"recoParentParams");
  	while(myvecReader6_0.Next()){
		oc->hrecoParentParams_param_0->Fill(myvecvar6_0->at(0));
	}
  	TTreeReader myvecReader6_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar6_1(myvecReader6_1,"recoParentParams");
  	while(myvecReader6_1.Next()){
		oc->hrecoParentParams_param_1->Fill(myvecvar6_1->at(1));
	}
  	TTreeReader myvecReader6_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar6_2(myvecReader6_2,"recoParentParams");
  	while(myvecReader6_2.Next()){
		oc->hrecoParentParams_param_2->Fill(myvecvar6_2->at(2));
	}
  	TTreeReader myvecReader6_3(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar6_3(myvecReader6_3,"recoParentParams");
  	while(myvecReader6_3.Next()){
		oc->hrecoParentParams_param_3->Fill(myvecvar6_3->at(3));
	}
  	TTreeReader myvecReader8_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_0(myvecReader8_0,"fitParentParams");
  	while(myvecReader8_0.Next()){
		oc->hfitParentParams_param_0->Fill(myvecvar8_0->at(0));
	}
  	TTreeReader myvecReader8_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_1(myvecReader8_1,"fitParentParams");
  	while(myvecReader8_1.Next()){
		oc->hfitParentParams_param_1->Fill(myvecvar8_1->at(1));
	}
  	TTreeReader myvecReader8_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_2(myvecReader8_2,"fitParentParams");
  	while(myvecReader8_2.Next()){
		oc->hfitParentParams_param_2->Fill(myvecvar8_2->at(2));
	}
  	TTreeReader myvecReader8_3(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_3(myvecReader8_3,"fitParentParams");
  	while(myvecReader8_3.Next()){
		oc->hfitParentParams_param_3->Fill(myvecvar8_3->at(3));
	}
  	TTreeReader myvecReader9_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar9_0(myvecReader9_0,"fitParentErrors");
  	while(myvecReader9_0.Next()){
		oc->hfitParentErrors_param_0->Fill(myvecvar9_0->at(0));
	}
  	TTreeReader myvecReader9_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar9_1(myvecReader9_1,"fitParentErrors");
  	while(myvecReader9_1.Next()){
		oc->hfitParentErrors_param_1->Fill(myvecvar9_1->at(1));
	}
  	TTreeReader myvecReader9_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar9_2(myvecReader9_2,"fitParentErrors");
  	while(myvecReader9_2.Next()){
		oc->hfitParentErrors_param_2->Fill(myvecvar9_2->at(2));
	}
  	TTreeReader myvecReader9_3(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar9_3(myvecReader9_3,"fitParentErrors");
  	while(myvecReader9_3.Next()){
		oc->hfitParentErrors_param_3->Fill(myvecvar9_3->at(3));
	}
	f2->Write();
}