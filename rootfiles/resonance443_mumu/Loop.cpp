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
  	TTreeReaderValue<std::vector<double> > myvecvar6_0(myvecReader6_0,"vertex");
  	while(myvecReader6_0.Next()){
		oc->hvertex_param_0->Fill(myvecvar6_0->at(0));
	}
  	TTreeReader myvecReader6_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar6_1(myvecReader6_1,"vertex");
  	while(myvecReader6_1.Next()){
		oc->hvertex_param_1->Fill(myvecvar6_1->at(1));
	}
  	TTreeReader myvecReader6_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar6_2(myvecReader6_2,"vertex");
  	while(myvecReader6_2.Next()){
		oc->hvertex_param_2->Fill(myvecvar6_2->at(2));
	}
  	TTreeReader myvecReader7_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar7_0(myvecReader7_0,"vertexErrors");
  	while(myvecReader7_0.Next()){
		oc->hvertexErrors_param_0->Fill(myvecvar7_0->at(0));
	}
  	TTreeReader myvecReader7_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar7_1(myvecReader7_1,"vertexErrors");
  	while(myvecReader7_1.Next()){
		oc->hvertexErrors_param_1->Fill(myvecvar7_1->at(1));
	}
  	TTreeReader myvecReader7_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar7_2(myvecReader7_2,"vertexErrors");
  	while(myvecReader7_2.Next()){
		oc->hvertexErrors_param_2->Fill(myvecvar7_2->at(2));
	}
  	TTreeReader myvecReader8_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_0(myvecReader8_0,"recoParentParams");
  	while(myvecReader8_0.Next()){
		oc->hrecoParentParams_param_0->Fill(myvecvar8_0->at(0));
	}
  	TTreeReader myvecReader8_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_1(myvecReader8_1,"recoParentParams");
  	while(myvecReader8_1.Next()){
		oc->hrecoParentParams_param_1->Fill(myvecvar8_1->at(1));
	}
  	TTreeReader myvecReader8_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_2(myvecReader8_2,"recoParentParams");
  	while(myvecReader8_2.Next()){
		oc->hrecoParentParams_param_2->Fill(myvecvar8_2->at(2));
	}
  	TTreeReader myvecReader8_3(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar8_3(myvecReader8_3,"recoParentParams");
  	while(myvecReader8_3.Next()){
		oc->hrecoParentParams_param_3->Fill(myvecvar8_3->at(3));
	}
  	TTreeReader myvecReader10_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar10_0(myvecReader10_0,"fitParentParams");
  	while(myvecReader10_0.Next()){
		oc->hfitParentParams_param_0->Fill(myvecvar10_0->at(0));
	}
  	TTreeReader myvecReader10_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar10_1(myvecReader10_1,"fitParentParams");
  	while(myvecReader10_1.Next()){
		oc->hfitParentParams_param_1->Fill(myvecvar10_1->at(1));
	}
  	TTreeReader myvecReader10_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar10_2(myvecReader10_2,"fitParentParams");
  	while(myvecReader10_2.Next()){
		oc->hfitParentParams_param_2->Fill(myvecvar10_2->at(2));
	}
  	TTreeReader myvecReader10_3(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar10_3(myvecReader10_3,"fitParentParams");
  	while(myvecReader10_3.Next()){
		oc->hfitParentParams_param_3->Fill(myvecvar10_3->at(3));
	}
  	TTreeReader myvecReader11_0(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar11_0(myvecReader11_0,"fitParentErrors");
  	while(myvecReader11_0.Next()){
		oc->hfitParentErrors_param_0->Fill(myvecvar11_0->at(0));
	}
  	TTreeReader myvecReader11_1(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar11_1(myvecReader11_1,"fitParentErrors");
  	while(myvecReader11_1.Next()){
		oc->hfitParentErrors_param_1->Fill(myvecvar11_1->at(1));
	}
  	TTreeReader myvecReader11_2(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar11_2(myvecReader11_2,"fitParentErrors");
  	while(myvecReader11_2.Next()){
		oc->hfitParentErrors_param_2->Fill(myvecvar11_2->at(2));
	}
  	TTreeReader myvecReader11_3(tree);
  	TTreeReaderValue<std::vector<double> > myvecvar11_3(myvecReader11_3,"fitParentErrors");
  	while(myvecReader11_3.Next()){
		oc->hfitParentErrors_param_3->Fill(myvecvar11_3->at(3));
	}


	/*pull distribution code chunk:
	*/
	TTreeReader fitParentRead(tree);
  	TTreeReaderValue<std::vector<double> > fitParentVec(fitParentRead,"fitParentParams");
	TTreeReader recoParentRead(tree);
	TTreeReaderValue<std::vector<double> > recoParentVec(recoParentRead,"recoParentParams");
	TTreeReader fitParentErrorRead(tree);
	TTreeReaderValue<std::vector<double> > fitParentErrorVec(fitParentErrorRead,"fitParentErrors");

  	while(fitParentRead.Next() && recoParentRead.Next() && fitParentErrorRead.Next() ){

		oc->hParentPull_px->Fill((recoParentVec->at(0)-fitParentVec->at(0))/fitParentErrorVec->at(0));
		oc->hParentPull_py->Fill((recoParentVec->at(1)-fitParentVec->at(1))/fitParentErrorVec->at(1));
		oc->hParentPull_pz->Fill((recoParentVec->at(2)-fitParentVec->at(2))/fitParentErrorVec->at(2));
		oc->hParentPull_E->Fill((recoParentVec->at(3)-fitParentVec->at(3))/fitParentErrorVec->at(3));
	}

	/*
	Vertex position pull
	*/
	TTreeReader VertexRead(tree);
  	TTreeReaderValue<std::vector<double> > VertexVec( VertexRead,"vertex");
	TTreeReader VertexErrorRead(tree);
	TTreeReaderValue<std::vector<double> > VertexErrorVec( VertexErrorRead,"vertexErrors");

	while(VertexRead.Next() && VertexErrorRead.Next()){
		oc->hVertexPull_x->Fill( VertexVec->at(0)/VertexErrorVec->at(0) );
		oc->hVertexPull_y->Fill( VertexVec->at(1)/ VertexErrorVec->at(1) );
		oc->hVertexPull_z->Fill( VertexVec->at(2)/ VertexErrorVec->at(2) );
	}
	
	/*end chunk*/

	/* track pulls split by charge */
	
  	TTreeReader pdgreader(tree);
  	TTreeReaderValue<std::vector<double> > pdgvec(pdgreader,"pdgs");
	
  	TTreeReader recoparamreader(tree);
  	TTreeReaderValue<std::vector< std::vector<double> > > recoparamvec(recoparamreader,"recoLocalParams");
	
  	TTreeReader fitparamreader(tree);
  	TTreeReaderValue<std::vector< std::vector<double> > > fitparamvec(fitparamreader,"fitLocalParams");
	
  	TTreeReader fiterrorreader(tree);
  	TTreeReaderValue<std::vector< std::vector<double> > > fiterrorvec(fiterrorreader,"fitLocalErrors");

	while(pdgreader.Next() && recoparamreader.Next() && fitparamreader.Next() && fiterrorreader.Next() ){
		//indexed by pdg
		//go and get the +13s then -13s
		for(int i=0; i<pdgvec->size(); i++){
			if( pdgvec->at(i) > 0 ){ 
			// this is + track
			 	oc->hPd0pull->Fill( (recoparamvec->at(i).at(0) - fitparamvec->at(i).at(0) )/ fiterrorvec->at(i).at(0) );
				oc->hPphipull->Fill( (recoparamvec->at(i).at(1) - fitparamvec->at(i).at(1) )/ fiterrorvec->at(i).at(1) );
				oc->hPomepull->Fill( (recoparamvec->at(i).at(2) - fitparamvec->at(i).at(2) )/ fiterrorvec->at(i).at(2) );
				oc->hPz0pull->Fill( (recoparamvec->at(i).at(3) - fitparamvec->at(i).at(3) )/ fiterrorvec->at(i).at(3) );
				oc->hPtlpull->Fill( (recoparamvec->at(i).at(4) - fitparamvec->at(i).at(4) )/ fiterrorvec->at(i).at(4) );
			}
			else{
			//this is - track
				oc->hMd0pull->Fill( (recoparamvec->at(i).at(0) - fitparamvec->at(i).at(0) )/ fiterrorvec->at(i).at(0) );
				oc->hMphipull->Fill( (recoparamvec->at(i).at(1) - fitparamvec->at(i).at(1) )/ fiterrorvec->at(i).at(1) );
				oc->hMomepull->Fill( (recoparamvec->at(i).at(2) - fitparamvec->at(i).at(2) )/ fiterrorvec->at(i).at(2) );
				oc->hMz0pull->Fill( (recoparamvec->at(i).at(3) - fitparamvec->at(i).at(3) )/ fiterrorvec->at(i).at(3) );
				oc->hMtlpull->Fill( (recoparamvec->at(i).at(4) - fitparamvec->at(i).at(4) )/ fiterrorvec->at(i).at(4) );
			}
		
		}	


	}
	
	/*end chunk */

	f2->Write();
}
