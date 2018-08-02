#include "TH1D.h"

class OutputClass{ 
     public: 
		TH1D* hRecoEnergy;
		TH1D* hFitEnergy;
		TH1D* hRecoMass;
		TH1D* hFitProbability;
		TH1D* hChisq;
		TH1D* hpdgs_param_0;
		TH1D* hpdgs_param_1;
		TH1D* hrecoParentParams_param_0;
		TH1D* hrecoParentParams_param_1;
		TH1D* hrecoParentParams_param_2;
		TH1D* hrecoParentParams_param_3;
		TH1D* hfitParentParams_param_0;
		TH1D* hfitParentParams_param_1;
		TH1D* hfitParentParams_param_2;
		TH1D* hfitParentParams_param_3;
		TH1D* hfitParentErrors_param_0;
		TH1D* hfitParentErrors_param_1;
		TH1D* hfitParentErrors_param_2;
		TH1D* hfitParentErrors_param_3;
		OutputClass();
};
OutputClass::OutputClass(){
	hRecoEnergy = new TH1D("hRecoEnergy", "default description", 64, 18.0246, 20.8695);
	hFitEnergy = new TH1D("hFitEnergy", "default description", 64, 19.4456, 20.5892);
	hRecoMass = new TH1D("hRecoMass", "default description", 64, 2.74273, 3.25865);
	hFitProbability = new TH1D("hFitProbability", "default description", 64, 0, 0.999763);
	hChisq = new TH1D("hChisq", "default description", 64, 8.79381e-08, 3126.86);
		hpdgs_param_0 = new TH1D("hpdgs_param_0","default description",64, -13, 13);
		hpdgs_param_1 = new TH1D("hpdgs_param_1","default description",64, -13, 13);
		hrecoParentParams_param_0 = new TH1D("hrecoParentParams_param_0","default description",64, -19.7644, 19.7728);
		hrecoParentParams_param_1 = new TH1D("hrecoParentParams_param_1","default description",64, -19.764, 19.8545);
		hrecoParentParams_param_2 = new TH1D("hrecoParentParams_param_2","default description",64, -20.1051, 19.8431);
		hrecoParentParams_param_3 = new TH1D("hrecoParentParams_param_3","default description",64, 18.0246, 20.8695);
		hfitParentParams_param_0 = new TH1D("hfitParentParams_param_0","default description",64, -19.7475, 19.7723);
		hfitParentParams_param_1 = new TH1D("hfitParentParams_param_1","default description",64, -19.7727, 19.7431);
		hfitParentParams_param_2 = new TH1D("hfitParentParams_param_2","default description",64, -19.882, 19.8537);
		hfitParentParams_param_3 = new TH1D("hfitParentParams_param_3","default description",64, 19.4456, 20.5892);
		hfitParentErrors_param_0 = new TH1D("hfitParentErrors_param_0","default description",64, 1.93622e-06, 0.00230038);
		hfitParentErrors_param_1 = new TH1D("hfitParentErrors_param_1","default description",64, 1.96541e-06, 0.00366521);
		hfitParentErrors_param_2 = new TH1D("hfitParentErrors_param_2","default description",64, 0.000139965, 0.153152);
		hfitParentErrors_param_3 = new TH1D("hfitParentErrors_param_3","default description",64, 1.42975e-08, 0.00866324);
}