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
		TH1D* hpdgs_param_2;
		TH1D* hrecoParentParams_param_0;
		TH1D* hrecoParentParams_param_1;
		TH1D* hrecoParentParams_param_2;
		TH1D* hrecoParentParams_param_3;
		TH1D* hfitParentParams_param_0;
		TH1D* hfitParentParams_param_1;
		TH1D* hfitParentParams_param_2;
		TH1D* hfitParentParams_param_3;
		OutputClass();
};
OutputClass::OutputClass(){
	hRecoEnergy = new TH1D("hRecoEnergy", "default description", 100, 10, 30);
	hFitEnergy = new TH1D("hFitEnergy", "default description", 100, 10, 30);
	hRecoMass = new TH1D("hRecoMass", "default description", 100, 0.25, 0.85);
	hFitProbability = new TH1D("hFitProbability", "default description", 100, 0, 1);
	hChisq = new TH1D("hChisq", "default description",100, 0, 20);
		hpdgs_param_0 = new TH1D("hpdgs_param_0","default description",84, -211, 211);
		hpdgs_param_1 = new TH1D("hpdgs_param_1","default description",84, -211, 211);
		hpdgs_param_2 = new TH1D("hpdgs_param_2","default description",84, -211, 211);
		hrecoParentParams_param_0 = new TH1D("hrecoParentParams_param_0","default description",84, -21.641, 21.2562);
		hrecoParentParams_param_1 = new TH1D("hrecoParentParams_param_1","default description",84, -21.1166, 21.289);
		hrecoParentParams_param_2 = new TH1D("hrecoParentParams_param_2","default description",84, -21.2454, 21.6354);
		hrecoParentParams_param_3 = new TH1D("hrecoParentParams_param_3","default description",84, 0.618566, 24.2318);
		hfitParentParams_param_0 = new TH1D("hfitParentParams_param_0","default description",84, -27.9297, 29.6517);
		hfitParentParams_param_1 = new TH1D("hfitParentParams_param_1","default description",84, -33.5383, 31.5549);
		hfitParentParams_param_2 = new TH1D("hfitParentParams_param_2","default description",84, -34.8836, 33.4752);
		hfitParentParams_param_3 = new TH1D("hfitParentParams_param_3","default description",84, 0.602035, 41.2006);
}
