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
		OutputClass();
};
OutputClass::OutputClass(){
	hRecoEnergy = new TH1D("hRecoEnergy", "default description", 87, 0.573775, 24.2318);
	hFitEnergy = new TH1D("hFitEnergy", "default description", 87, 0.600451, 31.9519);
	hRecoMass = new TH1D("hRecoMass", "default description", 87, 0.281587, 4.68418);
	hFitProbability = new TH1D("hFitProbability", "default description", 87, 0, 0.999976);
	hChisq = new TH1D("hChisq", "default description", 87, 9.00506e-10, 54974.1);
		hpdgs_param_0 = new TH1D("hpdgs_param_0","default description",87, -211, 211);
		hpdgs_param_1 = new TH1D("hpdgs_param_1","default description",87, -211, 211);
		hpdgs_param_2 = new TH1D("hpdgs_param_2","default description",87, -211, 211);
}