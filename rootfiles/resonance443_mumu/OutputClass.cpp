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
		TH1D* hvertex_param_0;
		TH1D* hvertex_param_1;
		TH1D* hvertex_param_2;
		TH1D* hvertexErrors_param_0;
		TH1D* hvertexErrors_param_1;
		TH1D* hvertexErrors_param_2;
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

		/*Parent pull distribution chunk*/
		TH1D* hParentPull_px;
		TH1D* hParentPull_py;
		TH1D* hParentPull_pz;
		TH1D* hParentPull_E;
		
		/* end extra stuff */
		/*vertex*/
		TH1D* hVertexPull_x;
		TH1D* hVertexPull_y;
		TH1D* hVertexPull_z;
		/* end vertex */

		/* track pulls */
		TH1D* hPd0pull;
		TH1D* hPphipull;
		TH1D* hPomepull;
		TH1D* hPz0pull; 
		TH1D* hPtlpull;
		TH1D* hMd0pull;
		TH1D* hMphipull;
		TH1D* hMomepull;
		TH1D* hMz0pull; 
		TH1D* hMtlpull;
		/*end track pulls */

		OutputClass();
};
OutputClass::OutputClass(){
	hRecoEnergy = new TH1D("hRecoEnergy", "default description", 50, 19.6, 20.4);
	hFitEnergy = new TH1D("hFitEnergy", "default description", 50, 19.6, 20.4);
	hRecoMass = new TH1D("hRecoMass", "default description", 50, 3.04, 3.16);
	hFitProbability = new TH1D("hFitProbability", "default description", 50, 0, 1);
	hChisq = new TH1D("hChisq", "default description", 50, 8.16837e-09, 15.0);
		hpdgs_param_0 = new TH1D("hpdgs_param_0","default description",50, -13, 13);
		hpdgs_param_1 = new TH1D("hpdgs_param_1","default description",50, -13, 13);
		hvertex_param_0 = new TH1D("hvertex_param_0","default description",50, -0.2, 0.2);
		hvertex_param_1 = new TH1D("hvertex_param_1","default description",50, -0.2, 0.2);
		hvertex_param_2 = new TH1D("hvertex_param_2","default description",50, -0.2, 0.2);
		hvertexErrors_param_0 = new TH1D("hvertexErrors_param_0","default description",50, 0.0009877, 4.63995);
		hvertexErrors_param_1 = new TH1D("hvertexErrors_param_1","default description",50,0.000968902, 4.76966);
		hvertexErrors_param_2 = new TH1D("hvertexErrors_param_2","default description",50, 0.00202624, 14.2397);
		hrecoParentParams_param_0 = new TH1D("hrecoParentParams_param_0","default description",19, -19.7638, 19.7028);
		hrecoParentParams_param_1 = new TH1D("hrecoParentParams_param_1","default description",19, -19.7938, 19.7507);
		hrecoParentParams_param_2 = new TH1D("hrecoParentParams_param_2","default description",19, -19.7073, 19.6615);
		hrecoParentParams_param_3 = new TH1D("hrecoParentParams_param_3","default description",19, 15.1534, 20.3613);
		hfitParentParams_param_0 = new TH1D("hfitParentParams_param_0","default description",19, -19.7642, 19.7027);
		hfitParentParams_param_1 = new TH1D("hfitParentParams_param_1","default description",19, -19.7939, 19.7507);
		hfitParentParams_param_2 = new TH1D("hfitParentParams_param_2","default description",19, -19.6982, 19.6482);
		hfitParentParams_param_3 = new TH1D("hfitParentParams_param_3","default description",19, 15.1534, 20.3387);
		hfitParentErrors_param_0 = new TH1D("hfitParentErrors_param_0","default description",19, 1.4232e-06, 0.00461959);
		hfitParentErrors_param_1 = new TH1D("hfitParentErrors_param_1","default description",19, 1.75409e-06, 0.00131926);
		hfitParentErrors_param_2 = new TH1D("hfitParentErrors_param_2","default description",19, 1.60245e-06, 0.119044);
		hfitParentErrors_param_3 = new TH1D("hfitParentErrors_param_3","default description",19, 0.00012489, 0.120484);

		/*The Pull dists*/
		hParentPull_px = new TH1D("hParentPull_px","px",70,-5,5);
		hParentPull_py= new TH1D("hParentPull_py","py",70,-5,5);
		hParentPull_pz= new TH1D("hParentPull_pz","pz",70,-5,5);
		hParentPull_E= new TH1D("hParentPull_E","E",70,-5,5);
		
		/* end pull chunk*/
		/* vertex chunk */
		hVertexPull_x = new TH1D("hVertexPull_x","x",70,-5,5);
		hVertexPull_y = new TH1D("hVertexPull_y","y",70,-5,5);
		hVertexPull_z = new TH1D("hVertexPull_z","z",70,-5,5);
		/* end vertex chunk */

		/*track pulls*/
		hPd0pull = new TH1D("hPd0pull","d0+",70,-5,5);
		hPphipull = new TH1D("hPphipull","phi+",70,-5,5);
		hPomepull = new TH1D("hPomepull","omega+",70,-5,5);
		hPz0pull = new TH1D("hPz0pull","z0+",70,-5,5);
		hPtlpull = new TH1D("hPtlpull","tanLambda+",70,-5,5);

		hMd0pull = new TH1D("hMd0pull","d0-",70,-5,5);
		hMphipull = new TH1D("hMphipull","phi-",70,-5,5);
		hMomepull = new TH1D("hMomepull","omega-",70,-5,5);
		hMz0pull = new TH1D("hMz0pull","z0-",70,-5,5);
		hMtlpull = new TH1D("hMtlpull","tanLambda-",70,-5,5);
		/*end track pulls*/

}
