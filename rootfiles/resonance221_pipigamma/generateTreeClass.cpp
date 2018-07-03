
#include "TTree.h"
#include "TFile.h"
#include "TObjArray.h"
#include "TBranch.h"
#include "TLeaf.h"
#include <vector>
#include <iostream>
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include <fstream>


TTree* tree=0;
std::string treename;

struct leafinfo{
	const char* name;
	const char* type;
	std::vector<std::vector<double> > min;
	std::vector<std::vector<double> > max;
	int nentries=0;
};
double getMin(std::vector<double> v){
	double min = v.at(0);
	for(int i=0; i<v.size(); i++){
		if(v.at(i) < min) min=v.at(i);
	}
	return min;
}
double getMax(std::vector<double> v){
	double max = v.at(0);
	for(int i=0; i<v.size(); i++){
		if(v.at(i) > max) max=v.at(i);
	}
	return max;
}
leafinfo* processDouble(const char* name, const char* type){
	leafinfo* L = new leafinfo();
	L->name = name;
	L->type = type;
	TTreeReader myReader(tree);
	TTreeReaderValue<Double_t> myvar(myReader, name);
	
	std::vector<double> data;
	while(myReader.Next()){
		data.push_back(*myvar);
		L->nentries++;
	}
	std::vector<double> minshell;
	std::vector<double> maxshell;
	minshell.push_back(getMin(data));
	maxshell.push_back(getMax(data));
	L->min.push_back(minshell);
	L->max.push_back(maxshell);
	return L;
}
leafinfo* processVector(const char* name, const char* type){
	leafinfo* L = new leafinfo();
	L->name = name;
	L->type = type;
	TTreeReader myReader(tree);
	TTreeReaderValue<std::vector<double> > myvar(myReader, name);

	std::vector<std::vector<double> > data;
	while(myReader.Next()){
		data.push_back(*myvar);
		L->nentries++;
	}
	
	//std::cout<<L->name<<" "<<L->type<< " "<< data.at(0).size() <<" "<<data.size()<<std::endl;
	//loop through one value at a time
	std::vector<double> subdata;
	int i=0;
	std::vector<double>minshell;
	std::vector<double>maxshell;
	for(int j=0; j<data.at(i).size() ; j++){
		for( i=0; i<data.size(); i++){
			subdata.push_back( data.at(i).at(j) );
		}
		
		minshell.push_back(getMin(subdata));
		maxshell.push_back(getMax(subdata));
		subdata.clear();
		i=0;
	}
	L->min.push_back(minshell);
	L->max.push_back(maxshell);
	
		
	return L;
}
leafinfo* processVectorInt(const char* name, const char* type){
	leafinfo* L = new leafinfo();
	L->name = name;
	L->type = type;
	TTreeReader myReader(tree);
	TTreeReaderValue<std::vector<int> > myvar(myReader, name);

	std::vector<std::vector<int> > data;
	while(myReader.Next()){
		data.push_back(*myvar);
		L->nentries++;
	}
	
	//std::cout<<L->name<<" "<<L->type<< " "<< data.at(0).size() <<" "<<data.size()<<std::endl;
	//loop through one value at a time
	std::vector<double> subdata;
	int i=0;
	std::vector<double>minshell;
	std::vector<double>maxshell;
	for(int j=0; j<data.at(i).size() ; j++){
		for( i=0; i<data.size(); i++){
			subdata.push_back( double(data.at(i).at(j)) );
		}
		
		minshell.push_back(getMin(subdata));
		maxshell.push_back(getMax(subdata));
		subdata.clear();
		i=0;
	}
	L->min.push_back(minshell);
	L->max.push_back(maxshell);
	
		
	return L;
}
 leafinfo* process2DVector(const char* name, const char* type){
	leafinfo* L = new leafinfo();
	L->name = name;
	L->type = type;
	TTreeReader myReader(tree);
	TTreeReaderValue<std::vector<std::vector<double> > > myvar(myReader, name);

	std::vector<std::vector<std::vector<double> > > data;
	while(myReader.Next()){
		data.push_back(*myvar);
		L->nentries++;
	}
	//std::cout<<name<<" "<<type<<" ";
	//std::cout<<data.size()<<" ";
	//std::cout<<data.at(0).size()<<" "<<data.at(0).at(0).size() <<" ";
	//std::cout<<data.at(0).at(0).at(0)<<std::endl;
	/*for(int i=0; i<data.size(); i++){
		
		for(int j=0; j<data.at(i).size(); j++){
			
			for(int k=0; k<data.at(i).at(j).size();k++){
				std::cout<<"event "<<i<<" ";
				std::cout<<"particle "<< j<<" ";
				std::cout<<"parameter "<<k<<" "<< data.at(i).at(j).at(k) <<std::endl;
			}
		}
	} */
	std::vector<double> subdata;  
	//i-event# j-particle# k-parameter#
	//fix j and k and collect all i of k
	std::vector<double> minshell;
	std::vector<double> maxshell;
	int i=0;
	int j=0;
	int k=0;
	for( j=0; j<data.at(i).size(); j++){
		//std::cout<<" i j k in j "<<i<<" "<<j<<" "<<k<<std::endl;
	//	std::cout<<"in j size "<<data.at(i).size()<<std::endl;
		for( k=0; k< data.at(i).at(j).size(); k++){
		//std::cout<<"i j k in k "<<i<<" "<<j<<" "<<k<<std::endl;
	//	std::cout<<"in k size "<<data.at(i).at(j).size()<<std::endl;		
			for( i=0; i<data.size(); i++){
	
	//			std::cout<<"event "<<i<<" ";
	//			std::cout<<"particle "<< j<<" ";
	//			std::cout<<"parameter "<<k<<" "<< data.at(i).at(j).at(k) <<std::endl;
				subdata.push_back( data.at(i).at(j).at(k) );
	//			std::cout<<"i j k pushed"<<std::endl;
				
			}
	//		std::cout<<"subdata size "<<subdata.size()<<std::endl;
			if(subdata.size() != 0){
			minshell.push_back(getMin(subdata));
			maxshell.push_back(getMax(subdata));
			}
			subdata.clear();
			i=0;			
		}
		L->min.push_back(minshell);
		L->max.push_back(maxshell);
		minshell.clear();
		maxshell.clear();
		k=0;
	}
	
	return L;
}
int getBins(double min, double max, int nentries ){
	int bins=0;
		//1% statistics per bin
		double binsize = (max-min)/( 0.01 * nentries );
		bins =(int)  nentries/100; //(max - min)/binsize;
	
	return bins;
}
void printLeafInfo(leafinfo* L){
	std::cout<<L->name<<" "<<L->type<<" "<<L->nentries<<std::endl;
//	for(int i=0; i<L->min.size(); i++){
	//	std::cout<< L->min.at(i)<< " "<<L->max.at(i)<<std::endl;
//	}
	//ith particle jth parameter
	for(int i=0; i<L->min.size(); i++){
		std::cout<< "Particle "<< i <<std::endl;
		for(int j=0; j<L->min.at(i).size(); j++){
			std::cout<<"Parameter "<< j <<" min: "<< L->min.at(i).at(j) << " max: "<< L->max.at(i).at(j) <<std::endl;
		}
	}
}
void writeClass(std::vector<leafinfo*> LV){

	ofstream outputFile("OutputClass.cpp");



	outputFile <<"#include \"TH1D.h\""<<std::endl;
	
	outputFile<<std::endl;

	outputFile << "class OutputClass{ "<< std::endl;
	outputFile << "     public: "<<std::endl;
	int bins;
	for(int i=0; i< LV.size(); i++){
		
		if(!strcmp(LV.at(i)->type,"Double_t")){
			outputFile<<"		TH1D* h"<<LV.at(i)->name<<";"<<std::endl;
		}
		if(!strcmp(LV.at(i)->type,"vector<double>" )){
			//jth particle kth parameter
			for(int j=0; j < LV.at(i)->min.size(); j++){
				for(int k=0; k < LV.at(i)->min.at(j).size(); k++){
					outputFile<<"		TH1D* h"<<LV.at(i)->name<<"_param_"<<k<<";"<<std::endl;
				}
			}
		}
		/*if(!strcmp(LV.at(i)->type,"vector<vector<double> >" )){			
			for(int j=0; j < LV.at(i)->min.size(); j++){
				for(int k=0; k < LV.at(i)->min.at(j).size(); k++){
					outputFile<<"		TH1D* h"<<LV.at(i)->name<<"_part_"<<j<<"_param_"<<k<<";"<<std::endl;
				}
			}
		}*/
		
	}
	outputFile <<"		OutputClass();"<<std::endl;	
	outputFile << "};"<<std::endl;
	outputFile << "OutputClass::OutputClass(){"<<std::endl;
	for(int i=0; i< LV.size(); i++){
			if(!strcmp(LV.at(i)->type,"Double_t")){
				bins = getBins(LV.at(i)->min.at(0).at(0), LV.at(i)->max.at(0).at(0), LV.at(i)->nentries);
				outputFile<<"	h"<<LV.at(i)->name<< " = new TH1D(\"h"<< LV.at(i)->name<<"\", \"default description\", "<<bins<<", "<<LV.at(i)->min.at(0).at(0)<<", "<<LV.at(i)->max.at(0).at(0)<<");"<<std::endl;
			}
			if(!strcmp(LV.at(i)->type,"vector<double>" )){
				for(int j=0; j < LV.at(i)->min.size(); j++){
					for(int k=0; k < LV.at(i)->min.at(j).size(); k++){
						bins = getBins(LV.at(i)->min.at(0).at(0), LV.at(i)->max.at(0).at(0), LV.at(i)->nentries);
						outputFile<<"		h"<<LV.at(i)->name<<"_param_"<<k<< " = new TH1D(\"h"<<LV.at(i)->name<<"_param_"<<k<<"\",\"default description\","<<bins<<", "<<LV.at(i)->min.at(j).at(k)<<", "<<LV.at(i)->max.at(j).at(k)<<");"<<std::endl;
					
				
					}
				}	
			}
	
			if(!strcmp(LV.at(i)->type,"vector<int>" )){
				for(int j=0; j < LV.at(i)->min.size(); j++){
					for(int k=0; k < LV.at(i)->min.at(j).size(); k++){
						bins = getBins(LV.at(i)->min.at(0).at(0), LV.at(i)->max.at(0).at(0), LV.at(i)->nentries);
						outputFile<<"		h"<<LV.at(i)->name<<"_param_"<<k<< " = new TH1D(\"h"<<LV.at(i)->name<<"_param_"<<k<<"\",\"default description\","<<bins<<", "<<LV.at(i)->min.at(j).at(k)<<", "<<LV.at(i)->max.at(j).at(k)<<");"<<std::endl;
					
				
					}
				}	
			}
			/*if(!strcmp(LV.at(i)->type,"vector<vector<double> >" )){
				for(int j=0; j < LV.at(i)->min.size(); j++){
					for(int k=0; k < LV.at(i)->min.at(j).size(); k++){
						bins = getBins(LV.at(i)->min.at(0).at(0), LV.at(i)->max.at(0).at(0), LV.at(i)->nentries);
						outputFile<<"		h"<<LV.at(i)->name<<"_part_"<<j<<"_param_"<<k<< " = new TH1D(\"h"<<LV.at(i)->name<<"_part_"<<j<<"_param_"<<k<<"\",\"default description\","<<bins<<", "<<LV.at(i)->min.at(j).at(k)<<", "<<LV.at(i)->max.at(j).at(k)<<");"<<std::endl;
					
				
					}
				}	
			}*/
	}
	outputFile<<"}";
	outputFile.close();
}
void writeLoop(std::vector<leafinfo*> LV){
	ofstream outputFile("Loop.cpp");
	outputFile <<"#include \"OutputClass.cpp\""<<std::endl;
	outputFile <<"#include \"TTree.h\""<<std::endl;
	outputFile <<"#include \"TTreeReader.h\""<<std::endl;
	outputFile <<"#include \"TTreeReaderValue.h\""<<std::endl;
	outputFile <<"#include \"TFile.h\""<<std::endl;
	outputFile << std::endl;
	outputFile <<"void Loop(){"<<std::endl;
	outputFile <<"	TTree* tree=0;"<<std::endl;
	outputFile <<"	TFile *f = TFile::Open(\"MassConstraint.root\");"<<std::endl;
	outputFile <<" 	TFile *f2 = new TFile(\"test.root\",\"RECREATE\");"<<std::endl;
	outputFile <<"	f->GetObject(\"resonance221_node0\",tree);"<<std::endl;
	outputFile <<"  OutputClass* oc = new OutputClass();"<<std::endl;;
	outputFile <<std::endl;
	
	outputFile<<std::endl;
	for(int i=0; i< LV.size(); i++){
			if(!strcmp(LV.at(i)->type,"Double_t")){
				outputFile<< "	TTreeReader mydoubleReader"<<i<< "(tree);"<<std::endl;
				outputFile<< "	TTreeReaderValue<Double_t> mydoublevar"<<i<<"(mydoubleReader"<<i<<",\""<<LV.at(i)->name<<"\");"<<std::endl;
				outputFile<< "	while(mydoubleReader"<<i<<".Next()){	"<<std::endl;
				outputFile<< "		oc->h"<<LV.at(i)->name<<"->Fill(*mydoublevar"<<i<<");"<<std::endl;
				outputFile<< "	}"<<std::endl;
				//outputFile<< " 	delete mydoubleReader"<<i<<";"<<std::endl;
			}
			if(!strcmp(LV.at(i)->type,"vector<double>" )){
		
				for(int j=0; j<LV.at(i)->min.size(); j++){
					for(int k=0; k<LV.at(i)->min.at(j).size(); k++){
						outputFile<<"  	TTreeReader myvecReader"<<i<<"_"<<k<<"(tree);"<<std::endl;
						outputFile<<"  	TTreeReaderValue<std::vector<double> > myvecvar"<<i<<"_"<<k<<"(myvecReader"<<i<<"_"<<k<<",\""<<LV.at(i)->name<<"\");"<<std::endl;
						outputFile<<"  	while(myvecReader"<<i<<"_"<<k<<".Next()){"<<std::endl;
						outputFile<<"		oc->h"<<LV.at(i)->name<<"_param_"<<k<<"->Fill(myvecvar"<<i<<"_"<<k<<"->at("<<k<<"));"<<std::endl;
						outputFile<<"	}"<<std::endl;;
					}
				}
			}

			if(!strcmp(LV.at(i)->type,"vector<int>" )){
		
				for(int j=0; j<LV.at(i)->min.size(); j++){
					for(int k=0; k<LV.at(i)->min.at(j).size(); k++){
						outputFile<<"  	TTreeReader myvecReader"<<i<<"_"<<k<<"(tree);"<<std::endl;
						outputFile<<"  	TTreeReaderValue<std::vector<int> > myvecvar"<<i<<"_"<<k<<"(myvecReader"<<i<<"_"<<k<<",\""<<LV.at(i)->name<<"\");"<<std::endl;
						outputFile<<"  	while(myvecReader"<<i<<"_"<<k<<".Next()){"<<std::endl;
						outputFile<<"		oc->h"<<LV.at(i)->name<<"_param_"<<k<<"->Fill(myvecvar"<<i<<"_"<<k<<"->at("<<k<<"));"<<std::endl;
						outputFile<<"	}"<<std::endl;;
					}
				}
			}
			/*if(!strcmp(LV.at(i)->type,"vector<vector<double> >" )){
				for(int j=0; j<LV.at(i)->min.size(); j++){
					for(int k=0; k<LV.at(i)->min.at(j).size(); k++){
						outputFile<<"  	TTreeReader my2dvecReader"<<i<<"_"<<j<<"_"<<k<<"(tree);"<<std::endl;
						outputFile<<"  	TTreeReaderValue<std::vector<std::vector<double> > > my2dvecvar"<<i<<"_"<<j<<"_"<<k<<"(my2dvecReader"<<i<<"_"<<j<<"_"<<k<<",\""<<LV.at(i)->name<<"\");"<<std::endl;
						outputFile<<"  	while(my2dvecReader"<<i<<"_"<<j<<"_"<<k<<".Next()){"<<std::endl;
						outputFile<<"		oc->h"<<LV.at(i)->name<<"_part_"<<j<<"_param_"<<k<<"->Fill(my2dvecvar"<<i<<"_"<<j<<"_"<<k<<"->at("<<j<<").at("<<k<<"));"<<std::endl;
						outputFile<<"	}"<<std::endl;;
					}
				}

			}*/
			
		
	}
	outputFile <<"	f2->Write();"<<std::endl;
	outputFile <<"}";
				
			
}

void generateTreeClass(){

	//TTree* tree=0;
	//mod 371 and 306 with treename hardcode

	TFile *f = TFile::Open("MassConstraint.root");
	f->GetObject("resonance221_node0",tree);

	TObjArray* listofbranches;
	TObjArray* listofleaves;

	listofbranches = tree->GetListOfBranches();

	std::vector<const char*> names;
	std::vector<const char*> types;
	
	TBranch* branch;
	TLeaf* leaf;

	//extract names and types from ttree
	for(unsigned int i=0; i<listofbranches->GetLast()+1; i++){

		branch = (TBranch*) listofbranches->At(i);
		listofleaves = branch->GetListOfLeaves();

		for(unsigned int j=0; j<listofleaves->GetLast()+1; j++){

			leaf = (TLeaf*) listofleaves->At(j);

			names.push_back(leaf->GetName());
			types.push_back(leaf->GetTypeName());

			std::cout<<"Found: "<<leaf->GetName()<<" "<<leaf->GetTypeName()<<std::endl;
		}		
	}
	
	//loop over named entries and extract information to generate code
	std::vector<leafinfo*> leafstructvec;

	for(unsigned int i=0; i<names.size(); i++){
	//loop over doubles
		if(!strcmp(types.at(i),"Double_t")){
			leafstructvec.push_back( processDouble(names.at(i), types.at(i)) );
		}
	//loop over vector<double>
		if(!strcmp(types.at(i),"vector<double>")){
			leafstructvec.push_back( processVector(names.at(i), types.at(i)) );
		}
		if(!strcmp(types.at(i),"vector<int>")){
			leafstructvec.push_back( processVectorInt(names.at(i), types.at(i)) );
		}
	//loop over vector<vector<double> >
		/*if(!strcmp(types.at(i),"vector<vector<double> >")){
			leafstructvec.push_back( process2DVector(names.at(i), types.at(i)) );
		}*/	
	}
	for(int i=0; i<leafstructvec.size(); i++){
		printLeafInfo(leafstructvec.at(i));
	}
	writeClass(leafstructvec);
	writeLoop(leafstructvec);
	
}
