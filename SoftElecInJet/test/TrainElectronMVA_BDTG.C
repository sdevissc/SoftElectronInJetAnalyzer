#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/MethodBase.h"
#include "TMVA/MethodCategory.h"
#endif

void TrainElectronMVA_BTDG_NonTrigV1_Cat1_2012() {
  
  TMVA::Tools::Instance();
  TFile* outputFile = TFile::Open("SoftElectronMVA__BDTG.root", "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory("MVA", outputFile, "!V:!Silent");
  


  cout<<"Adding the variables"<<endl;
  factory->AddVariable("fbrem", 'F');
  factory->AddVariable("EtotOvePin",'F');
  factory->AddVariable("EClusOverPout", 'F');
  factory->AddVariable("EBremOverDeltaP", 'F');
  factory->AddVariable("logSigmaEtaEta", 'F');
  factory->AddVariable("DeltaEtaTrackEcalSeed", 'F');
  factory->AddVariable("HOverE", 'F');
  factory->AddVariable("Chi2GSF", 'F');
  factory->AddVariable("CHi2KF", 'F');
  factory->AddVariable("nHits", 'F');
  factory->AddVariable("SigmaPtOverPt", 'F');
  factory->AddVariable("deta", 'F');
  factory->AddVariable("dphi", 'F');
  factory->AddVariable("detacalo", 'F');
  factory->AddVariable("see", 'F');
  factory->AddVariable("etawidth", 'F');
  factory->AddVariable("phiwidth", 'F');
  factory->AddVariable("e1x5e5x5", 'F');


  //  factory->AddVariable("", 'F');
  cout<<"Adding the spectators"<<endl;
  factory->AddSpectator( "gedgsfElecPt");  
  factory->AddSpectator( "gedgsfElecEta");
  //  factory->AddSpectator( "matchConv");


  std::vector<Double_t> varf( 20 );
  Float_t  treevars[20];
  TFile* input = TFile::Open("./ZZZ");	
  cout<<"input file defined"<<endl;
  TTree *tree   = (TTree*)input->Get("tree_purity");	
  cout<<"tree defined"<<endl;	
  tree->SetBranchAddress( "fbrem", &(treevars[0]) );
  tree->SetBranchAddress( "EtotOvePin", &(treevars[1]) );
  tree->SetBranchAddress( "EClusOverPout", &(treevars[2]) );
  tree->SetBranchAddress( "EBremOverDeltaP", &(treevars[3]) );
  tree->SetBranchAddress( "logSigmaEtaEta", &(treevars[4]) );
  tree->SetBranchAddress( "DeltaEtaTrackEcalSeed", &(treevars[5]) );
  tree->SetBranchAddress( "HOverE", &(treevars[6]) );
  tree->SetBranchAddress( "Chi2GSF", &(treevars[7]) );
  tree->SetBranchAddress( "CHi2KF", &(treevars[8]) );
  tree->SetBranchAddress( "nHits", &(treevars[9]) );
  tree->SetBranchAddress( "SigmaPtOverPt", &(treevars[10]) );	
  tree->SetBranchAddress( "deta", &(treevars[11]) );
  tree->SetBranchAddress( "dphi", &(treevars[12]) );
  tree->SetBranchAddress( "detacalo", &(treevars[13]) );
  tree->SetBranchAddress( "see", &(treevars[14]) );
  tree->SetBranchAddress( "etawidth", &(treevars[15]) );
  tree->SetBranchAddress( "phiwidth", &(treevars[16]) );
  tree->SetBranchAddress( "e1x5e5x5", &(treevars[17]) );


  tree->SetBranchAddress( "gedgsfElecPt",&(treevars[18]));
  tree->SetBranchAddress( "gedgsfElecEta",&(treevars[19])); 

 
  cout<<"setbranches done"<<endl;
  Int_t classif[3];
  tree->SetBranchAddress( "pdgId", &(classif[0]) );
  tree->SetBranchAddress( "tppdgId", &(classif[1]));
  tree->SetBranchAddress( "origin", &(classif[2]) );	 

  cout<<"Now will run on the events"<<endl;
  for (UInt_t i=0; i<tree->GetEntries(); i++) {
	tree->GetEntry(i);
	for (UInt_t ivar=0; ivar<18; ivar++) varf[ivar] = treevars[ivar];	
	cout<<"event --->";
	for(int u=0;u<20;u++)cout<<treevars[u]<<" ";
	cout<<classif[0]<<" "<<classif[1]<<" "<<classif[2]<<endl;
	if(fabs(classif[0])==11 && (classif[2]=4 || classif[2]==6)){
		if (i < tree->GetEntries()/2.0){
			cout<<"-> To Signal Training"<<endl;
			factory->AddSignalTrainingEvent( varf, 1 );	
		}
		else{
			cout<<"-> To Signal Test"<<endl;
			factory->AddSignalTestEvent( varf, 1 );
		}
		
	}
	else if((fabs(classif[1])==211/* && (classif[2]==1)*/)){
//        else if((fabs(classif[1])==11 &&  (classif[2]<0))){
		if (i < tree->GetEntries()/2.0){
			cout<<"-> To BackGround Training"<<endl;
                        factory->AddBackgroundTrainingEvent( varf, 1 );
                }
		else{
			cout<<"-> To BackGround Test"<<endl;	
                        factory->AddBackgroundTestEvent( varf, 1 );
                
		
		}
  	}
  }	 
  cout<<"finished"<<endl;
  factory->PrepareTrainingAndTestTree("","","nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
  

  
  factory->BookMethod(TMVA::Types::kBDT, "BDT", "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );



  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  outputFile->Close();

  delete factory;
}
