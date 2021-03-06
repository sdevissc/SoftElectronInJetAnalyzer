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

void TrainSEtagMVA_BTDG_NonTrigV1_Cat1_2012() {
  
  TMVA::Tools::Instance();
  TFile* outputFile = TFile::Open("SEtagMVA__BDTG_weight.root", "RECREATE");
  TMVA::Factory *factory = new TMVA::Factory("MVA", outputFile, "!V:!Silent:Transformations=I");
//  TFile* input = TFile::Open("./RecoToGenFill_GEDGSF_80_170.root");  
  TFile* input = TFile::Open("./RecoToGenFill_PFinJet_QCD80-170_13022014_ForTraining.root");
  TTree *GenElecB   = (TTree*)input->Get("tree_purity_GenElecB");
  TTree *GenPi      = (TTree*)input->Get("tree_purity_GenPi");
  TTree *GenX       = (TTree*)input->Get("tree_purity_GenX");
  TTree *NonGenElec = (TTree*)input->Get("tree_purity_NonGenElec");
  TTree *NonGenPi   = (TTree*)input->Get("tree_purity_NonGenPi");
  TTree *NonGenX    = (TTree*)input->Get("tree_purity_NonGenX");        

  int nGenElecB=GenElecB->GetEntries();
  int nGenPi=GenPi->GetEntries();
  int nGenX=GenX->GetEntries();
  int nNonGenElec=NonGenElec->GetEntries();
  int nNonGenPi=NonGenPi->GetEntries();
  int nNonGenX=NonGenX->GetEntries();

  double sigWeight = 1.0;
  double bkgWeight = 1.0;
  factory->SetInputTrees(GenElecB, GenPi, sigWeight, bkgWeight);

  cout<<"Adding the variables"<<endl;
  factory->AddVariable("sip2d", 'F');
  factory->AddVariable("sip3d",'F');
  factory->AddVariable("ptRel", 'F');
  factory->AddVariable("ratioRel", 'F');
//  factory->AddVariable("mva_e_pi", 'F');
//  factory->AddVariable( "pt");
//  factory->AddVariable( "eta");


  //  factory->AddVariable("", 'F');
  cout<<"Adding the spectators"<<endl;
  factory->AddSpectator( "pt");  
  factory->AddSpectator( "eta");
  factory->AddSpectator( "sip3d");
  factory->AddSpectator( "tppdgId");
  cout<<"finished"<<endl;

  factory->PrepareTrainingAndTestTree("pt>2 && sip3d>-200 && sip3d<200 && abs(eta)<2.4","abs(tppdgId)!=11 &&pt>2 && abs(eta)<2.4 && sip3d>-200 && sip3d<200","nTrain_Signal=0:nTrain_Background=0:SplitMode=Block:!V" );
//  factory->SetWeightExpression( "weight" );
  

  
  factory->BookMethod(TMVA::Types::kBDT, "SE_BDT_weight", "!H:!V:NTrees=2000:BoostType=Grad:Shrinkage=0.10:!UseBaggedGrad:nCuts=2000:nEventsMin=100:NNodesMax=5:UseNvars=4:PruneStrength=5:PruneMethod=CostComplexity:MaxDepth=6" );



  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();
  outputFile->Close();

  delete factory;
}
