//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Jan 27 12:14:28 2013 by ROOT version 5.32/00
// from TTree tree_efficiency/Reconst ntuple
// found on file: ElectronTree.root
//////////////////////////////////////////////////////////

#ifndef EventPurity_h
#define EventPurity_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class EventPurity {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   TTree *tree;
   // Declaration of leaf types
	
   Float_t EtotOvePin;
   Float_t EClusOverPout;
   Float_t fbrem;
   Float_t EBremOverDeltaP;
   Float_t  logSigmaEtaEta;
   Float_t DeltaEtaTrackEcalSeed;
   Float_t HOverE;
   Float_t gsfchi2;
   Float_t kfchi2;
   Float_t kfhits;
   Float_t SigmaPtOverPt;
   Float_t lnPt;	
   Float_t deta;
   Float_t dphi;
   Float_t detacalo;
   Float_t see;
   Float_t etawidth;
   Float_t phiwidth;
   Float_t e1x5e5x5;
   Float_t R9;
   Float_t HoE;
   Float_t spp;
   Float_t PreShowerOverRaw;	 	
   Float_t EoP;
   Float_t IoEmIoP;
   Float_t eleEoPout;
   Float_t d0;
   Float_t ip3d;
   Float_t ip3dSig;
   Float_t mva_e_pi; 	


   Int_t fromConversion;

   Int_t tppdgId;
   Int_t Vtx;
   Int_t pdgId;
   Int_t origin;
   Float_t pt;
   Float_t eta;
   Float_t phi;
   Float_t nPV;



   // List of branches
   TBranch        *b_EtotOvePin;
   TBranch        *b_EClusOverPout;
   TBranch        *b_fbrem;
   TBranch        *b_EBremOverDeltaP;
   TBranch        *b_logSigmaEtaEta;
   TBranch        *b_DeltaEtaTrackEcalSeed;
   TBranch        *b_HOverE;
   TBranch        *b_gsfchi2;
   TBranch        *b_kfchi2;
   TBranch        *b_kfhits;
   TBranch        *b_SigmaPtOverPt;
   TBranch        *b_lnPt;
   TBranch        *b_deta;
   TBranch        *b_dphi;
   TBranch        *b_detacalo;
   TBranch        *b_see;
   TBranch        *b_etawidth;
   TBranch        *b_phiwidth;
   TBranch        *b_e1x5e5x5;
   TBranch        *b_R9;
   TBranch        *b_HoE;
   TBranch        *b_EoP;
   TBranch        *b_IoEmIoP;
   TBranch        *b_eleEoPout;
   TBranch        *b_d0;
   TBranch        *b_ip3d;
   TBranch        *b_ip3dSig;
   TBranch        *b_mva_e_pi;

   TBranch        *b_spp;
   TBranch        *b_PreShowerOverRaw;	
   TBranch        *b_tppdgId;
   TBranch        *b_Vtx;
   TBranch        *b_pdgId;
   TBranch        *b_origin;
   TBranch        *b_pt;
   TBranch        *b_eta;
   TBranch        *b_phi;
   TBranch        *b_nPV;
   TBranch 	  *b_fromConversion; 	

   EventPurity( TString filename );
   virtual ~EventPurity();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
//   virtual Float_t StringToFloat(TString);	
};

#endif

#ifdef EventPurity_cxx

EventPurity::EventPurity( TString filename )
{

  tree = 0;

  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
      if (!f) {
         f = new TFile(filename);
      }
      tree = (TTree*)gDirectory->Get("tree_purity");

   }
   Init(tree);
}


/*EventPurity::EventPurity(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("ElectronTree.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("ElectronTree.root");
      }
      f->GetObject("tree_efficiency",tree);

   }
   Init(tree);
}*/

EventPurity::~EventPurity()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventPurity::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventPurity::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void EventPurity::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("EtotOvePin",&EtotOvePin,&b_EtotOvePin);
   fChain->SetBranchAddress("EClusOverPout",&EClusOverPout,&b_EClusOverPout);
   fChain->SetBranchAddress("fbrem",&fbrem,&b_fbrem);
   fChain->SetBranchAddress("EBremOverDeltaP",&EBremOverDeltaP,&b_EBremOverDeltaP);
   fChain->SetBranchAddress("logSigmaEtaEta",&logSigmaEtaEta,&b_logSigmaEtaEta);
   fChain->SetBranchAddress("DeltaEtaTrackEcalSeed",&DeltaEtaTrackEcalSeed,&b_DeltaEtaTrackEcalSeed);
   fChain->SetBranchAddress("HOverE",&HOverE,&b_HOverE);
   fChain->SetBranchAddress("gsfchi2",&gsfchi2,&b_gsfchi2);
   fChain->SetBranchAddress("kfchi2",&kfchi2,&b_kfchi2);
   fChain->SetBranchAddress("kfhits",&kfhits,&b_kfhits);
   fChain->SetBranchAddress("SigmaPtOverPt",&SigmaPtOverPt,&b_SigmaPtOverPt);
   fChain->SetBranchAddress("lnPt",&lnPt,&b_lnPt);
   fChain->SetBranchAddress("deta",&deta,&b_deta);
   fChain->SetBranchAddress("dphi",&dphi,&b_dphi);
   fChain->SetBranchAddress("detacalo",&detacalo,&b_detacalo);
   fChain->SetBranchAddress("see",&see,&b_see);
   fChain->SetBranchAddress("etawidth",&etawidth,&b_etawidth);
   fChain->SetBranchAddress("phiwidth",&phiwidth,&b_phiwidth);
   fChain->SetBranchAddress("e1x5e5x5",&e1x5e5x5,&b_e1x5e5x5);
   fChain->SetBranchAddress("R9",&R9,&b_R9);
   fChain->SetBranchAddress("HoE",&HoE,&b_HoE);
   fChain->SetBranchAddress("EoP",&EoP,&b_EoP);
   fChain->SetBranchAddress("IoEmIoP",&IoEmIoP,&b_IoEmIoP);
   fChain->SetBranchAddress("eleEoPout",&eleEoPout,&b_eleEoPout);
   fChain->SetBranchAddress("d0",&d0,&b_d0);
   fChain->SetBranchAddress("ip3d",&ip3d,&b_ip3d);
   fChain->SetBranchAddress("ip3dSig",&ip3dSig,&b_ip3dSig);
   fChain->SetBranchAddress("mva_e_pi",&mva_e_pi,&b_mva_e_pi);
   fChain->SetBranchAddress("spp",&spp,&b_spp);
   fChain->SetBranchAddress("PreShowerOverRaw",&PreShowerOverRaw,&b_PreShowerOverRaw);



   fChain->SetBranchAddress("tppdgId",&tppdgId,&b_tppdgId);
   fChain->SetBranchAddress("Vtx",&Vtx,&b_Vtx);
   fChain->SetBranchAddress("pdgId",&pdgId,&b_pdgId);
   fChain->SetBranchAddress("origin",&origin,&b_origin);
   fChain->SetBranchAddress("pt",&pt,&b_pt);
   fChain->SetBranchAddress("eta",&eta,&b_eta);
   fChain->SetBranchAddress("phi",&phi,&b_phi);
   fChain->SetBranchAddress("nPV",&nPV,&b_nPV);
   fChain->SetBranchAddress("fromConversion",&fromConversion,&b_fromConversion);
   Notify();
}

Bool_t EventPurity::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventPurity::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventPurity::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

#endif 

