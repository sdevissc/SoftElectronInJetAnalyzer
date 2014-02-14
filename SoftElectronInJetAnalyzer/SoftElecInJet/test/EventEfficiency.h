//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Sep 24 15:46:51 2013 by ROOT version 5.34/09
// from TTree tree_efficiency/Reconst ntuple
// found on file: test.root
//////////////////////////////////////////////////////////

#ifndef EventEfficiency_h
#define EventEfficiency_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.
const Int_t kMaxgen_match = 1;
const Int_t kMaxtrk_pt = 1;
const Int_t kMaxtrk_pTOB = 1;
const Int_t kMaxtrk_eta = 1;
const Int_t kMaxtrk_phi = 1;
const Int_t kMaxtrk_chi2 = 1;
const Int_t kMaxtrk_ndof = 1;
const Int_t kMaxtrk_nchi2 = 1;
const Int_t kMaxtrk_dpt = 1;
const Int_t kMaxtrk_nhits = 1;
const Int_t kMaxtrk_nmatched = 1;
const Int_t kMaxtrk_quality = 1;
const Int_t kMaxtrk_ecalDist = 1;
const Int_t kMaxtrk_ecalDeta = 1;
const Int_t kMaxtrk_ecalDphi = 1;
const Int_t kMaxecal_e = 1;
const Int_t kMaxtrk_ep = 1;
const Int_t kMaxtrk_dptGSF = 1;
const Int_t kMaxtrk_chiRatio = 1;
const Int_t kMaxtrk_chiReduced = 1;
const Int_t kMaxtrk_ecalChi = 1;
const Int_t kMaxecal_ps = 1;
const Int_t kMaxecal_ps1e = 1;
const Int_t kMaxecal_ps2e = 1;
const Int_t kMaxtrk_epCorr = 1;

class EventEfficiency {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           origin;
   Int_t           nPV;
   Float_t         ecalEnergy;
   Float_t         pIn;
   Float_t         pOut;
   Float_t         EemPinRatio;
   Float_t         EemPoutRatio;
   Float_t         GGE_p;
   Float_t         GGE_pt;
   Float_t         GGE_eta;
   Float_t         GGE_phi;
   Float_t         GGE_energy;
   Float_t         fbrem;
   Float_t         dRGsfTrackElectron;
   Float_t         inversedRFirstLastHit;
   Float_t         radiusFirstHit;
   Float_t         zFirstHit;
   Float_t         insidePFJet;
   Float_t         ptGen;
   Float_t         etaGen;
   Float_t         phiGen;
   Float_t         pGen;
   Float_t         dRJetReco;
   Float_t         ptRel;
   Float_t         mindistJR;
   Int_t           gen_match_;
   Float_t         trk_pt_;
   Float_t         trk_pTOB_;
   Float_t         trk_eta_;
   Float_t         trk_phi_;
   Float_t         trk_chi2_;
   Int_t           trk_ndof_;
   Float_t         trk_nchi2_;
   Float_t         trk_dpt_;
   Int_t           trk_nhits_;
   Int_t           trk_nmatched_;
   Float_t         trk_quality_;
   Float_t         trk_ecalDist_;
   Float_t         trk_ecalDeta_;
   Float_t         trk_ecalDphi_;
   Float_t         ecal_e_;
   Float_t         trk_ep_;
   Float_t         trk_dptGSF_;
   Float_t         trk_chiRatio_;
   Float_t         trk_chiReduced_;
   Float_t         trk_ecalChi_;
   Float_t         ecal_ps_;
   Float_t         ecal_ps1e_;
   Float_t         ecal_ps2e_;
   Float_t         trk_epCorr_;
   Float_t         gsftrkSeedPt;
   Float_t         gsftrkSeedEta;
   Float_t         gsftrkSeedPhi;
   Int_t           inRecoJet;
   Int_t           isMatchedWithASeed;
   Float_t         gedgsfElecPt;
   Float_t         gedgsfElecEta;
   Float_t         gedgsfElecPhi;
   Int_t           isMatchedWithAGedGsfElec;

   // List of branches
   TBranch        *b_origin;   //!
   TBranch        *b_nPV;   //!
   TBranch        *b_ecalEnergy;   //!
   TBranch        *b_pIn;   //!
   TBranch        *b_pOut;   //!
   TBranch        *b_EemPinRatio;   //!
   TBranch        *b_EemPoutRatio;   //!
   TBranch        *b_GGE_p;   //!
   TBranch        *b_GGE_pt;   //!
   TBranch        *b_GGE_eta;   //!
   TBranch        *b_GGE_phi;   //!
   TBranch        *b_GGE_energy;   //!
   TBranch        *b_fbrem;   //!
   TBranch        *b_dRGsfTrackElectron;   //!
   TBranch        *b_inversedRFirstLastHit;   //!
   TBranch        *b_radiusFirstHit;   //!
   TBranch        *b_zFirstHit;   //!
   TBranch        *b_insidePFJet;   //!
   TBranch        *b_ptGen;   //!
   TBranch        *b_etaGen;   //!
   TBranch        *b_phiGen;   //!
   TBranch        *b_pGen;   //!
   TBranch        *b_dRJetReco;   //!
   TBranch        *b_ptRel;   //!
   TBranch        *b_mindistJR;   //!
   TBranch        *b_gen_match_;   //!
   TBranch        *b_trk_pt_;   //!
   TBranch        *b_trk_pTOB_;   //!
   TBranch        *b_trk_eta_;   //!
   TBranch        *b_trk_phi_;   //!
   TBranch        *b_trk_chi2_;   //!
   TBranch        *b_trk_ndof_;   //!
   TBranch        *b_trk_nchi2_;   //!
   TBranch        *b_trk_dpt_;   //!
   TBranch        *b_trk_nhits_;   //!
   TBranch        *b_trk_nmatched_;   //!
   TBranch        *b_trk_quality_;   //!
   TBranch        *b_trk_ecalDist_;   //!
   TBranch        *b_trk_ecalDeta_;   //!
   TBranch        *b_trk_ecalDphi_;   //!
   TBranch        *b_ecal_e_;   //!
   TBranch        *b_trk_ep_;   //!
   TBranch        *b_trk_dptGSF_;   //!
   TBranch        *b_trk_chiRatio_;   //!
   TBranch        *b_trk_chiReduced_;   //!
   TBranch        *b_trk_ecalChi_;   //!
   TBranch        *b_ecal_ps_;   //!
   TBranch        *b_ecal_ps1e_;   //!
   TBranch        *b_ecal_ps2e_;   //!
   TBranch        *b_trk_epCorr_;   //!
   TBranch        *b_gsftrkSeedPt;   //!
   TBranch        *b_gsftrkSeedEta;   //!
   TBranch        *b_gsftrkSeedPhi;   //!
   TBranch        *b_inRecoJet;   //!
   TBranch        *b_isMatchedWithASeed;   //!
   TBranch        *b_gedgsfElecPt;   //!
   TBranch        *b_gedgsfElecEta;   //!
   TBranch        *b_gedgsfElecPhi;   //!
   TBranch        *b_isMatchedWithAGedGsfElec;   //!

   EventEfficiency(TTree *tree=0);
   virtual ~EventEfficiency();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef EventEfficiency_cxx
EventEfficiency::EventEfficiency(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test.root");
      }
      f->GetObject("tree_efficiency",tree);

   }
   Init(tree);
}

EventEfficiency::~EventEfficiency()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t EventEfficiency::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t EventEfficiency::LoadTree(Long64_t entry)
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

void EventEfficiency::Init(TTree *tree)
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

   fChain->SetBranchAddress("origin", &origin, &b_origin);
   fChain->SetBranchAddress("nPV", &nPV, &b_nPV);
   fChain->SetBranchAddress("ecalEnergy", &ecalEnergy, &b_ecalEnergy);
   fChain->SetBranchAddress("pIn", &pIn, &b_pIn);
   fChain->SetBranchAddress("pOut", &pOut, &b_pOut);
   fChain->SetBranchAddress("EemPinRatio", &EemPinRatio, &b_EemPinRatio);
   fChain->SetBranchAddress("EemPoutRatio", &EemPoutRatio, &b_EemPoutRatio);
   fChain->SetBranchAddress("GGE_p", &GGE_p, &b_GGE_p);
   fChain->SetBranchAddress("GGE_pt", &GGE_pt, &b_GGE_pt);
   fChain->SetBranchAddress("GGE_eta", &GGE_eta, &b_GGE_eta);
   fChain->SetBranchAddress("GGE_phi", &GGE_phi, &b_GGE_phi);
   fChain->SetBranchAddress("GGE_energy", &GGE_energy, &b_GGE_energy);
   fChain->SetBranchAddress("fbrem", &fbrem, &b_fbrem);
   fChain->SetBranchAddress("dRGsfTrackElectron", &dRGsfTrackElectron, &b_dRGsfTrackElectron);
   fChain->SetBranchAddress("inversedRFirstLastHit", &inversedRFirstLastHit, &b_inversedRFirstLastHit);
   fChain->SetBranchAddress("radiusFirstHit", &radiusFirstHit, &b_radiusFirstHit);
   fChain->SetBranchAddress("zFirstHit", &zFirstHit, &b_zFirstHit);
   fChain->SetBranchAddress("insidePFJet", &insidePFJet, &b_insidePFJet);
   fChain->SetBranchAddress("ptGen", &ptGen, &b_ptGen);
   fChain->SetBranchAddress("etaGen", &etaGen, &b_etaGen);
   fChain->SetBranchAddress("phiGen", &phiGen, &b_phiGen);
   fChain->SetBranchAddress("pGen", &pGen, &b_pGen);
   fChain->SetBranchAddress("dRJetReco", &dRJetReco, &b_dRJetReco);
   fChain->SetBranchAddress("ptRel", &ptRel, &b_ptRel);
   fChain->SetBranchAddress("mindistJR", &mindistJR, &b_mindistJR);
   fChain->SetBranchAddress("gen_match_", &gen_match_, &b_gen_match_);
   fChain->SetBranchAddress("trk_pt_", &trk_pt_, &b_trk_pt_);
   fChain->SetBranchAddress("trk_pTOB_", &trk_pTOB_, &b_trk_pTOB_);
   fChain->SetBranchAddress("trk_eta_", &trk_eta_, &b_trk_eta_);
   fChain->SetBranchAddress("trk_phi_", &trk_phi_, &b_trk_phi_);
   fChain->SetBranchAddress("trk_chi2_", &trk_chi2_, &b_trk_chi2_);
   fChain->SetBranchAddress("trk_ndof_", &trk_ndof_, &b_trk_ndof_);
   fChain->SetBranchAddress("trk_nchi2_", &trk_nchi2_, &b_trk_nchi2_);
   fChain->SetBranchAddress("trk_dpt_", &trk_dpt_, &b_trk_dpt_);
   fChain->SetBranchAddress("trk_nhits_", &trk_nhits_, &b_trk_nhits_);
   fChain->SetBranchAddress("trk_nmatched_", &trk_nmatched_, &b_trk_nmatched_);
   fChain->SetBranchAddress("trk_quality_", &trk_quality_, &b_trk_quality_);
   fChain->SetBranchAddress("trk_ecalDist_", &trk_ecalDist_, &b_trk_ecalDist_);
   fChain->SetBranchAddress("trk_ecalDeta_", &trk_ecalDeta_, &b_trk_ecalDeta_);
   fChain->SetBranchAddress("trk_ecalDphi_", &trk_ecalDphi_, &b_trk_ecalDphi_);
   fChain->SetBranchAddress("ecal_e_", &ecal_e_, &b_ecal_e_);
   fChain->SetBranchAddress("trk_ep_", &trk_ep_, &b_trk_ep_);
   fChain->SetBranchAddress("trk_dptGSF_", &trk_dptGSF_, &b_trk_dptGSF_);
   fChain->SetBranchAddress("trk_chiRatio_", &trk_chiRatio_, &b_trk_chiRatio_);
   fChain->SetBranchAddress("trk_chiReduced_", &trk_chiReduced_, &b_trk_chiReduced_);
   fChain->SetBranchAddress("trk_ecalChi_", &trk_ecalChi_, &b_trk_ecalChi_);
   fChain->SetBranchAddress("ecal_ps_", &ecal_ps_, &b_ecal_ps_);
   fChain->SetBranchAddress("ecal_ps1e_", &ecal_ps1e_, &b_ecal_ps1e_);
   fChain->SetBranchAddress("ecal_ps2e_", &ecal_ps2e_, &b_ecal_ps2e_);
   fChain->SetBranchAddress("trk_epCorr_", &trk_epCorr_, &b_trk_epCorr_);
   fChain->SetBranchAddress("gsftrkSeedPt", &gsftrkSeedPt, &b_gsftrkSeedPt);
   fChain->SetBranchAddress("gsftrkSeedEta", &gsftrkSeedEta, &b_gsftrkSeedEta);
   fChain->SetBranchAddress("gsftrkSeedPhi", &gsftrkSeedPhi, &b_gsftrkSeedPhi);
   fChain->SetBranchAddress("inRecoJet", &inRecoJet, &b_inRecoJet);
   fChain->SetBranchAddress("isMatchedWithASeed", &isMatchedWithASeed, &b_isMatchedWithASeed);
   fChain->SetBranchAddress("gedgsfElecPt", &gedgsfElecPt, &b_gedgsfElecPt);
   fChain->SetBranchAddress("gedgsfElecEta", &gedgsfElecEta, &b_gedgsfElecEta);
   fChain->SetBranchAddress("gedgsfElecPhi", &gedgsfElecPhi, &b_gedgsfElecPhi);
   fChain->SetBranchAddress("isMatchedWithAGedGsfElec", &isMatchedWithAGedGsfElec, &b_isMatchedWithAGedGsfElec);
   Notify();
}

Bool_t EventEfficiency::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void EventEfficiency::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t EventEfficiency::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef EventEfficiency_cxx
