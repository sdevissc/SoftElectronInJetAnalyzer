#ifndef jprtag
#define jprtag

#include <TROOT.h>
#include <TTree.h>
#include <TFile.h>


class JPsiReco{
	public:
	TFile *g;
	TTree *tree_JpsiUpsilon;
	inline void initJPsiRecoFillerObject();
	inline void WriteInFileAndCloseIt();
	float Mass;
	inline JPsiReco();
	inline ~JPsiReco();
};

#endif

JPsiReco::JPsiReco()
{
	g=new TFile("JpsiUpsilonFill.root","recreate");
	tree_JpsiUpsilon= new TTree("tree_JpsiUpsilon","JU ntuple");
	tree_JpsiUpsilon->Branch("Mass",&Mass,"Mass/F");
}

void JPsiReco::initJPsiRecoFillerObject(){
	Mass=-777.0;	
}

void JPsiReco::WriteInFileAndCloseIt(){ 
	g->cd();
        tree_JpsiUpsilon->Write();
        g->Close();
}

JPsiReco::~JPsiReco(){}	

#ifndef rtgtag
#define rtgtag

class RecoToGenFiller{
	public:
	inline RecoToGenFiller(const TString &);
	inline ~RecoToGenFiller();
	inline void initRecoToGenFillerObject();
	inline void WriteInFileAndCloseIt();
	TTree *tree_purity;
	TTree *tree_purity_GenElecB;
	TTree *tree_purity_GenPi;
	TTree *tree_purity_GenX;
	TTree *tree_purity_NonGenElec;
        TTree *tree_purity_NonGenPi;
        TTree *tree_purity_NonGenX;
	TFile *f;
	int origin;
        int nPV;
        float ecalEnergy;
        float pIn;
        float pOut;
        float EemPinRatio;
        float EemPoutRatio;
        float GGE_p;
        float GGE_pt;
        float GGE_eta;
        float GGE_phi;
        float GGE_energy;
        float fbrem;
        float dRGsfTrackElectron;
        float inversedRFirstLastHit;
        float radiusFirstHit;
        float zFirstHit;
        float insidePFJet;
        float ptGen;
        float etaGen;
        float phiGen;
        float ptRel;
        float mindistJR;
        float pGen;
	float dRJetReco;
	int gen_match_;
        float trk_pt_;
        float trk_pTOB_;
        float trk_eta_;
        float trk_phi_;
        float trk_chi2_;
        int trk_ndof_;
        float trk_nchi2_;
        float trk_dpt_;
        int trk_nhits_;
        int trk_nmatched_;
        float trk_quality_;
        float trk_ecalDist_;
        float trk_ecalDeta_;
        float trk_ecalDphi_;
        float ecal_e_;
        float trk_ep_;
        float trk_dptGSF_;
        float trk_chiRatio_;
        float trk_chiReduced_;
	float trk_ecalChi_;
	float trk_epCorr_;
	float ecal_ps_;
        float ecal_ps1e_;
        float ecal_ps2e_;
	float gsftrkSeedPt;
	float gsftrkSeedEta;
	float gsftrkSeedPhi;
	float gedgsfElecPt;
        float gedgsfElecEta;
        float gedgsfElecPhi;
	int inRecoJet;
	int isMatchedWithASeed;
	int isMatchedWithAGedGsfElec;
	int Vtx;
        int pdgId;
	int tppdgId;
	float EtotOvePin;
	float EClusOverPout;
        float EBremOverDeltaP;
        float logSigmaEtaEta;
        float DeltaEtaTrackEcalSeed;
        float HOverE;
        float Chi2GSF;
        float CHi2KF;
        float nHits;
        float SigmaPtOverPt;
        float lnPt;
	float deta;
        float dphi;
        float detacalo;
        float see;
        float etawidth;
        float phiwidth;
        float e1x5e5x5;
        float R9;
        float HoE;
        float EoP;
        float IoEmIoP;
        float eleEoPout;
        float d0;
        float ip3d;
        float ip3dSig;
	float mva_e_pi;
        float sip2d;
        float sip3d;
        float deltaR;
        float ptRel;
        float etaRel;
        float ratio;
        float ratioRel;

};

RecoToGenFiller::RecoToGenFiller(const TString & tag){
	f=new TFile("RecoToGenFill.root","recreate");
	tree_purity_GenElecB = new TTree("tree_purity_GenElecB","Reconst ntuple");
        tree_purity_GenElecB->Branch("origin",&origin,"origin/I");
        tree_purity_GenElecB->Branch("nPV",&nPV,"nPV/I");
        tree_purity_GenElecB->Branch("ecalEnergy",&ecalEnergy,"ecalEnergy/F");
        tree_purity_GenElecB->Branch("pIn",&pIn,"pIn/F");
        tree_purity_GenElecB->Branch("pOut",&pOut,"pOut/F");
        tree_purity_GenElecB->Branch("EemPinRatio",&EemPinRatio,"EemPinRatio/F");
        tree_purity_GenElecB->Branch("EemPoutRatio",&EemPoutRatio,"EemPoutRatio/F");
        tree_purity_GenElecB->Branch("GGE_p",&GGE_p,"GGE_p/F");
        tree_purity_GenElecB->Branch("GGE_pt",&GGE_pt,"GGE_pt/F");
        tree_purity_GenElecB->Branch("GGE_eta",&GGE_eta,"GGE_eta/F");
        tree_purity_GenElecB->Branch("GGE_phi",&GGE_phi,"GGE_phi/F");
        tree_purity_GenElecB->Branch("GGE_energy",&GGE_energy,"GGE_energy/F");
        tree_purity_GenElecB->Branch("fbrem",&fbrem,"fbrem/F");
        tree_purity_GenElecB->Branch("dRGsfTrackElectron",&dRGsfTrackElectron,"dRGsfTrackElectron/F");
        tree_purity_GenElecB->Branch("inversedRFirstLastHit",&inversedRFirstLastHit,"inversedRFirstLastHit/F");
        tree_purity_GenElecB->Branch("radiusFirstHit",&radiusFirstHit,"radiusFirstHit/F");
        tree_purity_GenElecB->Branch("zFirstHit",&zFirstHit,"zFirstHit/F");
        tree_purity_GenElecB->Branch("insidePFJet",&insidePFJet,"insidePFJet/F");
        tree_purity_GenElecB->Branch("ptGen",&ptGen,"ptGen/F");
        tree_purity_GenElecB->Branch("etaGen",&etaGen,"etaGen/F");
        tree_purity_GenElecB->Branch("phiGen",&phiGen,"phiGen/F");
        tree_purity_GenElecB->Branch("pGen",&pGen,"pGen/F");
        tree_purity_GenElecB->Branch("dRJetReco",&dRJetReco,"dRJetReco/F");
        tree_purity_GenElecB->Branch("ptRel",&ptRel,"ptRel/F");
        tree_purity_GenElecB->Branch("mindistJR",&mindistJR,"mindistJR/F");	
	tree_purity_GenElecB->Branch("gen_match_",&gen_match_,"gen_match_/I");
	tree_purity_GenElecB->Branch("trk_pt_",&trk_pt_,"trk_pt_/F");
	tree_purity_GenElecB->Branch("trk_pTOB_",&trk_pTOB_,"trk_pTOB_/F");
	tree_purity_GenElecB->Branch("trk_eta_",&trk_eta_,"trk_eta_/F");
	tree_purity_GenElecB->Branch("trk_phi_",&trk_phi_,"trk_phi_/F");
	tree_purity_GenElecB->Branch("trk_chi2_",&trk_chi2_,"trk_chi2_/F");
	tree_purity_GenElecB->Branch("trk_ndof_",&trk_ndof_,"trk_ndof_/I");
	tree_purity_GenElecB->Branch("trk_nchi2_",&trk_nchi2_,"trk_nchi2_/F");
	tree_purity_GenElecB->Branch("trk_dpt_",&trk_dpt_,"trk_dpt_/F");
	tree_purity_GenElecB->Branch("trk_nhits_",&trk_nhits_,"trk_nhits_/I");
	tree_purity_GenElecB->Branch("trk_nmatched_",&trk_nmatched_,"trk_nmatched_/I");
	tree_purity_GenElecB->Branch("trk_quality_",&trk_quality_,"trk_quality_/F");
	tree_purity_GenElecB->Branch("trk_ecalDist_",&trk_ecalDist_,"trk_ecalDist_/F");
	tree_purity_GenElecB->Branch("trk_ecalDeta_",&trk_ecalDeta_,"trk_ecalDeta_/F");
	tree_purity_GenElecB->Branch("trk_ecalDphi_",&trk_ecalDphi_,"trk_ecalDphi_/F");
	tree_purity_GenElecB->Branch("ecal_e_",&ecal_e_,"ecal_e_/F");
	tree_purity_GenElecB->Branch("trk_ep_",&trk_ep_,"trk_ep_/F");
	tree_purity_GenElecB->Branch("trk_dptGSF_",&trk_dptGSF_,"trk_dptGSF_/F");
	tree_purity_GenElecB->Branch("trk_chiRatio_",&trk_chiRatio_,"trk_chiRatio_/F");
	tree_purity_GenElecB->Branch("trk_chiReduced_",&trk_chiReduced_,"trk_chiReduced_/F");
	tree_purity_GenElecB->Branch("trk_ecalChi_",&trk_ecalChi_,"trk_ecalChi_/F");
        tree_purity_GenElecB->Branch("ecal_ps_",&ecal_ps_,"ecal_ps_/F");
        tree_purity_GenElecB->Branch("ecal_ps1e_",&ecal_ps1e_,"ecal_ps1e_/F");
        tree_purity_GenElecB->Branch("ecal_ps2e_",&ecal_ps2e_,"ecal_ps2e_/F");
	tree_purity_GenElecB->Branch("trk_epCorr_",&trk_epCorr_,"trk_epCorr_/F");
	tree_purity_GenElecB->Branch("gsftrkSeedPt",&gsftrkSeedPt,"gsftrkSeedPt/F");
        tree_purity_GenElecB->Branch("gsftrkSeedEta",&gsftrkSeedEta,"gsftrkSeedEta/F");
        tree_purity_GenElecB->Branch("gsftrkSeedPhi",&gsftrkSeedPhi,"gsftrkSeedPhi/F");
	tree_purity_GenElecB->Branch("inRecoJet",&inRecoJet,"inRecoJet/I");
	tree_purity_GenElecB->Branch("isMatchedWithASeed",&isMatchedWithASeed,"isMatchedWithASeed/I");
	tree_purity_GenElecB->Branch("pt",&gedgsfElecPt,"pt/F");
        tree_purity_GenElecB->Branch("eta",&gedgsfElecEta,"eta/F");
        tree_purity_GenElecB->Branch("gedgsfElecPhi",&gedgsfElecPhi,"gedgsfElecPhi/F");
	tree_purity_GenElecB->Branch("isMatchedWithAGedGsfElec",&isMatchedWithAGedGsfElec,"isMatchedWithAGedGsfElec/I");
	tree_purity_GenElecB->Branch("Vtx",&Vtx,"Vtx/I");
        tree_purity_GenElecB->Branch("pdgId",&pdgId,"pdgId/I");
	tree_purity_GenElecB->Branch("tppdgId",&tppdgId,"tppdgId/I");
	tree_purity_GenElecB->Branch("EtotOvePin",&EtotOvePin,"EtotOvePin/F");
        tree_purity_GenElecB->Branch("EClusOverPout",&EClusOverPout,"EClusOverPout/F");
        tree_purity_GenElecB->Branch("fbrem",&fbrem,"fbrem/F");
        tree_purity_GenElecB->Branch("EBremOverDeltaP",&EBremOverDeltaP,"EBremOverDeltaP/F");
        tree_purity_GenElecB->Branch("logSigmaEtaEta",&logSigmaEtaEta,"logSigmaEtaEta/F");
        tree_purity_GenElecB->Branch("DeltaEtaTrackEcalSeed",&DeltaEtaTrackEcalSeed,"DeltaEtaTrackEcalSeed/F");
        tree_purity_GenElecB->Branch("HOverE",&HOverE,"HOverE/F");
        tree_purity_GenElecB->Branch("gsfchi2",&Chi2GSF,"gsfchi2/F");
        tree_purity_GenElecB->Branch("kfchi2",&CHi2KF,"kfchi2/F");
        tree_purity_GenElecB->Branch("kfhits",&nHits,"kfhits/F");
        tree_purity_GenElecB->Branch("SigmaPtOverPt",&SigmaPtOverPt,"SigmaPtOverPt/F");
        tree_purity_GenElecB->Branch("lnPt",&lnPt,"lnPt/F");
	tree_purity_GenElecB->Branch("deta",&deta,"deta/F");
        tree_purity_GenElecB->Branch("dphi",&dphi,"dphi/F");
        tree_purity_GenElecB->Branch("detacalo",&detacalo,"detacalo/F");
        tree_purity_GenElecB->Branch("see",&see,"see/F");
        tree_purity_GenElecB->Branch("etawidth",&etawidth,"etawidth/F");
        tree_purity_GenElecB->Branch("phiwidth",&phiwidth,"phiwidth/F");
        tree_purity_GenElecB->Branch("e1x5e5x5",&e1x5e5x5,"e1x5e5x5/F");
        tree_purity_GenElecB->Branch("R9",&R9,"R9/F");
        tree_purity_GenElecB->Branch("HoE",&HoE,"HoE/F");
        tree_purity_GenElecB->Branch("EoP",&EoP,"EoP/F");
        tree_purity_GenElecB->Branch("IoEmIoP",&IoEmIoP,"IoEmIoP/F");
        tree_purity_GenElecB->Branch("eleEoPout",&eleEoPout,"eleEoPout/F");
        tree_purity_GenElecB->Branch("d0",&d0,"d0/F");
        tree_purity_GenElecB->Branch("ip3d",&ip3d,"ip3d/F");
        tree_purity_GenElecB->Branch("ip3dSig",&ip3dSig,"ip3dSig/F");
	tree_purity_GenElecB->Branch("mva_e_pi",&mva_e_pi,"mva_e_pi/F");
        tree_purity_GenElecB->Branch("sip2d",&sip2d,"sip2d/F");
        tree_purity_GenElecB->Branch("sip3d",&sip3d,"sip3d/F");
        tree_purity_GenElecB->Branch("deltaR",&deltaR,"deltaR/F");
        tree_purity_GenElecB->Branch("ptRel",&ptRel,"ptRel/F");
        tree_purity_GenElecB->Branch("etaRel",&etaRel,"etaRel/F");
        tree_purity_GenElecB->Branch("ratio",&ratio,"ratio/F");
        tree_purity_GenElecB->Branch("ratioRel",&ratioRel,"ratioRel/F");
	tree_purity_GenPi=(TTree*)tree_purity_GenElecB->CloneTree();
	tree_purity_GenX=(TTree*)tree_purity_GenElecB->CloneTree();
	tree_purity_NonGenElec=(TTree*)tree_purity_GenElecB->CloneTree();
	tree_purity_NonGenPi=(TTree*)tree_purity_GenElecB->CloneTree();
	tree_purity_NonGenX=(TTree*)tree_purity_GenElecB->CloneTree();
        tree_purity=(TTree*)tree_purity_GenElecB->CloneTree();
	tree_purity_GenPi->SetName("tree_purity_GenPi");
        tree_purity_GenX->SetName("tree_purity_GenX");
        tree_purity_NonGenElec->SetName("tree_purity_NonGenElec");
        tree_purity_NonGenPi->SetName("tree_purity_NonGenPi");
        tree_purity_NonGenX->SetName("tree_purity_NonGenX");	
        tree_purity->SetName("tree_purity");
	
}

void RecoToGenFiller::initRecoToGenFillerObject(){
	tppdgId=-777;	
	origin=				-777;
	nPV                         =-1;
        ecalEnergy                  = -777.0;
        pIn                         = -777.0;
        pOut                        = -777.0;
        EemPinRatio                 = -777.0;
        EemPoutRatio                = -777.0;
        GGE_p                      = -777.0;
        GGE_pt                      = -777.0;
        GGE_eta                     = -777.0;
        GGE_phi                     = -777.0;
        GGE_energy                  = -777.0;
        fbrem                       = -777.0;
        dRGsfTrackElectron          = -777.0;
        insidePFJet                 = 0;
        ptGen                       =-777.0;
        etaGen                      =-777.0;
        phiGen                      =-777.0;
        pGen                        =-777.0;

	gen_match_ = -1;
        trk_pt_       = -777.0;
        trk_pTOB_     = -777.0;
        trk_eta_      = -777.0;
        trk_phi_      = -777.0;
        trk_chi2_     = -777.0;
        trk_ndof_     = -777.0;
        trk_nchi2_    = -777.0;
        trk_dpt_      = -777.0;
        trk_nhits_    = -777.0;
        trk_nmatched_ = -777.0;
        trk_quality_  = -777.0;
	trk_ecalDist_ = -777.0;
        trk_ecalDeta_ = -777.0;
        trk_ecalDphi_ = -777.0;	
	ecal_e_ = -777.0;
        trk_ep_ = -777.0;
	trk_dptGSF_ = -777.0;
	trk_chiRatio_ = -777.0;
        trk_chiReduced_ = -777.0;
	ecal_ps_ = -777.0;
        ecal_ps1e_ = -777.0;
        ecal_ps2e_ =-777.0;
	trk_epCorr_= -777.0;	
	gsftrkSeedPt = -777.0;
	gsftrkSeedEta = -777.0;
	gsftrkSeedPhi = -777.0;
	gedgsfElecPt = -777.0;
        gedgsfElecEta = -777.0;
        gedgsfElecPhi = -777.0;
	inRecoJet = 0;
	isMatchedWithASeed=0;
	isMatchedWithAGedGsfElec=0;	
	Vtx=-777;
        pdgId=-777;
	EtotOvePin=-777;
        EClusOverPout=-777;
        fbrem=-777;
        EBremOverDeltaP=-777;;
        logSigmaEtaEta=-777;;
        DeltaEtaTrackEcalSeed=-777;;
        HOverE=-777;;
        Chi2GSF=-777;;
        CHi2KF=-777;;
        nHits=-777;;
        SigmaPtOverPt=-777;;
        lnPt=-777;;


	deta=-777;
        dphi=-777;
        detacalo=-777;
        see=-777;
	etawidth=-777;
        phiwidth=-777;
        e1x5e5x5=-777;
        R9=-777;
        HoE=-777;
        EoP=-777;
        IoEmIoP=-777;
        eleEoPout=-777;
        d0=-777;
	ip3d=-777;
	ip3dSig=-777;
	mva_e_pi = -777.0;
	sip2d = -777.0;
        sip3d = -777.0;
        deltaR = -777.0;
        ptRel = -777.0;
        etaRel = -777.0;
        ratio = -777.0;
        ratioRel = -777.0
}


void RecoToGenFiller::WriteInFileAndCloseIt(){
	f->cd();
	tree_purity->Write();
	tree_purity_GenElecB->Write();
	tree_purity_GenPi->Write();
        tree_purity_GenX->Write();
        tree_purity_NonGenElec->Write();
        tree_purity_NonGenPi->Write();
        tree_purity_NonGenX->Write();
	f->Close();
}
RecoToGenFiller::~RecoToGenFiller(){};

#endif
