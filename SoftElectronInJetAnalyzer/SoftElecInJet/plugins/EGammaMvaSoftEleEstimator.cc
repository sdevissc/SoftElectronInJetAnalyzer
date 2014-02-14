#include "SoftElectronInJetAnalyzer/SoftElecInJet/interface/EGammaMvaSoftEleEstimator.h"
#include <cmath>
#include <vector>
using namespace std;

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "DataFormats/Common/interface/RefToPtr.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include "DataFormats/Common/interface/RefToPtr.h"
using namespace reco;

//--------------------------------------------------------------------------------------------------
EGammaMvaSoftEleEstimator::EGammaMvaSoftEleEstimator() :
fMethodname("BDTG method"),
fisInitialized(kFALSE),
fMVAType(kTrig),
fUseBinnedVersion(kTRUE),
fNMVABins(0)
{
  // Constructor.  
}

//--------------------------------------------------------------------------------------------------
EGammaMvaSoftEleEstimator::~EGammaMvaSoftEleEstimator()
{
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
}

//--------------------------------------------------------------------------------------------------
void EGammaMvaSoftEleEstimator::initialize( std::string methodName,
                                       std::string weightsfile
                                       )
{
  
  std::vector<std::string> tempWeightFileVector;
  tempWeightFileVector.push_back(weightsfile);
  initialize(methodName,kFALSE,tempWeightFileVector);
}


//--------------------------------------------------------------------------------------------------
void EGammaMvaSoftEleEstimator::initialize( std::string methodName,
                                       Bool_t useBinnedVersion,
				       std::vector<std::string> weightsfiles
  ) {

  //clean up first
  for (unsigned int i=0;i<fTMVAReader.size(); ++i) {
    if (fTMVAReader[i]) delete fTMVAReader[i];
  }
  fTMVAReader.clear();

  //initialize
  fisInitialized = kTRUE;
  fMethodname = methodName;
  fUseBinnedVersion = useBinnedVersion;

  //Define expected number of bins
  // By defaultt, let's run on 3 bins in Pt and 3 bins in Eta, that means 9 bins in total
  UInt_t ExpectedNBins = 1;
  fNMVABins = ExpectedNBins;
  
  //Check number of weight files given
  if (fNMVABins != weightsfiles.size() ) {
    std::cout << "Error: Expected Number of bins = " << fNMVABins << " does not equal to weightsfiles.size() = " 
              << weightsfiles.size() << std::endl; 
 
    assert(fNMVABins == weightsfiles.size());
  }

  //Loop over all bins
  for (unsigned int i=0;i<fNMVABins; ++i) {


  
    	TMVA::Reader *tmpTMVAReader = new TMVA::Reader( "!Color:!Silent:Error" );  
    	tmpTMVAReader->SetVerbose(kTRUE);
	tmpTMVAReader->AddVariable("fbrem", 		&fMVAVar_fbrem);
	tmpTMVAReader->AddVariable("EtotOvePin",		&fMVAVar_EtotOvePin);
	tmpTMVAReader->AddVariable("EClusOverPout", 		&fMVAVar_eleEoPout);
	tmpTMVAReader->AddVariable("EBremOverDeltaP", 		&fMVAVar_EBremOverDeltaP);
	tmpTMVAReader->AddVariable("logSigmaEtaEta", 		&fMVAVar_logSigmaEtaEta);
	tmpTMVAReader->AddVariable("DeltaEtaTrackEcalSeed", 	&fMVAVar_DeltaEtaTrackEcalSeed);
	tmpTMVAReader->AddVariable("HoE", 			&fMVAVar_HoE);
	tmpTMVAReader->AddVariable("gsfchi2", 			&fMVAVar_gsfchi2);
	tmpTMVAReader->AddVariable("kfchi2", 			&fMVAVar_kfchi2);
	tmpTMVAReader->AddVariable("kfhits", 			&fMVAVar_kfhits);
	tmpTMVAReader->AddVariable("SigmaPtOverPt", 		&fMVAVar_SigmaPtOverPt);
	tmpTMVAReader->AddVariable("deta", 			&fMVAVar_deta);
	tmpTMVAReader->AddVariable("dphi", 			&fMVAVar_dphi);
	tmpTMVAReader->AddVariable("detacalo", 			&fMVAVar_detacalo);
	tmpTMVAReader->AddVariable("see", 			&fMVAVar_see);
        tmpTMVAReader->AddVariable("spp",                       &fMVAVar_spp);
        tmpTMVAReader->AddVariable("R9",                        &fMVAVar_R9);
	tmpTMVAReader->AddVariable("etawidth", 			&fMVAVar_etawidth);
	tmpTMVAReader->AddVariable("phiwidth", 			&fMVAVar_phiwidth);
	tmpTMVAReader->AddVariable("e1x5e5x5", 			&fMVAVar_OneMinusE1x5E5x5);
    	tmpTMVAReader->AddVariable("IoEmIoP",                   &fMVAVar_IoEmIoP);
    	tmpTMVAReader->AddVariable("PreShowerOverRaw",          &fMVAVar_PreShowerOverRaw);
	tmpTMVAReader->AddVariable("nPV",          		&fMVAVar_nPV);

	tmpTMVAReader->AddVariable( "pt",		&fMVAVar_pt);
	tmpTMVAReader->AddVariable( "eta",		&fMVAVar_eta); 

        tmpTMVAReader->AddSpectator( "pt",               &fMVAVar_pt);
        tmpTMVAReader->AddSpectator( "eta",              &fMVAVar_eta);

    	tmpTMVAReader->BookMVA(fMethodname , weightsfiles[i]);
//    	std::cout << "MVABin " << i << " : MethodName = " << fMethodname 
//          	    << "Load weights file : " << weightsfiles[i] 
//          	    << std::endl;
    	fTMVAReader.push_back(tmpTMVAReader);
  }
//  std::cout << "Electron ID MVA Completed\n";

}


//--------------------------------------------------------------------------------------------------
UInt_t EGammaMvaSoftEleEstimator::GetMVABin(double eta, double pt) const {
  
    //Default is to return the first bin
    unsigned int bin = 0;

	bool ptrange[3],etarange[3];
	ptrange[0]=pt > 2 && pt < 5;
	ptrange[1]=pt > 5 && pt < 10;
	ptrange[2]=pt > 10; 
	etarange[0]=fabs(eta) < 0.8;
	etarange[1]=fabs(eta) > 0.8 && fabs(eta) <1.4;
	etarange[2]=fabs(eta) > 1.4;

	int index=0;
	for(int kETA=0;kETA<3;kETA++)
	for(int kPT=0;kPT<3;kPT++){
		if (ptrange[kPT] && etarange[kETA]) bin=index;
		index++;
	}	
/*
      	if (ptrange[0] && etarange[0])	bin = 0;
	if (ptrange[1] && etarange[0])  bin = 1;
	if (ptrange[2] && etarange[0])  bin = 2;
	if (ptrange[0] && etarange[1])  bin = 3;
	if (ptrange[1] && etarange[1])  bin = 4;
	if (ptrange[2] && etarange[1])  bin = 5;
	if (ptrange[0] && etarange[2])  bin = 6;
	if (ptrange[1] && etarange[2])  bin = 7;
	if (ptrange[2] && etarange[2])  bin = 8;
  */  return bin;
}


//--------------------------------------------------------------------------------------------------
Double_t EGammaMvaSoftEleEstimator::mvaValue(const reco::GsfElectron& electron,const edm::Event & evt,Bool_t printDebug) {
  
  if (!fisInitialized) { 
    std::cout << "Error: EGammaMvaSoftEleEstimator not properly initialized.\n"; 
    return -9999;
  }


  edm::Handle<reco::VertexCollection> FullprimaryVertexCollection;
  evt.getByLabel("offlinePrimaryVertices", FullprimaryVertexCollection);
  const reco::VertexCollection pvc = *(FullprimaryVertexCollection.product());

  fMVAVar_fbrem			=electron.fbrem();
  fMVAVar_EtotOvePin		=electron.eSuperClusterOverP();
  fMVAVar_eleEoPout		=electron.eEleClusterOverPout();
  float etot			=electron.eSuperClusterOverP()*electron.trackMomentumAtVtx().R();
  float eEcal			=electron.eEleClusterOverPout()*electron.trackMomentumAtEleClus().R();
  float dP			=electron.trackMomentumAtVtx().R()-electron.trackMomentumAtEleClus().R();
  fMVAVar_EBremOverDeltaP	=(etot-eEcal)/dP;	
  fMVAVar_logSigmaEtaEta	=log(electron.sigmaEtaEta());
  fMVAVar_DeltaEtaTrackEcalSeed	=electron.deltaEtaEleClusterTrackAtCalo();
  fMVAVar_HoE			=electron.hcalOverEcalBc();
  bool validKF= false;
  reco::TrackRef myTrackRef 	= electron.closestCtfTrackRef();
  validKF 			= (myTrackRef.isAvailable());
  validKF 			= (myTrackRef.isNonnull());	
  fMVAVar_kfchi2		=(validKF) ? myTrackRef->normalizedChi2() : 0 ;
  fMVAVar_kfhits		=(validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
  fMVAVar_gsfchi2		=electron.gsfTrack()->normalizedChi2();
  fMVAVar_SigmaPtOverPt		=electron.gsfTrack().get()->ptModeError()/electron.gsfTrack().get()->ptMode() ;
  fMVAVar_deta			=electron.deltaEtaSuperClusterTrackAtVtx();
  fMVAVar_dphi			=electron.deltaPhiSuperClusterTrackAtVtx();
  fMVAVar_detacalo		=electron.deltaEtaSeedClusterTrackAtCalo();
  fMVAVar_see			=electron.sigmaIetaIeta();
  fMVAVar_spp                   =electron.sigmaIphiIphi();
  fMVAVar_R9                    =electron.r9();	
  fMVAVar_etawidth		=electron.superCluster()->etaWidth();
  fMVAVar_phiwidth		=electron.superCluster()->phiWidth();
  fMVAVar_OneMinusE1x5E5x5	=(electron.e5x5()) !=0. ? 1.-(electron.e1x5()/electron.e5x5()) : -1. ;
  fMVAVar_IoEmIoP               =  (1.0/electron.ecalEnergy()) - (1.0 / electron.p());
  fMVAVar_pt			=electron.pt();
  fMVAVar_eta			=electron.eta();
  fMVAVar_nPV			=pvc.size();
 fMVAVar_PreShowerOverRaw=electron.superCluster()->preshowerEnergy() / electron.superCluster()->rawEnergy();

  if(printDebug) {
    cout << " *** Inside the class fMethodname " << fMethodname << endl;
    cout << " fbrem " <<  fMVAVar_fbrem<< endl;
    cout        << "fbrem "   <<fMVAVar_fbrem<< endl;
    cout        << "EtotOvePin "<<                fMVAVar_EtotOvePin<< endl;
     cout       << "EClusOverPout "<<             fMVAVar_eleEoPout<< endl;
    cout        << "EBremOverDeltaP "<<           fMVAVar_EBremOverDeltaP<< endl;
    cout        << "logSigmaEtaEta "<<            fMVAVar_logSigmaEtaEta<< endl;
    cout        << "DeltaEtaTrackEcalSeed "<<     fMVAVar_DeltaEtaTrackEcalSeed<< endl;
    cout        << "HoE "      <<                 fMVAVar_HoE<< endl;
    cout        << "gsfchi2 "    <<               fMVAVar_gsfchi2<< endl;
    cout        << "kfchi2 "       <<             fMVAVar_kfchi2<< endl;
    cout        << "kfhits "        <<            fMVAVar_kfhits<< endl;
    cout        << "SigmaPtOverPt "   <<          fMVAVar_SigmaPtOverPt<< endl;
    cout        << "deta "          <<            fMVAVar_deta<< endl;
    cout        << "dphi "            <<          fMVAVar_dphi<< endl;
    cout        << "detacalo "       <<           fMVAVar_detacalo<< endl;
    cout        << "see "            <<           fMVAVar_see<< endl;
    cout        << "spp "             <<          fMVAVar_spp<< endl;
    cout        << "R9 "             <<           fMVAVar_R9<< endl;
    cout        << "etawidth "        <<          fMVAVar_etawidth<< endl;
    cout        << "phiwidth "       <<           fMVAVar_phiwidth<< endl;
    cout        << "e1x5e5x5 "       <<           fMVAVar_OneMinusE1x5E5x5<< endl;
    cout        << "IoEmIoP "        <<           fMVAVar_IoEmIoP<< endl;
     cout       << "PreShowerOverRaw " << fMVAVar_PreShowerOverRaw <<endl;

  }



  bindVariables();
  Double_t mva = -9999;  
  if (fUseBinnedVersion) {
  //  mva = fTMVAReader[GetMVABin(fMVAVar_eta,fMVAVar_pt)]->EvaluateMVA(fMethodname);
    mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
  } else {
    mva = fTMVAReader[0]->EvaluateMVA(fMethodname);
  }
//cout << " ### MVA " << mva << endl;


  return mva;
}
//--------------------------------------------------------------------------------------------------------



void EGammaMvaSoftEleEstimator::bindVariables() {

  // this binding is needed for variables that sometime diverge. 


  if(fMVAVar_fbrem < -1.)
    fMVAVar_fbrem = -1.;	
  
  fMVAVar_deta = fabs(fMVAVar_deta);
  if(fMVAVar_deta > 0.06)
    fMVAVar_deta = 0.06;
  
  
  fMVAVar_dphi = fabs(fMVAVar_dphi);
  if(fMVAVar_dphi > 0.6)
    fMVAVar_dphi = 0.6;
  

  if(fMVAVar_EoP > 20.)
    fMVAVar_EoP = 20.;
  
  if(fMVAVar_eleEoPout > 20.)
    fMVAVar_eleEoPout = 20.;
  
  
  fMVAVar_detacalo = fabs(fMVAVar_detacalo);
  if(fMVAVar_detacalo > 0.2)
    fMVAVar_detacalo = 0.2;
  
  if(fMVAVar_OneMinusE1x5E5x5 < -1.)
    fMVAVar_OneMinusE1x5E5x5 = -1;
  
  if(fMVAVar_OneMinusE1x5E5x5 > 2.)
    fMVAVar_OneMinusE1x5E5x5 = 2.; 
  
  
  
  if(fMVAVar_gsfchi2 > 200.)
    fMVAVar_gsfchi2 = 200;
  
  
  if(fMVAVar_kfchi2 > 10.)
    fMVAVar_kfchi2 = 10.;
  
  
  
  return;
}








