//--------------------------------------------------------------------------------------------------
// $Id $
//
// EGammaMvaSoftEleEstimator
//
// Helper Class for applying MVA electron ID selection
//
// Authors: D.Benedetti, E.DiMaro, S.Xie
//--------------------------------------------------------------------------------------------------


/// --> NOTE if you want to use this class as standalone without the CMSSW part 
///  you need to uncomment the below line and compile normally with scramv1 b 
///  Then you need just to load it in your root macro the lib with the correct path, eg:
///  gSystem->Load("/data/benedet/CMSSW_5_2_2/lib/slc5_amd64_gcc462/pluginEGammaEGammaAnalysisTools.so");

//#define STANDALONE   // <---- this line

#ifndef EGammaMvaSoftEleEstimator_H
#define EGammaMvaSoftEleEstimator_H

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidateFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "EgammaAnalysis/ElectronTools/interface/ElectronEffectiveArea.h"
#include <vector>
#include <TROOT.h>
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

class EGammaMvaSoftEleEstimator{
  public:
    EGammaMvaSoftEleEstimator();
    ~EGammaMvaSoftEleEstimator(); 
  
    enum MVAType {
      kTrig = 0,                     // MVA for triggering electrons     
      kTrigNoIP = 1,                     // MVA for triggering electrons without IP info
      kNonTrig = 2,                      // MVA for non-triggering electrons     
      kIsoRings,                     // Isolation MVA for non-trigger electrons
      kTrigIDIsoCombined,            // ID+Iso Combined MVA for triggering electrons
      kTrigIDIsoCombinedPUCorrected  // ID+Iso Combined MVA for triggering electrons
    };
  
    void     initialize( std::string methodName,
                         std::string weightsfile);
    void     initialize( std::string methodName,
                         Bool_t useBinnedVersion,
                         std::vector<std::string> weightsfiles );
    
    Bool_t   isInitialized() const { return fisInitialized; }
    UInt_t   GetMVABin(double eta,double pt ) const;
    
    void bindVariables();
    
    Double_t mvaValue(const reco::GsfElectron& electron,const edm::Event & evt,bool printDebug = false);


  private:

    std::vector<TMVA::Reader*> fTMVAReader;
    std::string                fMethodname;
    Bool_t                     fisInitialized;
    MVAType                    fMVAType;
    Bool_t                     fUseBinnedVersion;
    UInt_t                     fNMVABins;

    Float_t                    fMVAVar_fbrem;
    Float_t 	               fMVAVar_EtotOvePin;
    Float_t                    fMVAVar_EBremOverDeltaP;	
    Float_t                    fMVAVar_logSigmaEtaEta;	
    Float_t                    fMVAVar_DeltaEtaTrackEcalSeed; 
    Float_t                    fMVAVar_kfchi2;
    Float_t                    fMVAVar_kfhits;    //number of layers
    Float_t                    fMVAVar_gsfchi2;
    Float_t	               fMVAVar_SigmaPtOverPt; 


    Float_t                    fMVAVar_deta;
    Float_t                    fMVAVar_dphi;
    Float_t                    fMVAVar_detacalo;

    Float_t                    fMVAVar_see;
    Float_t                    fMVAVar_etawidth;
    Float_t                    fMVAVar_phiwidth;
    Float_t                    fMVAVar_OneMinusE1x5E5x5;

    Float_t                    fMVAVar_HoE;
    Float_t                    fMVAVar_EoP;
    Float_t                    fMVAVar_eleEoPout;

    Float_t 			fMVAVar_spp;
    Float_t 			fMVAVar_R9;
    Float_t    			fMVAVar_IoEmIoP;
    Float_t        		fMVAVar_PreShowerOverRaw;	

    Float_t                    fMVAVar_eta;
    Float_t                    fMVAVar_pt;
    Float_t 		       fMVAVar_nPV;		  
 
};

#endif
