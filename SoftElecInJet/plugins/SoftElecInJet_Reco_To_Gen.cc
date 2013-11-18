// -*- C++ -*-
//
// Package:    SoftElectronInJetAnalyzer/SoftElecInJetRecoToGen
// Class:      SoftElecInJetRecoToGen
// 
/**\class SoftElecInJetRecoToGen SoftElecInJetRecoToGen.cc SoftElectronInJetAnalyzer/SoftElecInJetRecoToGen/plugins/SoftElecInJetRecoToGen.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Simon De Visscher
//         Created:  Fri, 13 Sep 2013 11:43:57 GMT
//
//



#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedTrackerVertex.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "TrackingTools/TransientTrack/interface/TransientTrackBuilder.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "RecoParticleFlow/PFTracking/interface/ElectronSeedMerger.h"
#include "DataFormats/TrajectorySeed/interface/PropagationDirection.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/HepMCCandidate/interface/GenParticleFwd.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "RecoParticleFlow/PFTracking/interface/GoodSeedProducer.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/ParticleFlowReco/interface/PreId.h"
#include "DataFormats/ParticleFlowReco/interface/PreIdFwd.h"
#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TSCBLBuilderNoMaterial.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "DataFormats/Math/interface/deltaR.h"


//add for the input variables
#include "RecoParticleFlow/PFTracking/interface/PFTrackTransformer.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFResolutionMap.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryFitter.h"
#include "TrackingTools/PatternTools/interface/TrajectorySmoother.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "FastSimulation/BaseParticlePropagator/interface/BaseParticlePropagator.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoParticleFlow/PFClusterTools/interface/LinkByRecHit.h"
#include "RecoParticleFlow/PFClusterTools/interface/PFEnergyCalibration.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecTrack.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateMode.h"
#include "TrackingTools/GsfTools/interface/MultiTrajectoryStateTransform.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "DataFormats/HeavyIonEvent/interface/EvtPlane.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Track/interface/CoreSimTrack.h"
#include "SimDataFormats/EncodedEventId/interface/EncodedEventId.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "SimDataFormats/TrackingHit/interface/UpdatablePSimHit.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticleFwd.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingVertexContainer.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorBase.h"
#include "SimTracker/TrackAssociation/interface/TrackAssociatorByChi2.h"
#include "SimTracker/TrackAssociation/interface/QuickTrackAssociatorByHits.h"
#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"
#include "SimTracker/Records/interface/TrackAssociatorRecord.h"
#include "SimTracker/TrackAssociation/plugins/ParametersDefinerForTPESProducer.h"
#include "SimDataFormats/TrackingAnalysis/interface/TrackingParticle.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "TrackingTools/PatternTools/interface/TwoTrackMinimumDistance.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
#include "RecoVertex/AdaptiveVertexFinder/interface/TracksClusteringFromDisplacedSeed.h"

#include "RecoVertex/ConfigurableVertexReco/interface/ConfigurableVertexReconstructor.h"

#include "DataFormats/EgammaReco/interface/ElectronSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"

#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h"

#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "DataFormats/TauReco/interface/PFTau.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "DataFormats/Common/interface/Association.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "RecoParticleFlow/PFProducer/interface/Utils.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexUpdator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexTrackCompatibilityEstimator.h"
#include "RecoVertex/KalmanVertexFit/interface/KalmanVertexSmoother.h"
#include "RecoVertex/MultiVertexFit/interface/MultiVertexFitter.h"

#include <TLorentzVector.h>

#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>
#include "TMath.h"
#include "TVector2.h"
#include "Math/VectorUtil.h"

#include "GenToRecoClass.h"
#include "RecoToGenClass.h"
//
// class declaration
//


using namespace edm;
using namespace std;
using namespace reco;

#define DEBUG2  0

class SoftElecInJetRecoToGen : public edm::EDAnalyzer {
   public:
      SoftElecInJetRecoToGen(const edm::ParameterSet&);
      ~SoftElecInJetRecoToGen();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      bool trackFilter(const reco::TrackRef &track) const;	
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      bool findMatch(const reco::Candidate& genCand, Handle<TrackingParticleCollection> TPCollectionH,reco::SimToRecoCollection q,reco::Track& track,int& indexTrack,float& dR,float& sharedHits,RefToBase<Track> &trRefToBase); 
      bool  isaBhadron( int );
      bool  isaDhadron( int );
      bool  isaV( int );
      int GetOrigin(reco::GenParticle&);
      std::vector<TransientVertex>	 GetVCandidates(std::vector<TransientTrack>);	
//        void GenToReco(GsfElectronCollection ,GenParticleCollection,PFJetCollection,VertexCollection,reco::SimToRecoCollection,edm::Handle<TrackingParticleCollection>,vector<Trajectory>);
//        void RecoToGen(GsfElectronCollection ,GenParticleCollection,PFJetCollection,const Vertex);
	  // checks if two of the tracks are innermost false if track1 is near, true if track2 is near, returns true for sameLayer if they are in the same layer
  	virtual bool isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer);
	double testPreshowerDistance(PFCluster eeclus,
                               PFCluster psclus);
	virtual void beginRun(const edm::Run & run, const edm::EventSetup &) override;

         edm::InputTag gsfElectronTag_,pfJetsTag_;
	GenToRecoFiller *gtrf;
	RecoToGenFiller *rtgf;
	JPsiReco *jpr;
	edm::InputTag gedgsfElectronTag_,genParticleTag_,PFJetTag_,PrimaryVerticesTag_,TrackingParticleTag_,tracksTag_,tracksTag;
	unsigned int                            minHits;
        unsigned int                            maxNTracks;
        double                                  maxLIP;
	double 					minPt;

	std::auto_ptr<VertexReconstructor>      vtxReco;
        std::auto_ptr<TracksClusteringFromDisplacedSeed>        clusterizer;

        double                                  vertexMinAngleCosine;
        double                                  vertexMinDLen2DSig;
	double 					vertexMinDLenSig;
	edm::InputTag pfCLusTagECLabel_;
  	edm::InputTag pfCLusTagHCLabel_;
  	edm::InputTag pfCLusTagPSLabel_;
  	std::vector<edm::InputTag> tracksContainers_;
	std::string   fitterName_;
  	std::string   smootherName_;
	//  edm::InputTag PreIdMapLabel_;
  	edm::InputTag gsfTrackLabel_;
  	edm::InputTag pfNuclear_;
  	edm::InputTag pfTrackLabel_;
  	edm::InputTag simtracksTag;
  	///TRACK QUALITY
  	bool useQuality_;
  	reco::TrackBase::TrackQuality trackQuality_;
  	Double_t      clusThreshold_;
  	Double_t      minEp_;
  	Double_t      maxEp_;
  	Double_t      minPt_;
  	Double_t      maxPt_;
  	Bool_t        usePreshower_;
	PFTrackTransformer* pfTkTransformer_;
	math::XYZVector B_;
	PFTrackTransformer *pfTransformer_;
  	edm::ParameterSet conf_;
  	PFResolutionMap* resMapEtaECAL_;
  	PFResolutionMap* resMapPhiECAL_;	
	 ///Vector of clusters of the PreShower
  	std::vector<reco::PFCluster> ps1Clus;
  	std::vector<reco::PFCluster> ps2Clus;
	
	  edm::ESHandle<TrajectoryFitter> fitter_;
  	///Smoother
  	edm::ESHandle<TrajectorySmoother> smoother_;
	float thr[150];
  	float thrPS[20];
	

};
SoftElecInJetRecoToGen::SoftElecInJetRecoToGen(const edm::ParameterSet& iConfig):
	gedgsfElectronTag_(iConfig.getParameter<edm::InputTag>("gedgsfElectronTag")),
	genParticleTag_(iConfig.getParameter<edm::InputTag>("genParticleTag")),
	PFJetTag_(iConfig.getParameter<edm::InputTag>("PFJetTag")),
	PrimaryVerticesTag_(iConfig.getParameter<edm::InputTag>("PrimaryVerticesTag")),
	TrackingParticleTag_(iConfig.getParameter<edm::InputTag>("TrackingParticleTag")),
	tracksTag_(iConfig.getParameter<edm::InputTag>("tracksTag")),
	minHits(iConfig.getParameter<unsigned int>("minHits")),
	maxNTracks(iConfig.getParameter<unsigned int>("maxNTracks")),
       	maxLIP(iConfig.getParameter<double>("maximumLongitudinalImpactParameter")),	
	vtxReco(new ConfigurableVertexReconstructor(iConfig.getParameter<edm::ParameterSet>("vertexReco"))),
        clusterizer(new TracksClusteringFromDisplacedSeed(iConfig.getParameter<edm::ParameterSet>("clusterizer"))),
	vertexMinAngleCosine(iConfig.getParameter<double>("vertexMinAngleCosine")), //0.98
        vertexMinDLen2DSig(iConfig.getParameter<double>("vertexMinDLen2DSig")), //2.5
        vertexMinDLenSig(iConfig.getParameter<double>("vertexMinDLenSig")), //0.5
	pfTkTransformer_(0),
	pfTransformer_(0),
  	conf_(iConfig),
  	resMapEtaECAL_(0),
  	resMapPhiECAL_(0)
{
	pfTrackLabel_     = iConfig.getParameter<edm::InputTag>( "PFRecTrackLabel" );
  	pfNuclear_        = iConfig.getParameter<edm::InputTag>( "PFNuclear" );
  	gsfTrackLabel_    = iConfig.getParameter<edm::InputTag>( "GsfTrackModuleLabel" );

  	//YM get the input tags for the collections to read from the event
  	tracksTag         = iConfig.getParameter<edm::InputTag>("tracksTag");
  	simtracksTag      = iConfig.getParameter<edm::InputTag>("simtracksTag");

  	tracksContainers_ = iConfig.getParameter< vector < InputTag > >("TkColList");
  	useQuality_       = iConfig.getParameter<bool>("UseQuality");
  	trackQuality_     = TrackBase::qualityByName(iConfig.getParameter<std::string>("TrackQuality"));

  	fitterName_       = iConfig.getParameter<string>("Fitter");
  	smootherName_     = iConfig.getParameter<string>("Smoother");
  	pfCLusTagECLabel_ = iConfig.getParameter<InputTag>("PFEcalClusterLabel");
  	pfCLusTagHCLabel_ = iConfig.getParameter<InputTag>("PFHcalClusterLabel");
  	pfCLusTagPSLabel_ = iConfig.getParameter<InputTag>("PFPSClusterLabel");
  	clusThreshold_    = iConfig.getParameter<double>("ClusterThreshold");
  	minEp_            = iConfig.getParameter<double>("MinEOverP");
  	maxEp_            = iConfig.getParameter<double>("MaxEOverP");
  	minPt_            = iConfig.getParameter<double>("MinPt");
  	maxPt_            = iConfig.getParameter<double>("MaxPt");
  	usePreshower_     = iConfig.getParameter<bool>("UsePreShower");

	gtrf= new GenToRecoFiller("GTRC");
	rtgf= new RecoToGenFiller("RTGC");	
	jpr= new JPsiReco();
}


SoftElecInJetRecoToGen::~SoftElecInJetRecoToGen(){}


//
// member functions
//

// ------------ method called for each event  ------------
void
SoftElecInJetRecoToGen::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

	edm::ESHandle<TrackAssociatorBase> theHitsAssociator;
  	iSetup.get<TrackAssociatorRecord>().get("quickTrackAssociatorByHits", theHitsAssociator);
  	TrackAssociatorBase* associatorByHits = (TrackAssociatorBase *)theHitsAssociator.product();


	//Loading the GED GSF Electron candidates
	edm::Handle<reco::GsfElectronCollection> gedgsfCandidates;
        iEvent.getByLabel(gedgsfElectronTag_, gedgsfCandidates);
        GsfElectronCollection gedgsfCollection = *(gedgsfCandidates.product());

	//Loading the GenParticles
         Handle<GenParticleCollection> GPC;
        iEvent.getByLabel(genParticleTag_, GPC);
        GenParticleCollection gpc = *(GPC.product());

	//Loading the jets
        Handle<PFJetCollection> pfjets;
        iEvent.getByLabel(PFJetTag_, pfjets);
        PFJetCollection pfjc=*(pfjets.product());

	//Loading the PV collection
        edm::Handle<reco::VertexCollection> FullprimaryVertexCollection;
        iEvent.getByLabel(PrimaryVerticesTag_, FullprimaryVertexCollection);
        const reco::VertexCollection pvc = *(FullprimaryVertexCollection.product());

	//Loading the trackingParticle
	edm::Handle<TrackingParticleCollection>  TPCollectionH ;
	iEvent.getByLabel(TrackingParticleTag_,TPCollectionH );
//	const TrackingParticleCollection *tPC   = TPCollectionH.product();


//	Handle<View<Track> > trackCollectionH;
//  	iEvent.getByLabel(tracksTag, trackCollectionH);
//  	const View<Track> tC = *(trackCollectionH.product());

	
	Handle<View<Track> > trackCollectionH;
        iEvent.getByLabel(tracksTag, trackCollectionH);
        const View<Track> tC = *(trackCollectionH.product());


// 	Handle<reco::TrackCollection> theTrackCollection;
//  	iEvent.getByLabel(tracksTag, theTrackCollection);
//  	reco::TrackCollection  Tk=*(theTrackCollection.product());

  	Handle<vector<Trajectory> > tjCollection;
  	iEvent.getByLabel(tracksTag, tjCollection);
  	vector<Trajectory> Tj=*(tjCollection.product());
	
	
	Handle<PFClusterCollection> theECPfClustCollection;
  	iEvent.getByLabel(pfCLusTagECLabel_,theECPfClustCollection);
  	vector<PFCluster> basClus;
  	vector<PFCluster>::const_iterator iklus;
  	for (iklus=theECPfClustCollection.product()->begin();
  	     iklus!=theECPfClustCollection.product()->end();
  	     iklus++){
  	  if((*iklus).energy()>clusThreshold_) basClus.push_back(*iklus);
  	}


	edm::Handle<reco::GsfPFRecTrackCollection> gsfPfRecTracksCollection ;
	iEvent.getByLabel("pfTrackElec",gsfPfRecTracksCollection);
	reco::GsfPFRecTrackCollection  gsfpfTk=*(gsfPfRecTracksCollection.product());

	cout<<"test"<<endl;	
	edm::Handle<reco::GsfTrackCollection> gTCol;
	iEvent.getByLabel( "electronGsfTracks", gTCol);
        GsfTrackCollection gsftc=*(gTCol.product());	
	cout<<"EGSFT ok"<<endl;
	iSetup.get<TrajectoryFitter::Record>().get(fitterName_, fitter_);
  	iSetup.get<TrajectoryFitter::Record>().get(smootherName_, smoother_);

	reco::SimToRecoCollection StRC =  associatorByHits->associateSimToReco(trackCollectionH, TPCollectionH, &iEvent,&iSetup);
	reco::RecoToSimCollection RtSC =  associatorByHits->associateRecoToSim(trackCollectionH, TPCollectionH, &iEvent,&iSetup);

	cout<<"Executing the analysis code"<<endl;


	edm::ESHandle<TransientTrackBuilder> trackBuilder;
	iSetup.get<TransientTrackRecord>().get("TransientTrackBuilder",trackBuilder);
	std::vector<TransientTrack> tts;

	edm::Handle<BeamSpot> beamSpot;
        iEvent.getByLabel("offlineBeamSpot", beamSpot);


	edm::Handle< edm::ValueMap<float> > mvasoftHandle;
  	iEvent.getByLabel("mvaSoft", mvasoftHandle);
	edm::ValueMap<float> mvasoft = *mvasoftHandle;

	for (int u = 0 ; u < (int)gedgsfCollection.size(); ++u)
        {
		rtgf->initRecoToGenFillerObject();
                 //To Estimate Gen->RECO efficiency
		RefToBase<Track> tkRef  = RefToBase<Track>(gedgsfCollection[u].gsfTrack() );
		std::vector<std::pair<TrackingParticleRef, double> > tp;
		TrackingParticleRef tpr;
		if(RtSC.find(tkRef) != RtSC.end()){
			tp = RtSC[tkRef];
			tpr = tp.begin()->first;
			rtgf->tppdgId=tpr->pdgId();
			TrackingParticle::genp_iterator j, b = tpr->genParticle_begin(), e = tpr->genParticle_end();	
			for( j = b; j != e; ++ j ) {
				const reco::GenParticle * p = j->get();
				reco::GenParticle *ncp=(reco::GenParticle*)p;
				rtgf->Vtx=1;
				rtgf->pdgId=ncp->pdgId();
				rtgf->origin=GetOrigin(*ncp);
			}
		}	
		GsfElectronRef elref( gedgsfCandidates, u);	
		cout<<"Testing the MVASoftOutput "<< (mvasoft)[elref]<<endl;
		rtgf->mva_e_pi=(mvasoft)[elref];	
                rtgf->gedgsfElecPt=gedgsfCollection[u].pt();
                rtgf->gedgsfElecEta=gedgsfCollection[u].eta();
                rtgf->gedgsfElecPhi=gedgsfCollection[u].phi();
		//For MVA filling ########################################################################
                rtgf->EtotOvePin=gedgsfCollection[u].eSuperClusterOverP();
                rtgf->EClusOverPout=gedgsfCollection[u].eEleClusterOverPout();
                rtgf->fbrem=gedgsfCollection[u].fbrem()<-1 ? gedgsfCollection[u].fbrem():-1;
		float etot=gedgsfCollection[u].eSuperClusterOverP()*gedgsfCollection[u].trackMomentumAtVtx().R();
		float eEcal=gedgsfCollection[u].eEleClusterOverPout()*gedgsfCollection[u].trackMomentumAtEleClus().R();
		float dP=gedgsfCollection[u].trackMomentumAtVtx().R()-gedgsfCollection[u].trackMomentumAtEleClus().R();
                rtgf->EBremOverDeltaP=(etot-eEcal)/dP ;
                rtgf->logSigmaEtaEta=log(gedgsfCollection[u].sigmaEtaEta());
                rtgf->DeltaEtaTrackEcalSeed=gedgsfCollection[u].deltaEtaEleClusterTrackAtCalo();
                rtgf->HOverE=gedgsfCollection[u].hadronicOverEm();
		rtgf->Chi2GSF=gedgsfCollection[u].gsfTrack()->normalizedChi2()<200?gedgsfCollection[u].gsfTrack()->normalizedChi2():200;
		bool validKF= false; 
  		reco::TrackRef myTrackRef = gedgsfCollection[u].closestCtfTrackRef();
 		validKF = (myTrackRef.isAvailable());
		validKF = (myTrackRef.isNonnull());
		float tempchi2kf=(validKF) ? myTrackRef->normalizedChi2() : 0 ;
		rtgf->CHi2KF=tempchi2kf<10?tempchi2kf:10 ;
		rtgf->nHits=(validKF) ? myTrackRef->hitPattern().trackerLayersWithMeasurement() : -1. ;
		rtgf->SigmaPtOverPt=gedgsfCollection[u].gsfTrack().get()->ptModeError()/gedgsfCollection[u].gsfTrack().get()->ptMode() ;

		rtgf->deta=gedgsfCollection[u].deltaEtaSuperClusterTrackAtVtx()<0.06?gedgsfCollection[u].deltaEtaSuperClusterTrackAtVtx():0.06;
		rtgf->dphi=gedgsfCollection[u].deltaPhiSuperClusterTrackAtVtx()<0.6 ?gedgsfCollection[u].deltaPhiSuperClusterTrackAtVtx():0.6;
		rtgf->detacalo=gedgsfCollection[u].deltaEtaSeedClusterTrackAtCalo();
		rtgf->see=gedgsfCollection[u].sigmaIetaIeta();
		rtgf->etawidth=gedgsfCollection[u].superCluster()->etaWidth();
  		rtgf->phiwidth =  gedgsfCollection[u].superCluster()->phiWidth();
		float temp=(gedgsfCollection[u].e5x5()) !=0. ? 1.-(gedgsfCollection[u].e1x5()/gedgsfCollection[u].e5x5()) : -1. ;
		rtgf->e1x5e5x5=  temp>-1?temp:-1;
		rtgf->e1x5e5x5=  temp<2?temp:2;
		rtgf->R9              =  0;
		rtgf->HoE             =  gedgsfCollection[u].hadronicOverEm();
  		rtgf->EoP             =  gedgsfCollection[u].eSuperClusterOverP()<20?gedgsfCollection[u].eSuperClusterOverP():20;
  		rtgf->IoEmIoP         =  (1.0/gedgsfCollection[u].ecalEnergy()) - (1.0 / gedgsfCollection[u].p());  // in the future to be changed with ele.gsfTrack()->p()
  		rtgf->eleEoPout       =  gedgsfCollection[u].eEleClusterOverPout()<20?gedgsfCollection[u].eEleClusterOverPout():20;
    		//d0
    		if (gedgsfCollection[u].gsfTrack().isNonnull()) {
      			rtgf->d0 = (-1.0)*gedgsfCollection[u].gsfTrack()->dxy(pvc[0].position()); 
    		} else if (gedgsfCollection[u].closestCtfTrackRef().isNonnull()) {
      			rtgf->d0 = (-1.0)*gedgsfCollection[u].closestCtfTrackRef()->dxy(pvc[0].position()); 
    		} else {
      			rtgf->d0 = -9999.0;
    		}
   

 
    		//default values for IP3D
    		if (gedgsfCollection[u].gsfTrack().isNonnull()) {
      			const double gsfsign   = ( (-gedgsfCollection[u].gsfTrack()->dxy(pvc[0].position()))   >=0 ) ? 1. : -1.;
      
      			const reco::TransientTrack tt = trackBuilder->build(gedgsfCollection[u].gsfTrack()); 
      			const std::pair<bool,Measurement1D> &ip3dpv =  IPTools::absoluteImpactParameter3D(tt,pvc[0]);
      			if (ip3dpv.first) {
				double ip3d = gsfsign*ip3dpv.second.value();
				double ip3derr = ip3dpv.second.error();  
				rtgf->ip3d = ip3d; 
        			rtgf->ip3dSig = ip3d/ip3derr;
      			}
    		}
  		
		//#############################################################################################
	
		rtgf->nPV=pvc.size();	
                rtgf->trk_pt_=0;
                rtgf->trk_pTOB_=0;
                rtgf->trk_eta_=0;
                rtgf->trk_phi_=0;
                rtgf->trk_chi2_=0;
                rtgf->trk_ndof_=0;
                rtgf->trk_nchi2_=0;
                rtgf->trk_dpt_=0;
                rtgf->trk_nhits_=0;
                rtgf->trk_nmatched_=0;
                rtgf->trk_quality_=0;
                rtgf->trk_ecalDist_=0;
                rtgf->trk_ecalDeta_=0;
                rtgf->trk_ecalDphi_=0;
                rtgf->ecal_e_=0;
                rtgf->trk_ep_=0;
                rtgf->trk_dptGSF_=0;
                rtgf->trk_chiRatio_=0;
                rtgf->trk_chiReduced_=0;
                rtgf->trk_ecalChi_=0;
                rtgf->trk_epCorr_=0;
                rtgf->ecal_ps_=0;
                rtgf->ecal_ps1e_=0;
                rtgf->ecal_ps2e_=0;
		
		//Filling the TTrees
		rtgf->tree_purity->Fill();
		if(fabs(rtgf->pdgId)==11 && (rtgf->origin==4 || rtgf->origin==6)){
			rtgf->tree_purity_GenElecB->Fill();
			cout<<"An gen elec from B"<<endl;
		}
		else if (fabs(rtgf->pdgId)==211 && (rtgf->origin>0)){
			rtgf->tree_purity_GenPi->Fill();
			cout<<"An gen pion from B"<<endl;
		}
		else if (rtgf->Vtx==1){
			rtgf->tree_purity_GenX->Fill();
			cout<<"An gen X"<<endl;
		}
		else if (fabs(rtgf->tppdgId)==11){
			rtgf->tree_purity_NonGenElec->Fill();
			cout<<"An non gen elec"<<endl;	
		}
		else if (fabs(rtgf->tppdgId)==211){
			rtgf->tree_purity_NonGenPi->Fill();
			cout<<"An non gen pion"<<endl;
		}
		else {
			rtgf->tree_purity_NonGenX->Fill();
			cout<<"An non gen X"<<endl;
		}
		
		if(!(myTrackRef.isAvailable() && myTrackRef.isNonnull()))
			continue;
                if( std::abs(myTrackRef->dz(pvc[0].position())) > maxLIP)
			continue;
		TransientTrack tt = trackBuilder->build(gedgsfCollection[u].closestCtfTrackRef());
		tt.setBeamSpot(*beamSpot);
		//cout<<"pushing back a track "<<tt.track().pt()<<" "<<tt.track().eta()<<" "<<tt.track().phi()<<endl;
		tts.push_back(tt);
        }

	if(tts.size()>=2){
		jpr->initJPsiRecoFillerObject();
//		cout<<"--->event with "<<tts.size()<<" transient tracks"<<endl;		
		std::vector<TransientVertex> VCandidates;
		VCandidates=GetVCandidates(tts);	
//		cout<<"--->Combining the transient tracks two by two, we have formed "<<VCandidates.size()<<" vertices"<<endl;
		for(std::vector<TransientVertex>::const_iterator v = VCandidates.begin();v != VCandidates.end(); ++v) {
			VertexDistance3D vdist;
                	VertexDistanceXY vdist2d;
			Measurement1D dlen= vdist.distance(pvc[0],*v);
                        Measurement1D dlen2= vdist2d.distance(pvc[0],*v);
			GlobalVector dir;  
			GlobalPoint ppv(pvc[0].position().x(),pvc[0].position().y(),pvc[0].position().z());
                        GlobalPoint sv((*v).position().x(),(*v).position().y(),(*v).position().z());
			std::vector<reco::TransientTrack> ts = v->originalTracks();
                //        std::cout<<"ok, ts.size="<<ts.size()<<std::endl;
                        for(std::vector<reco::TransientTrack>::iterator i = ts.begin();i != ts.end(); ++i) {
                        	reco::TrackRef t = i->trackBaseRef().castTo<reco::TrackRef>();
                                float w = v->trackWeight(*i);
                                if (w > 0.5) dir+=i->impactPointState().globalDirection();
                        }
			float vscal = dir.unit().dot((sv-ppv).unit()) ;
		//	cout<<"-------> "<<dlen.significance()<<" "<<vscal<<" "<<v->normalisedChiSquared()<<" "<<dlen2.significance()<<": Passing the cuts???"<<endl;
			if(dlen.significance() > vertexMinDLenSig  && vscal > vertexMinAngleCosine &&  v->normalisedChiSquared() < 10 && dlen2.significance() > vertexMinDLen2DSig){
				reco::Vertex vv(*v);
                                TLorentzVector TLVtrack[2];
                                int TLVindex=0;
		//		cout<<"----------> Yes"<<endl;
                                if(vv.nTracks()==2){
                  //              	cout<<"---------------> Two tracks"<<endl;
                                        int charge=0;
                                        for (reco::Vertex::trackRef_iterator itTr = vv.tracks_begin(); itTr != vv.tracks_end(); ++itTr) {
                                        	charge+=(*itTr)->charge();
                    //                            cout<<"Using a track "<<(*itTr)->pt()<<" "<<(*itTr)->eta()<<" "<<(*itTr)->phi()<<endl;
                                                TLVtrack[TLVindex].SetPtEtaPhiM((*itTr)->pt(),(*itTr)->eta(),(*itTr)->phi(),0.0);       
                                                TLVindex++;
                                        }               
                                        if(charge==0){
                      //                  	cout<<"-----------------> A JPsi/Upsilon candidate, recording for later...(M="<<(TLVtrack[0]+TLVtrack[1]).M()<<")"<<endl; 
                                                jpr->Mass=(TLVtrack[0]+TLVtrack[1]).M();
                                                jpr->tree_JpsiUpsilon->Fill();
                                        }
                                                
                               }
			}
		}	
	}	
	for (int j = 0 ; j < (int)gpc.size(); ++j)
        {
                gtrf->initGenToRecoFillerObject();
        //Consider generated electrons, pions and Kaons 
                bool isFinal=gpc[j].status()==1;
                bool isElecOrPionOrKaon=fabs(gpc[j].pdgId())==11 || fabs(gpc[j].pdgId())==211 || fabs(gpc[j].pdgId())==321;
                //bool isElec=fabs(gpc[j].pdgId())==11;
                bool isKineOk=gpc[j].pt()>0.8 && fabs(gpc[j].eta())<2.4;
                if(isFinal && isElecOrPionOrKaon && isKineOk){
                        gtrf->pdgId=gpc[j].pdgId();
                        gtrf->origin=GetOrigin(gpc[j]);
                        gtrf->ptGen=gpc[j].pt();
                        gtrf->etaGen=gpc[j].eta();
                        gtrf->phiGen=gpc[j].phi();
                        gtrf->pGen=gpc[j].p();
                        //Running over jets
                        bool inRecoJet=false;
                        for (int k = 0 ; k < (int)pfjc.size(); ++k)
                        {
                                //Ask for reasonable kinematics
                                float dR=deltaR(gpc[j].eta(), gpc[j].phi(), pfjc[k].eta(), pfjc[k].phi());
                                if(pfjc[k].pt()>20 && fabs(pfjc[k].eta())<2.4 && dR<0.5)inRecoJet=true;
                        }
                        gtrf->inRecoJet=inRecoJet;
                                        //Preparing information about track and trackingParticle        
                        Track track;
                        float sharedHits = 0;
                        float dR = 100;
                        int indexTrack;
                        RefToBase<Track> tr;
                        bool result = findMatch(gpc[j], TPCollectionH, StRC,track, indexTrack, dR, sharedHits,tr);
                        if ( result ) gtrf->gen_match_ = 1;
                        else          gtrf->gen_match_ = 0;

                        // Trying to find the proper seed and gedgsfelectron
                        for (int u = 0 ; u < (int)gsfpfTk.size(); ++u)
                        {
                                if(tr.get()==gsfpfTk[u].gsfTrackRef().get()){
                                        gtrf->isMatchedWithASeed=1;
                                        gtrf->gsftrkSeedPt=gsfpfTk[u].gsfTrackRef().get()->pt();
                                        gtrf->gsftrkSeedEta=gsfpfTk[u].gsfTrackRef().get()->eta();
                                        gtrf->gsftrkSeedPhi=gsfpfTk[u].gsfTrackRef().get()->phi();
                                }
                        }
                        for (int u = 0 ; u < (int)gedgsfCollection.size(); ++u)
                        {
                                if(tr.get()==gedgsfCollection[u].gsfTrack().get()){
                                        gtrf->isMatchedWithAGedGsfElec=1;
                                        gtrf->gedgsfElecPt=gedgsfCollection[u].pt();
                                        gtrf->gedgsfElecEta=gedgsfCollection[u].eta();
                                        gtrf->gedgsfElecPhi=gedgsfCollection[u].phi();
					GsfElectronRef elref( gedgsfCandidates, u);
                			cout<<"Testing the MVASoftOutput "<< (mvasoft)[elref]<<endl;
                			gtrf->mva_e_pi=(mvasoft)[elref];

                                }
                        }
                        gtrf->nPV=pvc.size();
                        gtrf->tree_efficiency->Fill();
                }
         }
}


std::vector<TransientVertex> SoftElecInJetRecoToGen::GetVCandidates(std::vector<TransientTrack> tts){
	std::vector<TransientVertex> result;
	int index[100];
	for(int u=0;u<100;u++)index[u]=0;
	for(int u=1;u<(int)tts.size();u++)
	for(int v=0;v<u;v++){
//		cout<<"----->checking track "<<u<<" and "<<v<<endl;	
		std::vector<reco::TransientTrack> tracksToVertex;
		tracksToVertex.push_back(tts[u]);
		tracksToVertex.push_back(tts[v]);
		double sigmacut = 3.0;
                double Tini = 256.;
                double ratio = 0.25;
		AdaptiveVertexFitter theAdaptiveFitter(
                	GeometricAnnealing(sigmacut, Tini, ratio),
                	DefaultLinearizationPointFinder(),
                	KalmanVertexUpdator<5>(),
                	KalmanVertexTrackCompatibilityEstimator<5>(),
                	KalmanVertexSmoother() 
                );
		TransientVertex theVertex=theAdaptiveFitter.vertex(tracksToVertex);
//		cout<<"------>"<<theVertex.isValid()<<" "<<theVertex.totalChiSquared()<<" "<<theVertex.degreesOfFreedom() <<endl;
		if(theVertex.isValid() && theVertex.totalChiSquared() >= 0. && theVertex.degreesOfFreedom() > 0){
			if (theVertex.totalChiSquared() / theVertex.degreesOfFreedom() < 3.) {
				result.push_back(theVertex);	
				index[u]++;
				index[v]++;
		//		if(index[u]>1 || index[v]>1)cout<<"################### CAREFUL, tracks already used for another vertex ################################"<<endl;
				
			}
		}
	}
	return result;
	
};

bool SoftElecInJetRecoToGen::isInnerMost(const reco::Track track1, const reco::Track track2, bool& sameLayer) {
  // copied from http://cmssw.cvs.cern.ch/cgi-bin/cmssw.cgi/CMSSW/RecoParticleFlow/PFTracking/plugins/PFElecTkProducer.cc?revision=1.23&view=markup
  reco::HitPattern hitPattern1 = track1.hitPattern();
  reco::HitPattern hitPattern2 = track2.hitPattern();

  // retrieve the first valid hit
  int hitCounter1 = 0;
  trackingRecHit_iterator hitsIt1;
  for(hitsIt1 = track1.recHitsBegin(); hitsIt1 != track1.recHitsEnd(); ++hitsIt1, ++hitCounter1)
    { if (((**hitsIt1).isValid())) break; }
  int hitCounter2 = 0;
  trackingRecHit_iterator hitsIt2;
  for(hitsIt2 = track2.recHitsBegin(); hitsIt2 != track2.recHitsEnd(); ++hitsIt2, ++hitCounter2)
    { if (((**hitsIt2).isValid())) break; }
  uint32_t hit1 = hitPattern1.getHitPattern(hitCounter1);
  uint32_t hit2 = hitPattern2.getHitPattern(hitCounter2);
  if ( hitPattern1.getSubStructure(hit1) != hitPattern2.getSubStructure(hit2) )
    return hitPattern2.getSubStructure(hit2) < hitPattern1.getSubStructure(hit1);
  else if ( hitPattern1.getLayer(hit1) != hitPattern2.getLayer(hit2) )
    return hitPattern2.getLayer(hit2) < hitPattern1.getLayer(hit1);
  else {
    sameLayer = true;
    return false;
  }
}


int SoftElecInJetRecoToGen::GetOrigin(reco::GenParticle& part){
        int isBanAncestor=0;
        int isDanAncestor=0;
        int isVanAncestor=0;
        if(part.numberOfMothers()!=0 ){
                const reco::Candidate *m1=part.mother(0);
                if(isaBhadron(fabs(m1->pdgId())))isBanAncestor=1;
                if(isaDhadron(fabs(m1->pdgId())))isDanAncestor=1;
                if(isaV(fabs(m1->pdgId())))isVanAncestor=1;
                if(m1->numberOfMothers()!=0 ){
                       const reco::Candidate * m2=m1->mother(0);
                        if(isaBhadron(fabs(m2->pdgId())))isBanAncestor=1;
                        if(isaDhadron(fabs(m2->pdgId())))isDanAncestor=1;
                        if(isaV(fabs(m2->pdgId())))isVanAncestor=1;
                        if(m2->numberOfMothers()!=0 ){
                                const reco::Candidate *m3=m2->mother(0);
                                if(isaBhadron(fabs(m3->pdgId())))isBanAncestor=1;
                                if(isaDhadron(fabs(m3->pdgId())))isDanAncestor=1;
                                if(isaV(fabs(m3->pdgId())))isVanAncestor=1;
                                if(m3->numberOfMothers()!=0 ){
                                        const reco::Candidate *m4=m3->mother(0);
                                        if(isaBhadron(fabs(m4->pdgId())))isBanAncestor=1;
                                        if(isaDhadron(fabs(m4->pdgId())))isDanAncestor=1;
                                        if(isaV(fabs(m4->pdgId())))isVanAncestor=1;
                                        if(m4->numberOfMothers()!=0 ){
                                                const reco::Candidate *m5=m4->mother(0);
                                                if(isaBhadron(fabs(m5->pdgId())))isBanAncestor=1;
                                                if(isaDhadron(fabs(m5->pdgId())))isDanAncestor=1;
                                                if(isaV(fabs(m5->pdgId())))isVanAncestor=1;
                                                if(m5->numberOfMothers()!=0 ){
                                                        const reco::Candidate *m6=m5->mother(0);
                                                        if(isaBhadron(fabs(m6->pdgId())))isBanAncestor=1;
                                                        if(isaDhadron(fabs(m6->pdgId())))isDanAncestor=1;
                                                        if(isaV(fabs(m6->pdgId())))isVanAncestor=1;
                                                }
                                        }

                                }
                        }

                }
        }
        int result=4*isBanAncestor+2*isDanAncestor+isVanAncestor;
        return result;
}

bool SoftElecInJetRecoToGen::isaV(int pidAbs){
        int res=false;
        if(pidAbs==24 || pidAbs==23 || pidAbs==22)res=true;
        return res;
}


bool SoftElecInJetRecoToGen::isaBhadron(int pidAbs){

        bool isB = false;
        if(pidAbs>500){

                pidAbs/= 100;

                if(pidAbs<60 && pidAbs>50)pidAbs/= 10;
                int mod10 = pidAbs % 5;

                if(mod10 == 0) {
                        isB = true;
                }
        }
        else  {
          //              std::cout<<"PID too low to be a B meson"<<std::endl;
        }
        return isB;
}

bool SoftElecInJetRecoToGen::isaDhadron(int pidAbs){

        bool isD = false;
        if(pidAbs>400){
                pidAbs/= 100;
                if(pidAbs<50 && pidAbs>40)pidAbs/= 10;
                int mod10 = pidAbs % 4;

                if(mod10 == 0 && pidAbs<10) {
                        isD = true;
                }
        }
        else{
                //              std::cout<<"PID too low to be a B meson"<<std::endl;
        }
        return isD;
}


bool SoftElecInJetRecoToGen::trackFilter(const reco::TrackRef &track) const
{
	if (track->hitPattern().numberOfValidHits() < (int)minHits)
//	if (track->hitPattern().trackerLayersWithMeasurement() < (int)minHits)
		return false;
	if (track->pt() < minPt )
		return false;
 
	return true;
}

bool SoftElecInJetRecoToGen::findMatch(const reco::Candidate& genCand,
                              Handle<TrackingParticleCollection> TPCollectionH,
                              reco::SimToRecoCollection q,
                              reco::Track& track,
                              int& indexTrack,
                              float& dR,
                              float& sharedHits,
			      RefToBase<Track> &trRefToBase) {

  const TrackingParticleCollection tPC = *(TPCollectionH.product());
  if ( tPC.size() == 0 ) return false;

  unsigned int nTPs = tPC.size();
  for(unsigned int iTP = 0; iTP < nTPs; ++iTP) {
    const TrackingParticle& lTP(tPC[iTP]);
    float deta = lTP.eta() - genCand.eta();
    float dphi = acos(cos(lTP.phi() - genCand.phi()));
    dR = sqrt(deta*deta + dphi*dphi);

    if ( /*fabs(lTP.pt() - genCand.pt()) > 0.05 || */lTP.pdgId() != genCand.pdgId() || dR > 0.05 ) continue;
    edm::Ref<TrackingParticleCollection> tp(TPCollectionH, iTP);
    try {
      	std::vector<std::pair<RefToBase<Track>, double> > trackV = q[tp];
      if ( trackV.size() == 1 ) {
        std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
        RefToBase<Track> tr = it->first;
        indexTrack = tr.key();
        track = *tr;
	trRefToBase=tr;
        sharedHits = it->second;
        return true;
      }
      if ( trackV.size() > 1 ) {
        std::vector<std::pair<RefToBase<Track>, double> >::const_iterator it = trackV.begin();
        RefToBase<Track> tr = it->first;
        Track track1 = *tr;
        int indexTrack1 = tr.key();
        float nHitsRatio1 = it->second;
        ++it; tr = it->first;
        int indexTrack2 = tr.key();
        Track track2 = *tr;
        float nHitsRatio2 = it->second;
        bool sameLayer = false;
        bool result = isInnerMost(track1, track2, sameLayer);
        Track finalTrack;
        if ( !sameLayer ) {
          if ( result ) {
            finalTrack = track2;
            sharedHits = nHitsRatio2;
            indexTrack = indexTrack2;
          }
          else {
            //cout << "track1 is inner most, with pT = " << track1.pt() << '\t' << nHitsRatio1 << endl;
            finalTrack = track1;
            sharedHits = nHitsRatio1;
            indexTrack = indexTrack1;
          }
        } else {
          //cout << "track1 and track2 are in the same layer!" << endl;
          if ( nHitsRatio2 > nHitsRatio1 ) {
            // use track2
            finalTrack = track2;
            sharedHits = nHitsRatio2;
            indexTrack = indexTrack2;
          }
          else {
            // use track1
            finalTrack = track1;
            sharedHits = nHitsRatio1;
            indexTrack = indexTrack1;
          }
        }
        track = finalTrack;
        return true;
      }
    } catch (Exception event) {
      // no track is found to be matched, do nothing
      return false;
    }
  }
  return false;
}



double SoftElecInJetRecoToGen::testPreshowerDistance(PFCluster eeclus,
                                            PFCluster psclus) {

  const reco::PFCluster::REPPoint& pspos = psclus.positionREP();
  const reco::PFCluster::REPPoint& eepos = eeclus.positionREP();
  // same side of the detector?
  if ( eeclus.z()*psclus.z() < 0 ) return -1.0;

  const double dphi = std::abs(TVector2::Phi_mpi_pi(eepos.phi() - pspos.phi()));
  if ( dphi > 0.6 ) return -1.0;
  const double deta = std::abs(eepos.eta() - pspos.eta());
  if ( deta > 0.3 ) return -1.0;
  return LinkByRecHit::testECALAndPSByRecHit(eeclus, psclus, false);
}

// ------------ method called once each job just before starting event loop  ------------
void 
SoftElecInJetRecoToGen::beginJob()
{
}
	

// ------------ method called once each job just after ending the event loop  ------------
void 
SoftElecInJetRecoToGen::endJob() 
{
	gtrf->WriteInFileAndCloseIt();
	rtgf->WriteInFileAndCloseIt();
	jpr->WriteInFileAndCloseIt();
}

// ------------ method called when starting to processes a run  ------------

void SoftElecInJetRecoToGen::beginRun(const edm::Run& run, const edm::EventSetup & es)
{

  //Magnetic Field
  ESHandle<MagneticField> magneticField;
  es.get<IdealMagneticFieldRecord>().get(magneticField);
  B_=magneticField->inTesla(GlobalPoint(0,0,0));

  pfTransformer_= new PFTrackTransformer(B_);
  pfTransformer_->OnlyProp();

  pfTkTransformer_= new PFTrackTransformer(math::XYZVector(magneticField->inTesla(GlobalPoint(0,0,0))));
  pfTkTransformer_->OnlyProp();

  //Resolution maps
  FileInPath ecalEtaMap(conf_.getParameter<string>("EtaMap"));
  FileInPath ecalPhiMap(conf_.getParameter<string>("PhiMap"));
  resMapEtaECAL_ = new PFResolutionMap("ECAL_eta",ecalEtaMap.fullPath().c_str());
  resMapPhiECAL_ = new PFResolutionMap("ECAL_phi",ecalPhiMap.fullPath().c_str());

  //read threshold
  FileInPath parFile(conf_.getParameter<string>("ThresholdFile"));
  ifstream ifs(parFile.fullPath().c_str());
  for (int iy=0;iy<72;iy++) ifs >> thr[iy];

  //read PS threshold
  FileInPath parPSFile(conf_.getParameter<string>("PSThresholdFile"));
  ifstream ifsPS(parPSFile.fullPath().c_str());
  for (int iy=0;iy<12;iy++) ifsPS >> thrPS[iy];

}


// ------------ method called when ending the processing of a run  ------------
/*
void 
SoftElecInJetRecoToGen::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
SoftElecInJetRecoToGen::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
SoftElecInJetRecoToGen::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
SoftElecInJetRecoToGen::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(SoftElecInJetRecoToGen);
