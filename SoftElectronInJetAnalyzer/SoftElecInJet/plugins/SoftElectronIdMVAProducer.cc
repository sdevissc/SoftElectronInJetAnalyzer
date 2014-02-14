// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SoftElectronInJetAnalyzer/SoftElecInJet/interface/EGammaMvaSoftEleEstimator.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"
//
// class declaration
//

using namespace std;
using namespace reco;
class SoftElectronIdMVAProducer : public edm::EDFilter {
	public:
		explicit SoftElectronIdMVAProducer(const edm::ParameterSet&);
		~SoftElectronIdMVAProducer();

	private:
		virtual bool filter(edm::Event&, const edm::EventSetup&);

		// ----------member data ---------------------------
                bool verbose_;
		edm::InputTag vertexTag_;
		edm::InputTag electronTag_;
                edm::InputTag reducedEBRecHitCollection_;
                edm::InputTag reducedEERecHitCollection_;
  
                double _Rho;
                string method_;
                vector<string> mvaWeightFiles_;
                bool Trig_;
                bool NoIP_;
 
                EGammaMvaSoftEleEstimator* mvaID_;

};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
SoftElectronIdMVAProducer::SoftElectronIdMVAProducer(const edm::ParameterSet& iConfig) {
        verbose_ = iConfig.getUntrackedParameter<bool>("verbose", false);
	electronTag_ = iConfig.getParameter<edm::InputTag>("electronTag");
	method_ = iConfig.getParameter<string>("method");
	std::vector<string> fpMvaWeightFiles = iConfig.getParameter<std::vector<std::string> >("mvaWeightFile");

        produces<edm::ValueMap<float> >("");

        mvaID_ = new EGammaMvaSoftEleEstimator();
 
        bool manualCat_ = true;

	string path_mvaWeightFileEleID;
	for(unsigned ifile=0 ; ifile < fpMvaWeightFiles.size() ; ++ifile) {
	  path_mvaWeightFileEleID = edm::FileInPath ( fpMvaWeightFiles[ifile].c_str() ).fullPath();
//	  path_mvaWeightFileEleID = fpMvaWeightFiles[ifile].c_str();
	  mvaWeightFiles_.push_back(path_mvaWeightFileEleID);
	}
	
        mvaID_->initialize(method_, manualCat_, mvaWeightFiles_);

}


SoftElectronIdMVAProducer::~SoftElectronIdMVAProducer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called on each new Event  ------------
bool SoftElectronIdMVAProducer::filter(edm::Event& iEvent, const edm::EventSetup& iSetup) {
	using namespace edm;

        std::auto_ptr<edm::ValueMap<float> > out(new edm::ValueMap<float>() );

	Handle<reco::GsfElectronCollection> egCollection;
	iEvent.getByLabel(electronTag_,egCollection);
        const reco::GsfElectronCollection egCandidates = (*egCollection.product());

	_Rho=0;
	edm::Handle<double> rhoPtr;
	const edm::InputTag eventrho("kt6PFJets", "rho");
	iEvent.getByLabel(eventrho,rhoPtr);
	_Rho=*rhoPtr;

        std::vector<float> values;
        values.reserve(egCollection->size());
   
        for ( reco::GsfElectronCollection::const_iterator egIter = egCandidates.begin(); egIter != egCandidates.end(); ++egIter) {

          double mvaVal = -999999;
          mvaVal = mvaID_->mvaValue( *egIter,iEvent, verbose_);
	  values.push_back( mvaVal ); 
	}

        edm::ValueMap<float>::Filler filler(*out);
        filler.insert(egCollection, values.begin(), values.end() );
	filler.fill();

	iEvent.put(out);

	return true;
}

//define this as a plug-in
DEFINE_FWK_MODULE(SoftElectronIdMVAProducer);
