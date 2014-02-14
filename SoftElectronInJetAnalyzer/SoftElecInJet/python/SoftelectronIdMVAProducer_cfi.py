import FWCore.ParameterSet.Config as cms

mvaSoft = cms.EDFilter("SoftElectronIdMVAProducer",
                         verbose = cms.untracked.bool(False),
                         electronTag = cms.InputTag('gedGsfElectrons'),
                         method = cms.string("BDT"),
			mvaWeightFile = cms.vstring(
                                "SoftElectronInJetAnalyzer/SoftElecInJet/data/TMVA_BDTSoftElectrons_27Jan2014.weights.xml"
				),
			)
				


