import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")


process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("FWCore.MessageService.MessageLogger_cfi")
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Reconstruction_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('TrackingTools.GsfTracking.FwdAnalyticalPropagator_cfi')
process.load('RecoParticleFlow.PFTracking.trackerDrivenElectronSeeds_cfi')
process.load('SimTracker.TrackAssociation.TrackAssociatorByHits_cfi')
#from SimTracker.TrackAssociation.TrackAssociatorByHits_cfi import *
process.TrackAssociatorByHits.Quality_SimToReco = cms.double(0.01)

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(300) )

# Other statements
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:startup', '')


import FWCore.Utilities.FileUtils as FileUtils
mylist = FileUtils.loadListFromFile ('list2') 
readFiles = cms.untracked.vstring( *mylist)
process.source = cms.Source('PoolSource', fileNames = readFiles)

process.load("PhysicsTools.JetMCAlgos.CaloJetsMCFlavour_cfi")
process.AK5PFbyRef=process.AK5byRef.clone(jets="ak5PFJets")
process.AK5PFbyValAlgo=process.AK5byValAlgo.clone(srcByReference="AK5PFbyRef")

process.demo = cms.EDAnalyzer('SoftElecInJetRecoToGen',

  JetFTag=cms.InputTag("AK5PFbyValAlgo"),
 gedgsfElectronTag= cms.InputTag("gedGsfElectrons"),
 genParticleTag=cms.InputTag("genParticles"),
 PFJetTag=cms.InputTag("ak5PFJets"), 
 PrimaryVerticesTag=cms.InputTag("offlinePrimaryVertices"),	 
 simtracksTag = cms.InputTag("g4SimHits"),
 tracksTag = cms.InputTag("electronGsfTracks"),
# tracksTag = cms.InputTag("gsfTracks"),
 TrackingParticleTag=cms.InputTag("mix","MergedTrackTruth"),
 EtaMap = cms.string('RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_eta.dat'),
 PhiMap = cms.string('RecoParticleFlow/PFBlockProducer/data/resmap_ECAL_phi.dat'),
 HcalWindow=cms.double(0.184),
 GsfTrackModuleLabel = cms.InputTag("gsfElectrons"),
 PFRecTrackLabel = cms.InputTag("particleFlowRecHitECAL"),
 PFNuclear = cms.InputTag("particleFlowBlock"),
 useNuclear = cms.bool(False),
 TkColList = cms.VInputTag(cms.InputTag("generalTracks")),
 UseQuality = cms.bool(True),
 TrackQuality = cms.string('highPurity'),
 Smoother = cms.string('GsfTrajectorySmoother_forPreId'),
 Fitter = cms.string('GsfTrajectoryFitter_forPreId'),
 PFEcalClusterLabel = cms.InputTag("particleFlowClusterECAL"),
 PFHcalClusterLabel = cms.InputTag("particleFlowClusterHCAL"),
 PSThresholdFile = cms.string('RecoParticleFlow/PFTracking/data/PSThreshold.dat'),
 ClusterThreshold = cms.double(0.5),
 MinEOverP = cms.double(0.3),
 MaxEOverP = cms.double(3.0),
 MinPt = cms.double(2.0),
 MaxPt = cms.double(50.0),
 UsePreShower =cms.bool(False),
 ApplyIsolation = cms.bool(False),
 EcalStripSumE_deltaPhiOverQ_minValue = cms.double(-0.1),
 EcalStripSumE_minClusEnergy = cms.double(0.1),
 EcalStripSumE_deltaEta = cms.double(0.03),
 EcalStripSumE_deltaPhiOverQ_maxValue = cms.double(0.5),
 PFPSClusterLabel = cms.InputTag("particleFlowClusterPS"),
 EOverPLead_minValue = cms.double(0.95),
 HOverPLead_maxValue = cms.double(0.05),
 ThresholdFile = cms.string('RecoParticleFlow/PFTracking/data/Threshold.dat'),

 clusterizer = cms.PSet(
           seedMin3DIPSignificance = cms.double(1.2),
           seedMin3DIPValue = cms.double(0.005),
           clusterMaxDistance = cms.double(0.05), #500um
           clusterMaxSignificance = cms.double(4.5), #4.5 sigma
           clusterScale = cms.double(1), 
           clusterMinAngleCosine = cms.double(0.5), # only forward decays
       ),

 minHits = cms.uint32(8),
 maximumLongitudinalImpactParameter = cms.double(0.3),
 minPt = cms.double(0.8),
 maxNTracks = cms.uint32(30),
 vertexMinAngleCosine = cms.double(0.95), # scalar prod direction of tracks and flight dir
 vertexMinDLen2DSig = cms.double(2.5), #2.5 sigma
 vertexMinDLenSig = cms.double(0.5), #0.5 sigma
 vertexReco = cms.PSet(
               finder = cms.string('avr'),
               primcut = cms.double(1.0),
               seccut = cms.double(3),
               smoothing = cms.bool(True)
       )


)

#process.load("SimGeneral.TrackingAnalysis.simHitTPAssociation_cfi");

process.reconstruction_step = cms.Path(process.reconstruction_fromRECO)
#process.endjob_step = cms.EndPath(process.endOfProcess)
#
process.load("RecoParticleFlow.PFTracking.particleFlowTrack_cff")
process.make_pftracks = cms.Path(process.pfTrackingGlobalReco)
#
#process.p = cms.Path(process.simHitTPAssocProducer*process.demo)

process.load('SoftElectronInJetAnalyzer.SoftElecInJet.SoftelectronIdMVAProducer_cfi')
process.softeidMVASequence = cms.Sequence(  process.mvaSoft )

########################

import SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi 
from SimTracker.TrackerHitAssociation.clusterTpAssociationProducer_cfi import *
process.load("SimTracker.TrackAssociation.quickTrackAssociatorByHits_cfi")
process.quickTrackAssociatorByHits.useClusterTPAssociation = cms.bool(True)
process.load("SimTracker.TrackerHitAssociation.clusterTpAssociationProducer_cfi")

process.test = cms.Sequence(
     process.tpClusterProducer #*
#     process.quickTrackAssociatorByHits

)

process.dump=cms.EDAnalyzer('EventContentAnalyzer')
process.p = cms.Path(  
  process.myPartons *
  process.AK5PFbyRef *
  process.AK5PFbyValAlgo *
  process.AK5Flavour *
  process.softeidMVASequence*
  process.demo
)

######################

#process.RECODEBUGoutput = cms.OutputModule("PoolOutputModule",
#    splitLevel = cms.untracked.int32(0),
#    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
#    outputCommands = process.RECODEBUGEventContent.outputCommands,
#    fileName = cms.untracked.string('output.root'),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string(''),
#        dataTier = cms.untracked.string('RECODEBUG')
#    ),
##    SelectEvents = cms.untracked.PSet(
##        SelectEvents = cms.vstring('generation_step')
##    )
#)

#process.RECODEBUGoutput_step = cms.EndPath(process.RECODEBUGoutput)


# Schedule definition
process.schedule = cms.Schedule(process.reconstruction_step,
#				process.tpClusterProducer,
                                process.make_pftracks,
                                 process.p#,
			#	 process.RECODEBUGoutput_step
                                )
process.gedGsfElectronsTmp.PreSelectMVA = -2


