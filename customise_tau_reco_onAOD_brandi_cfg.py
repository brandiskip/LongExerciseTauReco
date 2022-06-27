import FWCore.ParameterSet.Config as cms

'''
# trick to make python know where to look for the imports
# import sys
# sys.path.append('..')

# for example: here
# from rerunTauRecoOnMiniAOD import process
'''
def addTauReReco(process):
        process.load('PhysicsTools.PatAlgos.producersLayer1.tauProducer_cff')
        process.load('PhysicsTools.PatAlgos.selectionLayer1.selectedPatCandidates_cff')
        process.load("RecoTauTag.Configuration.RecoPFTauTag_cff")
        # process.ptau = cms.Path(process.PFTau)
        process.PATTauSequence = cms.Sequence(process.PFTau+process.makePatTaus+process.selectedPatTaus)
        process.TauReco = cms.Path(process.PATTauSequence)
        # process.task.add(process.TauReco)
process = cms.Process("TAURECO")
process.load("Configuration.StandardSequences.MagneticField_cff") # for CH reco
process.load("Configuration.Geometry.GeometryRecoDB_cff")




# runSignal = False # Set to False to read in QCD file instead of ZTT

maxEvents = -1

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource", fileNames=readFiles, secondaryFileNames=secFiles)


addTauReReco(process)

## to work on stau sample
process.tauGenJets.GenParticles = cms.InputTag("genParticlePlusGeant")
process.tauMatch. matched = cms.InputTag("genParticlePlusGeant")
process.tauMatch.mcStatus = cms.vint32(8)



# process.GlobalTag.globaltag = 'auto:phase1_2021_realistic'
process.GlobalTag = cms.ESSource("PoolDBESSource",
    DBParameters = cms.PSet(
        authenticationPath = cms.untracked.string('.'),
        connectionRetrialPeriod = cms.untracked.int32(10),
        connectionRetrialTimeOut = cms.untracked.int32(60),
        connectionTimeOut = cms.untracked.int32(0),
        enableConnectionSharing = cms.untracked.bool(True),
        enablePoolAutomaticCleanUp = cms.untracked.bool(False),
        enableReadOnlySessionOnUpdateConnection = cms.untracked.bool(False),
        idleConnectionCleanupPeriod = cms.untracked.int32(10),
        messageLevel = cms.untracked.int32(0)
    ),
    DumpStat = cms.untracked.bool(False),
    ReconnectEachRun = cms.untracked.bool(False),
    RefreshAlways = cms.untracked.bool(False),
    RefreshEachRun = cms.untracked.bool(False),
    RefreshOpenIOVs = cms.untracked.bool(False),
    connect = cms.string('frontier://FrontierProd/CMS_CONDITIONS'),
    globaltag = cms.string('auto:phase1_2021_realistic'),
    pfnPostfix = cms.untracked.string('None'),
    snapshotTime = cms.string(''),
    toGet = cms.VPSet()
)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:phase1_2021_realistic', '')


readFiles.extend([
    'root://cms-xrd-global.cern.ch//store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC/crab_gmsb_m200_100mm_1230pre3_v33_SingleTauL1Seed_TauForPixelPt30_L2Filter30_eta2p2_dEta03_pt1p2_onAOD_v4/220616_081141/0000/outputHLT_7.root',
    'root://cms-xrd-global.cern.ch//store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC/crab_gmsb_m200_100mm_1230pre3_v33_SingleTauL1Seed_TauForPixelPt30_L2Filter30_eta2p2_dEta03_pt1p2_onAOD_v4/220616_081141/0000/outputHLT_8.root',
    'root://cms-xrd-global.cern.ch//store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC/crab_gmsb_m200_100mm_1230pre3_v33_SingleTauL1Seed_TauForPixelPt30_L2Filter30_eta2p2_dEta03_pt1p2_onAOD_v4/220616_081141/0000/outputHLT_9.root',
    'root://cms-xrd-global.cern.ch//store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC/crab_gmsb_m200_100mm_1230pre3_v33_SingleTauL1Seed_TauForPixelPt30_L2Filter30_eta2p2_dEta03_pt1p2_onAOD_v4/220616_081141/0000/outputHLT_10.root',
    ])

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( maxEvents )
)

process.combinatoricRecoTaus.builders[0].qualityCuts.signalQualityCuts.maxDeltaZ = cms.double(100.)
process.combinatoricRecoTaus.builders[0].qualityCuts.signalQualityCuts.maxTrackChi2 = cms.double(1000.)
process.combinatoricRecoTaus.builders[0].qualityCuts.signalQualityCuts.maxTransverseImpactParameter = cms.double(100.)
process.combinatoricRecoTaus.builders[0].qualityCuts.pvFindingAlgo = cms.string('highestPtInEvent')

# chargedPFCandidates from PFChargedHadrons ... well this is obvious
process.ak4PFJetsRecoTauChargedHadrons.builders[0].qualityCuts.signalQualityCuts.maxDeltaZ = cms.double(100.)
process.ak4PFJetsRecoTauChargedHadrons.builders[0].qualityCuts.signalQualityCuts.maxTrackChi2 = cms.double(1000.)
process.ak4PFJetsRecoTauChargedHadrons.builders[0].qualityCuts.signalQualityCuts.maxTransverseImpactParameter = cms.double(100.)
process.ak4PFJetsRecoTauChargedHadrons.builders[0].qualityCuts.pvFindingAlgo = cms.string('highestPtInEvent') # closest in dz makes no sense for displaced stuff

# chargedPFCandidates from lostTracks ... this collections exists in miniAODs
# the aim is to use as many track candidates as possible to build taus
# in order to maximise the efficiency
process.ak4PFJetsRecoTauChargedHadrons.builders[1].qualityCuts.signalQualityCuts.maxDeltaZ = cms.double(100.)
process.ak4PFJetsRecoTauChargedHadrons.builders[1].qualityCuts.signalQualityCuts.maxTrackChi2 = cms.double(1000.)
process.ak4PFJetsRecoTauChargedHadrons.builders[1].qualityCuts.signalQualityCuts.maxTransverseImpactParameter = cms.double(100.)
process.ak4PFJetsRecoTauChargedHadrons.builders[1].qualityCuts.pvFindingAlgo = cms.string('highestPtInEvent')

# chargedPFCandidates from PFNeutralHadrons ... yes, from neutrals too, nothing is thrown away
process.ak4PFJetsRecoTauChargedHadrons.builders[2].qualityCuts.signalQualityCuts.maxDeltaZ = cms.double(100.)
process.ak4PFJetsRecoTauChargedHadrons.builders[2].qualityCuts.signalQualityCuts.maxTrackChi2 = cms.double(1000.)
process.ak4PFJetsRecoTauChargedHadrons.builders[2].qualityCuts.signalQualityCuts.maxTransverseImpactParameter = cms.double(100.)
process.ak4PFJetsRecoTauChargedHadrons.builders[2].qualityCuts.pvFindingAlgo = cms.string('highestPtInEvent')


# change the output file name, don't overwrite the original file!
# process.output.fileName = cms.untracked.string('{}_miniAOD_rerunTauRECO.root'.format("ZTT" if runSignal else "QCD"))

process.output = cms.OutputModule('PoolOutputModule',
                                         fileName=cms.untracked.string('outputFULL.root'),
                                         fastCloning=cms.untracked.bool(False),
                                         dataset=cms.untracked.PSet(
                                             dataTier=cms.untracked.string('RECO'),
                                             filterName=cms.untracked.string('')
                                         ),
                                         outputCommands=cms.untracked.vstring('keep *'),
                                         SelectEvents=cms.untracked.PSet(
                                             SelectEvents=cms.vstring('*',)
                                         )
                                 )

process.out = cms.EndPath(process.output)


process.output.fileName = cms.untracked.string('outputHLT_3.root')
# process.output.fileName = cms.untracked.string('/eos/user/b/bskipwor/outputHLT/1_3.root')
process.output.outputCommands.append('keep *_genParticlePlusGeant_*_*')
# process.output.outputComands = cms.untracked.vstring(
#     'drop *',
#     ''
# )
