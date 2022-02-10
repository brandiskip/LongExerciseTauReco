import FWCore.ParameterSet.Config as cms

# trick to make python know where to look for the imports
import sys
sys.path.append('..')

# for example: here
from rerunTauRecoOnMiniAOD import process

runSignal = False # Set to False to read in QCD file instead of ZTT

maxEvents = -1

readFiles = cms.untracked.vstring()
secFiles = cms.untracked.vstring()
process.source = cms.Source(
    "PoolSource", fileNames=readFiles, secondaryFileNames=secFiles)
# add , eventsToProcess = cms.untracked.VEventRange('1:958444-1:958444','2:100-2:101') right after secondaryFileNames=secFiles above to run on specific event(s)

# print('\t Max events:', process.maxEvents.input.value())


readFiles.extend([
    'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/HeavyNeutrino_trilepton_M-5_V-0_00639_tau_massiveAndCKM_LO/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/FCE8BD0F-2C74-EA40-A4FC-12C3EE4B0285.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/HeavyNeutrino_trilepton_M-5_V-0_00639_tau_massiveAndCKM_LO/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/F0C6A50F-B346-1B47-9359-5005A3129424.root',
    'root://cms-xrd-global.cern.ch//store/mc/RunIIAutumn18MiniAOD/HeavyNeutrino_trilepton_M-5_V-0_00639_tau_massiveAndCKM_LO/MINIAODSIM/102X_upgrade2018_realistic_v15-v1/270000/CFC82CF4-FE57-7840-ADC5-4B81AA85899B.root',
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
process.output.fileName = cms.untracked.string('/eos/user/b/bskipwor/HNL_5_Samples/HNL_miniAOD_17_3.root')
# process.output.outputComands = cms.untracked.vstring(
#     'drop *',
#     ''
# )
