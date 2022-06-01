'''
Loops on the events and operates the matching between reconstructed and generated taus.
It produces two flat ntuples:
    - one with an entry for each gen tau (useful for efficiencies)
    - one with an entry for each reconstructed tau (useful for fake studies)
'''
import numpy as np
import ROOT
import sys, os, pdb, argparse
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from treeVariables import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer, genDecayModeGEANT, isAncestor # utility functions
from math import sqrt, pow
from copy import deepcopy as dc
from itertools import product as itertools_product


parser = argparse.ArgumentParser(
		    description="Convert MiniAOD to flat ntuples!")
parser.add_argument(
	"--sample",
# 	choices=['ZTT','SMS', 'test'],
	required=True,
	help='Specify the sample you want to flatten')

parser.add_argument(
	"--file",
	required=False,
	default=0,
	type=int,
	help='Specify the file number in the list from das')

args = parser.parse_args()
ifile = args.file
sample = args.sample

#########################################################################################
feat_list = ['lxy', 'dxy', 'visdxy', 'cosxy', 'momct', 'momct2d', 'mom_mass', 'pi_lxy', 'pi_cosxy' ]
mom_pdgId = [1000015]
dR_cone = 0.1
if 'HNL' in sample:  mom_pdgId = [9900012]
if 'HNL' in sample and 'Dirac' in sample:  mom_pdgId = [9990012]

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('tau_gentau_tuple_{}_{}.root'.format(sample,ifile), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(list(zip(branches, [-99.]*len(branches)))) # initialise all branches to unphysical -99

# outfile_jet = ROOT.TFile('tau_jet_tuple_{}_{}.root'.format(sample,ifile), 'recreate')
# ntuple_jet = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
# tofill_jet = OrderedDict(zip(branches, [-99.]*len(branches))) # initialise all branches to unphysical -99

##########################################################################################
# Get ahold of the events
f = open('%s_filelist.txt'%args.sample)
infile = f.readlines()[ifile]

#events = Events('root://cms-xrd-global.cern.ch/'+infile) # make sure this corresponds to your file name!
events = Events(infile.strip()) # make sure this corresponds to your file name!
#events = Events('/eos/user/b/bskipwor/HNL_M_5/HNL_miniAOD_1_3.root') # make sure this corresponds to your file name!
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files

def isAncestor(a, p):
    if a == p :
        return True
    for i in range(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
    return False

def findMatchToGen(gen_taus, hlt_taus, hlt_tau):
    gen_taus_match = gen_taus
    for gg in gen_taus_match:
        setattr(gg, hlt_tau, None)
        bestcom = bestMatch(gg.visp4, hlt_taus)
        print(gg.visp4)
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
            setattr(gg,hlt_tau,bestcom[0])
#             if hlt_tau=='hlt_pftau_displ':  pdb.set_trace()
#            print("found one", hlt_tau)
    return gen_taus_match


##########################################################################################
# instantiate the handles to the relevant collections.
# Do this *outside* the event loop
# Reminder: you can see the names of the collections stored in your input file by doing:
# edmDumpEventContent outputFULL.root

# PAT taus
#label_taus = ('slimmedTaus', '', 'PAT')
label_taus = ('selectedPatTaus', '', 'TAURECO')
handle_taus = Handle('std::vector<pat::Tau>')

# PAT jets
#label_jets = ('slimmedJets', '', 'PAT')
#handle_jets = Handle('std::vector<pat::Jet>')

# gen particles
label_gen  = ('genParticlePlusGeant', '', 'SIM')
handle_gen = Handle('std::vector<reco::GenParticle>')

# AOD
# vertices
handle_vtx = Handle('std::vector<reco::Vertex>')
label_vtx  = ('offlinePrimaryVertices','','RECO')

# L1 taus
handle_l1 = Handle('BXVector<l1t::Tau>')
label_l1  = ('caloStage2Digis','Tau','RECO')

# L1 jet
handle_l1j = Handle('BXVector<l1t::Jet>')
label_l1j  = ('caloStage2Digis','Jet','RECO')

# general tracks
label_tracks = ('generalTracks', '', 'RECO')
handle_tracks = Handle('std::vector<reco::Track>')

# hlt pftaus
label_hlt_pftaus_displ = ('hltHpsPFTauProducerDispl', '',  'MYHLT')
handle_hlt_pftaus_displ = Handle('std::vector<reco::PFTau>')

# MINIAOD
#########################################################################################
# gen particles
#label_gen  = ('genParticle', '', 'SIM')
#handle_gen = Handle('std::vector<reco::GenParticle>')

# vertices
#handle_vtx = Handle('std::vector<reco::Vertex>')
#label_vtx  = ('offlineSlimmedPrimaryVertices','','PAT')

# lost tracks
#label_lost = ('lostTracks', '', 'PAT')
#handle_lost = Handle('std::vector<pat::PackedCandidate>')

# packed PFCandidates
#label_packed = ('packedPFCandidates')
#handle_packed = Handle('std::vector<pat::PackedCandidate')
#########################################################################################

# instantiate the handles to the relevant collections.
# handles = OrderedDict()
# handles['taus'          ] = ('slimmedTaus'                     , Handle('std::vector<pat::Tau>'))
# handles['muon'          ] = ('slimmedJets'                     , Handle('std::vector<pat::Jet>'))
# handles['gen_particles' ] = ('prunedGenParticles'              , Handle('std::vector<reco::GenParticle>'))
# handles['vertices'      ] = ('offlineSlimmedPrimaryVertices'   , Handle('std::vector<reco::Vertex>'))
# handles['l1_taus'       ] = (('caloStage2Digis', 'Tau', 'RECO'), Handle('BXVector<l1t::Tau>'))
# handles['l1_jets'       ] = (('caloStage2Digis', 'Jet', 'RECO'), Handle('BXVector<l1t::Jet>'))
#
##########################################################################################
# example histogram
histos = OrderedDict()
histos['pull_pt'] = ROOT.TH1F('pull_pt', 'pull_pt', 50, -2, 2)
histos['pull_pt'].GetXaxis().SetTitle('(p_{T}^{off} - p_{T}^{gen})/p_{T}^{gen}')
histos['pull_pt'].GetYaxis().SetTitle('counts')

##########################################################################################
# start looping on the events
for i, ev in enumerate(events):

    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break

    if i%100==0:
        print('===> processing %d / %d event' %(i, totevents))

#     for k, v in handles.iteritems():
#         ev.getByLabel(v[0], v[1])
#         setattr(ev, k, v[1].product())

    ######################################################################################
    # access the taus
    ev.getByLabel(label_taus, handle_taus)
    taus = handle_taus.product()

    # loosely filter the reco taus
    taus = [tau for tau in taus if tau.pt()>18.]

    ######################################################################################
    # access the jets
#    ev.getByLabel(label_jets, handle_jets)
#    jets = handle_jets.product()
#
     # loosely filter jets
#     jets = [jet for jet in jets if jet.pt()>18. and abs(jet.eta())<2.5]

    ######################################################################################
    # access the vertices
    ev.getByLabel(label_vtx, handle_vtx)
    vertices = handle_vtx.product()

    ######################################################################################
    # access the gen taus
    ev.getByLabel (label_gen, handle_gen)
    gen_particles = handle_gen.product()

    # select only hadronically decaying taus
    gen_taus = [pp for pp in gen_particles if abs(pp.pdgId())==15 and \
                pp.status()==8 and isGenHadTau(pp) and\
                pp.pt()>10 and abs(pp.eta()) < 2.1]

    if len(gen_taus) == 0:  continue

    # determine gen decaymode and find mother
    for gg in gen_taus:

        ## reset other attributes
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])
        if 'gmsb' in sample:
            dm_string = genDecayModeGEANT( [d for d in finalDaughters(gg) if abs(d.pdgId()) not in [12, 14, 16]])
            gg.decayMode = tauDecayModes.nameToInt(dm_string)

        gg.bestmom = None
        if gg.numberOfMothers() > 1:
            ## print 'more than 1 one, taking first one'
            tau_moms = [imom for imom in ev.gen_particles if isAncestor(imom, gg) and abs(imom.pdgId()) in mom_pdgId]
            gg.bestmom = tau_moms[0]
        elif gg.numberOfMothers() == 1 and abs(gg.mother(0).pdgId()) in mom_pdgId:
            gg.bestmom = gg.mother(0)

        ### find first dau to be used for vertices
        gg.dau = None
        if gg.numberOfDaughters() > 0:  gg.dau = gg.daughter(0)
        if gg.dau == None:  print('is none')
#         else: pdb.set_trace()

        if gg.bestmom == None or gg.dau == None:
            continue


#     for gg in gen_taus:
#         print
#         pdb.set_trace()
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()

        if 'taugun' in sample:
            gg.lxy = sqrt(pow(gg.dau.vx(),2)+pow(gg.dau.vy(),2))
#             gg.visdxy = (gg.dau.vy()*gg.vispx() - gg.dau.vx()*gg.vispy())/gg.vispt()
        else:
            gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))
#             gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()

        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))

        ## gg = tau   bestmom = stau    vx = production point
        ## now find tau decay point to stau production point distance
        gg.pi_lxy = sqrt(pow(gg.dau.vx()-gg.bestmom.vx(),2) + pow(gg.dau.vy()-gg.bestmom.vy(),2) )

#             stau_vectorP = np.array([bestmom.px(), bestmom.py(), 0])
#             stau_vectorL = np.array([gg.vx()-bestmom.vx(), gg.vy()-bestmom.vy(), 0])
#             gg.tau_cosxy  = tau_vectorL.dot(tau_vectorP)/((np.linalg.norm(tau_vectorL) * np.linalg.norm(tau_vectorP)))
#             l3d = sqrt(pow(gg.vx()-bestmom.vx(),2) + pow(gg.vy()-bestmom.vy(),2) + pow(gg.vz()-bestmom.vz(),2))
#             if bestmom.p() != 0: gg.momct    = c_const*l3d*bestmom.mass()/bestmom.p()
        if gg.bestmom.pt() != 0:  gg.momct2d = c_const*gg.lxy*gg.bestmom.mass()/gg.bestmom.pt()

    ######################################################################################
    # match reco taus to gen taus
    for tt in taus    : tt.gen_tau  = None # first initialise the matching to None
    for gg in gen_taus: gg.reco_tau = None # first initialise the matching to None

    gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched

    for tt in taus:
        matches = [gg for gg in gen_taus_copy if deltaR(tt.p4(), gg.visp4)<0.3]
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(tt.p4(), gg.visp4))
        bestmatch = matches[0]
        tt.gen_tau = bestmatch
        bestmatch.reco_tau = tt
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

        for recotau in taus:
            if recotau.leadChargedHadrCand().isNonnull():
                mytrk = recotau.leadChargedHadrCand().get()
#                pdb.set_trace()
                print("collection", mytrk.bestTrack().originalAlgo(), mytrk.bestTrack().algo())

    ######################################################################################
    # match reco taus to reco jets
#     for jj in jets : jj.tau = None # first initialise the matching to None
#
#     taus_copy = taus # we'll cyclically remove any tau that gets matched
#
#     for jj in jets:
#         matches = [tt for tt in taus_copy if deltaR(jj.p4(), tt.p4())<0.3]
#         if not len(matches):
#             continue
#         matches.sort(key = lambda tt : deltaR(jj.p4(), tt.p4()))
#         bestmatch = matches[0]
#         jj.tau = bestmatch
#         taus_copy = [tt for tt in taus_copy if tt != bestmatch]

    ######################################################################################
    try:
        ev.getByLabel(label_hlt_pftaus_displ, handle_hlt_pftaus_displ)
        hlt_pftau = handle_hlt_pftaus_displ.product()

        gen_taus = findMatchToGen(gen_taus, hlt_pftau, 'hlt_pftau_displ')

        for gg in gen_taus:
            if hasattr(gg, 'hlt_pftau_displ') and gg.hlt_pftau_displ:

                theLeadChargedCand = gg.hlt_pftau_displ.leadChargedHadrCand()
                theLeadPFCand      = gg.hlt_pftau_displ.leadCand()
                gg.hlt_pftau_displ.leadChargedCandPt     = theLeadChargedCand.pt()
                gg.hlt_pftau_displ.leadChargedCandPdgId  = theLeadChargedCand.pdgId()
                gg.hlt_pftau_displ.leadCandPt            = theLeadPFCand.pt()
                gg.hlt_pftau_displ.leadCandPdgId         = theLeadPFCand.pdgId()

                gg.hlt_pftau_displ.maxHCALPFClusterEt    = gg.hlt_pftau_displ.maximumHCALPFClusterEt()
                gg.hlt_pftau_displ.nChargedHad           = gg.hlt_pftau_displ.signalChargedHadrCands().size()
                gg.hlt_pftau_displ.nGamma                = gg.hlt_pftau_displ.signalGammaCands().size()

                sum_pt_charged = 0
                for ich in range(gg.hlt_pftau_displ.nChargedHad  ) :
                    sum_pt_charged += gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pt()
#                   print '\t ch: ', gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pdgId(), '\t', gg.hlt_pftau_displ.signalChargedHadrCands()[ich].pt()

                sum_pt_neutral = 0
                for ineu in range(gg.hlt_pftau_displ.nGamma  ) :
                    sum_pt_neutral += gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()
#                   print '\t neu: ', gg.hlt_pftau_displ.signalGammaCands()[ineu].pdgId(), '\t', gg.hlt_pftau_displ.signalGammaCands()[ineu].pt()

                gg.hlt_pftau_displ.sum_pt_charged =  sum_pt_charged
                gg.hlt_pftau_displ.sum_pt_neutral =  sum_pt_neutral
                print("hlt_pftau_displ collection found")

    except:
        print("collection not found")

    ######################################################################################
    # fill histograms
    for gg in gen_taus:
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            pull = (gg.reco_tau.pt() - gg.vispt())/gg.vispt()
            histos['pull_pt'].Fill(pull)

    ######################################################################################
    # find ancestor
#     print 'event', i
    c_const = 299.792458
    for gg in gen_taus :
        gg.lxy   = -9999. # first initialise to None
        gg.myvx  = -9999. # first initialise to None
        gg.myvy  = -9999. # first initialise to None
        gg.myvz  = -9999. # first initialise to None
        gg.cosxy = -9999. # first initialise to None
        gg.momct = -9999. # first initialise to None
        gg.momct2d  = -9999. # first initialise to None
        gg.mom_mass = -9999. # first initialise to None

    for gg in gen_taus:
        #tau_moms_tmp = [imom for imom in gen_particles if isAncestor(imom, gg)]
        #for imom in tau_moms_tmp:  print imom.pdgId()
        tau_moms = [imom for imom in gen_particles if isAncestor(imom, gg) and abs(imom.pdgId())==mom_pdgId]
        if len(tau_moms)>0 and tau_moms[0]!= None:
            bestmom = tau_moms[0]

            gg.lxy = sqrt(pow(gg.vx()-bestmom.vx(),2)+pow(gg.vy()-bestmom.vy(),2))
            gg.myvx = bestmom.vx()
            gg.myvy = bestmom.vy()
            gg.myvz = bestmom.vz()

            vectorP = np.array([gg.px(), gg.py(), 0])
            vectorL = np.array([gg.vx()-bestmom.vx(), gg.vy()-bestmom.vy(), 0])
            gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))

            l3d = sqrt(pow(gg.vx()-bestmom.vx(),2) + pow(gg.vy()-bestmom.vy(),2) + pow(gg.vz()-bestmom.vz(),2))
#             pdb.set_trace()
            gg.momct    = c_const*l3d*bestmom.mass()/bestmom.p()
            gg.momct2d  = c_const*gg.lxy*bestmom.mass()/bestmom.pt()
            gg.mom_mass = bestmom.mass()

        else:
            gg.lxy    = -9999.
            gg.myvx   = -9999.
            gg.myvy   = -9999.
            gg.myvz   = -9999.
            gg.cosxy  = -9999.
            gg.momct  = -9999.
            gg.momct2d  = -9999.
            gg.mom_mass = -9999.



#       math::XYZVector pperp(refMuP->px() + refMuM->px(), refMuP->py() + refMuM->py(), 0.);
#
#       GlobalPoint displacementFromBeamspot(-1*((bs.x0() - mumuVertex.position().x()) + (mumuVertex.position().z() - bs.z0()) * bs.dxdz()),
#                                            -1*((bs.y0() - mumuVertex.position().y()) + (mumuVertex.position().z() - bs.z0()) * bs.dydz()), 0);
#
#       double LxyJpsi        = displacementFromBeamspot.perp();
#       reco::Vertex::Point vperp(displacementFromBeamspot.x(),displacementFromBeamspot.y(),0.);
#       double cosJpsiXY      = vperp.Dot(pperp)/(vperp.R()*pperp.R());
#
#

    ######################################################################################
    # access the l1 taus
    ev.getByLabel (label_l1, handle_l1)
    l1_taus = handle_l1.product()

    l1_taus = [tau for tau in l1_taus if l1_taus.getFirstBX()==0]
    l1_taus = [tau for tau in l1_taus if tau.pt()>15.]

    # match reco taus to gen taus
    for l1 in l1_taus : l1.gen_tau  = None # first initialise the matching to None
    for gg in gen_taus: gg.l1_tau   = None # first initialise the matching to None

    gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched

    for l1 in l1_taus:
        matches = [gg for gg in gen_taus_copy if deltaR(l1.p4(), gg.visp4)<0.3]
        ave_matches = len(matches)
        if ave_matches > 1:  print(ave_matches)
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(l1.p4(), gg.visp4))
        bestmatch = matches[0]
        l1.gen_tau = bestmatch
        bestmatch.l1_tau = l1
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

    ####################################################################
    # access the l1 jets
    ev.getByLabel (label_l1j, handle_l1j)
    l1_jets = handle_l1j.product()

    l1_jets = [tau for tau in l1_jets if l1_jets.getFirstBX()==0]
    l1_jets = [tau for tau in l1_jets if tau.pt()>15.]

    # match reco jets to gen taus
    for l1 in l1_jets : l1.gen_tau  = None # first initialise the matching to None
    for gg in gen_taus: gg.l1_jet   = None # first initialise the matching to None

    gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched

    for l1 in l1_jets:
        matches = [gg for gg in gen_taus_copy if deltaR(l1.p4(), gg.visp4)<0.3]
        ave_matches = len(matches)
        if ave_matches > 1:  print(ave_matches)
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(l1.p4(), gg.visp4))
        bestmatch = matches[0]
        l1.gen_tau = bestmatch
        bestmatch.l1_jet = l1
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

    ######################################################################################
    #access general tracks
    ev.getByLabel(label_tracks, handle_tracks)
    track = handle_tracks.product()

    track = [tau for tau in track if tau.pt()>15]
#    print("number of tracks",len(track))

    # match general tracks to gen taus
    for rr in track : rr.gen_tau = None # first initialise the matching to None
    for gg in gen_taus: gg.tracks  = None # first initialise the matching to None

    gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched

    for rr in track:
        matches = [gg for gg in gen_taus_copy if deltaR(rr.momentum(), gg.visp4)<0.7 and abs(rr.pt() - gg.vispt())<0.6*gg.vispt()]
#        print("number of matches",len(matches))
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(rr.momentum(), gg.visp4))
        bestmatch = matches[0]
        rr.gen_tau = bestmatch
        bestmatch.tracks = rr
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]

    for gg in gen_taus:
        if hasattr(gg, 'tracks') and gg.tracks:  
            print("attribute")  


    ######################################################################################
    # access the lost tracks
    #ev.getByLabel(label_lost, handle_lost)
    #lost = handle_lost.product()

    # only keep pion candidates
    #lost_tracks = [ll for ll in lost if abs(ll.pdgId())==211]

    ######################################################################################
    # access packed PFCandidates
    #ev.getByLabel(label_packed, handle_packed)
    #packed = handle_packed.product()

    # only keep pion candidates
    #packed_tracks = [ff for ff in packed if abs(ll.pdgId())==211]

    ######################################################################################
    # add together lost tracks and packed PFCandidates
    #comtracks = lost_tracks + packed_tracks

    '''
    # takes event number as input if matched reco and not matched comtracks to gen tau for 1-prong
    if ev.eventAuxiliary().event() == 416082:
        al_pt  = []
        al_eta = []
        al_phi = []
        for jj in comtracks:
            for gg in gen_taus:
                apt  = abs(jj.pt() - gg.pt())
                al_pt.append(apt)
                aeta = abs(jj.eta() - gg.eta())
                al_eta.append(aeta)
                aphi = abs(jj.phi() - gg.phi())
                al_phi.append(aphi)
        print(("min_diff_pt", min(al_pt)))
        print(("min_diff_eta", min(al_eta)))
        print(("min_diff_phi", min(al_phi)))


    ## skim comtracks collection to keep only PFcands that are "close" in dR and dPt from at least a gen tau
    comtracks_matchable = []
    for cc in comtracks:
        found = False
        for gg in gen_taus :
            if found == False and (abs(cc.pt() - gg.vispt())<0.2*gg.vispt()) and ( deltaR(cc.p4(), gg.visp4)<0.3) :
                comtracks_matchable.append(cc)
                found = True

    ## for each gen tau, find the PFCand that best matches in dR
    for gg in gen_taus  : gg.up_com_tau = None # first initialise the matching to None
    for gg in gen_taus:
      bestcom = bestMatch(gg, comtracks_matchable )
      if bestcom[0] != None:
          #print deltaR(bestcom[0],gg)
          gg.up_com_tau = bestcom[0]



    ## original version, using same approach as for reco_taus
    # match combined lost and PFCandidate tracks to gen taus
    for cc in comtracks : cc.gen_tau = None # first initialise the matching to None
    for gg in gen_taus  : gg.com_tau = None # first initialise the matching to None

    gen_taus_copy = gen_taus # we'll cyclically remove any gen taus that gets matched

    for cc in comtracks:
        matches = [gg for gg in gen_taus_copy if deltaR(cc.p4(), gg.visp4)<0.3 and abs(cc.pt() - gg.vispt())<0.2*gg.vispt()]
        if not len(matches):
            continue
        matches.sort(key = lambda gg : deltaR(cc.p4(), gg.visp4))
        bestmatch = matches[0]
        cc.gen_tau = bestmatch
        bestmatch.com_tau = cc
        gen_taus_copy = [gg for gg in gen_taus_copy if gg != bestmatch]
        '''

    ######################################################################################
    # fill the ntuple: each gen tau makes an entry

    for gg in gen_taus:
        for k, v in tofill_gen.items(): tofill_gen[k] = -99. # initialise before filling
        tofill_gen['run'               ] = ev.eventAuxiliary().run()
        tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
        tofill_gen['event'             ] = ev.eventAuxiliary().event()
        tofill_gen['nvtx'              ] = vertices.size()
        tofill_gen['PV_x'              ] = vertices[0].x()
        tofill_gen['PV_y'              ] = vertices[0].y()
        tofill_gen['PV_z'              ] = vertices[0].z()
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            tofill_gen['tau_reco_mass'     ] = gg.reco_tau.mass()
            tofill_gen['tau_reco_pt'       ] = gg.reco_tau.pt()
            tofill_gen['tau_reco_eta'      ] = gg.reco_tau.eta()
            tofill_gen['tau_reco_phi'      ] = gg.reco_tau.phi()
            tofill_gen['tau_reco_charge'   ] = gg.reco_tau.charge()
            tofill_gen['tau_reco_decaymode'] = gg.reco_tau.decayMode()
            tofill_gen['tau_reco_ip3d'     ] = gg.reco_tau.ip3d()
            tofill_gen['tau_reco_dxy'      ] = gg.reco_tau.dxy()
#            tofill_gen['tau_reco_pixel'    ] = gg.reco_tau.leadChargedHadrCand().numberOfPixelHits()

        if hasattr(gg, 'l1_tau') and gg.l1_tau:
#            tofill_gen['tau_l1_mass'     ] = gg.reco_tau.mass()
            tofill_gen['tau_l1_pt'       ] = gg.l1_tau.pt()
            tofill_gen['tau_l1_eta'      ] = gg.l1_tau.eta()
            tofill_gen['tau_l1_phi'      ] = gg.l1_tau.phi()
            tofill_gen['tau_l1_charge'   ] = gg.l1_tau.charge()
            tofill_gen['tau_l1_iso'      ] = gg.l1_tau.hwIso()

        if hasattr(gg, 'l1_jet') and gg.l1_jet:
#            tofill_gen['tau_l1_mass'     ] = gg.reco_tau.mass()
            tofill_gen['jet_l1_pt'       ] = gg.l1_jet.pt()
            tofill_gen['jet_l1_eta'      ] = gg.l1_jet.eta()
            tofill_gen['jet_l1_phi'      ] = gg.l1_jet.phi()
            tofill_gen['jet_l1_charge'   ] = gg.l1_jet.charge()

        if hasattr(gg, 'tracks') and gg.tracks:
            tofill_gen['tau_track_pt'      ] = gg.tracks.pt()
            tofill_gen['tau_track_eta'     ] = gg.tracks.eta()
            tofill_gen['tau_track_phi'     ] = gg.tracks.phi()
            tofill_gen['tau_track_charge'  ] = gg.tracks.charge()

#        if hasattr(gg, 'com_tau') and gg.com_tau:
#            tofill_gen['tau_com_pt'       ] = gg.com_tau.pt()
#            tofill_gen['tau_com_eta'      ] = gg.com_tau.eta()
#            tofill_gen['tau_com_phi'      ] = gg.com_tau.phi()
#            tofill_gen['tau_com_charge'   ] = gg.com_tau.charge()
#            tofill_gen['tau_com_pixel'    ] = gg.com_tau.numberOfPixelHits()

#        if hasattr(gg, 'up_com_tau') and gg.up_com_tau:
#            tofill_gen['tau_up_com_pt'       ] = gg.up_com_tau.pt()
#            tofill_gen['tau_up_com_eta'      ] = gg.up_com_tau.eta()
#            tofill_gen['tau_up_com_phi'      ] = gg.up_com_tau.phi()
#            tofill_gen['tau_up_com_charge'   ] = gg.up_com_tau.charge()

        if hasattr(gg, 'hlt_pftau_displ') and gg.hlt_pftau_displ:
            tofill_gen['tau_hltPFdispltau_mass'     ] = gg.hlt_pftau_displ.mass()
            tofill_gen['tau_hltPFdispltau_pt'       ] = gg.hlt_pftau_displ.pt()
            tofill_gen['tau_hltPFdispltau_eta'      ] = gg.hlt_pftau_displ.eta()
            tofill_gen['tau_hltPFdispltau_phi'      ] = gg.hlt_pftau_displ.phi()
            tofill_gen['tau_hltPFdispltau_charge'   ] = gg.hlt_pftau_displ.charge()
            tofill_gen['tau_hltPFdispltau_decaymode'] = gg.hlt_pftau_displ.decayMode()
            ## new vars
            tofill_gen['tau_hltPFdispltau_leadChargedCandPdgId'] = gg.hlt_pftau_displ.leadChargedCandPdgId
            tofill_gen['tau_hltPFdispltau_leadChargedCandPt'   ] = gg.hlt_pftau_displ.leadChargedCandPt
            tofill_gen['tau_hltPFdispltau_leadPFCandPdgId'     ] = gg.hlt_pftau_displ.leadCandPdgId
            tofill_gen['tau_hltPFdispltau_leadPFCandPt'        ] = gg.hlt_pftau_displ.leadCandPt
            tofill_gen['tau_hltPFdispltau_maxHCALPFClusterEt'  ] = gg.hlt_pftau_displ.maxHCALPFClusterEt
            tofill_gen['tau_hltPFdispltau_nChargedHad'         ] = gg.hlt_pftau_displ.nChargedHad
            tofill_gen['tau_hltPFdispltau_nGamma'              ] = gg.hlt_pftau_displ.nGamma
            tofill_gen['tau_hltPFdispltau_sumPtCharged'        ] = gg.hlt_pftau_displ.sum_pt_charged
            tofill_gen['tau_hltPFdispltau_sumPtNeutral'        ] = gg.hlt_pftau_displ.sum_pt_neutral

#            tofill_gen['tau_hltPFdispltau_dxy'              ] = gg.hlt_pftau_displ.dxy
#            tofill_gen['tau_hltPFdispltau_dxyerr'           ] = gg.hlt_pftau_displ.dxyerr
#            tofill_gen['tau_hltPFdispltau_ip3d'             ] = gg.hlt_pftau_displ.ip3d
#            tofill_gen['tau_hltPFdispltau_ip3derr'          ] = gg.hlt_pftau_displ.ip3derr
#            tofill_gen['tau_hltPFdispltau_passChargedIso'   ] = gg.hlt_pftau_displ.passChargedIso
#            tofill_gen['tau_hltPFdispltau_passRelChargedIso'] = gg.hlt_pftau_displ.passRelChargedIso
#            tofill_gen['tau_hltPFdispltau_passAbsChargedIso'] = gg.hlt_pftau_displ.passAbsChargedIso
#            tofill_gen['tau_hltPFdispltau_isoval'           ] = gg.hlt_pftau_displ.isoVal
#            tofill_gen['tau_hltPFdispltau_passFilters'      ] = gg.hlt_pftau_displ.passFilters


        tofill_gen['tau_gen_pt'        ] = gg.pt()
        tofill_gen['tau_gen_eta'       ] = gg.eta()
        tofill_gen['tau_gen_phi'       ] = gg.phi()
        tofill_gen['tau_gen_charge'    ] = gg.charge()
        tofill_gen['tau_gen_decaymode' ] = gg.decayMode
        tofill_gen['tau_gen_lxy'       ] = gg.lxy
        tofill_gen['tau_gen_vx'        ] = gg.myvx  ## user defined ones (mother prod. vertex)
        tofill_gen['tau_gen_vy'        ] = gg.myvy
        tofill_gen['tau_gen_vz'        ] = gg.myvz
        tofill_gen['tau_gen_cosxy'     ] = gg.cosxy
        tofill_gen['tau_gen_momct'     ] = gg.momct
        tofill_gen['tau_gen_momct2d'   ] = gg.momct2d
        tofill_gen['tau_gen_mom_mass'  ] = gg.mom_mass
        tofill_gen['tau_gen_vis_mass'  ] = gg.vismass()
        tofill_gen['tau_gen_vis_pt'    ] = gg.vispt()
        tofill_gen['tau_gen_vis_eta'   ] = gg.viseta()
        tofill_gen['tau_gen_vis_phi'   ] = gg.visphi()
        ntuple_gen.Fill(array('f',list(tofill_gen.values())))

        # use for miniAOD
        # if reco matches gen but comtracks does not match gen, prints event ID and pt, eta, phi for gen and reco
#        if hasattr(gg, 'reco_tau') and gg.reco_tau and gg.decayMode==0:
#            if hasattr(gg, 'up_com_tau') and gg.up_com_tau:
#                print("both")
#                print(("event_id", ev.eventAuxiliary().event()))
#                print(("gen_pt", gg.pt()))
#                print(("reco_pt",gg.reco_tau.pt()))
#            else:
#                print(("event_id", ev.eventAuxiliary().event()))
#                print(("gen_pt", gg.pt()))
#                print(("reco_pt",gg.reco_tau.pt()))
#                print(("gen_eta", gg.eta()))
#                print(("reco_eta", gg.reco_tau.eta()))
#                print(("gen_phi", gg.phi()))
#                print(("reco_phi", gg.reco_tau.phi()))

    # fill the ntuple: each jet makes an entry
#     for jj in jets:
#         for k, v in tofill_jet.iteritems(): tofill_jet[k] = -99. # initialise before filling
#         tofill_jet['run'        ] = ev.eventAuxiliary().run()
#         tofill_jet['lumi'       ] = ev.eventAuxiliary().luminosityBlock()
#         tofill_jet['event'      ] = ev.eventAuxiliary().event()
#         tofill_jet['nvtx'       ] = vertices.size()
#         tofill_jet['jet_mass'   ] = jj.mass()
#         tofill_jet['jet_pt'     ] = jj.pt()
#         tofill_jet['jet_eta'    ] = jj.eta()
#         tofill_jet['jet_phi'    ] = jj.phi()
#         tofill_jet['jet_charge' ] = jj.charge()
#         if hasattr(jj, 'tau') and jj.tau:
#             tofill_jet['tau_reco_mass'     ] = jj.tau.mass()
#             tofill_jet['tau_reco_pt'       ] = jj.tau.pt()
#             tofill_jet['tau_reco_eta'      ] = jj.tau.eta()
#             tofill_jet['tau_reco_phi'      ] = jj.tau.phi()
#             tofill_jet['tau_reco_charge'   ] = jj.tau.charge()
#             tofill_jet['tau_reco_decaymode'] = jj.tau.decayMode()
#             if hasattr(jj.tau, 'gen_tau') and jj.tau.gen_tau:
#                 tofill_gen['tau_gen_pt'        ] = jj.tau.gen_tau.pt()
#                 tofill_gen['tau_gen_eta'       ] = jj.tau.gen_tau.eta()
#                 tofill_gen['tau_gen_phi'       ] = jj.tau.gen_tau.phi()
#                 tofill_gen['tau_gen_charge'    ] = jj.tau.gen_tau.charge()
#                 tofill_gen['tau_gen_decaymode' ] = jj.tau.gen_tau.decayMode
#                 tofill_gen['tau_gen_vis_mass'  ] = jj.tau.gen_tau.vismass()
#                 tofill_gen['tau_gen_vis_pt'    ] = jj.tau.gen_tau.vispt()
#                 tofill_gen['tau_gen_vis_eta'   ] = jj.tau.gen_tau.viseta()
#                 tofill_gen['tau_gen_vis_phi'   ] = jj.tau.gen_tau.visphi()
#         ntuple_jet.Fill(array('f',tofill_jet.values()))

    ######################################################################################
    # printout some info, if you want
    # printer(taus, gen_taus)

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()

# outfile_jet.cd()
# ntuple_jet.Write()
# outfile_jet.Close()
