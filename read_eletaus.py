'''
Loops on the events and operates the matching between reconstructed and generated taus.
It produces two flat ntuples:
    - one with an entry for each gen tau (useful for efficiencies)
    - one with an entry for each reconstructed tau (useful for fake studies)
'''
import ROOT, sys, os, pdb, argparse
import numpy as np
from array import array
from collections import OrderedDict
from DataFormats.FWLite import Events, Handle
from PhysicsTools.HeppyCore.utils.deltar import deltaR, deltaPhi, bestMatch
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from treeVariables import branches # here the ntuple branches are defined
from utils import isGenHadTau, finalDaughters, printer , genDecayModeGEANT, isAncestor, isGenLepTau# utility functions
from math import sqrt, pow
from copy import deepcopy as dc
from itertools import product as itertools_product


parser = argparse.ArgumentParser(description="Convert MiniAOD to flat ntuples!")
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

##########################################################################################
feat_list = ['lxy', 'dxy', 'visdxy', 'cosxy', 'momct', 'momct2d', 'mom_mass', 'pi_lxy', 'pi_cosxy' ]
pftau_feat_list = ['leadChargedCandPt', 'leadChargedCandPdgId', 'leadCandPt', 'leadCandPdgId', 'maxHCALPFClusterEt', 'nChargedHad', 'nGamma', 'sum_pt_charged', 'sum_pt_neutral', 'dxy', 'dxyerr', 'ip3d', 'ip3derr', 'isoVal', 'passChargedIso', 'passRelChargedIso', 'passAbsChargedIso', 'passFilters']
ele_feat_list = ['dxy']
c_const = 299.792458
dR_cone = 0.1
dPt_cone = 0.2
good_gen_status = 2
mom_pdgId = [1000015]
if 'HNL' in sample:  mom_pdgId = [9900012]
if 'HNL' in sample and 'Dirac' in sample:  mom_pdgId = [9990012]
if 'gmsb' in sample:  
    mom_pdgId = [2000015, 1000015]
    good_gen_status = 8
if 'DY' in sample:  mom_pdgId = [23]

##########################################################################################
# initialise output files to save the flat ntuples
outfile_gen = ROOT.TFile('tau_ele_tuple_{}_{}.root'.format(sample,ifile), 'recreate')
ntuple_gen = ROOT.TNtuple('tree', 'tree', ':'.join(branches))
tofill_gen = OrderedDict(list(zip(branches, [-9999.]*len(branches)))) # initialise all branches to unphysical -99       
print((outfile_gen.GetName())) 

##########################################################################################
# Get ahold of the events
f = open('%s_filelist.txt'%args.sample)
infile = f.readlines()[ifile]

os.environ['X509_USER_PROXY'] = '/afs/cern.ch/user/b/bskipwor/x509up_u58808' 
redirector = 'root://cms-xrd-global.cern.ch//'
if 'fnal' in sample:
    redirector = 'root://cmseos.fnal.gov//'
if 'eos' in sample:
    redirector = ''

events = Events(redirector+infile.strip())
#events = Events([
#                  '/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/CMSSW_12_3_0_pre6/src/HLTrigger/Configuration/test/outputHLT.root'
#                  '/eos/cms/store/group/phys_bphys/fiorendi/p5prime/displTaus/Staus_M_200_100mm_14TeV_Run3MC_Summer21/crab_ntuples_mutau_singlephoton_v18_gmsb_M200_100mm_summer21/220627_083854/0000/outputHLT_111.root',
#                 '/afs/cern.ch/work/f/fiorendi/private/displacedTaus/hlt/CMSSW_12_4_0/src/HLTrigger/Configuration/test/outputHLT.root'
#                ])

print('using a local file!!!!!!!!!!!!!!!!!!!!!')


print(infile)
maxevents = -1 # max events to process
totevents = events.size() # total number of events in the files


##########################################################################################
def findMatchToGen(gen_taus, hlt_taus, hlt_tau): 
 
    gen_taus_match = gen_taus 
    for gg in gen_taus_match: 
        setattr(gg, hlt_tau, None)
        bestcom = bestMatch(gg.visp4, hlt_taus )
        if bestcom[0] != None and sqrt(bestcom[1]) < dR_cone :
            setattr(gg,hlt_tau,bestcom[0])
    return gen_taus_match    


def printForDebug(cc, gg):
    print(('    pt: %.2f'%( cc.pt()), \
          '\t eta: %.2f'%( cc.eta()), \
          '\t phi: %.2f'%( cc.phi()), \
          '\t pdgID: %d'%( cc.pdgId()), \
          '\t dpT/pT = %.2f,  dR = %.2f '%(abs(cc.pt() - gg.vispt())/gg.vispt(), deltaR(cc.eta(),cc.phi(),gg.viseta(), gg.visphi()))))

def printGEN(gg):
    print(('GEN pt: %.2f'%(gg.vispt(),), \
          '\t eta: %.2f'%(gg.viseta()), \
          '\t phi: %.2f'%(gg.visphi()), '\n'))

   
##########################################################################################
# instantiate the handles to the relevant collections.

process_name = 'RECO'
if 'UL' or 'HNL' in sample:  process_name = 'PAT'

## gen particles
handle_gen              = Handle('std::vector<reco::GenParticle>')
## offline objects
handle_reco_taus        = Handle('std::vector<pat::Tau>')
handle_vtx              = Handle('std::vector<reco::Vertex>')
# L1 taus
handle_l1_tau = Handle('BXVector<l1t::Tau>')
handle_l1_ele = Handle('BXVector<l1t::Muon>')
## HLT objects 
handle_hlt_taus = Handle('std::vector<reco::PFTau>')
handle_hlt_ele  = Handle('std::vector<reco::RecoEcalCandidate>')
handle_beamspot  = Handle('reco::BeamSpot') 
## HLT ancillary variables
handle_IP_displ     = Handle('edm::AssociationVector<reco::PFTauRefProd,std::vector<reco::PFTauTransverseImpactParameterRef>>')
handle_iso_displ    = Handle('reco::PFTauDiscriminator')

## HLT filters
handle_hlt_filter     = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l1ele_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_l3ele_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 
handle_isoele_filter  = Handle('trigger::TriggerFilterObjectWithRefs') 

# handle_glbd_filter   = Handle('trigger::TriggerFilterObjectWithRefs') 



handles = {}
handles['gen_particles']     = [('prunedGenParticles', '', process_name) , handle_gen, False]
if 'gmsb' in sample:
    handles['gen_particles'] = [('genParticlePlusGeant', '', 'SIM') , handle_gen, False]


## L1 objects
handles['l1_eles'] = [('caloStage2Digis','EGamma','RECO'), handle_l1_ele, False]
handles['l1_taus'] = [('caloStage2Digis','Tau','RECO'), handle_l1_tau, False]
## HLT objects 
handles['hlt_taus'] = [('hltHpsPFTauProducerDispl', '', 'MYHLT'), handle_hlt_taus, False]
handles['hlt_eles'] = [('hltEgammaCandidates', '', 'MYHLT'), handle_hlt_ele, False]
handles['beamspot'] = [('hltOnlineBeamSpot', '', 'MYHLT'), handle_beamspot, False]
## HLT ancillary variables
handles['hlt_IP_displ']   = [('hltHpsPFTauTransverseImpactParameters', '', 'MYHLT'), handle_IP_displ, False]
handles['hlt_iso_displ']  = [('hltHpsDisplPFTauMediumAbsOrRelChargedIsolationDiscriminator', '', 'MYHLT'), handle_iso_displ, False]

## HLT filters
handles['hlt_filter']  = [('hltHpsDisplacedMuMediumChargedIsoDisplPFTau20TrackPt1L1HLTMatchedGlob', '', 'MYHLT'), handle_hlt_filter, False]

handles['hlt_l1ele_filter']  = [('hltL1sSingleEGNonIsoOrWithJetAndTauNoPS', '', 'MYHLT'), handle_l1ele_filter, False]
handles['hlt_eg_filter']  = [('hltEGL1SingleEGNonIsoOrWithJetAndTauNoPSFilter', '', 'MYHLT'), handle_l3ele_filter, False]
handles['hlt_isoele_filter'] = [('hltL3crIsoL1sSingleMu22L1f0L2f10QL3f24QL3trkIsoFiltered0p07', '', 'MYHLT'), handle_isoele_filter, False]
# handles['hlt_glbdisp_filter'] = [('hltL3fSingleMuL1f0L2NVf7L3GlbDispl10', '', 'MYHLT'), handle_glbd_filter, False]

## offline objects
handles['reco_taus']         = [('slimmedTaus', '', process_name),        handle_reco_taus, False]
handles['vtx']               = [('offlineSlimmedPrimaryVertices','','PAT'), handle_vtx, False]


##########################################################################################
# start looping on the events
for i, ev in enumerate(events):
    ######################################################################################
    # controls on the events being processed
    if maxevents>0 and i>maxevents:
        break
        
    if i%100==0:
        print(('===> processing %d / %d event' %(i, totevents)))
    
    for k, v in list(handles.items()):
        setattr(ev, k, None)
        v[2] = False
        try:
            ev.getByLabel(v[0], v[1])
            setattr(ev, k, v[1].product())
            v[2] = True
        except:    
            v[2] = False
    
    # select only hadronically decaying taus
    gen_had_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
                    pp.status()==good_gen_status and isGenHadTau(pp)]

    # select only taus decaying to electrons
    gen_ele_taus = [pp for pp in ev.gen_particles if abs(pp.pdgId())==15 and \
                    pp.status()==good_gen_status and isGenLepTau(pp, 11)]

    # skip events where we do not have a tau_h + tau_ele
    if len(gen_had_taus) < 1 or  len(gen_ele_taus) < 1 :  continue

    # info about tau_h
    for gg in gen_had_taus:
        gg.bestmom = None
        if gg.numberOfMothers() > 1:
            ## print 'more than 1 one, taking first one'
            tau_moms = [imom for imom in ev.gen_particles if isAncestor(imom, gg) and abs(imom.pdgId()) in mom_pdgId]
            gg.bestmom = tau_moms[0]
        elif gg.numberOfMothers() == 1 and abs(gg.mother(0).pdgId()) in mom_pdgId:
            gg.bestmom = gg.mother(0)

    # info about tau_ele
    for gg in gen_ele_taus:
        gg.bestmom = None
        if gg.numberOfMothers() > 1:
            ## print 'more than 1 one, taking first one'
            tau_moms = [imom for imom in ev.gen_particles if isAncestor(imom, gg) and abs(imom.pdgId()) in mom_pdgId]
            gg.bestmom = tau_moms[0]
        elif gg.numberOfMothers() == 1 and abs(gg.mother(0).pdgId()) in mom_pdgId:
            gg.bestmom = gg.mother(0)


    for gg in gen_ele_taus: print(('ele mom: ', gg.bestmom.pdgId()))
    for gg in gen_had_taus: print(('had mom: ', gg.bestmom.pdgId()))

    if len(gen_ele_taus) == 0 or  len(gen_had_taus) == 0 : 
      continue
    
    ## temporary!
    if len(gen_ele_taus) > 1 : 
      print('more than one gen ele tau')
      continue
    if len(gen_had_taus) > 1 : 
      print('more than one gen had tau')
      continue
      
    
    if gen_ele_taus[0].bestmom.pdgId() != -gen_had_taus[0].bestmom.pdgId() :
      print('tau_ele and tau_had coming from different moms')

    ######################################################################
    ## save variables for tau->ele
    for gg in gen_ele_taus:
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

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
        if gg.dau == None:  print ('is none')

        if gg.bestmom == None or gg.dau == None:
            continue
        
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()

        if 'taugun' in sample:
            gg.lxy = sqrt(pow(gg.dau.vx(),2)+pow(gg.dau.vy(),2))
        else:    
            gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))

        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))


    ######################################################################
    ## save variables for tau->had
    for gg in gen_had_taus:

        ## reset other attributes
        for ifeat in feat_list:
            setattr(gg, ifeat, -9999.)

        gg.decayMode = tauDecayModes.genDecayModeInt([d for d in finalDaughters(gg) \
                                                      if abs(d.pdgId()) not in [12, 14, 16]])
        if 'gmsb' in sample:
            dm_string = genDecayModeGEANT( [d for d in finalDaughters(gg) if abs(d.pdgId()) not in [12, 14, 16]])
            gg.decayMode = tauDecayModes.nameToInt(dm_string)
        ### find first dau to be used for vertices 
        gg.dau = None
        if gg.numberOfDaughters() > 0:  gg.dau = gg.daughter(0)
        if gg.dau == None:  print ('is none')
        gg.dxy = (gg.dau.vy()*gg.px() - gg.dau.vx()*gg.py())/gg.pt()
        gg.lxy = sqrt(pow(gg.vx()-gg.bestmom.vx(),2)+pow(gg.vy()-gg.bestmom.vy(),2))
        vectorP = np.array([gg.px(), gg.py(), 0])
        vectorL = np.array([gg.vx()-gg.bestmom.vx(), gg.vy()-gg.bestmom.vy(), 0])
        gg.cosxy = vectorL.dot(vectorP)/((np.linalg.norm(vectorL) * np.linalg.norm(vectorP)))


#     ######################################################################

    if handles['hlt_l1ele_filter'][2]:  
#         pdb.set_trace()
        ele_filter_product = ev.hlt_l1ele_filter.l1tegammaRefs()
        gen_ele_taus = findMatchToGen(gen_ele_taus, ele_filter_product, 'l1_ele')

    if handles['hlt_eg_filter'][2]:  
        ele_filter_product = ev.hlt_eg_filter.electronRefs()
        gen_ele_taus = findMatchToGen(gen_ele_taus, ele_filter_product, 'egamma')
# 
#     if handles['hlt_isoele_filter'][2]:  
#         ele_filter_product = ev.hlt_isoele_filter.eleonRefs()
#         gen_ele_taus = findMatchToGen(gen_ele_taus, ele_filter_product, 'iso_ele')
# 
#         for gg in gen_ele_taus:
#           if hasattr(gg, 'glb_ele') and gg.glb_ele:
#             gg.glb_ele.dxy = gg.glb_ele.track().dxy(ev.beamspot)


    # fill the ntuple: each gen tau makes an entry
    for gg in gen_ele_taus:
        if gg.bestmom == None or gg.dau == None:  continue

        for k, v in tofill_gen.items(): tofill_gen[k] = -9999. # initialise before filling
        tofill_gen['run'               ] = ev.eventAuxiliary().run()
        tofill_gen['lumi'              ] = ev.eventAuxiliary().luminosityBlock()
        tofill_gen['event'             ] = ev.eventAuxiliary().event()

        tofill_gen['gen_ele_ndau'      ] = gg.numberOfDaughters()
        tofill_gen['gen_ele_pt'        ] = gg.pt()
        tofill_gen['gen_ele_eta'       ] = gg.eta()
        tofill_gen['gen_ele_phi'       ] = gg.phi()
        tofill_gen['gen_ele_charge'    ] = gg.charge()
        tofill_gen['gen_ele_lxy'       ] = gg.lxy
        tofill_gen['gen_ele_dxy'       ] = gg.dxy
        tofill_gen['gen_ele_vx'        ] = gg.bestmom.vx()  ## user defined ones (mother prod. vertex)
        tofill_gen['gen_ele_vy'        ] = gg.bestmom.vy()
        tofill_gen['gen_ele_vz'        ] = gg.bestmom.vz()
        tofill_gen['gen_ele_cosxy'     ] = gg.cosxy
        tofill_gen['gen_ele_mom_mass'  ] = gg.bestmom.mass()
        tofill_gen['gen_ele_mom_pt'    ] = gg.bestmom.pt()
        tofill_gen['gen_ele_mom_eta'   ] = gg.bestmom.eta()
        tofill_gen['gen_ele_mom_phi'   ] = gg.bestmom.phi()
        tofill_gen['gen_ele_momct2d'   ] = gg.momct2d
        tofill_gen['gen_ele_vis_mass'  ] = gg.vismass()
        tofill_gen['gen_ele_vis_pt'    ] = gg.vispt()
        tofill_gen['gen_ele_vis_eta'   ] = gg.viseta()
        tofill_gen['gen_ele_vis_phi'   ] = gg.visphi()

        if hasattr(gg, 'l1_ele') and gg.l1_ele:
            tofill_gen['l1_ele_pt'       ] = gg.l1_ele.pt()
            tofill_gen['l1_ele_eta'      ] = gg.l1_ele.eta()
            tofill_gen['l1_ele_phi'      ] = gg.l1_ele.phi()
            tofill_gen['l1_ele_charge'   ] = gg.l1_ele.charge()

        if hasattr(gg, 'egamma') and gg.egamma:
            tofill_gen['hlt_ele_pt'       ] = gg.egamma.pt()
            tofill_gen['hlt_ele_eta'      ] = gg.egamma.eta()
            tofill_gen['hlt_ele_phi'      ] = gg.egamma.phi()
            tofill_gen['hlt_ele_charge'   ] = gg.egamma.charge()

#         if hasattr(gg, 'iso_ele') and gg.iso_ele:
#             tofill_gen['filter_ele_pt'       ] = gg.iso_ele.pt()
#             tofill_gen['filter_ele_eta'      ] = gg.iso_ele.eta()
#             tofill_gen['filter_ele_phi'      ] = gg.iso_ele.phi()
#             tofill_gen['filter_ele_charge'   ] = gg.iso_ele.charge()
# 
        ntuple_gen.Fill(array('f',tofill_gen.values()))

    ######################################################################################
    # printout some info, if you want
    # printer(taus, gen_taus)

##########################################################################################
# write the ntuples and close the files
outfile_gen.cd()
ntuple_gen.Write()
outfile_gen.Close()