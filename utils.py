'''
Define here utility functions to be imported in the main FWLite analyser.
'''

# identify hadronically decaying generator level taus
# def isGenHadTau(gen_tau):
#     ndau = gen_tau.numberOfDaughters()
#     # identify as had decaying if no electron/muon appears as daughter
#     isgentau = sum([abs(gen_tau.daughter(dd).pdgId()) in (11,13) for dd in range(ndau)]) == 0
#     # assign the gen tau a new attribut, visp4, 
#     # which is the visible 4-momentum of the generated tau.
#     # This is done by subtracting the neutrino p4's from the total gen tau p4.
#     gen_tau.visp4 = gen_tau.p4()
#     for neu in [gen_tau.daughter(dd) for dd in range(ndau) if abs(gen_tau.daughter(dd).pdgId()) in (12,14,16)]:
#         gen_tau.visp4 -= neu.p4()
#     gen_tau.vispt   = gen_tau.visp4.pt
#     gen_tau.viseta  = gen_tau.visp4.eta
#     gen_tau.visphi  = gen_tau.visp4.phi
#     gen_tau.vismass = gen_tau.visp4.mass
#     return isgentau


def isGenHadTau(gen_tau):
    ndau = gen_tau.numberOfDaughters()
    # identify as had decaying if no electron/muon appears as daughter
    isgentau = sum([abs(gen_tau.daughter(dd).pdgId()) in (11,13) for dd in range(ndau)]) == 0
    # assign the gen tau a new attribut, visp4, 
    # which is the visible 4-momentum of the generated tau.
    # This is done by adding the p4's of all charged hadrons/pi0/muons/ele
    for i,cha in enumerate([gen_tau.daughter(dd) for dd in range(ndau) if abs(gen_tau.daughter(dd).pdgId()) not in (12,14,16,22)]):
        if i==0:
            gen_tau.visp4 = cha.p4()
        else:
            gen_tau.visp4 += cha.p4()
    gen_tau.vispt   = gen_tau.visp4.pt
    gen_tau.viseta  = gen_tau.visp4.eta
    gen_tau.visphi  = gen_tau.visp4.phi
    gen_tau.vismass = gen_tau.visp4.mass
    return isgentau



# grab the final daughters (after any possible intermediate decay and radiative process)
# def finalDaughters(gen, daughters=None):
#     if daughters is None:
#         daughters = []
#     for i in range(gen.numberOfDaughters()):
#         daughter = gen.daughter(i)
#         if daughter.numberOfDaughters() == 0:
#             daughters.append(daughter)
#         else:
#             finalDaughters(daughter, daughters)
#     return daughters

# grab the final daughters (after any possible intermediate decay and radiative process)
def finalDaughters(gen, daughters=None):
    if daughters is None:
        daughters = []
    for i in range(gen.numberOfDaughters()):
        daughters.append(gen.daughter(i))
    return daughters


def genDecayModeGEANT(daughters):
    ''' Returns the generated tau decay mode based on a passed list of all
    final daughters before further decay (as contained in miniAOD).
    '''
    numElectrons = 0
    numMuons = 0
    numChargedHadrons = 0
    numNeutralHadrons = 0
    numPhotons = 0
    numNeutralPi = 0
    
    for daughter in daughters:
        pdg_id = abs(daughter.pdgId())
        if pdg_id == 22:
            numPhotons += 1
        elif pdg_id == 11:
            numElectrons +=1
        elif pdg_id == 13:
            numMuons += 1
        elif pdg_id == 111:
            numNeutralPi +=1
        else:
            if daughter.charge() != 0:
                numChargedHadrons += 1
            elif pdg_id not in [12, 14, 16]:
                numNeutralHadrons += 1
    
#     print 'numPhotons = ', numPhotons
#     print 'numElectrons = ', numElectrons
#     print 'numMuons = ', numMuons
#     print 'numChargedHadrons = ', numChargedHadrons
#     print 'numNeutralHadrons = ', numNeutralHadrons
#     print 'numNeutralPi = ', numNeutralPi
    
    if numElectrons == 1:
        return "electron"
    if numMuons == 1:
        return "muon"
    
    if numChargedHadrons == 1:
        if numNeutralHadrons != 0:
            return "kOneProngOther"
        if numNeutralPi == 0:
            return "kOneProng0PiZero"
        elif numNeutralPi == 1:
            return "kOneProng1PiZero"
        elif numNeutralPi == 2:
            return "kOneProng2PiZero"
        elif numNeutralPi == 3:
            return "kOneProng3PiZero"
        else:
            return "kOneProngNPiZero"
    elif numChargedHadrons == 3:
        if numNeutralHadrons != 0:
            return "kThreeProngOther"
        if numNeutralPi == 0:
            return "kThreeProng0PiZero"
        elif numNeutralPi == 1:
            return "kThreeProng1PiZero"
        elif numNeutralPi == 2:
            return "kThreeProng2PiZero"
        elif numNeutralPi == 3:
            return "kThreeProng3PiZero"
        else:
            return "kThreeProngNPiZero"
    
    return "kRareDecayMode"


def isAncestor(a, p):
    if a == p :
        return True
    for i in range(0,p.numberOfMothers()):
        if isAncestor(a,p.mother(i)):
            return True
    return False

# print tau reco/gen information in a nicely laid out format
def printer(taus, gen_taus):
    print('\n\n===========> RECO TAUS')
    for tt in taus: 
        print('reco tau         pt %.3f eta %.3f phi %.3f' %(tt.pt(), tt.eta(), tt.phi()))
        if hasattr(tt, 'gen_tau') and tt.gen_tau:
            print('matched gen tau  pt %.3f eta %.3f phi %.3f' %(tt.gen_tau.vispt(), tt.gen_tau.viseta(), tt.gen_tau.visphi()))
        else:
            print('matched gen tau  MISSING')
    print('===========> GEN TAUS')
    for gg in gen_taus: 
        print('gen tau          pt %.3f eta %.3f phi %.3f' %(gg.vispt(), gg.viseta(), gg.visphi()))
        if hasattr(gg, 'reco_tau') and gg.reco_tau:
            print('matched reco tau pt %.3f eta %.3f phi %.3f' %(gg.reco_tau.pt(), gg.reco_tau.eta(), gg.reco_tau.phi()))
        else:
            print('matched reco tau MISSING')