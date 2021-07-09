'''
Example to read the ntuples and produce efficiency plots
'''

import ROOT
import numpy as np
import sys
import argparse
import cms_style
import time
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes

execfile("basic_plotting.py")

#Command "True" prevents plots from popping up while running. Make "false" to turn off.
ROOT.gROOT.SetBatch(True)
cms_style.setTDRStyle()

parser = argparse.ArgumentParser(description="Plot quantities from ntuples")
parser.add_argument("--file",
        choices=['ZTT','QCD'],
        required=False,
        help='Specify the sample you want to use for plotting')
args = parser.parse_args()
sample = args.file

'''
# Location of files in your directory
# * pulls all files with that name
floc = 'ntuples/'
fnames = [
	"taus_hnlSample_unmod/tau_gentau_tuple_HNL_M_10_*.root",
	"taus_hnlSample_mod/tau_gentau_tuple_HNL_M_10_ReMINI_*.root",
    ]

# definition to pull in files, "getTree" defined in basic_plotting.py
trees = {}
for fname in fnames:
    fullname = floc+fname
    tname = "tree"
    displaced = fname.split("_")[2].split("/")[0]
    trees[displaced] = getTree(fullname, tname)

# variable interested in plotting
variables = {
        "lxy": {"varname": "tau_gen_lxy", "nbins": 20, "xmin": 0, "xmax": 20, "label": "displacement"},
            }

# To get index of decay mode from PhysicsTools.Heppy.physicsutils.TauDecayModes
#prong = tauDecayModes.nameToInt('kOneProng0PiZero')
#print(prong)


hists = {}
for v in variables:

    hists[v] = {}
    for displaced in trees:
        # "getHist" defined in basic_plotting.py
        hists[v][displaced] = {}

        # Plotting with no prong cut
        #hists[v][displaced]["den"] = getHist(trees[displaced], "h_%s_%s_den"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1')
        #hists[v][displaced]["num"] = getHist(trees[displaced], "h_%s_%s_num"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_reco_pt>0')

        # Plotting 1-prong decay with reco_decaymodes
        #hists[v][displaced]["den"] = getHist(trees[displaced], "h_%s_%s_den"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode >= 0 && tau_gen_decaymode <= 4')
        #hists[v][displaced]["num"] = getHist(trees[displaced], "h_%s_%s_num"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_reco_pt>=0 && tau_gen_decaymode >= 0 && tau_gen_decaymode <= 4 && tau_reco_decaymode >= 0 && tau_reco_decaymode <= 4')

        #Plotting 3-prong decay
        #hists[v][displaced]["den"] = getHist(trees[displaced], "h_%s_%s_den"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode >= 10 && tau_gen_decaymode <= 14')
        #hists[v][displaced]["num"] = getHist(trees[displaced], "h_%s_%s_num"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_reco_pt>0 && tau_gen_decaymode >= 10 && tau_gen_decaymode <= 14 && tau_reco_decaymode >= 10 && tau_reco_decaymode <= 14')



eff = {}
for v in variables:

    eff[v] = {}
    for displaced in trees:
        eff[v][displaced] = ROOT.TEfficiency(hists[v][displaced]["num"], hists[v][displaced]["den"])
            
    # "plotEfficiencies defined in basic_plotting.py"
    plotEfficiencies(eff[v], "h_%s.pdf"%v, xlabel="transverse plane displacement (cm)", ylabel="reconstruction efficiency")
'''

##################################################################################################################################    
# Plotting combined lost tracks and packed PFCharged

# Location of files in your directory
# * pulls all files with that name
direct = 'ntuples/'
locations = [
    "taus_hnlSample_combined/tau_gentau_tuple_HNL_M_10_ReMINI_*.root",
    ]

# definition to pull in files, "getTree" defined in basic_plotting.py
treez = {}
for location in locations:
    completename = direct + location
    tname = "tree"
    displace = location.split("_")[2].split("/")[0]
    treez[displace] = getTree(completename, tname)

# variable interested in plotting
variable = {
        "lxy": {"varname": "tau_gen_lxy", "nbins": 20, "xmin": 0, "xmax": 20, "label": "displacement"},
            }

# To get index of decay mode from PhysicsTools.Heppy.physicsutils.TauDecayModes
#prong = tauDecayModes.nameToInt('kOneProng0PiZero')
#print(prong)


hist = {}
for v in variable:

    hist[v] = {}
    for displace in treez:
        # "getHist" defined in basic_plotting.py
        hist[v][displace] = {}

        # Plotting lost tracks and packed PFcandidates
        #hist[v][displace]["den"] = getHist(trees[displace], "c_%s_%s_den"%(v,displace), variable[v]["varname"], getBinStr(variable[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1')
        #hist[v][displace]["num"] = getHist(trees[displace], "c_%s_%s_num"%(v,displace), variable[v]["varname"], getBinStr(variable[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_com_pt>=0')

        # Plotting 1-prong decay with combined tracks
        hist[v][displace]["den"] = getHist(treez[displace], "c_%s_%s_den"%(v,displace), variable[v]["varname"], getBinStr(variable[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode >= 0 && tau_gen_decaymode <= 4')
        hist[v][displace]["num"] = getHist(treez[displace], "c_%s_%s_num"%(v,displace), variable[v]["varname"], getBinStr(variable[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_com_pt>=0 && tau_gen_decaymode >= 0 && tau_gen_decaymode <= 4')

effic = {}
for v in variable:

    effic[v] = {}
    for displace in treez:
        effic[v][displace] = ROOT.TEfficiency(hist[v][displace]["num"], hist[v][displace]["den"])
            
    # "plotEfficiencies defined in basic_plotting.py"
    plotEfficiencies(effic[v], "c_%s.pdf"%v, xlabel="transverse plane displacement (cm)", ylabel="reconstruction efficiency")