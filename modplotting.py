'''
Example to read the ntuples and produce plots
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

floc = 'ntuples/'
fnames = [
	"taus_hnlSample_unmod/tau_gentau_tuple_HNL_M_10_*.root",
	"taus_hnlSample_mod/tau_gentau_tuple_HNL_M_10_ReMINI_*.root",
]

trees = {}
for fname in fnames:
    fullname = floc+fname
    tname = "tree"
    displaced = fname.split("_")[2].split("/")[0]
    trees[displaced] = getTree(fullname, tname)

nameToInt('kOneProng0PiZero')
 
print("here")

variables = {
        "lxy":       {"varname": "tau_gen_lxy", "nbins": 20, "xmin": 0, "xmax": 20, "label": "displacement"},
        
            }
   
hists = {}
for v in variables:

    hists[v] = {}
    for displaced in trees:

        hists[v][displaced] = {}
        hists[v][displaced]["den"] = getHist(trees[displaced], "h_%s_%s_den"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1')         
        hists[v][displaced]["num"] = getHist(trees[displaced], "h_%s_%s_num"%(v,displaced), variables[v]["varname"], getBinStr(variables[v]), sel_name='tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_reco_pt>0')
                
        print(type(hists[v][displaced]["num"]))        


eff = {}        
for v in variables:
    
    eff[v] = {}
    for displaced in trees:
        eff[v][displaced] = ROOT.TEfficiency(hists[v][displaced]["num"], hists[v][displaced]["den"])
  
    plotEfficiencies(eff[v], "h_%s.pdf"%v, xlabel="transverse plane displacement (cm)", ylabel="reconstruction efficiency")

