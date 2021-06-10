'''
Example to read the ntuples and produce plots
'''

import ROOT
import numpy as np
import sys
import argparse
import cms_style
import time

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
    displaced = fname.split("_")[2]
    trees[displaced] = getTree(fullname, tname)

 
variables = {
        "tau_gen_lxy": {"varname": "tau_gen_lxy", "nbins": 10, "xmin": 0, "xmax": 10, "label": "displacement"},
            }
   
print("test")

hists = {}
for v in variables:

    hists[v] = {}
    for displaced in trees:
        if ('tau_gen_vis_pt>20 & abs(tau_gen_vis_eta)<2.1'):
            if ('tau_reco_pt>0'):
                print("cat")
                hists[v][displaced] = getHist(trees[displaced], "h_%s_%s"%(v,displaced), v, getBinStr(variables[v]), sel_name="num")
                trees[displaced].Draw('tau_gen_lxy >> num', 'tau_gen_vis_pt>20 & abs(tau_gen_vis_eta)<2.1 & tau_reco_pt>0')
            else :  
                hists[v][displaced] = getHist(trees[displaced], "h_%s_%s"%(v,displaced), v, getBinStr(variables[v]), sel_name="den")
                trees[displaced].Draw('tau_gen_lxy >> den', 'tau_gen_vis_pt>20 & abs(tau_gen_vis_eta)<2.1')
print("test 2")        

eff = {}        
for e in trees:
    eff[e] = ROOT.TEfficiency(num, den)
  
    plotEfficiencies(eff[e], "h_%s.pdf"%e, xlabel="transverse plane displacement (cm)", ylabel="reconstruction efficiency")

print("test 3")

