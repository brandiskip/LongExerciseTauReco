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
        hists[v][displaced] = getHist(trees[displaced], "h_%s_%s"%(v,displaced), v, getBinStr(variables[v]))

  
print("test 2")


# Everything below from old code
'''
histo_den = ROOT.TH1F('den', 'den', 10, 0, 10)
histo_num = ROOT.TH1F('num', 'num', 10, 0, 10)

trees[displaced].Draw('tau_gen_lxy >> den', 'tau_gen_vis_pt>20 & abs(tau_gen_vis_eta)<2.1')
trees[displaced].Draw('tau_gen_lxy >> num', 'tau_gen_vis_pt>20 & abs(tau_gen_vis_eta)<2.1 & tau_reco_pt>0')

c = ROOT.TCanvas("num")
c = ROOT.TCanvas("den")
c.cd()
eff = ROOT.TEfficiency(histo_num, histo_den)
eff.SetTitle('unmodified displaced #tau_{h}^{gen} ;transverse plane displacement (cm) ; reconstruction efficiency')
eff.SetMarkerStyle(8)
eff.Draw('AB')

c.SaveAs('reco_efficiency_{}.pdf'.format(sample))
c.SaveAs('reco_efficiency_{}.png'.format(sample))
c.Clear()
'''