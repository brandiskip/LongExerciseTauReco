# Read ntuples and create efficiency plots for displaced electrons coming from tau decay

import ROOT
import numpy as np
import sys
import argparse
import cms_style
import time
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes
from argparse import ArgumentParser

exec(compile(open("basic_plotting.py", "rb").read(), "basic_plotting.py", 'exec'))

#Command "True" prevents plots from popping up while running. Make "false" to turn off.
ROOT.gROOT.SetBatch(True)
cms_style.setTDRStyle()

parser = ArgumentParser()
parser.add_argument("--branch", help = "which branch to plot", default = 'pt', required = False)
parser.add_argument("--level", help = "which level electron to plot e.g. l1 or hlt", default = 'l1', required = False)
args = parser.parse_args()

branch = args.branch
level = args.level

# open input files
floc = 'ntuples/'
fnames = [
	"tau_eleSample/tau_ele_tuple_HLTstau200ele_gmsb_*.root",
    ]

# definition to pull in files, "getTree" defined in basic_plotting.py
trees = {}
for fname in fnames:
    fullname = floc+fname
    tname = "tree"
    displaced = fname.split("_")[0].split("/")[0]
    trees[displaced] = getTree(fullname, tname)
'''
if branch==eta:
	cut = 'gen_ele_pt>20 && abs(gen_ele_eta)<2.5 && abs(l1_ele_eta)<2.5 && l1_ele_pt>30'
else:
	cut = 'gen_ele_pt>20 && abs(gen_ele_eta)<2.5 && l1_ele_pt>30'

if level==hlt:
	if branch==eta:
		cut = 'gen_ele_pt>20 && abs(gen_ele_eta)<2.5 && abs(hlt_ele_eta)<2.5 && hlt_ele_pt>30'
	else:
		cut = 'gen_ele_pt>20 && abs(gen_ele_eta)<2.5 && hlt_ele_pt>30'
'''

# variable interested in plotting
variable = {
        "%s"%(level): {"varname": "gen_ele_%s"%(branch), "nbins": 100, "xmin": 0, "xmax": 500, "label": "%s reconstruction"%(branch)},
        }

hist = {}
for v in variable:

    hist[v] = {}
    for displaced in trees:
    	# "getHist" defined in basic_plotting.py
        hist[v][displaced] = {}

        hist[v][displaced]["den"] = getHist(trees[displaced], "c_%s_%s_den"%(v,displaced), variable[v]["varname"], getBinStr(variable[v]), sel_name='gen_ele_pt>20 && abs(gen_ele_eta)<2.5')
        hist[v][displaced]["num"] = getHist(trees[displaced], "c_%s_%s_num"%(v,displaced), variable[v]["varname"], getBinStr(variable[v]), sel_name='gen_ele_pt>20 && abs(gen_ele_eta)<2.5 && %s_ele_pt>30'%(level))

effic = {}
for v in variable:

    effic[v] = {}
    for displaces in trees:
        effic[v][displaced] = ROOT.TEfficiency(hist[v][displaced]["num"], hist[v][displaced]["den"])      

    # "plotEfficiencies defined in basic_plotting.py"
    plotEfficiencies(effic[v], "c_%s_%s.png"%(v,branch), xlabel="%s_ele_%s"%(level,branch), ylabel="efficiency")