import ROOT
import numpy as np
import sys
import argparse
import cms_style
import time
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes

#Command "True" prevents plots from popping up while running. Make "false" to turn off.
ROOT.gROOT.SetBatch(True)
cms_style.setTDRStyle()

#open and read file
tch = ROOT.TChain("tree")
tch.Add("outputHLT/HLT_*.root")

# define selections for generated and l1 taus
pass_gen  = 'tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode == 0'

# selections for l1
pass_l1 = 'tau_l1_pt>20 && abs(tau_l1_eta)<2.1'
# selections for reco tau
pass_reco = 'tau_reco_pt>20'
# selections for hlt
pass_hlt  = 'tau_hltPFdispltau_pt>20'

# create selections that gives events fully reconstructed offline but fail hlt
pass_all_cut = pass_gen + "&&" + pass_l1 + "&&" + pass_reco + "&&" + pass_hlt
pass_all_but_hlt = pass_gen + "&&" + pass_l1 + "&&" + pass_reco + "&& !(" + pass_hlt + ")"
# print(pass_all_but_hlt)

c = ROOT.TCanvas("c","c",600,500)
# den is the number of taus passing only generated tau selection
den = ROOT.TH1F("den", "den", 40, 0, 40)
den.SetBinError(2, 2)
den.SetFillColor(ROOT.kBlue)

# add axis labels
den.GetXaxis().SetTitle("algorithm index")

# fill histogram
tch.Draw("tau_reco_origAlgo >> den", pass_all_but_hlt)

# create box for legend
legend = ROOT.TLegend(.75, .75, .9, .95)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.03)
legend.AddEntry(tch,"no HLT tau","l")
legend.Draw()

# save the current canvas
c.SaveAs('track_origAlgo.pdf')
c.SaveAs('track_origAlgo.png')
c.Clear()
