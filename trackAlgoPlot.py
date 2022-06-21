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
tch.Add("tau_gentau_tuple_HLTstau200_0.root")

# define selections for generated taus
gen_sel  = 'tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode == 0'

# selections for reco tau and hlt
reco_sel = 'tau_reco_pt>0 && tau_reco_decaymode == 0'
hlt_sel  = 'tau_hltPFdispltau_pt>0 && tau_hltPFdispltau_decaymode == 0'

c = ROOT.TCanvas("c","c",600,500)
# den is the number of taus passing only generated tau selection
den = ROOT.TH1F("den", "den", 40, 0, 40)
den.SetFillColor(ROOT.kBlue)

# add axis labels
den.GetXaxis().SetTitle("algorithm index")

# fill histogram
tch.Draw("tau_reco_origAlgo >> den", gen_sel + "&&" + reco_sel + "&&" + hlt_sel)

# create box for legend
legend = ROOT.TLegend(.75, .75, .9, .95)
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.SetTextSize(0.03)
legend.AddEntry(tch,"L1 selection","l")
legend.Draw()

# save the current canvas
c.SaveAs('track_origAlgo.pdf')
c.SaveAs('track_origAlgo.png')
c.Clear()