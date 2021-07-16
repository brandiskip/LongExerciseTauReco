'''
Plotting efficiency of combined tracks with reco tracks
'''

import ROOT
import numpy as np
import sys
import argparse
import cms_style
import glob
from array import array
from ROOT import TLatex
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes

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

# open the input file and get the tree
tch = ROOT.TChain("tree")
tch.Add("ntuples/taus_hnlSample_combined/tau_gentau_tuple_HNL_M_10_ReMINI_*.root")

histo_den = ROOT.TH1F('den', 'den', 20, 0, 20)
histo_num = ROOT.TH1F('num', 'num', 20, 0, 20)

################################################################################################################################

tch.Draw('tau_gen_lxy >> den', 'tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode == 0')
num_sel = ['tau_gen_lxy >> num', 'tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode == 0 && tau_up_com_pt>0', 'tau_gen_lxy >> num', 'tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode == 0 && tau_reco_pt>0']
for nn in num_sel:
	tch.Draw(nn)
	

################################################################################################################################
c = ROOT.TCanvas("num")
c = ROOT.TCanvas("den")
c.cd()
eff = ROOT.TEfficiency(histo_num, histo_den)

eff.SetTitle('displaced #tau_{h}^{gen} ;transverse plane displacement (cm) ; efficiency')
eff.Draw('alpe')

legend = ROOT.TLegend(.66, .64, .8, .88)
legend.SetHeader("displaced #tau_{h}^{gen}")
legend.SetBorderSize(0)
legend.SetFillStyle(0)
legend.Draw()


# save the current canvas
c.SaveAs('com_efficiency_{}.pdf'.format(sample))
c.SaveAs('com_efficiency_{}.png'.format(sample))
c.Clear()