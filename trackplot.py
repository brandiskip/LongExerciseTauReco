import ROOT
import cms_style
from copy import deepcopy as dc

# creates list of colors
colors = [ROOT.kBlue, ROOT.kGreen, ROOT.kRed]

#Command "True" prevents plots from popping up while running. Make "false" to turn off.
ROOT.gROOT.SetBatch(True)
cms_style.setTDRStyle()

#open and read file
tch = ROOT.TChain("tree")
tch.Add("ntuples/taus_hnlSample_combined/tau_gentau_tuple_HNL_M_10_ReMINI_*.root")

#define fiducial selection (types of taus we care about)
fid_sel = 'tau_gen_vis_pt>20 && abs(tau_gen_vis_eta)<2.1 && tau_gen_decaymode == 0'

#define cut selections (plot multiple cut selections)
cut_sels = {
	"reco taus": "tau_reco_pt>0",
	"com taus" : "tau_up_com_pt>0",
	}

# den is the number of taus passing only fiducial selection
den = ROOT.TH1F("den", "den", 20, 0, 20)
# fill histogram
tch.Draw("tau_gen_lxy >> den", fid_sel)

# define empty dictionary to fill with efficiency plots
effs_con = {}

#loop over cut selections
for c in cut_sels:
	cut_sel = cut_sels[c]
	# num is the number of taus passing the cut selection plus fiducial selection
	# make object
	num = ROOT.TH1F("num " + c, "num " + c, 20, 0, 20)
	# fill histogram
	tch.Draw("tau_gen_lxy >> num " + c, cut_sel + " && " + fid_sel)
	#	efficiency is num/den
	eff = ROOT.TEfficiency(num, den)
	effs_con[c] = dc(eff)

#plot efficiencies on canvas
# inside () is just the name you give the canvas
c = ROOT.TCanvas()
# create box for legend
legend = ROOT.TLegend(.66, .64, .8, .88)
legend.SetBorderSize(0)
legend.SetFillStyle(0)

# loop on efficiency dictionary to create plots
# e is the string of the key
for i,e in enumerate(effs_con):
	eff_con = effs_con[e]
	eff_con.SetLineColor(colors[i])
	eff_con.SetMarkerColor(colors[i])
	if i==0:
		eff_con.Draw('alpe') # a creates the axis
	else:
		eff_con.Draw('lpe same') # if use a again, will remove axis
	legend.AddEntry(eff_con,e)

legend.Draw()

# save the current canvas
c.SaveAs('com_efficiency.pdf')
c.SaveAs('com_efficiency.png')
c.Clear()

