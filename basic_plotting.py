import glob
import ROOT
from array import array
from ROOT import TLatex
from PhysicsTools.Heppy.physicsutils.TauDecayModes import tauDecayModes

colors = [ROOT.TColor.GetColor('#e23f59'), ROOT.TColor.GetColor('#eab508'), ROOT.TColor.GetColor('#7fb00c'), ROOT.TColor.GetColor('#8a15b1'),
                ROOT.TColor.GetColor('#57a1e8'), ROOT.TColor.GetColor('#e88000'), ROOT.TColor.GetColor('#1c587e'), ROOT.TColor.GetColor('#b8a91e'),
                ROOT.TColor.GetColor('#19cbaf'), ROOT.TColor.GetColor('#322a2d'),]

more_colors = [ROOT.kAzure+2, ROOT.kGreen+2, ROOT.kPink+4, ROOT.kOrange+10, ROOT.kOrange, ROOT.kSpring+7]

theme_colors = [ROOT.TColor.GetColor("#263EA8"),    # Dark Blue
                ROOT.TColor.GetColor("#60BC94"),    # Seafoam
                ROOT.TColor.GetColor("#F7AA02"),    # Saffron
                ROOT.TColor.GetColor("#E86C6D"),    # Coral
                ROOT.TColor.GetColor("#2D84C0"),    # ATLAS Blue
                ROOT.TColor.GetColor("#6A5ACD"),    # Slate blue
                ROOT.TColor.GetColor("#666666"),    # Slate blue
                ROOT.TColor.GetColor("#000000")]    # Slate blue

theme_colors2 = [ROOT.TColor.GetColor("#FFA900"),   # Sunny Orange
                ROOT.TColor.GetColor("#3629AC"),    # Dark Blue
                ROOT.TColor.GetColor("#2FC494"),    # Seafoam
                ROOT.TColor.GetColor("#F65866"),    # Pink
                ROOT.TColor.GetColor("#0E81C4"),    # Light Blue
                ]

# Get one of our trees -- takes an input with wildcards
def getTree(file_name, tree_name="trees_SR_highd0_", verbose=False):
    files = glob.glob(file_name)
    if len(files)==0:
        print "Couldn't find any files matching string", file_name
        return -1

    print "Found %s files"%len(files)
    t = ROOT.TChain(tree_name)
    for f in files:
        if verbose: print "Adding", f
        t.Add(f)
    return t

# Get one of our trees -- takes a list of files
def getTreeFromList(files, tree_name="trees_SR_highd0_"):
    if len(files)==0:
        print "Couldn't find any files matching string", file_name
        return -1

    print "Found %s files"%len(files)
    t = ROOT.TChain(tree_name)
    for f in files: t.Add(f)
    return t

# Get a histogram with TTree Draw
def getHist(t, hname, value_name, binstr="", sel_name="", binning="", verbose=False):
    #hname = hname.replace("(","-").replace(")","-")
    hname_nobin = hname
    if binning != "": 
        hname = hname + "(" + binning +")"
    print hname
    draw_cmd = '%s>>%s' % (value_name, hname)
    if binstr != "":
        draw_cmd += "(%s)" % (binstr)
    
    print draw_cmd
    print sel_name
    if sel_name != "": t.Draw( draw_cmd, '(%s)' % (sel_name) ,"COLZ")
    else: t.Draw( draw_cmd )
    if verbose: print "\"%s\", \"%s\""%(draw_cmd, sel_name)
    h = ROOT.gROOT.FindObject(hname_nobin)
    print h
    return h

# From a map of histograms (key as label name), plot them all on a canvas and save
def plotHistograms(h_map, save_name, xlabel="", ylabel="", interactive=False, logy=False, atltext=""):

    can = ROOT.TCanvas("can", "can")
    h_keys = h_map.keys()
    h_values = h_map.values()

    if len(h_keys)<1:
        print "No histograms found in h_map. Drawing blank canvas."
        return

    # Get maxes/mins of all hists
    maxy = 1.5*getMaximum(h_values)
    miny = 0
    if logy:
        maxy*=1e4
        miny = 1e-1

    # Draw histograms
    colorHists(h_values)
    h_map[h_keys[0]].GetXaxis().SetTitle(xlabel)
    h_map[h_keys[0]].GetYaxis().SetTitle(ylabel)
    h_values[0].SetMinimum(miny)
    h_values[0].SetMaximum(maxy)
    h_map[h_keys[0]].Draw("hist")
    if logy: ROOT.gPad.SetLogy(1)
    for k in h_keys[1:]:
        h_map[k].Draw("hist same")


    leg = ROOT.TLegend(.66, .64, .8, .88)
    for k in h_keys:
        leg.AddEntry(h_map[k], k, "l")
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.Draw()
    if interactive: raw_input("...")

    if atltext != "":
        ROOT.ATLASLabel(0.2,0.85, atltext[0])
        text = ROOT.TLatex()
        text.SetNDC()
        text.SetTextSize(0.04)
        for i, t in enumerate(atltext[1:]):
            text.DrawLatex(0.2,0.85-0.07*(i+1), t)

    can.SaveAs(save_name)
    return

# From a map of efficiencies (key as label name), plot them all on a canvas and save
def plotEfficiencies(eff_map, save_name, xlabel="", ylabel=""):

    can = ROOT.TCanvas("can", "can")
    if len(eff_map)<1:
        print "No efficiency plots found in eff_map. Drawing blank canvas."
        return

    colorHists(eff_map.values())

    for i, k in enumerate(eff_map):
        print k
        if i==0:
            eff_map[k].Draw("alpe")
            eff_map[k].SetTitle(";%s;%s"%(xlabel,ylabel))
            ROOT.gPad.Update()

            eff_map[k].GetPaintedGraph().GetXaxis().SetTitle(xlabel)
            eff_map[k].GetPaintedGraph().GetYaxis().SetTitle(ylabel)
            eff_map[k].GetPaintedGraph().SetMinimum(0)
            eff_map[k].GetPaintedGraph().SetMaximum(1)
            eff_map[k].Draw("alpe")
            ROOT.gPad.Update()

        else: eff_map[k].Draw("lpe same")
        ROOT.gPad.Update()

    
    leg = ROOT.TLegend(.66, .74, .8, .88)
    for k in eff_map:
        leg.AddEntry(eff_map[k], k, "lp")
    leg.Draw()

    can.SaveAs(save_name)
    return
    

def colorHists(hists):
    i=0
    for h in hists:
        h.SetLineColor(more_colors[i])
        h.SetMarkerColor(more_colors[i])
        i += 1
        i %= len(more_colors)
    return


def getBinStr(entry):
    binstr = "%d,%f,%f"%(entry["nbins"],entry["xmin"],entry["xmax"])
    return binstr


def getHistFromEff(eff):
    
    h_out = eff.GetPassedHistogram().Clone()
    h_out.Clear()

    for i in xrange(h_out.GetNbinsX()+1):
        h_out.SetBinContent(i,eff.GetEfficiency(i))
        h_out.SetBinError(i,eff.GetEfficiencyErrorUp(i))
    return h_out

