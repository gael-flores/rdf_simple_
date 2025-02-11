import ROOT, math, os, itertools, glob, sys, ast
sys.path.append("/uscms/home/gfavila/nobackup/rdf_simple_/")
ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
import common.CMS_lumi as CMS_lumi
#from startVH import *
import startVH
from python.VHTools.VHcuts import *
from analysis.VH import *
from python.VHTools.plotEditor import plotter as editor
plotEditor = editor()
ROOT.gROOT.SetBatch(True)

import optparse
parser = optparse.OptionParser()
parser.add_option("-y","--year",dest="year",default="2018",help="2016 or 2017 or 2018")
parser.add_option("-d", "--date", dest="date", default="04_05_24",help="current date (used to name output directory)")
parser.add_option("-p", "--prod", dest="prod", default = "03_26_24", help="Date of production")
parser.add_option("-c", "--datacard", dest="datacard", default = "08_16_23", help="date for datacard directory")
parser.add_option("-m", "--mass", dest="mass" , type=int)
(options,args) = parser.parse_args()
year = options.year
date = options.date
prod = options.prod
mass = options.mass

# Defining output directories
path = "AN_{}".format(date)
if not os.path.isdir(path):
    os.mkdir(path)
os.chdir(path)
if not os.path.isdir("figs"):
    os.mkdir("figs")
os.chdir("/uscms/home/gfavila/nobackup/rdf_simple_/")

masses = [options.mass]

plotters = startVH.getPlotters(year, prod, "/uscms/home/gfavila/nobackup/rdf_simple_/DDP/",masses) #prod 03_26_24
dataEMU  = plotters['dataEMU']
signal   = plotters['signal']
dyPlotters = plotters['DYJets']
wjPlotters = plotters['WJets']



mcStack = {} # Combined plotter of mc bkg+sig
VHStack = {} # Combined plotter of mc bkg + data
mc = {} # merged plotter of all MC bkg as 1 contribution
for cat in ['zmm2g', 'zee2g', 'wen2g', 'wmn2g']:

    VHStack[cat] = startVH.combined_plotter()
    VHStack[cat].add_plotter(dyPlotters[cat], "DYJets", "DY+Jets", "background")
    VHStack[cat].add_plotter(wjPlotters[cat], "WJets", "W+Jets", "background")
    VHStack[cat].add_plotter(dataEMU[cat], "data_obs", "Data", "data")

    mcStack[cat] = {}
    mc[cat] = startVH.merged_plotter(dyPlotters[cat].plotters + wjPlotters[cat].plotters)
    mc[cat].setFillProperties(1001, ROOT.kAzure+5)
    mc[cat].setLineProperties(1, ROOT.kAzure+5, 3)

    # mc stack without signal
    mcStack[cat][0] = startVH.combined_plotter()
    mcStack[cat][0].add_plotter(dyPlotters[cat], "DYJets", "DY+Jets", "background")
    mcStack[cat][0].add_plotter(wjPlotters[cat], "WJets", "W+Jets", "background")

    
    for mass in masses:
        mcStack[cat][mass] = {}
        #mcStack[cat][channel][mass] = {}
        for ct in startVH.ctaus:
            signal['ZH'][cat][mass][ct]['2G2Q'].linecolor = ROOT.kBlue
            signal['ggZH'][cat][mass][ct]['2G2Q'].linecolor = ROOT.kMagenta
            signal['WH'][cat][mass][ct]['2G2Q'].linecolor = ROOT.kRed
            mcStack[cat][mass][ct] = startVH.combined_plotter()
            mcStack[cat][mass][ct].add_plotter(dyPlotters[cat], "DYJets", "DY+Jets", "background")
            mcStack[cat][mass][ct].add_plotter(wjPlotters[cat], "WJets", "W+Jets", "background")
            mcStack[cat][mass][ct].add_plotter(signal['ZH'][cat][mass][ct]['2G2Q'], 'ZH', "ZH Signal", "ZH signal")
            mcStack[cat][mass][ct].add_plotter(signal['ggZH'][cat][mass][ct]['2G2Q'], 'ggZH', "ggZH Signal", "ggZH signal")
            mcStack[cat][mass][ct].add_plotter(signal['WH'][cat][mass][ct]['2G2Q'], 'WH', "WH Signal", "WH signal")
            #mcStack[cat][mass][ct].add_plotter(signal['ttH'][cat][mass][ct]['2G2Q'], 'VH', "Signal", "signal")
sfs = {}
for v in ['W', 'Z']:
    sfs[v] = {}
    for l in ['ELE', 'MU']:
        sfs[v][l] = {}
        for m in [mass]:
            sfs[v][l][m] = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m])

# Make vertex related plots for section 3 of the AN
def plotSec3(year = "2018"):
    if not os.path.isdir("AN_{}/figs/vertex".format(date)):
        os.mkdir("AN_{}/figs/vertex".format(date))
    os.chdir("/uscms/home/gfavila/nobackup/rdf_simple_/DDP/")
    dir_out = "../AN_{}/figs/vertex/".format(date)

    # Plot Lxy distributions for various signal lifetimes
    colors = {0:ROOT.kBlack, 50:ROOT.kRed, 100:ROOT.kBlue, 1000:ROOT.kGreen+2}
    hists = {}
    for ct in [0, 50, 100, 1000]:
        signal_wmn2g = signal['WH']['wmn2g'][options.mass][ct]['2G2Q'].hist1d("best_2g_dxy_m30", cuts['photons'][options.mass]+"&&GenPart_nSignal==2", startVH.lumi[year], ("sig_dxy_ct{}".format(ct), "", 20, -10, 100))
        signal_wen2g = signal['WH']['wen2g'][options.mass][ct]['2G2Q'].hist1d("best_2g_dxy_m30", cuts['photons'][options.mass]+"&&GenPart_nSignal==2", startVH.lumi[year], ("sig_dxy_ct{}".format(ct), "", 20, -10, 100))
        #for channel in ['ZH','ggZH','ttH']:
        #    signal_wmn2g += signal[channel]['wmn2g'][options.mass][ct]['2G2Q'].hist1d("best_2g_dxy_m30", cuts['photons'][options.mass]+"&&GenPart_nSignal==2", startVH.lumi[year], ("sig_dxy_ct{}".format(ct), "", 20, -10, 100))
        #    signal_wen2g += signal[channel]['wen2g'][options.mass][ct]['2G2Q'].hist1d("best_2g_dxy_m30", cuts['photons'][options.mass]+"&&GenPart_nSignal==2", startVH.lumi[year], ("sig_dxy_ct{}".format(ct), "", 20, -10, 100))

        hists[ct] = signal_wmn2g
        wen2g = signal_wen2g
        hists[ct].Add(wen2g.GetValue())
        hists[ct].Scale(1.0/hists[ct].Integral())
        hists[ct].SetMarkerColor(colors[ct])
        hists[ct].SetLineColor(colors[ct])

    c = plotEditor.getCanvas("c")
    hists[0].Draw()
    hists[0].SetStats(0)
    hists[0].SetMinimum(0)
    hists[50].Draw("same")
    hists[100].Draw("same")
    hists[1000].Draw("same")

    hists[0].GetXaxis().SetTitle("L_{xy} [cm]")
    hists[0].GetYaxis().SetTitle("a.u.")

    leg = ROOT.TLegend(.65, .65, .9, .9)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for ct in [0, 50, 100, 1000]:
        leg.AddEntry(hists[ct].GetValue(), "c#tau = {} mm".format(ct), "l")

    leg.Draw("same")

    CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Simulation Work in progress")
    c.SaveAs(dir_out+"{}_signalLxy_comparison.pdf".format(year))

    ROOT.gInterpreter.Declare(
    '''
    float genLxy(RVecF vx, RVecF vy){
       float out = 0;
       if (vx.size() > 0)
          out = std::sqrt(vx[0]*vx[0]+vy[0]*vy[0]);
       return out;
    }
    ''')

    # 2d plot of signal resolution vs gen Lxy and fitSlicesY
    for v in ['W', 'Z']:
        means = {}
        sigmas = {}
        for m in masses:
            h = ROOT.TH2D("h{}".format(m), "", 10, 0, 100, 11, -11, 11)
            for ct in startVH.ctaus:
                for l in ['ELE', 'MU']:
                    signal['WH'][startVH.ana[v][l]][m][ct]['2G2Q'].define("GenSignal_lxy_m{}".format(m), "genLxy(GenPart_vx[GenPart_isSignal], GenPart_vy[GenPart_isSignal])")
                    signal['WH'][startVH.ana[v][l]][m][ct]['2G2Q'].define("deltaLxy_m{}".format(m), "best_2g_dxy_m{}-GenSignal_lxy_m{}".format(m,m))
                    if ct in [100, 1000]:
                        h_tmp = signal['WH'][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("GenSignal_lxy_m{}".format(m), "deltaLxy_m{}".format(m), cuts['photons'][m]+"&&GenPart_nSignal==2", startVH.lumi[year], ("gen_dxy_{}_{}_m{}_ct{}".format(v,l,m,ct), "", 10, 0, 100, 22, -11, 11))
                        h.Add(h_tmp.GetValue())
            h.SetStats(0)
            c = plotEditor.getCanvas2D("c")
            h.Scale(1.0/h.Integral())
            h.Draw("colz")
            h.GetXaxis().SetTitle("L_{xy}^{Gen} [cm]")
            h.GetYaxis().SetTitle("L_{xy}^{Gen} - L_{xy}^{reco} [cm]")
            h.GetZaxis().SetTitle("a.u.")
            h.GetXaxis().SetLabelSize(0.04)  # Set the label size for the X-axis
            h.GetXaxis().SetTitleSize(0.05)  # Set the title size for the X-axis
            h.GetYaxis().SetLabelSize(0.04)  # Set the label size for the Y-axis
            h.GetYaxis().SetTitleSize(0.05)  # Set the title size for the Y-axis
            CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Simulation Work in progress")
            c.SaveAs(dir_out+"{}_deltaLxyVsLxy_{}_m{}.pdf".format(year,v,m))
            h.FitSlicesY()
            mean = ROOT.gDirectory.Get("h{}_1".format(m))
            mean.SetMarkerStyle(22)
            mean.SetLineColor(ROOT.kBlack)
            mean.SetTitle("")
            mean.SetStats(0)
            c = plotEditor.getCanvas("c")
            mean.Draw("p")
            mean.GetXaxis().SetTitle("L_{xy}^{Gen} [cm]")
            mean.GetYaxis().SetTitle("Mean #Delta L_{xy} [cm]")
            mean.GetYaxis().SetRangeUser(-10, 10)
            CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Simulation Work in progress")
            c.SaveAs(dir_out+"{}_meanLxyVsLxy_{}_m{}.pdf".format(year,v,m))
            res = ROOT.gDirectory.Get("h{}_2".format(m))
            res.SetMarkerStyle(22)
            res.SetLineColor(ROOT.kBlack)
            res.SetMarkerColor(ROOT.kBlack)
            c = plotEditor.getCanvas("c")
            res.SetStats(0)
            res.Draw("p")
            res.SetTitle("")
            res.GetXaxis().SetTitle("L_{xy}^{Gen} [cm]")
            res.GetYaxis().SetTitle("#sigma_{Lxy} [cm]")
            res.GetYaxis().SetRangeUser(0, 14)
            CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Simulation Work in progress")
            ###plotEditor.drawPrelim(c)
            ###plotEditor.drawCOM(c, intLumi[year])
            c.SaveAs(dir_out+"{}_sigmaLxyVsLxy_{}_m{}.pdf".format(year,v,m))
            means[m] = mean.Clone("mean{}".format(m))
            sigmas[m] = res.Clone("sigma{}".format(m))

        # Plot means of deltaLxy vs Lxy
        colors = {7: ROOT.kOrange-3, 15: ROOT.kRed, 20: ROOT.kGreen+2, 30: ROOT.kCyan, 40: ROOT.kBlue, 50: ROOT.kViolet, 55: ROOT.kBlack}
        c = plotEditor.getCanvas("c")
        leg = ROOT.TLegend(.181, .627, .512, .884)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        first = True
        for m in masses:
            means[m].SetLineColor(colors[m])
            means[m].SetMarkerColor(colors[m])
            leg.AddEntry(means[m], "m_{#Phi}="+"{} GeV".format(m), "lp")
            if first:
                means[m].Draw("p")
                means[m].GetXaxis().SetTitle("L_{xy}^{Gen} [cm]")
                means[m].GetYaxis().SetTitle("Mean #Delta L_{xy} [cm]")
                means[m].GetYaxis().SetRangeUser(-10, 10)
                first = False
            else:
                means[m].Draw("same,p")
        leg.Draw("same")
        CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Simulation Work in progress")
        c.SaveAs(dir_out+"{}_meanLxyVsLxy_{}_all.pdf".format(year,v))
        c.SaveAs(dir_out+"{}_meanLxyVsLxy_{}_all.root".format(year,v))
        c.SaveAs(dir_out+"{}_meanLxyVsLxy_{}_all.png".format(year,v))

        # Plot sigmas of deltalxy vs lxy
        c = plotEditor.getCanvas("c")
        first = True
        for m in masses:
            sigmas[m].SetLineColor(colors[m])
            sigmas[m].SetMarkerColor(colors[m])
            if first:
                sigmas[m].Draw("p")
                sigmas[m].GetXaxis().SetTitle("L_{xy}^{Gen} [cm]")
                sigmas[m].GetYaxis().SetTitle("#sigma_{Lxy} [cm]")
                sigmas[m].GetYaxis().SetRangeUser(0, 10)
                first = False
            else:
                sigmas[m].Draw("same,p")
        CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Simulation Work in progress")

        leg.Draw("same")
        c.SaveAs(dir_out+"{}_sigmaLxyVsLxy_{}_all.pdf".format(year,v))
        c.SaveAs(dir_out+"{}_sigmaLxyVsLxy_{}_all.png".format(year,v))
        c.SaveAs(dir_out+"{}_sigmaLxyVsLxy_{}_all.root".format(year,v))

    # 1D masses under incorrect vertex assumptions
    colors = {7: ROOT.kOrange-3, 15: ROOT.kRed, 20: ROOT.kGreen+2, 30: ROOT.kCyan, 40: ROOT.kBlue, 50: ROOT.kBlue, 55: ROOT.kViolet}
    for v in ['W', 'Z']:
        for m1 in masses:
            stk = ROOT.THStack()
            leg = ROOT.TLegend(.191, .651, .493, .838)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            for m2 in [15, 20, 40, 55]:
                h = ROOT.TH1D("h{}".format(m), "", 25, -100, 100)
                for ct in startVH.ctaus:
                    for l in ['ELE', 'MU']:
                        if (m1 != m2):
                            signal['WH'][startVH.ana[v][l]][m1][ct]['2G2Q'].define("deltaLxy_m{}".format(m2), "best_2g_dxy_m{} - GenSignal_lxy_m{}".format(m2,m1))
                        tmp = signal['WH'][startVH.ana[v][l]][m1][ct]['2G2Q'].hist1d("deltaLxy_m{}".format(m2), cuts['photons'][m2]+"&&GenPart_nSignal==2", lumi[year], ("h", "", 25, -100, 100))
                        h.Add(tmp.GetValue())
                h.SetStats(0)
                c = plotEditor.getCanvas2D("c")
                if h.Integral() > 0:
                    h.Scale(1.0/h.Integral())
                h.Draw("colz")
                h.GetXaxis().SetTitle("L_{xy}^{Gen} [cm]")
                h.GetYaxis().SetTitle("L_{xy}^{Gen} - L_{xy}^{reco} [cm]")
                h.GetXaxis().SetLabelSize(0.04)  # Set the label size for the X-axis
                h.GetXaxis().SetTitleSize(0.05)  # Set the title size for the X-axis
                h.GetYaxis().SetLabelSize(0.04)  # Set the label size for the Y-axis
                h.GetYaxis().SetTitleSize(0.05)  # Set the title size for the Y-axis
                h.GetZaxis().SetTitle("a.u.")
                h.SetStats(0)
                ###plotEditor.drawPrelim(c)
                ###plotEditor.drawCOM(c, intLumi[year])
                CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Work in progress")
                h.SetMarkerStyle(22)
                h.SetLineColor(colors[m2])
                h.SetMarkerColor(colors[m2])
                h.SetLineWidth(3)
                leg.AddEntry(h, "Vertex mass = {} GeV".format(m2), "l")
                stk.Add(h)

            c = plotEditor.getCanvas("c")
            stk.Draw("hist,nostack")
            leg.Draw("same")
            stk.GetXaxis().SetTitle("L_{xy}^{reco}-L_{xy}^{Gen} [cm]")
            stk.GetYaxis().SetTitle("a.u.")
            ###plotEditor.drawPrelim(c)
            ###plotEditor.drawCOM(c, intLumi[year])
            CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Work in progress")
    
            c.SaveAs(dir_out+"{}_deltaLxy_{}_m{}.pdf".format(year, v, m1))


def plotSec4(year = "2018"):
    if not os.path.isdir("AN_{}/figs/selection".format(date)):
        os.mkdir("AN_{}/figs/selection".format(date))
    os.chdir("/uscms/home/gfavila/nobackup/rdf_simple_/DDP/")
    dir_out = "../AN_{}/figs/selection/".format(date)
    
    # Preselection data/MC agreement
    #for i in [1, 2]:
    """
    for i in [1]:
        sfx = "_med" if i==1 else "_tight"
        for l in ['ELE', 'MU']:
            lep = "#mu" if l == "MU" else "e"
            z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass", cuts['cr']['Z'][l][options.mass]+"&&best_2g_raw_mass_m{}>65".format(options.mass), startVH.lumi[year], ("z_mass_{}_{}".format(l,sfx), "", 20, 70, 110), SFs = sfs['Z'][l][options.mass], titlex = "m({}{}) [GeV]".format(lep,lep))
            #z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass",cuts['Z'][l], startVH.lumi[year], ("z_mass_{}_preselection{}".format(l,sfx), "", 20, 70, 110), SFs = sfs['Z'][l][options.mass], titlex = "m({}{}) [GeV]".format(lep,lep))
            
            #z['canvas'].SaveAs(year + "_ZX_Z_mass_{}_preselection{}.pdf".format(l,sfx))
            #z['canvas'].SaveAs(dir_out+year + "_ZX_Z_mass_{}_preselection{}_noCuts.pdf".format(l,sfx))
            z['canvas'].SaveAs(dir_out+year + "_ZX_Z_mass_{}{}_antiID_m{}.png".format(l,sfx,options.mass))
            
            w = VHStack[startVH.ana['W'][l]].draw_stack("W_mt", cuts['cr']['W'][l][options.mass]+"&&best_2g_raw_mass_m{}>65".format(options.mass), startVH.lumi[year], ("w_mt_{}_preselection{}".format(l,sfx), "", 10, 0, 200), SFs = sfs['W'][l][options.mass], titlex = "m_{T}"+"({}+MET) [GeV]".format(lep))
            #w = VHStack[startVH.ana['W'][l]].draw_stack("W_mt",cuts['W'][l] ,startVH.lumi[year], ("w_mt_{}_preselection{}".format(l,sfx), "", 10, 0, 200), SFs = sfs['W'][l][options.mass], titlex = "m_{T}"+"({}+MET) [GeV]".format(lep))
            #w['canvas'].SaveAs(dir_out+year + "_WX_W_mt_{}_preselection{}.pdf".format(l,sfx))
            w['canvas'].SaveAs(dir_out+year + "_WX_W_mt_{}{}_antiID_m{}.png".format(l,sfx,options.mass))
    

    for i in [1]:
        sfx = "_med" if i==1 else "_tight"
        for l in ['ELE', 'MU']:
            lep = "#mu" if l == "MU" else "e"
            z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass", cuts['Z'][l]+"&&"+cuts['photons'][options.mass]+"&&best_2g_sumID_m30=={}".format(i), startVH.lumi[year], ("z_mass_{}_preselection{}".format(l,sfx), "", 20, 70, 110), SFs = sfs['Z'][l][options.mass], titlex = "m({}{}) [GeV]".format(lep,lep))
            #z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass",cuts['Z'][l], startVH.lumi[year], ("z_mass_{}_preselection{}".format(l,sfx), "", 20, 70, 110), SFs = sfs['Z'][l][options.mass], titlex = "m({}{}) [GeV]".format(lep,lep))
            
            #z['canvas'].SaveAs(year + "_ZX_Z_mass_{}_preselection{}.pdf".format(l,sfx))
            #z['canvas'].SaveAs(dir_out+year + "_ZX_Z_mass_{}_preselection{}_noCuts.pdf".format(l,sfx))
            z['canvas'].SaveAs(dir_out+year + "_ZX_Z_mass_{}_preselection{}_antiID.png".format(l,sfx))
            
            w = VHStack[startVH.ana['W'][l]].draw_stack("W_mt", cuts['W'][l]+"&&"+cuts['photons'][options.mass]+"&&best_2g_sumID_m30=={}".format(i), startVH.lumi[year], ("w_mt_{}_preselection{}".format(l,sfx), "", 10, 0, 200), SFs = sfs['W'][l][options.mass], titlex = "m_{T}"+"({}+MET) [GeV]".format(lep))
            #w = VHStack[startVH.ana['W'][l]].draw_stack("W_mt",cuts['W'][l] ,startVH.lumi[year], ("w_mt_{}_preselection{}".format(l,sfx), "", 10, 0, 200), SFs = sfs['W'][l][options.mass], titlex = "m_{T}"+"({}+MET) [GeV]".format(lep))
            #w['canvas'].SaveAs(dir_out+year + "_WX_W_mt_{}_preselection{}.pdf".format(l,sfx))
            w['canvas'].SaveAs(dir_out+year + "_WX_W_mt_{}_preselection{}_antiID.png".format(l,sfx))
    
    # FSR plots
    for l in ['ELE', 'MU']:
        lep = "#mu" if l == 'MU' else 'e'
        
        # MC comparisons
        for m in [mass]:
            
            # Z mass
            z = mcStack[startVH.ana['Z'][l]][0].draw_comp("Z_mass", cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==1".format(m), ("z_mass_{}_preFSR_comp".format(m), "", 20, 70, 110), SFs=sfs['Z'][l][m], titlex="m({}{}) [GeV]".format(lep,lep), prelim="Simulation Work in progress")
            
            sig_ggZH = signal['ggZH'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_z_mass_{}_preFSR_comp".format(m), "", 20, 70, 110))
            
            sig_ZH = signal['ZH'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_z_mass_{}_preFSR_comp".format(m), "", 20, 70, 110))
            
            sig_WH = signal['WH'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_z_mass_{}_preFSR_comp".format(m), "", 20, 70, 110))
            
            sig_ttH = signal['ttH'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_z_mass_{}_preFSR_comp".format(m), "", 20, 70, 110))
    
            sig_ggZH.Scale(1.0/sig_ggZH.Integral() if sig_ggZH.Integral() > 0 else 1)
            sig_ZH.Scale(1.0/sig_ZH.Integral() if sig_ZH.Integral() > 0 else 1)
            sig_WH.Scale(1.0/sig_WH.Integral() if sig_WH.Integral() > 0 else 1)
            sig_ttH.Scale(1.0/sig_ttH.Integral() if sig_ttH.Integral() > 0 else 1)

            sig_ggZH.SetLineColor(ROOT.kMagenta)
            sig_WH.SetLineColor(ROOT.kRed)
            sig_ZH.SetLineColor(ROOT.kBlue)

            sig_ggZH.Draw("hist,same")
            sig_ZH.Draw("hist,same")
            sig_WH.Draw("hist,same")
            
            z['legend'].AddEntry(sig_ggZH.GetValue(), "ggZH Signal", "l")
            z['legend'].AddEntry(sig_ZH.GetValue(), "ZH Signal", "l")
            z['legend'].AddEntry(sig_WH.GetValue(), "WH Signal", "l")
            z['legend'].AddEntry(sig_ttH.GetValue(), "ttH Signal", "l")

            z['stack'].SetMaximum(1.1*max(z['stack'].GetMaximum(), sig_ZH.GetMaximum()))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_preFSR_mx{}_comp.pdf".format(l,m))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_preFSR_mx{}_comp.png".format(l,m))
            
            
            # Z+g1 mass
            z = mcStack[startVH.ana['Z'][l]][0].draw_comp("best_2g_fsr1_m{}".format(m), cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==1".format(m), ("zg1_mass_{}_preFSR_comp".format(l), "", 20, 70, 110), SFs=sfs['Z'][l][m], titlex="m({}{}".format(lep,lep)+"#gamma_{1}) [GeV]", prelim="Simulation Work in progress")
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("best_2g_fsr1_m{}".format(m), "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_zg1_mass_{}_preFSR_comp".format(m), "", 20, 70, 110)) 
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("best_2g_fsr1_m{}".format(m), "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_zg1_mass_{}_preFSR_comp".format(m), "", 20, 70, 110)) 
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("best_2g_fsr1_m{}".format(m), "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_zg1_mass_{}_preFSR_comp".format(m), "", 20, 70, 110)) 

            sig.Scale(1.0/sig.Integral() if sig.Integral() > 0 else 1)
            sig.Scale(1.0/sig.Integral() if sig.Integral() > 0 else 1)
            sig.Scale(1.0/sig.Integral() if sig.Integral() > 0 else 1)
            sig.Draw("hist,same")
            sig.Draw("hist,same")
            sig.Draw("hist,same")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            
            z['stack'].SetMaximum(1.1*max(z['stack'].GetMaximum(), sig.GetMaximum()))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Zg1_mass_{}_preFSR_mx{}_comp.pdf".format(l,m))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Zg1_mass_{}_preFSR_mx{}_comp.png".format(l,m))
            
            
            # Z+g2 mass
            z = mcStack[startVH.ana['Z'][l]][0].draw_comp("best_2g_fsr2_m{}".format(m), cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==1".format(m), ("zg2_mass_{}_preFSR_comp".format(l), "", 20, 70, 110), SFs=sfs['Z'][l][m], titlex="m({}{}".format(lep,lep)+"#gamma_{2}) [GeV]", prelim="Simulation Work in progress")
            
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("best_2g_fsr2_m{}".format(m), "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_zg2_mass_{}_preFSR_comp".format(m), "", 20, 70, 110))
            
            sig.Scale(1.0/sig.Integral() if sig.Integral() > 0 else 1)
            sig.Draw("hist,same")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            z['stack'].SetMaximum(1.1*max(z['stack'].GetMaximum(), sig.GetMaximum()))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Zg2_mass_{}_preFSR_mx{}_comp.pdf".format(l,m))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Zg2_mass_{}_preFSR_mx{}_comp.png".format(l,m))
            # Z+g1+g2 mass
            z = mcStack[startVH.ana['Z'][l]][0].draw_comp("best_2g_fsr3_m{}".format(m), cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==1".format(m), ("zg3_mass_{}_preFSR_comp".format(l), "", 20, 70, 110), SFs=sfs['Z'][l][m], titlex="m({}{}".format(lep,lep)+"#gamma_{1}#gamma_{2}) [GeV]", prelim="Simulation Work in progress")
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("best_2g_fsr3_m{}".format(m), "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("sig_zg3_mass_{}_preFSR_comp".format(m), "", 20, 70, 110))
            sig.Scale(1.0/sig.Integral() if sig.Integral() > 0 else 1)
            sig.Draw("hist,same")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            z['stack'].SetMaximum(1.1*max(z['stack'].GetMaximum(), sig.GetMaximum()))
            z['canvas'].SaveAs(dir_out+year+"_ZX_mass_{}_preFSR_mx{}_comp.pdf".format(l,m))
            z['canvas'].SaveAs(dir_out+year+"_ZX_mass_{}_preFSR_mx{}_comp.png".format(l,m))
            
            # FSR Tag distribution
            signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].define("best_2g_g1_isFSR_m{}".format(m), "Photon_isFSR[best_2g_idx1_m{}]".format(m))
            signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].define("best_2g_g2_isFSR_m{}".format(m), "Photon_isFSR[best_2g_idx2_m{}]".format(m))
            signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].define("best_2g_fsrTag_m{}".format(m), "3*(best_2g_g1_isFSR_m{m}&&best_2g_g2_isFSR_m{m}) + 2*(!best_2g_g1_isFSR_m{m}&&best_2g_g2_isFSR_m{m}) + 1*(best_2g_g1_isFSR_m{m}&&!best_2g_g2_isFSR_m{m})".format(m=m))
            mcStack[startVH.ana['Z'][l]][0].define("best_2g_g1_isFSR_m{}".format(m), "Photon_isFSR[best_2g_idx1_m{}]".format(m))
            mcStack[startVH.ana['Z'][l]][0].define("best_2g_g2_isFSR_m{}".format(m), "Photon_isFSR[best_2g_idx2_m{}]".format(m))
            mcStack[startVH.ana['Z'][l]][0].define("best_2g_fsrTag_m{}".format(m), "3*(best_2g_g1_isFSR_m{m}&&best_2g_g2_isFSR_m{m}) + 2*(!best_2g_g1_isFSR_m{m}&&best_2g_g2_isFSR_m{m}) + 1*(best_2g_g1_isFSR_m{m}&&!best_2g_g2_isFSR_m{m})".format(m=m))
            z = mcStack[startVH.ana['Z'][l]][0].draw_comp("best_2g_fsrTag_m{}".format(m), cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==1".format(m), ("fsr_tag_{}_m{}_comp".format(l,m), "", 4, 0, 4))
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("best_2g_fsrTag_m{}".format(m), cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}==2".format(m), "1", ("fsr_tag_{}_m{}_comp".format(l,m), "", 4, 0, 4))
            sig.Scale(1.0/sig.Integral())
            sig.Draw("hist,same")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            z['canvas'].SaveAs(dir_out+"2018_ZX_FSR_tag_{}_mx{}.pdf".format(l,m))
            z['canvas'].SaveAs(dir_out+"2018_ZX_FSR_tag_{}_mx{}.png".format(l,m))
            
            # Z mass post FSR, bkg/sig MC comparison
            z = mcStack[startVH.ana['Z'][l]][0].draw_comp("Z_mass", cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&"+cuts['fsr']['Z'][m]+"&&best_2g_sumID_m{}==1".format(m), ("z_mass_{}_postFSR_comp".format(l), "", 20, 70, 110), titlex="m({}{}) [GeV]".format(lep,lep), SFs = sfs['Z'][l][m])
            sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&"+cuts['fsr']['Z'][m]+"&&best_2g_sumID_m{}==2".format(m)+")*("+sfs['Z'][l][m]+")", "1", ("z_mass_{}_postFSR_comp".format(l), "", 20, 70, 110))
            sig.Scale(1.0/sig.Integral())
            sig.Draw("hist,same")
            z['legend'].AddEntry(sig.GetValue(), "Signal", "l")
            z['stack'].SetMaximum(max(z['stack'].GetMaximum(), sig.GetMaximum()))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_postFSR_mx{}_comp.pdf".format(l,m))
            z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_postFSR_mx{}_comp.png".format(l,m))
        
        # Data vs MC post FSR
        z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass", cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&"+cuts['fsr']['Z'][m]+"&&best_2g_sumID_m{}==2".format(m), startVH.lumi[year], ("z_mass_{}_postFSR".format(l), "", 20, 70, 110), SFs = sfs['Z'][l][m], titlex="m({}{}) [GeV]".format(lep, lep))
        z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_postFSR_tight.pdf".format(l))
        z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_postFSR_tight.png".format(l))

        z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass", cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&"+cuts['fsr']['Z'][m]+"&&best_2g_sumID_m{}==1".format(m), startVH.lumi[year], ("z_mass_{}_postFSR".format(l), "", 20, 70, 110), SFs = sfs['Z'][l][m], titlex="m({}{}) [GeV]".format(lep, lep))
        z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_postFSR_med.pdf".format(l))
        z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_postFSR_tight.png".format(l))
    
    # FSR printouts (exact numbers for AN)
    for l in ['ELE', 'MU']:
        lep = "e" if l=='ELE' else '#mu'
        for i in [1, 2]:
            mc_all = mc[startVH.ana['Z'][l]].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][options.mass]+"&&best_2g_sumID_m{}=={})*(".format(options.mass,i)+sfs['Z'][l][options.mass]+")", startVH.lumi[year], ("z_{}_sumid{}_all".format(l,i), "", 20, 70, 110))
            mc_fsr = mc[startVH.ana['Z'][l]].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][options.mass]+"&&"+cuts['fsr']['Z'][options.mass]+"&&best_2g_sumID_m{}=={})*(".format(options.mass,i)+sfs['Z'][l][options.mass]+")", startVH.lumi[year], ("z_{}_sumid{}_all".format(l,i), "", 20, 70, 110))
            dy_all = dyPlotters[startVH.ana['Z'][l]].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][options.mass]+"&&best_2g_sumID_m{}=={})*(".format(options.mass,i)+sfs['Z'][l][options.mass]+")", startVH.lumi[year], ("z_{}_sumid{}_all".format(l,i), "", 20, 70, 110))
            dy_fsr = dyPlotters[startVH.ana['Z'][l]].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][options.mass]+"&&"+cuts['fsr']['Z'][options.mass]+"&&best_2g_sumID_m{}=={})*(".format(options.mass,i)+sfs['Z'][l][options.mass]+")", startVH.lumi[year], ("z_{}_sumid{}_all".format(l,i), "", 20, 70, 110))
            print("SumID {} Z->{} FSR Cut Efficiency:".format(i, l))
            print("   DYJ: {}->{} %rejection={}".format(dy_all.Integral(), dy_fsr.Integral(), (dy_all.Integral()-dy_fsr.Integral())/dy_all.Integral()))
            print("   MC bkg: {}->{} %rejection={}".format(mc_all.Integral(), mc_fsr.Integral(), (mc_all.Integral()-mc_fsr.Integral())/mc_all.Integral()))

            for m in [mass]:
                sig = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{}=={})*(".format(m, i)+sfs['Z'][l][m]+")", startVH.lumi[year], ("sig_{}_sumid{}_m{}_all".format(l,i,m), "", 20, 70, 110))
                sig_fsr = signal['total'][startVH.ana['Z'][l]][m][0]['2G2Q'].hist1d("Z_mass", "("+cuts['Z'][l]+"&&"+cuts['photons'][m]+"&&"+cuts['fsr']['Z'][m]+"&&best_2g_sumID_m{}=={})*(".format(m, i)+sfs['Z'][l][m]+")", startVH.lumi[year], ("sig_{}_sumid{}_m{}_fsr".format(l,i,m), "", 20, 70, 110))
                print("   signal m={}: {}->{} %rejection={}".format(m, sig.Integral(), sig_fsr.Integral(), (sig.Integral()-sig_fsr.Integral())/sig.Integral()))

    
    
    # W misID cuts:
    dyj = dyPlotters['wen2g'].hist1d("W_mt", "(" + cuts['W']['ELE'] + ")*(" + cuts['photons'][options.mass]+")*("+sfs['W']['ELE'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("dy_misID", "", 20, 0, 300))
    dyj_fsr = dyPlotters['wen2g'].hist1d("W_mt", "(" + cuts['W']['ELE'] + ")*(" + cuts['photons'][options.mass] + ")*(" + sfs['W']['ELE'][options.mass]+")*(" + cuts['misID']['W']['ELE'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("dy_misID_cut", "", 20, 0, 300))

    mc_all = mc['wen2g'].hist1d("W_mt", "(" + cuts['W']['ELE'] + ")*(" + cuts['photons'][options.mass]+")*("+sfs['W']['ELE'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("mc_misID", "", 20, 0, 300))
    mc_fsr = mc['wen2g'].hist1d("W_mt", "(" + cuts['W']['ELE'] + ")*(" + cuts['photons'][options.mass] + ")*(" + sfs['W']['ELE'][options.mass]+")*(" + cuts['misID']['W']['ELE'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("mc_misID_cut", "", 20, 0, 300))

    # printouts for exact cut efficiencies
    print ("m {} GeV {} W->enu misID Cut efficiency:".format(mass, year))
    print("   DYJ: {}->{} %rejection={}".format(dyj.Integral(), dyj_fsr.Integral(), (dyj.Integral()-dyj_fsr.Integral())/dyj.Integral()))
    print("   MC Bkg: {}->{} %rejection={}".format(mc_all.Integral(), mc_fsr.Integral(), (mc_all.Integral()-mc_fsr.Integral())/mc_all.Integral()))

    for m in [mass]:
        sig = signal['total']['wen2g'][m][0]['2G2Q'].hist1d("W_mt", "(" + cuts['W']['ELE'] + ")*(" + cuts['photons'][options.mass]+")*("+sfs['W']['ELE'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("sig_misID", "", 20, 0, 300))
        sig_fsr = signal['total']['wen2g'][m][0]['2G2Q'].hist1d("W_mt", "(" + cuts['W']['ELE'] + ")*(" + cuts['photons'][options.mass] + ")*(" + sfs['W']['ELE'][options.mass]+")*(" + cuts['misID']['W']['ELE'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("sig_misID_cut", "", 20, 0, 300))
        print("   signal m={}: {}->{} %rejection={}".format(m, sig.Integral(), sig_fsr.Integral(), (sig.Integral()-sig_fsr.Integral())/sig.Integral()))

    """
    # W misID signal/MC comparison
    for m in [mass]:
        mc_misid1 = mcStack['wen2g'][m][0].draw_comp("best_2g_misID1_m{}".format(m), "("+cuts['W']['ELE']+")*("+cuts['photons'][m]+")*("+sfs['W']['ELE'][m]+")*(best_2g_sumID_m{}==2)".format(m), ("wenu_mass_l1x", "", 50, 0, 200), prelim="Simulation Work in progress", titlex="m(e#gamma_{1}) [GeV]")
        mc_misid1['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g1_mx{}_comp.pdf".format(m))
        mc_misid1['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g1_mx{}_comp.png".format(m))

        mc_misid2 = mcStack['wen2g'][m][0].draw_comp("best_2g_misID2_m{}".format(m), "("+cuts['W']['ELE']+")*("+cuts['photons'][m]+")*("+sfs['W']['ELE'][m]+")*(best_2g_sumID_m{}==2)".format(m), ("wenu_mass_l1x", "", 50, 0, 200), prelim="Simulation Work in progress", titlex="m(e#gamma_{2}) [GeV]")
        mc_misid2['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g2_mx{}_comp.pdf".format(m))
        mc_misid2['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g2_mx{}_comp.png".format(m))

        mc_misid3 = mcStack['wen2g'][m][0].draw_comp("best_2g_misID3_m{}".format(m), "("+cuts['W']['ELE']+")*("+cuts['photons'][m]+")*("+sfs['W']['ELE'][m]+")*(best_2g_sumID_m{}==2)".format(m), ("wenu_mass_l1x", "", 50, 0, 200), prelim="Simulation Work in progress", titlex="m(e#gamma_{1}#gamma_{2}) [GeV]")
        mc_misid3['canvas'].SaveAs(dir_out+year+"_WX_mass_l1X_mx{}_comp.pdf".format(m))
        mc_misid3['canvas'].SaveAs(dir_out+year+"_WX_mass_l1X_mx{}_comp.png".format(m))
    return
    """
    # W misID data/MC comparison
    canv = VHStack['wen2g'].draw_stack("best_2g_misID3_m30", "("+cuts['W']['ELE']+")*("+cuts['photons'][options.mass]+")*(best_2g_sumID_m{}==2)".format(options.mass), startVH.lumi[year], ("WX_mass_l1X", "", 50, 0, 200), titlex="m(e#gamma_{1}#gamma_{2}) [GeV]", SFs = sfs['W']['ELE'][options.mass], verbose= False)
    canv['canvas'].SaveAs(dir_out+year+"_WX_mass_l1X.pdf")
    canv['canvas'].SaveAs(dir_out+year+"_WX_mass_l1X.png")
    
    canv = VHStack['wen2g'].draw_stack("best_2g_misID1_m30", "("+cuts['W']['ELE']+")*("+cuts['photons'][options.mass]+")*(best_2g_sumID_m{}==2)*(abs(best_2g_misID3_m{}-91)>15)".format(m,m), startVH.lumi[year], ("WX_mass_l1g1", "", 50, 0, 200), titlex="m(e#gamma_{1}) [GeV]", SFs = sfs['W']['ELE'][options.mass], verbose = False)
    canv['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g1.pdf")
    canv['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g1.png")
    
    canv = VHStack['wen2g'].draw_stack("best_2g_misID2_m30", "("+cuts['W']['ELE']+")*("+cuts['photons'][options.mass]+")*(best_2g_sumID_m{}==2)*(abs(best_2g_misID3_m{}-91)>15)*(abs(best_2g_misID1_m{}-91)>10)".format(m,m,m), startVH.lumi[year], ("WX_mass_l1g2", "", 50, 0, 200), titlex="m(e#gamma_{2}) [GeV]", SFs = sfs['W']['ELE'][options.mass], verbose = False)
    canv['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g2.pdf")
    canv['canvas'].SaveAs(dir_out+year+"_WX_mass_l1g2.png")
    
    # Photon pt cuts: Signal vs Background
    for v in ['W', 'Z']:
        for m in [mass]:
            for l in ['ELE', 'MU']:
                for i in [1, 2]:
                    dataEMU[startVH.ana[v][l]].define("best_2g_pt{}_m{}".format(i, m), "Photon_pt[best_2g_idx{}_m{}]".format(i,m))
                    signal['total'][startVH.ana[v][l]][m][0]['2G2Q'].define("best_2g_pt{}_m{}".format(i,m), "Photon_pt[best_2g_idx{}_m{}]".format(i,m))
                    bkg = dataEMU[startVH.ana[v][l]].hist1d("best_2g_pt{}_m{}".format(i,m), cuts[v][l]+"&&"+cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{m}==1&&best_2g_raw_mass_m{m}<{mp5}".format(m=m, mp5 = m+5), "1", ("{}H_{}_m{}_g{}_pt_bkg".format(v,l,m,i), "", 20, 20, 100))
                    bkg.SetFillStyle(0)
                    sig = signal['total'][startVH.ana[v][l]][m][0]['2G2Q']
                    sig=sig.hist1d("best_2g_pt{}_m{}".format(i,m), "("+cuts[v][l]+"&&"+cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{m}==2&&best_2g_raw_mass_m{m}<{mp5}".format(m=m, mp5=m+5)+")*("+sfs[v][l][m]+")", "1", ("{}H_{}_m{}_g{}_pt_sig".format(v,l,m,i), "", 20, 20, 100))
                    bkg.Scale(1.0/bkg.Integral() if bkg.Integral() > 0 else 1.0)
                    sig.Scale(1.0/sig.Integral() if sig.Integral() > 0 else 1.0)
                    c = plotEditor.getCanvas("c")
                    bkg.Draw("hist")
                    sig.Draw("hist,same")
                    bkg.SetStats(0)
                    sig.SetStats(0)
                    bkg.GetYaxis().SetRangeUser(0, 1.05 * max(bkg.GetMaximum(), sig.GetMaximum()))
                    bkg.GetXaxis().SetTitle("#gamma_{{{}}}".format(i)+" p_{T} [GeV]")
                    bkg.GetYaxis().SetTitle("a.u.")
                    CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(startVH.lumi[year])/1000.), extraText = "Work in progress")
                    leg = ROOT.TLegend(0.5, 0.6, 0.92, 0.9)
                    leg.SetBorderSize(0)
                    leg.SetFillStyle(0)
                    leg.AddEntry(bkg.GetValue(), "Control Region Background", "l")
                    leg.AddEntry(sig.GetValue(), "Signal", "l")
                    leg.Draw("same")
                    c.SaveAs(dir_out+year+"_{}X_g{}_pt_mx{}_{}_comp.pdf".format(v,i,m,l))
                    c.SaveAs(dir_out+year+"_{}X_g{}_pt_mx{}_{}_comp.png".format(v,i,m,l))
                
                # Photon pt cut efficiency:
                bkg1 = dataEMU[startVH.ana[v][l]].hist1d("{}_phi".format(v), cuts[v][l]+"&&"+cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{m}==1&&best_2g_raw_mass_m{m}<{mp5}".format(m=m, mp5 = m+5), "1", ("{}H_{}_m{}_prePt_bkg".format(v,l,m), "", 1, -4, 4))
                bkg2 = dataEMU[startVH.ana[v][l]].hist1d("{}_phi".format(v), cuts[v][l]+"&&"+cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m]+"&&"+cuts['photons'][m]+"&&"+cuts['pt'][m]+"&&best_2g_sumID_m{m}==1&&best_2g_raw_mass_m{m}<{mp5}".format(m=m, mp5 = m+5), "1", ("{}H_{}_m{}_postPt_bkg".format(v,l,m), "", 1, -4, 4))

                sig1 = signal['total'][startVH.ana[v][l]][m][0]['2G2Q'].hist1d("{}_phi".format(v), "("+cuts[v][l]+"&&"+cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m]+"&&"+cuts['photons'][m]+"&&best_2g_sumID_m{m}==2&&best_2g_raw_mass_m{m}<{mp5}".format(m=m, mp5=m+5)+")*("+sfs[v][l][m]+")", "1", ("{}H_{}_m{}_prePt_sig".format(v,l,m), "", 1, -4, 4))
                sig2 = signal['total'][startVH.ana[v][l]][m][0]['2G2Q'].hist1d("{}_phi".format(v), "("+cuts[v][l]+"&&"+cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m]+"&&"+cuts['photons'][m]+"&&"+cuts['pt'][m]+"&&best_2g_sumID_m{m}==2&&best_2g_raw_mass_m{m}<{mp5}".format(m=m, mp5=m+5)+")*("+sfs[v][l][m]+")", "1", ("{}H_{}_m{}_prePt_sig".format(v,l,m), "", 1, -4, 4))

                b1 = bkg1.Integral()
                b2 = bkg2.Integral()
                s1 = sig1.Integral()
                s2 = sig2.Integral()

                print("Photon PT Cut Efficiency {} {} m={}:\n   bkg={}/{}={:.4f} sig={:.0f}/{:.0f}={:.4f}".format(v,l,m,b2, b1, 1.0*b2/b1 if b1>0 else 0, s2, s1, 1.0*s2/s1 if s1 > 0 else 0))
    
    # Final agreement plots overlaid with signal
    
    for l in ['ELE', 'MU']:
        lep = "e" if l == 'ELE' else '#mu'

        # Signal Region Plots
        z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass", cuts['sr']['Z'][l][options.mass], startVH.lumi[year], ("Z_{}_mass_final".format(l), "", 20, 70, 110), titlex="m({}{}) [GeV]".format(lep,lep), SFs = sfs['Z'][l][options.mass], verbose = False)
        #sig_z = signal['total'][startVH.ana['Z'][l]][options.mass][0]['2G2Q']
        #sig_z = sig_z.hist1d("Z_mass", "("+cuts['sr']['Z'][l][options.mass]+")*("+sfs['Z'][l][options.mass]+")", startVH.lumi[year], ("Z_{}_mass_sig_final".format(l), "", 20, 70, 110))
        #sig_z.Scale(0.2)
        #sig_z.Draw("hist,same")
        #z['legend'].AddEntry(sig_z.GetValue(), "Signal, BR(H#rightarrow#Phi#Phi)#times BR(#Phi#rightarrow#gamma#gamma)=0.1", "l")
        #z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_final_tight.pdf".format(l))
        z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_final_tight.png".format(l))

        w = VHStack[startVH.ana['W'][l]].draw_stack("W_mt", cuts['sr']['W'][l][options.mass], startVH.lumi[year], ("W_{}_mt_final".format(l), "", 10, 0, 200), titlex = "m_{T}"+"({}+MET) [GeV]".format(lep), SFs = sfs['W'][l][options.mass], verbose = False)
        #sig_w  = signal['total'][startVH.ana['W'][l]][options.mass][0]['2G2Q']
        #sig_w  = sig_w.hist1d("W_mt", "("+cuts['sr']['W'][l][options.mass]+")*("+sfs['W'][l][options.mass]+")", startVH.lumi[year], ("W_{}_mt_sig_final".format(l), "", 10, 0, 200))
        #sig_w.Scale(0.3)
        #sig_w.Draw("hist,same")
        #w['legend'].AddEntry(sig_w.GetValue(), "Signal, BR(H#rightarrow#Phi#Phi)#times BR(#Phi#rightarrow#gamma#gamma)=0.15", "l")
        #w['canvas'].SaveAs(dir_out+year+"_WX_W_mt_{}_final_tight.pdf".format(l))
        w['canvas'].SaveAs(dir_out+year+"_WX_W_mt_{}_final_tight.png".format(l))

        # Control Region Plots
        z = VHStack[startVH.ana['Z'][l]].draw_stack("Z_mass", cuts['cr']['Z'][l][options.mass], startVH.lumi[year], ("Z_{}_mass_final_loose".format(l), "", 20, 70, 110), titlex="m({}{}) [GeV]".format(lep,lep), SFs = sfs['Z'][l][options.mass], verbose = False)
        #sig_z.Scale(10)
        #sig_z.Draw("hist,same")
        #z['legend'].AddEntry(sig_z.GetValue(), "Signal, BR(H#rightarrow#Phi#Phi)#times BR(#Phi#rightarrow#gamma#gamma)=1", "l")
        #z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_final_loose.pdf".format(l))
        z['canvas'].SaveAs(dir_out+year+"_ZX_Z_mass_{}_final_loose.png".format(l))

        w = VHStack[startVH.ana['W'][l]].draw_stack("W_mt", cuts['cr']['W'][l][options.mass], startVH.lumi[year], ("W_{}_mt_final_loose".format(l), "", 10, 0, 200), titlex = "m_{T}"+"({}+MET) [GeV]".format(lep), SFs = sfs['W'][l][options.mass], verbose = False)
        #sig_w.Scale(8)
        #sig_w.Draw("hist,same")
        #w['legend'].AddEntry(sig_w.GetValue(), "Signal, 4#timesBR(H#rightarrow#Phi#Phi)#times BR(#Phi#rightarrow#gamma#gamma)=1", "l")
        #w['canvas'].SaveAs(dir_out+year+"_WX_W_mt_{}_final_loose.pdf".format(l))
        w['canvas'].SaveAs(dir_out+year+"_WX_W_mt_{}_final_loose.png".format(l))
    
    """
    """
    print("MC Background expected yields")
    # Final mc background cutflow:
    for v in ['W', 'Z']:
        for l in ['ELE', 'MU']:
            print("   {}->{}".format(v,l))
            preselect = mc[startVH.ana[v][l]].hist1d("{}_phi".format(v), "("+cuts[v][l]+"&&"+cuts['photons'][options.mass]+"&&best_2g_sumID_m{}==2".format(30)+")*("+sfs[v][l][options.mass]+")", startVH.lumi[year], ("preselect", "", 1, -10, 10))
            print("      preselection = {}".format(preselect.Integral()))
            misid_fsr = mc[startVH.ana[v][l]].hist1d("{}_phi".format(v), "("+cuts[v][l]+"&&"+cuts['photons'][options.mass]+"&&"+cuts['fsr'][v][options.mass]+"&&"+cuts['misID'][v][l][options.mass]+"&&best_2g_sumID_m{}==2".format(30)+")*("+sfs[v][l][options.mass]+")", startVH.lumi[year], ("misid_fsr", "", 1, -10, 10))
            print("      misid/fsr = {}".format(misid_fsr.Integral()))
            ptcuts = mc[startVH.ana[v][l]].hist1d("{}_phi".format(v), "("+cuts[v][l]+"&&"+cuts['photons'][options.mass]+"&&"+cuts['fsr'][v][options.mass]+"&&"+cuts['misID'][v][l][options.mass]+"&&"+cuts['pt'][options.mass]+"&&best_2g_sumID_m{}==2".format(30)+")*("+sfs[v][l][options.mass]+")", startVH.lumi[year], ("ptcuts", "", 1, -10, 10))
            print("      pt cuts = {}".format(ptcuts.Integral()))
            print("      Total background rejection: {}".format(1 - ptcuts.Integral()/preselect.Integral()))
    return
    # Only calculate for 2018 simulation
    if year != '2018':
        return

    # Section 4 signal efficiencies
    
    report_steps = {}
    report_steps['wen2g'] = {0: 'passed_HLT',
                             1: 'exactly_1_tight_electron',
                             2: 'at_least_2_preselection_photons'}
    report_steps['wmn2g'] = {0: 'passed_HLT',
                             1: 'exactly_1_tight_muon',
                             2: 'at_least_2_preselection_photons'}
    report_steps['zee2g'] = {0: 'passed_HLT',
                             1: 'leading_electron_ptOver35',
                             2: 'at_least_2_preselection_photons'}
    report_steps['zmm2g'] = {0: 'passed_HLT',
                             1: 'leading_mu_pt_over28' if year=='2017' else 'leading_mu_pt_over25',
                             2: 'at_least_2_preselection_photons'}                                

    hlt = {}
    hlt['2018'] = {'MU': 'HLT_IsoMu24',
                   'ELE': 'HLT_Ele32_WPTight_Gsf'}
    hlt['2017'] =  {'MU': 'HLT_IsoMu27',
                   'ELE': 'HLT_Ele32_WPTight_Gsf_L1DoubleEG'}
    hlt['2016'] = {'MU': 'HLT_IsoMu24||HLT_IsoTkMu24',
                   'ELE': 'HLT_Ele27_WPTight_Gsf'}
    # Get efficiencies for each cut, save output so this doesn't have to be rerun to remake plots (unless changes are made to efficiencies)
    '''
    yields = {}
    denoms = {}
    for nsig in [2, 4]:
        yields[nsig] = {}
        denoms[nsig] = {}
        for v in ['W', 'Z']:
            yields[nsig][v] = {}
            denoms[nsig][v] = {}
            for l in ['ELE', 'MU']:
                yields[nsig][v][l] = {}
                denoms[nsig][v][l] = {}
                for m in masses:
                    yields[nsig][v][l][m] = {}
                    denoms[nsig][v][l][m] = {}
                    for ct in ctaus:
                        print ("Doing signal cutflow for {}->{}, m={} ct={}".format(v,l,m,ct))
                        yields[nsig][v][l][m][ct] = {}
                        denoms[nsig][v][l][m][ct] = 0
                        raw_samples = glob.glob("/home/tyler/samples/DDP/rdf_samples/signalSamples/{year}/*{v}*H_HTo2LongLivedTo*MFF-{m}_ctau-{ct}cm*.root".format(year=year, prod=prod, v=v, m=m, ct=ct))
                        rdf_raw = ROOT.RDataFrame("Events", raw_samples)
                        rdf_raw = rdf_raw.Define("sample_isMC", "1")
                        rdf_raw = genAna(rdf_raw)
                        rdf_raw = muonAna(rdf_raw,  era = '2016preVFP' if year == '2016' else year)
                        rdf_raw = electronAna(rdf_raw,  era = '2016preVFP' if year == '2016' else year)

                        cutflow = rdf_raw.Filter("GenPart_nSignal=={}".format(nsig))
                        cutflow = cutflow.Filter(hlt[year][l], "passed_HLT")

                        if v=='Z' and l == 'MU':

                            cutflow = cutflow.Filter("Sum(Muon_charge[tight_muon]==1)>0 && Sum(Muon_charge[tight_muon]==-1)>0","Tight_muplus_muminus")
                            cutflow = makeZ(cutflow, "Muon")  
                            cutflow = cutflow.Filter("Z_mass>70&&Z_mass<110", "dimuon_mass_70to110")

                            #Apply Thresholds to the muon pts and cut on muon pf iso
                            ptThresh = 28 if year == '2017' else 25
                            cutflow = cutflow.Filter("(Muon_pt[Z_idx[0]]>{p}||Muon_pt[Z_idx[1]]>{p})".format(p=ptThresh), "leading_mu_pt_over{}".format(ptThresh))

                            #require just a super cluster
                            cutflow = cutflow.Filter("nPhoton>0","at_least_1_photon")

                            #Apply photon ID (no ISO)
                            cutflow = photonAna(cutflow, year)

                            # FSR Recovery with loose muons and preselection photons
                            cutflow=cutflow.Define("Photon_isFSR","fsr_recovery(Z_idx,Muon_pt, Muon_eta, Muon_phi, Muon_mass,Photon_pt,Photon_eta,Photon_phi,Photon_preselection)");
                            # Correct isolation for FSR
                            cutflow = cutflow.Define("Photon_pfRelIso03_fsrCorr", "correct_gammaIso_for_muons(Z_idx, Muon_pt, Muon_eta, Muon_phi, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_isFSR)")

                            ##### gamma gamma +X analysis ######
                            ####################################

                            #At least Two good photons
                            cutflow=cutflow.Filter('Sum(Photon_preselection==1)>1','at_least_2_preselection_photons')

                        elif v == 'Z' and l == 'ELE':
                            cutflow = cutflow.Filter("Sum(Electron_charge[tight_electron]==1)>0 && Sum(Electron_charge[tight_electron]==-1)>0","Tight_eplus_eminus")

                            #create the best Zmumu candidate and filter
                            cutflow = makeZ(cutflow, "Electron")  
                            cutflow = cutflow.Filter("Z_mass>70&&Z_mass<110", "dielectron_mass_70to110")  
                            #Apply Thresholds to the muon pts and cut on muon pf iso
                            ptThresh = 35
                            cutflow = cutflow.Filter("(Electron_pt[Z_idx[0]]>{p}||Electron_pt[Z_idx[1]]>{p})".format(p=ptThresh), "leading_electron_ptOver{}".format(int(ptThresh)))

                            #require just a super cluster
                            cutflow = cutflow.Filter("nPhoton>0","at_least_1_photon")

                            #Apply photon ID (no ISO)
                            cutflow = photonAna(cutflow, year)

                            # FSR Recovery with loose muons and preselection photons
                            cutflow=cutflow.Define("Photon_isFSR","fsr_recovery(Z_idx,Electron_pt, Electron_eta, Electron_phi, Electron_mass,Photon_pt,Photon_eta,Photon_phi,Photon_preselection)");
                            # Correct isolation for FSR
                            cutflow = cutflow.Define("Photon_pfRelIso03_fsrCorr", "correct_gammaIso_for_muons(Z_idx, Electron_pt, Electron_eta, Electron_phi, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_isFSR)")

                            ##### gamma gamma +X analysis ######
                            ####################################

                            #Al least Two good photons
                            cutflow = cutflow.Filter("Sum(Photon_preselection==1)>1", "at_least_2_preselection_photons")
                        elif v=='W' and l=='MU':
                            cutflow = cutflow.Filter("Muon_ntight==1", "exactly_1_tight_muon")
                            ptThresh = 28 if year=='2017' else 25
                            cutflow = cutflow.Filter("Sum(Muon_pt[tight_muon]>{})>0".format(ptThresh), "muon_pt_over{}".format(ptThresh))
                            cutflow = makeW(cutflow, "Muon")

                            cutflow = cutflow.Filter("nPhoton>0", "at_least_1_photon")

                            cutflow = photonAna(cutflow, year)

                            cutflow = cutflow.Filter('Sum(Photon_preselection==1)>1', "at_least_2_preselection_photons")

                        elif v=='W' and l=='ELE':
                            cutflow = cutflow.Filter("Electron_ntight==1", "exactly_1_tight_electron")
                            ptThresh = 35
                            cutflow = cutflow.Filter("Sum(Electron_pt[tight_electron]>{})>0".format(ptThresh), "electron_pt_over{}".format(ptThresh))
                            cutflow = makeW(cutflow, "Electron")

                            cutflow = cutflow.Filter("nPhoton>0", "at_least_1_photon")

                            cutflow = photonAna(cutflow, year)

                            cutflow = cutflow.Filter('Sum(Photon_preselection==1)>1', "at_least_2_preselection_photons")

                        r = cutflow.Report()
                        for i in report_steps[ana[v][l]]:
                            for cut in r:
                                if cut.GetName() == report_steps[ana[v][l]][i]:
                                    yields[nsig][v][l][m][ct][i] = cut.GetPass()
                                    if i == 0:
                                        denoms[nsig][v][l][m][ct] = cut.GetAll()

                        samples = glob.glob("/home/tyler/samples/DDP/rdf_samples/MC{}_{}/*{}*H*_M{}_ctau{}_*.root".format(year, prod, v, m, ct))
                        rdf = ROOT.RDataFrame(ana[v][l], samples)
                        rdf = rdf.Filter("GenPart_nSignal=={}".format(nsig))

                        rdf = rdf.Filter(cuts['photons'][m], "photon_kinematics")

                        rdf = rdf.Filter("best_2g_sumID_m{}==2".format(m), "photon_id")

                        rdf = rdf.Filter(cuts['misID'][v][l][m]+"&&"+cuts['fsr'][v][m], "fsr")

                        rdf = rdf.Filter(cuts['pt'][m], "photon_pt")

                        rdf = rdf.Filter("best_2g_raw_mass_m{}<{}".format(m, m+5), "lxy")

                        report = rdf.Report()

                        for i, cut in enumerate(report):
                            yields[nsig][v][l][m][ct][i+3] = cut.GetPass()
   
    with open('yields.txt', 'w') as data:
        data.write(str(yields))
    with open('denoms.txt', 'w') as data:
        data.write(str(denoms))
    '''
    # get data and yields from txt output files instead of redoing them
    with open('yields.txt') as f:
        data = f.read()
    yields = ast.literal_eval(data)

    with open('denoms.txt') as f:
        data = f.read()
    denoms = ast.literal_eval(data)
    
    lines = ROOT.TH1D("line", "", 6, 0, 36)
    lines.Fill(0, 1000000)
    lines.Fill(13, 1000000)
    lines.Fill(25, 1000000)
    lines.SetFillStyle(0)
    lines.SetMinimum(1)
    lines.SetLineStyle(2)
    lines.SetLineColor(ROOT.kBlack)
    colors = {0: ROOT.kBlack,
              1: ROOT.kViolet-8,
              2: ROOT.kMagenta-7,
              3: ROOT.kBlue,
              4: ROOT.kAzure+8,
              5: ROOT.kGreen+2,
              6: ROOT.kOrange-3,
              7: ROOT.kRed}

    for nsig in [2, 4]:
        for v in ['W', 'Z']:
            for l in ['ELE', 'MU']:
                stack = ROOT.THStack()
                graphs = {}
                axis = ROOT.TH1D("h", "", 6*6, 0, 6*6)
                for i in range(8):
                    num = ROOT.TH1D("num_{}_{}_step{}_yields".format(v,l,i), "", 6*6, 0, 6*6)
                    denom = ROOT.TH1D("denom_{}_{}_step{}_yields".format(v,l,i), "", 6*6, 0, 6*6)
                    for j, ct in enumerate(ctaus):
                        for k,m in enumerate(masses):
                            nBin = 6*j+k
                            num.SetBinContent(num.GetXaxis().FindBin(nBin), yields[nsig][v][l][m][ct][i])
                            denom.SetBinContent(denom.GetXaxis().FindBin(nBin), denoms[nsig][v][l][m][ct])
                    graph = ROOT.TGraphAsymmErrors(num, denom)
                    graph.SetLineColor(colors[i])
                    graph.SetMarkerColor(colors[i])
                    graph.SetLineWidth(3)
                    for n in range(1, axis.GetNbinsX()+1):
                        axis.GetXaxis().SetBinLabel(n, "{}".format(masses[(n-1)%6]))
                    graphs[i] = graph
                c = ROOT.TCanvas("c", "", 850, 600)
                axis.Draw()
                axis.SetStats(0)
                axis.GetXaxis().SetTitle("m_{#Phi} [GeV]")
                if nsig == 2:
                    axis.GetYaxis().SetTitle("N(events passing cuts) / N(#Phi#Phi#rightarrow2#gamma2q)")
                elif nsig == 4:
                    axis.GetYaxis().SetTitle("N(events passing cuts) / N(#Phi#Phi#rightarrow4#gamma)")
                for g in graphs:
                    graphs[g].Draw("pe,same")
                axis.GetYaxis().SetRangeUser(0.0001, 1)
                lines.Draw("hist,same")
                c.Update()
                a2 = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), 0, 6, 6, "-")
                a2.SetLabelOffset(-0.06)
                a2.ChangeLabel(1, -1, 0.03, -1, -1, -1, "c#tau=0 mm")
                a2.ChangeLabel(2, -1, 0.03, -1, -1, -1, "c#tau=10 mm")
                a2.ChangeLabel(3, -1, 0.03, -1, -1, -1, "c#tau=20 mm")
                a2.ChangeLabel(4, -1, 0.03, -1, -1, -1, "c#tau=50 mm")
                a2.ChangeLabel(5, -1, 0.03, -1, -1, -1, "c#tau=100 mm")
                a2.ChangeLabel(6, -1, 0.03, -1, -1, -1, "c#tau=1000 mm")
                a2.CenterLabels()
                a2.Draw()
                legend = ROOT.TLegend(0.1504, 0.19, 0.58, 0.46)
                legend.SetHeader(headers[v][l])
                legend.SetBorderSize(0)
                legend.SetFillStyle(0)
                legend.AddEntry(graphs[0], "HLT", "l")
                legend.AddEntry(graphs[1], "+Lepton preselection", "l")
                legend.AddEntry(graphs[2], "+Photon preselection", "l")
                legend.AddEntry(graphs[3], "+Kinematic cuts", "l")
                legend.AddEntry(graphs[4], "+Photon ID", "l")
                legend.AddEntry(graphs[5], "+FSR/misID", "l")
                legend.AddEntry(graphs[6], "+Tight p_{T}(#gamma) cuts", "l")
                legend.AddEntry(graphs[7], "+m_{#gamma#gamma} < m_{#Phi} + 5 GeV", "l")
                legend.Draw("same")
                c.SetLogy()
                c.SetGridy()
                CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(lumi[year])/1000.), extraText = "Simulation Work in progress")
                if nsig == 2:
                    c.SaveAs(dir_out+year+"_signal_2G2Q_{}_{}_efficiency_raw.pdf".format(v,l))
                    c.SaveAs(dir_out+year+"_signal_2G2Q_{}_{}_efficiency_raw.png".format(v,l))
                elif nsig == 4:
                    c.SaveAs(dir_out+year+"_signal_4G_{}_{}_efficiency_raw.pdf".format(v,l))
                    c.SaveAs(dir_out+year+"_signal_4G_{}_{}_efficiency_raw.png".format(v,l))

                    
    # Branching ratio limit scale factors
    if not os.path.isdir("AN_{}/figs/signal".format(date)):
        os.mkdir("AN_{}/figs/signal".format(date))
    os.chdir("/home/tyler/DDP/rdf_simple/")
    def getParamAndError(num1, denom1, num2, denom2):
        n1 = ROOT.TH1D("n1", "", 1, 0, 1)
        n2 = ROOT.TH1D("n2", "", 1, 0, 1)
        d1 = ROOT.TH1D("d1", "", 1, 0, 1)
        d2 = ROOT.TH1D("d2", "", 1, 0, 1)
        for i in range(num1):
            n1.Fill(.5)
        for i in range(num2):
            n2.Fill(.5)
        for i in range(denom1):
            d1.Fill(.5)
        for i in range(denom2):
            d2.Fill(.5)

        err1 = ROOT.TGraphAsymmErrors(n1, d1)
        err2 = ROOT.TGraphAsymmErrors(n2, d2)
        e1 = 1.0*num1/denom1
        e2 = 1.0*num2/denom2
        param = e2/(2.0*e1)
        err1_up = err1.GetErrorYhigh(0)
        err1_down = err1.GetErrorYlow(0)
        err2_up = err2.GetErrorYhigh(0)
        err2_down = err2.GetErrorYlow(0)

        err1_avg = (err1_up+err1_down)/2./e1
        err2_avg = (err2_up+err2_down)/2./e2

        err = math.sqrt(err1_avg**2+err2_avg**2)*param

        return [param, err]

        
    for v in ['Z', 'W']:
        for l in ['ELE', 'MU']:
            for m in masses:
                f = {}
                c = ROOT.TCanvas("c", "", 600, 600)
                c.SetHighLightColor(2)
                c.Range(-0.17119, .305227, 1.077244, .392534)
                c.SetFillColor(0)
                c.SetBorderMode(0)
                c.SetBorderSize(2)
                c.SetLeftMargin(0.1371237)
                c.SetRightMargin(0.06187291)
                c.SetTopMargin(0.08333333)
                c.SetBottomMargin(0.1163194)
                c.SetFrameBorderMode(0)
                c.SetFrameBorderMode(0)
                h = ROOT.TH1D("h", "", 1, 0, 1)
                h.Draw()
                h.SetTitle("")
                h.SetStats(0)
                h.GetYaxis().SetRangeUser(.4, 1.2)
                h.GetYaxis().SetTitle("f(#beta)")
                h.GetXaxis().SetTitle("#beta")
                leg = ROOT.TLegend(.2, .15, .58, .39)
                leg.SetBorderSize(0)
                leg.SetFillStyle(0)
                leg.SetHeader(headers[v][l]+", m_{#Phi} = "+"{} GeV".format(m))
                colors = {0: ROOT.kBlack,
                          10: ROOT.kRed,
                          20: ROOT.kOrange-3,
                          50: ROOT.kGreen,
                          100: ROOT.kBlue,
                          1000: ROOT.kMagenta}

                minY = 1.0
                maxY = 0.0
                for ct in ctaus:
                    e1 = 1.0*yields[2][v][l][m][ct][6]/yields[2][v][l][m][ct][1]
                    e2 = 1.0*yields[4][v][l][m][ct][6]/yields[4][v][l][m][ct][1]
                    param, err = 0, 0
                    if e1 == 0:
                        continue
                    if e2 != 0:
                        param, err = getParamAndError(yields[2][v][l][m][ct][6], yields[2][v][l][m][ct][1], yields[4][v][l][m][ct][6], yields[4][v][l][m][ct][1])
                    f[ct] = ROOT.TF1("{}_{}_{}_{}_2g".format(v,l,m,ct), "1/(1+([0]-1)*x)", 0, 1)
                    f[ct].SetParameter(0, param)
                    f[ct].SetLineColor(colors[ct])
                    f[ct].SetLineWidth(3)
                    f[ct].Draw("same")
                    leg.AddEntry(f[ct], "c#tau = {} mm".format(ct), "l")
                    if e2 > 0:
                        if 2.0*e1/e2 < minY:
                            minY = 2.0*e1/e2
                        if 2.0*e1/e2 > maxY:
                            maxY = 2.0*e1/e2
                h.SetMaximum(max(1.1*maxY, 1.2))
                h.SetMinimum(min(.6*minY, .35))

                leg.Draw("same")
                header = leg.GetListOfPrimitives().First()
                header.SetTextFont(62)
                c.Modified()
                c.SaveAs("AN_{}/figs/signal/".format(date)+"BR_{}_{}_{}.pdf".format(v,l,m))
                c.SaveAs("AN_{}/figs/signal/".format(date)+"BR_{}_{}_{}.png".format(v,l,m))
                
    
    # Cut based ID efficiency
    
    fCuts = ROOT.TFile("photonID.root")
    colors = {'preselection': ROOT.kBlack,
              'idNoIso': ROOT.kRed,
              'looseId': ROOT.kRed,
              'id_corr': ROOT.kBlue,
              'HOE': ROOT.kViolet,
              'Sieie': ROOT.kGreen+2,
              'ChIso': ROOT.kOrange-3,
              'NeuIso': ROOT.kBlue,
              'PhIso': ROOT.kRed}

    legend = {'preselection': "Preselection",
              'idNoIso': "Preselection + H/E + #sigma_{i#etai#eta}",
              'looseId': "Preselection + loose cut based ID",
              'id_corr': "Preselection + Nearby photon removed ID",
              'HOE': "H/E",
              'Sieie': "#sigma_{i#etai#eta}",
              'ChIso': 'Charged Hadron Isolation',
              'NeuIso': 'Neutral Hadron Isolation',
              'PhIso': "Photon Isolation"}
    for v in ['W', 'Z']:
        for m in masses:
            denom = fCuts.Get("denom_{}_m{}".format(v,m))
            nums = {}
            effs = {}
            for cat in ['preselection', 'looseId', 'id_corr']:
                nums[cat] = fCuts.Get("num_{}_m{}_{}".format(v, m, cat))
                effs[cat] = ROOT.TGraphAsymmErrors(nums[cat], denom)
                effs[cat].SetLineWidth(3)
                effs[cat].SetLineColor(colors[cat])
                effs[cat].SetMarkerColor(colors[cat])
            c = plotEditor.getCanvas("c")
            effs['preselection'].Draw()
            effs['looseId'].Draw("same")
            effs['id_corr'].Draw("same")

            effs['preselection'].GetXaxis().SetLimits(0, 100)
            effs['preselection'].GetYaxis().SetRangeUser(0, 1)
            effs['preselection'].GetXaxis().SetTitle("Gen L_{xy} [cm]")
            effs['preselection'].GetYaxis().SetTitle("Efficiency")

            leg = ROOT.TLegend(.11, .13, .64, .33)
            leg.SetHeader(headers[v][l]+", m_{#Phi} = "+"{} GeV".format(m))
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.032)
            for cat in ['preselection', 'looseId', 'id_corr']:
                leg.AddEntry(effs[cat], legend[cat], 'l')
            leg.Draw("same")
            #plotEditor.drawCOM(c, 59830)
            #plotEditor.drawPrelim(c)
            CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(lumi[year])/1000.), extraText = "Simulation Work in progress")
            c.SaveAs(dir_out+"cutBasedID_effVsLxy_{}_m{}_{}.pdf".format(v,m,options.year))
            c.SaveAs(dir_out+"cutBasedID_effVsLxy_{}_m{}_{}.png".format(v,m,options.year))
            #c.SaveAs("cutBasedID_effVsLxy_{}_m{}_{}.root".format(v,m,options.year))

            leg = ROOT.TLegend(.11, .13, .63, .37)
            leg.SetHeader(headers[v][l]+", m_{#Phi} = "+"{} GeV".format(m))
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.SetTextSize(0.028)
            leg.AddEntry(effs['preselection'], 'Preselection', 'l')

            c = plotEditor.getCanvas("c")
            effs['preselection'].Draw()
            
            for var in ['HOE', 'Sieie', 'ChIso', 'NeuIso', 'PhIso']:
                nums[var] = fCuts.Get("num_{}_m{}_pass{}".format(v, m, var))
                effs[var] = ROOT.TGraphAsymmErrors(nums[var], denom)
                effs[var].SetLineWidth(3)
                effs[var].SetLineColor(colors[var])
                effs[var].SetMarkerColor(colors[var])
                leg.AddEntry(effs[var], "Preselection + "+legend[var], 'l')
                effs[var].Draw("same")
            leg.Draw("same")
            CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = str(float(lumi[year])/1000.), extraText = "Simulation Work in progress")
            c.SaveAs(dir_out+"cutBasedID_effVsLxy_{}_m{}_cats_{}.pdf".format(v,m,options.year))
            c.SaveAs(dir_out+"cutBasedID_effVsLxy_{}_m{}_cats_{}.png".format(v,m,options.year))

    fCuts.Close()
    """             
    
# This one is just copying files made from the python/VHTools/drawClosurePlots.py script
def plotSec5(year, datacard_date = "04_05_24"):
    if not os.path.isdir("AN_{}/figs/background".format(date)):
        os.mkdir("AN_{}/figs/background".format(date))
    os.chdir("/home/tyler/DDP/rdf_simple/")
    dir_out = "AN_{}/figs/background/".format(date)

    closure = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/validation_{}.root".format(datacard_date, year))
    
    for v in ['Z']:
        for l in ['ELE', 'MU']:
            for m in masses:
                os.system("cp datacards_{}/closure_{}H_{}_m{}_sideband_{}.pdf {}".format(datacard_date, v, l, m, year, dir_out))
                os.system("cp datacards_{}/closure_{}H_{}_m{}_signal_{}.pdf {}".format(datacard_date, v, l, m, year, dir_out))
                os.system("cp datacards_{}/{}H_{}_m{}_data_{}.pdf {}".format(datacard_date, v, l, m, year, dir_out))


# Section 6 limits need the results from running combine and plots using combineTools
def plotSec6(year, datacard_date = "04_05_24"):
    if not os.path.isdir("AN_{}/figs/signal".format(date)):
        os.mkdir("AN_{}/figs/signal".format(date))
    os.chdir("/uscms/home/gfavila/nobackup/rdf_simple_/DDP/")
    dir_out = "../AN_{}/figs/signal/".format(date)

    # 2d model dependent limits
    g = ROOT.TGraph2D()
    n = 0
    masses = [15,20,30,40,50,55]
    ctaus = [0,10,20,50,100,1000]
    for ct in ctaus:
        for m in masses:

            f = ROOT.TFile("/uscms/home/gfavila/nobackup/rdf_simple_/datacards_12_20_24/higgsCombine.VH_m{}_ctau{}_{}.AsymptoticLimits.mH125.root".format(m,ct,year))
            limit = f.Get("limit")
            limit.GetEntry(2)
            #if limit.limit > 1.0:
            #    continue
            g.SetPoint(n, m, max(.01,ct), limit.limit * 0.5 * 0.1)
            n+=1
    c = plotEditor.getCanvas2D("c")
    c.SetLogy()
    c.SetLogz()
    g.Draw("colz")
    g.SetMaximum(1.0)
    g.SetTitle(";m_{#Phi} [GeV];c#tau [mm];BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)")
    plotEditor.drawPrelim(c, "Preliminary")
    plotEditor.drawCOM(c, startVH.intLumi[year])
    
    output_file = ROOT.TFile(dir_out+"limit_2D_mVsCtau_{}.root".format(year), "RECREATE")
    g.Write("limit_2D_graph")  # Save with a recognizable name
    c.Write("canvas")  # Save the canvas too
    output_file.Close()
    c.SaveAs(dir_out+"limit_2D_mVsCtau_{}.pdf".format(year))
    c.SaveAs(dir_out+"limit_2D_mVsCtau_{}.png".format(year))
    #c.SaveAs(dir_out+"limit_2D_mVsCtau_{}.root".format(year))
    
    legendDict = {('W', 'MU'): "W#rightarrow#mu#nu", ('W', 'ELE'): "W#rightarrow"+"e#nu",
                  ('Z', 'MU'): "Z#rightarrow#mu#mu", ('Z', 'ELE'): "Z#rightarrow"+"ee"}

    """
    # 1d limits
    for m in masses:
        f = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/limits_2g_{}.root".format(datacard_date, year))
        limit = f.Get("limitVsCtau_m{}_{}".format(m, year))
        limit.SaveAs(dir_out+"limitVsCtau_m{}_{}.pdf".format(m, year))
        limit.SaveAs(dir_out+"limitVsCtau_m{}_{}.png".format(m, year))
        
        # 1d limits categorized
        c = plotEditor.getCanvas("c")
        c.SetLogy()
        c.SetLogx()
        c.SetGrid()
        central = f.Get("central_m{}_{}".format(m,year))
        central.SetTitle("")
        central.SetMarkerSize(0)
        central.Draw()
        central.GetXaxis().SetTitle("c#tau [mm]")
        central.GetYaxis().SetTitle("BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)")
        central.GetYaxis().SetRangeUser(0.001, 1)
        central.GetXaxis().SetLimits(0.1, 1000)
        central.GetYaxis().SetNdivisions(10)
        leg = ROOT.TLegend(.12, .66, .6, .86)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetHeader("m_{#Phi} = "+"{} GeV, BR(#Phi#rightarrow#gamma#gamma) = 0.5".format(m))
        leg.AddEntry(central, "Combined expected limit", "l")
        for v in ['Z', 'W']:
            for l in ['ELE', 'MU']:
                g = f.Get("{}_{}_m{}_{}".format(v,l,m,year))
                g.Draw("same")
                leg.AddEntry(g, legendDict[(v,l)], "l")
        leg.Draw("same")
        plotEditor.drawPrelim(c)
        plotEditor.drawCOM(c, intLumi[year])
        c.SaveAs(dir_out+"limitVsCtau_m{}_{}_categorized.pdf".format(m,year))
        c.SaveAs(dir_out+"limitVsCtau_m{}_{}_categorized.png".format(m,year))
        f.Close()

        f = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/limits_2g_modelIndependent_{}.root".format(datacard_date, year))
        for v in ['Z', 'W']:
            limit = f.Get("limitVsCtau_{}H_m{}_{}".format(v,m,year))
            limit.SaveAs(dir_out+"limitVsCtau_{}H_m{}_{}_modelIndependent.pdf".format(v,m,year))
            limit.SaveAs(dir_out+"limitVsCtau_{}H_m{}_{}_modelIndependent.png".format(v,m,year))

    
    # 2d model independent limits
    for v in ['W', 'Z']:
        g = ROOT.TGraph2D()
        n = 0
        for ct in ctaus:
            for m in masses:
                f = ROOT.TFile("/uscms/home/gfavila/nobackup/rdf_simple_/datacards_11_13_24/higgsCombine.{}H_m{}_ctau{}_{}_modelIndependent.AsymptoticLimits.mH125.root".format(v,m,ct,year))
                limit = f.Get("limit")
                limit.GetEntry(2)
                g.SetPoint(n, m, max(.01,ct), 0.5*0.1*limit.limit)
                n+=1
        c = plotEditor.getCanvas2D("c")
        c.SetLogy()
        c.SetLogz()
        g.Draw("colz")
        g.SetTitle(";m_{#Phi} [GeV];c#tau [mm];#sigma"+"({}+#Phi#Phi)".format(v)+"#timesBR(#Phi#rightarrow#gamma#gamma) [pb]")
        plotEditor.drawPrelim(c, "Preliminary")
        plotEditor.drawCOM(c, startVH.intLumi[year])
        c.SaveAs(dir_out+"limit_2D_{}H_mVsCtau_{}_modelIndependent.pdf".format(v,year))
        c.SaveAs(dir_out+"limit_2D_{}H_mVsCtau_{}_modelIndependent.png".format(v,year))
        
    """
def plotSec6_ZH(year, datacard_date = "04_05_24"):

    if not os.path.isdir("AN_{}/figs/signal".format(date)):
        os.mkdir("AN_{}/figs/signal".format(date))
    os.chdir("/home/tyler/DDP/rdf_simple/")
    dir_out = "AN_{}/figs/signal/".format(date)

    # 2d model dependent limits
    g = ROOT.TGraph2D()
    n = 0
    for ct in ctaus:
        for m in masses:

            f = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/higgsCombine.ZH_m{}_ctau{}_{}.AsymptoticLimits.mH125.root".format(datacard_date,m,ct,year))
            limit = f.Get("limit")
            limit.GetEntry(2)
            #if limit.limit > 1.0:
            #    continue
            g.SetPoint(n, m, max(.01,ct), limit.limit * 0.5 * 0.1)
            n+=1
    c = plotEditor.getCanvas2D("c")
    c.SetLogy()
    c.SetLogz()
    g.Draw("colz")
    g.SetMaximum(1.0)
    g.SetTitle(";m_{#Phi} [GeV];c#tau [mm];BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)")
    plotEditor.drawPrelim(c, "Work in progress")
    plotEditor.drawCOM(c, intLumi[year])
    c.SaveAs(dir_out+"limit_2D_mVsCtau_ZH_{}.pdf".format(year))
    c.SaveAs(dir_out+"limit_2D_mVsCtau_ZH_{}.png".format(year))
    c.SaveAs(dir_out+"limit_2D_mVsCtau_ZH_{}.root".format(year))

    legendDict = {('W', 'MU'): "W#rightarrow#mu#nu", ('W', 'ELE'): "W#rightarrow"+"e#nu",
                  ('Z', 'MU'): "Z#rightarrow#mu#mu", ('Z', 'ELE'): "Z#rightarrow"+"ee"}

    # 1d limits
    for m in masses:
        f = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/limits_ZH_{}.root".format(datacard_date, year))
        limit = f.Get("limitVsCtau_m{}_{}".format(m, year))
        limit.SaveAs(dir_out+"limitVsCtau_ZH_m{}_{}.pdf".format(m, year))
        limit.SaveAs(dir_out+"limitVsCtau_ZH_m{}_{}.png".format(m, year))

        # 1d limits categorized
        c = plotEditor.getCanvas("c")
        c.SetLogy()
        c.SetLogx()
        c.SetGrid()
        central = f.Get("central_m{}_{}".format(m,year))
        central.SetTitle("")
        central.SetMarkerSize(0)
        central.Draw()
        central.GetXaxis().SetTitle("c#tau [mm]")
        central.GetYaxis().SetTitle("BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)")
        central.GetYaxis().SetRangeUser(0.001, 1)
        central.GetXaxis().SetLimits(0.1, 1000)
        central.GetYaxis().SetNdivisions(10)
        leg = ROOT.TLegend(.12, .66, .6, .86)
        leg.SetFillStyle(0)
        leg.SetBorderSize(0)
        leg.SetHeader("m_{#Phi} = "+"{} GeV, BR(#Phi#rightarrow#gamma#gamma) = 0.5".format(m))
        leg.AddEntry(central, "Combined expected limit", "l")
        for v in ['Z']:
            for l in ['ELE', 'MU']:
                g = f.Get("{}_{}_m{}_{}".format(v,l,m,year))
                g.Draw("same")
                leg.AddEntry(g, legendDict[(v,l)], "l")
        leg.Draw("same")
        plotEditor.drawPrelim(c, "Work in progress")
        plotEditor.drawCOM(c, intLumi[year])
        c.SaveAs(dir_out+"limitVsCtau_ZH_m{}_{}_categorized.pdf".format(m,year))
        c.SaveAs(dir_out+"limitVsCtau_ZH_m{}_{}_categorized.png".format(m,year))
        f.Close()

        f = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/limits_2g_modelIndependent_{}.root".format(datacard_date, year))
        for v in ['Z']:
            limit = f.Get("limitVsCtau_{}H_m{}_{}".format(v,m,year))
            limit.SaveAs(dir_out+"limitVsCtau_{}H_m{}_{}_modelIndependent.pdf".format(v,m,year))
            limit.SaveAs(dir_out+"limitVsCtau_{}H_m{}_{}_modelIndependent.png".format(v,m,year))

        
    # 2d model independent limits
    for v in ['Z']:
        g = ROOT.TGraph2D()
        n = 0
        for ct in ctaus:
            for m in masses:
                f = ROOT.TFile("/home/tyler/DDP/rdf_simple/datacards_{}/higgsCombine.{}H_m{}_ctau{}_{}_modelIndependent.AsymptoticLimits.mH125.root".format(datacard_date,v,m,ct,year))
                limit = f.Get("limit")
                limit.GetEntry(2)
                g.SetPoint(n, m, max(.01,ct), 0.5*0.1*limit.limit)
                n+=1
        c = plotEditor.getCanvas2D("c")
        c.SetLogy()
        c.SetLogz()
        g.Draw("colz")
        g.SetTitle(";m_{#Phi} [GeV];c#tau [mm];#sigma"+"({}+#Phi#Phi)".format(v)+"#timesBR(#Phi#rightarrow#gamma#gamma) [pb]")
        plotEditor.drawPrelim(c, "Work in progress")
        plotEditor.drawCOM(c, intLumi[year])
        c.SaveAs(dir_out+"limit_2D_{}H_mVsCtau_{}_modelIndependent.pdf".format(v,year))
        c.SaveAs(dir_out+"limit_2D_{}H_mVsCtau_{}_modelIndependent.png".format(v,year))

 
#plotSec6("Run2",year)   
plotSec4(year)
#plotSec5(year, options.datacard)
#plotSec6_ZH("Run2", options.datacard)
