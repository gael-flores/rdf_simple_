# Draw closure plots for a given mass hypothesis
import sys
sys.path.append("/uscms/home/gfavila/nobackup/rdf_simple_/")
import common.tdrstyle as tdrstyle
import common.CMS_lumi as CMS_lumi
from python.VHTools.config import *
from python.VHTools.py_helpers import *
sty = tdrstyle.setTDRStyle()
import ROOT, math, os
#from python.plotEditor import plotter as plotEditor
#pe = plotEditor()
ROOT.gROOT.SetBatch(True)
import numpy as np

import optparse
parser = optparse.OptionParser()
parser.add_option("-e","--era",dest="era",default="2018",help="2016 or 2017 or 2018")
parser.add_option("-d", "--date", dest="date", default="03_26_24",help="current date")
parser.add_option("-b", "--blind", dest="blind", action="store_true",help="current date")
parser.add_option("-m", "--mass", dest = "mass", default = '15')
drawprelim = True

(options,args) = parser.parse_args()

date = options.date
if not os.path.isdir("datacards_{}".format(date)):
    os.mkdir("datacards_{}".format(date))
os.chdir("datacards_{}".format(date))

if not os.path.isdir("plots/closure"):
    os.mkdir("plots/closure")

dir_to_save = f"/uscms/home/gfavila/nobackup/rdf_simple_/datacards_{date}/plots/closure/"

# Poisson points for histogram that will be scaled (needs x errors for shaded error bands)
def getPoisson2(h, scale):
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        #if N == 0:
        #    continue
        gRate.SetPoint(n, thresh, scale*N)
        error = fetchError(q, N)
        gRate.SetPointError(n, h.GetBinWidth(i)/2.0, h.GetBinWidth(i)/2.0, scale*(N-error[0]), scale*(error[1]-N))
        n+=1
    return gRate

# Comparing Lxy distribution in negative Lxy sideband
def compareSideband(loose, tight, m, name, header):
    loose.SetLineColor(ROOT.kAzure+5)
    loose.SetLineWidth(3)
    loose.SetFillColor(ROOT.kAzure+5)
    loose.SetFillStyle(1001)
    c = ROOT.TCanvas("c", "", 600, 600)
    c.SetLeftMargin(0.2)
    c.SetRightMargin(.05)
    c.SetTopMargin(0.08506945)
    c.SetBottomMargin(0.1145833)

    tight_pois = getPoisson(tight)
    tight_pois.SetLineColor(ROOT.kBlack)
    tight_pois.SetMarkerColor(ROOT.kBlack)
    tight_pois.SetMarkerStyle(20)

    errors = getPoisson(loose, tight.Integral()/loose.Integral())
    errors.SetFillStyle(3244)
    errors.SetFillColor(ROOT.kAzure-6)
    errors.SetMarkerSize(0)
    errors.SetLineWidth(0)

    loose.Scale(tight.Integral()/loose.Integral())
    loose.SetLineColor(ROOT.kAzure+5)
    loose.SetFillColor(ROOT.kAzure+5)
    loose.SetFillStyle(1001)

    
    upper = ROOT.TPad("upper", "", 0.0025, 0.3, 0.9975, 0.9975)
    upper.Draw()
    upper.cd()
    upper.Range(-912.5,-5.428389,212.5,103.1394)
    upper.SetFillColor(0)
    upper.SetBorderMode(0)
    upper.SetBorderSize(2)
    upper.SetBottomMargin(0.05)
    upper.SetLeftMargin(.13)
    upper.SetRightMargin(0.07)
    upper.SetFrameBorderMode(0)

    loose.SetStats(0)
    loose.GetYaxis().SetLabelSize(0.044)
    loose.GetXaxis().SetLabelSize(0.044)
    loose.GetYaxis().SetTitleSize(0.044)

    tight.SetStats(0)
    loose.Draw("hist")
    tight_pois.Draw("p,same")
    errors.Draw("e2 same")
    cat = ROOT.TLatex()
    cat.SetTextSize(0.04)
    cat.DrawLatexNDC(0.17, 0.84, header)
    cat.Draw("same")

    leg = ROOT.TLegend(.16, .68, .46, .82)
    leg.SetBorderSize(0)
    leg.AddEntry(loose, "ID + antiID data", "lf")
    leg.AddEntry(tight_pois, "ID + ID data", "lp")
    leg.Draw("same")
    
    loose.GetYaxis().SetRangeUser(0, 1.4*max(loose.GetMaximum(), tight.GetMaximum()))
    loose.GetYaxis().SetTitle("Events")
    upper.cd()
    if drawprelim:
        CMS_lumi.CMS_lumi(upper, 4, 0, relPosX=0.077, lumi_13TeV = lumifb[options.era], extraText = "Preliminary")
    c.cd()

    lower = ROOT.TPad("lower", "", 0.0025, 0.0025, .9975, .3)
    lower.Draw()
    lower.cd()
    lower.SetGridy()
    lower.Range(-912.5,0.0384615,212.5,1.576923)
    lower.SetFillColor(0)
    lower.SetBorderMode(0)
    lower.SetBorderSize(2)
    lower.SetTopMargin(0.05)
    lower.SetLeftMargin(.13)
    lower.SetRightMargin(0.07)
    lower.SetBottomMargin(0.3)
    lower.SetFrameBorderMode(0)
    tight.GetXaxis().SetLabelSize(0.1)
    tight.GetYaxis().SetLabelSize(0.08)
    tight.GetXaxis().SetTitleSize(0.1)
    tight.GetYaxis().SetTitleOffset(0.5)
    tight.GetYaxis().SetTitleSize(0.1)
    tight.Divide(loose)
    tight.SetLineColor(ROOT.kBlack)
    tight.SetMarkerColor(ROOT.kBlack)
    tight.SetMarkerStyle(20)

    tight.Draw("p")
    tight.GetYaxis().SetRangeUser(0, 2)
    tight.GetYaxis().SetNdivisions(5)
    tight.GetYaxis().SetTitle("Obs/Exp")
    tight.GetXaxis().SetTitle("L_{xy} [cm]")
    c.Write(dir_to_save+name, ROOT.TObject.kOverwrite)
    c.SaveAs(dir_to_save+name+".png")
    c.SaveAs(dir_to_save+name+".pdf")
    del tight
    del loose
    del c

def compareSignal(background, signal0, signal100, m, name, header,v, blind= False):
    
    background.SetLineColor(ROOT.kAzure+5)
    background.SetLineWidth(3)
    background.SetFillColor(ROOT.kAzure+5)
    background.SetFillStyle(1001)
    signal0.SetLineColor(ROOT.kRed)
    signal0.SetLineWidth(3)
    signal0.SetFillStyle(0)
    signal100.SetLineColor(ROOT.kRed)
    signal100.SetLineWidth(3)
    signal100.SetFillStyle(0)
    signal100.Add(background)
    signal100.SetLineStyle(2)
    #Stack both signals on top of background
    stk = ROOT.THStack()
    stk.Add(background)
    stk.Add(signal0)
    #stk.SetMaximum(1.2*max(stk.GetMaximum(), signal100.GetMaximum()))
    stk.SetMaximum(75*max(stk.GetMaximum(), signal100.GetMaximum()))
    #will be on same canvas
    c = ROOT.TCanvas("c", '', 800, 600)
    c.SetRightMargin(0.05)
    c.SetLogy()
    stk.Draw("hist")
    signal100.Draw("hist,same")

    if drawprelim:
        CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = lumifb[options.era], extraText = "Preliminary")
    errors = ROOT.TGraphErrors(background.GetNbinsX())
    
    #Stack everything together
    for i in range(background.GetNbinsX()):
        errors.SetPoint(i, background.GetBinCenter(i+1), background.GetBinContent(i+1))
        errors.SetPointError(i, background.GetBinWidth(i+1)/2, background.GetBinContent(i+1) * 0.4)
    errors.SetFillStyle(3244)
    errors.SetFillColor(ROOT.kAzure-6)
    errors.SetMarkerSize(0)
    errors.SetLineWidth(0)
    errors.Draw("e2, same")

    cat = ROOT.TLatex()
    cat.SetTextSize(0.035)
    cat.DrawLatexNDC(0.16, 0.8, header)
    cat.Draw("same")
    leg = ROOT.TLegend(0.15, 0.56, 0.56, 0.78)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(background, "Predicted background", "lf")
    leg.AddEntry(signal0, "Signal, c#tau = 0 mm", "lf")
    leg.AddEntry(signal100, "Signal, c#tau = 100 mm", "lf")
    br = ROOT.TH1D()
    br.SetLineWidth(0)
    br.SetMarkerSize(0)
    br.SetFillStyle(0)
    #leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)=0.05", "l")
    leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)=0.00933", "l")
    leg.Draw("same")
    stk.GetYaxis().SetTitle("Events")
    stk.GetXaxis().SetLabelSize(0)
    lines = ROOT.TH1D("lines", "", binsM, 0, binsLxy*binsM)
    for i in range(1, binsM//2 + 1):
        lines.SetBinContent(2*i, 10000)
    lines.SetFillStyle(0)
    lines.SetLineColor(ROOT.kBlack)
    lines.SetLineWidth(2)
    lines.SetLineStyle(2)
    lines.Draw("hist,same")
    c.Update()
    a1 = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), 0, binsM*11, binsM*11, "+")
    for n in range(0, 11*binsM):
        if (n%11)%2 == 0:
            a1.ChangeLabel(n+1, -1, 0.025, -1, -1, -1, "{}".format((n*10)%110-10))
        else:
            a1.ChangeLabel(n+1, -1, 0.0, -1, -1, -1, "")
    a1.ChangeLabel(11*binsM+1, -1, 0.0, -1, -1, -1, "")
    a1.SetTitle("L_{xy} [cm]")
    a1.Draw()
    a2 = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), 0, binsM, binsM, "N-")
    for n in range(binsM):
        a2.ChangeLabel(n+1, -1, 0.022, -1, -1, -1, "{:.2f} #leq m(#gamma,#gamma) < {:.2f} GeV".format(4+(m+1.)/binsM*n, 4+(m+1.)/binsM*(n+1)))
    #a2.SetLabelOffset(-0.065)
    a2.SetLabelOffset(0.23)
    a2.SetLineWidth(0)
    a2.SetLineColor(0)
    a2.CenterLabels()
    a2.Draw()


    c.cd()
    
    c.Write(dir_to_save+name, ROOT.TObject.kOverwrite)
    c.SaveAs(dir_to_save+name+".png")
    c.SaveAs(dir_to_save+name+".pdf")
    del background
    del signal0
    del signal100
    del c

# Compare background, obs data, and signal distributions
def compareAll(background, data, signal0, signal100, m,l, name, header,v,scale,sigScale=0.1,blind= False):

    signal0.Scale(sigScale)
    signal100.Scale(sigScale)
    
    data.Write(f"{v}H_{l}_m{m}_data", ROOT.TObject.kOverwrite)
    
    #this is a Tgraph of the errors
    bkg_pois = getPoisson(background,scale)

    #Now the background can be rescaled
    background.Scale(scale)
    background.SetLineColor(ROOT.kAzure+5)
    #background.SetLineColor(ROOT.kBlack)
    background.SetLineWidth(3)
    background.SetFillColor(ROOT.kAzure+5)
    background.SetFillStyle(1001)
    background.Write(f"{v}H_{l}_m{m}_bkg", ROOT.TObject.kOverwrite)

    

    Nbins = signal0.GetNbinsX() 
    for bin in range(1,Nbins+1):
        print(f"bin {bin} sig0    ---", signal0.GetBinContent(bin))
        print(f"bin {bin} sig100  ---", signal100.GetBinContent(bin))
        print(f"bin {bin} bkg     ---", background.GetBinContent(bin))
        print(f"bin {bin} data    ---", data.GetBinContent(bin))
        print("")

    signal0.SetLineColor(ROOT.kRed)
    signal0.SetLineWidth(3)
    signal0.SetFillStyle(0)
    signal0.Write(f"{v}H_{l}_m{m}_sig0", ROOT.TObject.kOverwrite)
    signal100.SetLineColor(ROOT.kRed)
    signal100.SetLineWidth(3)
    signal100.SetFillStyle(0)
    signal100.SetLineStyle(2)
    signal100.Write(f"{v}H_{l}_m{m}_sig100", ROOT.TObject.kOverwrite)
    signal100.Add(background)
    
    data_pois = getPoisson(data)
    data_pois.SetLineColor(ROOT.kBlack)
    data_pois.SetMarkerColor(ROOT.kBlack)
    data_pois.SetMarkerStyle(20)

    stk = ROOT.THStack()
    stk.Add(background)
    if draw_sig:
        stk.Add(signal0)
        stk.SetMaximum(1.75*max(stk.GetMaximum(), signal100.GetMaximum(), (data.GetBinContent(data.GetMaximumBin())+data.GetBinError(data.GetMaximumBin()))))
    else:
        stk.SetMaximum(1.5*max(stk.GetMaximum(), get_max_y_with_error(data_pois)))

    

    #stk.SetMaximum(80)
    c = ROOT.TCanvas("c", '', 800, 600)
    c.SetRightMargin(0.05)
    stk.Draw("hist")
    
    if draw_sig:
        signal100.Draw("hist,same")
    data_pois.Draw("p,same")
    if drawprelim:
        CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = lumifb[options.era], extraText = "Preliminary")
    errors = ROOT.TGraphErrors(background.GetNbinsX())
    #TODO: change these to poisson errors
    #for i in range(background.GetNbinsX()):
    #    errors.SetPoint(i, background.GetBinCenter(i+1), background.GetBinContent(i+1))
    #    errors.SetPointError(i, background.GetBinWidth(i+1)/2, background.GetBinContent(i+1) * 0.4)
    errors = bkg_pois
    errors.SetFillStyle(3244)
    errors.SetFillColor(ROOT.kAzure-6)
    errors.SetMarkerSize(0)
    errors.SetLineWidth(0)
    errors.Draw("e2, same")

    cat = ROOT.TLatex()
    cat.SetTextSize(0.035)
    cat.DrawLatexNDC(0.16, 0.8, header)
    cat.Draw("same")
    

    if show_chi_sq:
        chi_sq, dof = chi_squared_graphs(data_pois,bkg_pois)
        cat2 = ROOT.TLatex()
        cat2.SetTextSize(0.035)
        cat2.DrawLatexNDC(0.73, 0.8, r"#chi^{2} = "+f"{chi_sq:.3f}" + f" dof: {dof}")
        cat2.Draw("same")
    
    leg = ROOT.TLegend(0.15, 0.56, 0.56, 0.78)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    if not closure:
        if postfit:
            leg.AddEntry(background, "Post-Fit Predicted Background", "f")
        else:
            leg.AddEntry(background, "Pre-Fit Predicted Background", "f")
        leg.AddEntry(data_pois, "Observed", "lp")
    else:
        leg.AddEntry(background, "antiID + antiID data", "f")
        leg.AddEntry(data_pois, "antiID + ID data", "lp")
    if draw_sig:
        if postfit:
            leg.AddEntry(signal0, "Post-Fit Signal, c#tau = 0 mm", "lf")
            leg.AddEntry(signal100, "Post-Fit Signal, c#tau = 100 mm", "lf")
        else:
            leg.AddEntry(signal0, "Pre-Fit Signal, c#tau = 0 mm", "lf")
            leg.AddEntry(signal100, "Pre-Fit Signal, c#tau = 100 mm", "lf")
    br = ROOT.TH1D()
    br.SetLineWidth(0)
    br.SetMarkerSize(0)
    br.SetFillStyle(0)
    if draw_sig:
        if postfit:
            leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma) = "+f"r*{0.5*sigScale:.2f}", "l")
        else:
            leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma) =  0.5", "l")   
    leg.Draw("same")
    stk.GetYaxis().SetTitle("Events")
    stk.GetXaxis().SetLabelSize(0)
    lines = ROOT.TH1D("lines", "", binsM, 0, binsLxy*binsM)
    for i in range(1, binsM//2 + 1):
        lines.SetBinContent(2*i, 10000)
    lines.SetFillStyle(0)
    lines.SetLineColor(ROOT.kBlack)
    lines.SetLineWidth(2)
    lines.SetLineStyle(2)
    lines.Draw("hist,same")
    c.Update()
    a1 = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), 0, binsM*11, binsM*11, "+")
    for n in range(0, 11*binsM):
        if (n%11)%2 == 0:
            a1.ChangeLabel(n+1, -1, 0.025, -1, -1, -1, "{}".format((n*10)%110-10))
        else:
            a1.ChangeLabel(n+1, -1, 0.0, -1, -1, -1, "")
    a1.ChangeLabel(11*binsM+1, -1, 0.0, -1, -1, -1, "")
    a1.SetTitle("L_{xy} [cm]")
    a1.Draw()
    a2 = ROOT.TGaxis(ROOT.gPad.GetUxmin(), ROOT.gPad.GetUymax(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), 0, binsM, binsM, "-")
    for n in range(binsM):
        a2.ChangeLabel(n+1, -1, 0.022, -1, -1, -1, "{:.2f} #leq m(#gamma,#gamma) < {:.2f} GeV".format(4+(m+1.)/binsM*n, 4+(m+1.)/binsM*(n+1)))
    a2.SetLabelOffset(-0.065)
    a2.CenterLabels()
    a2.Draw()

    c.cd()
    
    c.Write(name, ROOT.TObject.kOverwrite)
    if closure:
        c.SaveAs(dir_to_save+name+"_closure.png")
        c.SaveAs(dir_to_save+name+"_closure.pdf")
    else:
        c.SaveAs(dir_to_save+name+".png")
        c.SaveAs(dir_to_save+name+".pdf")

    c.Close()


fOut = ROOT.TFile("validation_{}.root".format(options.era), "UPDATE")
binsLxy = 110

postfit = False
show_chi_sq = False
draw_sig = False
closure = True
if closure == True:
    loose_name = "closure_sideband"
    tight_name = "loose_sideband"   
    bkg = "closure_loose_scaled"
    data= "closure_tight"
else:    
    loose_name = "loose_sideband"
    tight_name = "tight_sideband"
    bkg = "background"
    data = "data_obs"


for v in ['W','Z']:
    for l in ['MU','ELE']:
        m = int(options.mass)
        if v == 'Z' and m == 15 and l == 'ELE':
            binsM = 2
        elif v == 'Z' and m ==20 and l == 'ELE':
            binsM = 2
        else:
            binsM =4
        if options.era == 'Run2':
            if not postfit:
                f16         = ROOT.TFile(f"datacardInputs_{v}H_{l}_m{m}_2016.root")
                f17         = ROOT.TFile(f"datacardInputs_{v}H_{l}_m{m}_2017.root")
                f18         = ROOT.TFile(f"datacardInputs_{v}H_{l}_m{m}_2018.root")

                loose_sb    = f16.Get(loose_name) + f17.Get(loose_name) + f18.Get(loose_name)
                tight_sb    = f16.Get(tight_name) + f17.Get(tight_name) + f18.Get(tight_name)


                loose       = f16.Get(bkg) + f17.Get(bkg) + f18.Get(bkg)
                tight       = f16.Get(data) + f17.Get(data) + f18.Get(data)
                
                sig0   =  f16.Get("ZH_HToPhiPhi_ctau0")    + f16.Get("ggZH_HToPhiPhi_ctau0") + \
                        f17.Get("ZH_HToPhiPhi_ctau0")    + f17.Get("ggZH_HToPhiPhi_ctau0") + \
                        f18.Get("ZH_HToPhiPhi_ctau0")    + f18.Get("ggZH_HToPhiPhi_ctau0")   
                
                sig100 =  f16.Get("ZH_HToPhiPhi_ctau100")  + f16.Get("ggZH_HToPhiPhi_ctau100") + \
                        f17.Get("ZH_HToPhiPhi_ctau100")  + f17.Get("ggZH_HToPhiPhi_ctau100") + \
                        f18.Get("ZH_HToPhiPhi_ctau100")  + f18.Get("ggZH_HToPhiPhi_ctau100") 
    
                try: 
                    sig0   += f16.Get("ttH_HToPhiPhi_ctau0")    + f17.Get("ttH_HToPhiPhi_ctau0")    + f18.Get("ttH_HToPhiPhi_ctau0")
                    sig100 += f16.Get("ttH_HToPhiPhi_ctau100")  + f17.Get("ttH_HToPhiPhi_ctau100")  + f18.Get("ttH_HToPhiPhi_ctau100")
                except NotImplementedError:
                    print("no ttH contribution")
                try:
                    sig0   += f16.Get("WH_HToPhiPhi_ctau0")    + f17.Get("WH_HToPhiPhi_ctau0")    + f18.Get("WH_HToPhiPhi_ctau0")  
                    sig100 += f16.Get("WH_HToPhiPhi_ctau100")  + f17.Get("WH_HToPhiPhi_ctau100")  + f18.Get("WH_HToPhiPhi_ctau100")  
                except NotImplementedError:
                    print('no WH contribution')
                sigScale = 0.1
            else:
                f         = ROOT.TFile(f"postFit_{v}H_{l}_m{m}_Run2.root")
                
                loose_sb    = f.Get(loose_name)
                tight_sb    = f.Get(tight_name)

                loose       = f.Get(bkg)
                tight       = f.Get(data)
                
                sig0 = f.Get("sig0")
                sig100 = f.Get("sig100")

                sigScale = 1
            

        else:

            if not postfit:
                f           = ROOT.TFile(f"datacardInputs_{v}H_{l}_m{m}_{options.era}.root")
            else:
                f           = ROOT.TFile(f"postFit_{v}H_{l}_m{m}_{options.era}.root")

            loose_sb    = f.Get(loose_name)
            tight_sb    = f.Get(tight_name)

            loose       = f.Get(bkg)
            tight       = f.Get(data)

            
            if not postfit:
                sig0   = f.Get("ZH_HToPhiPhi_ctau0")    + f.Get("ggZH_HToPhiPhi_ctau0")
                
                sig100 = f.Get("ZH_HToPhiPhi_ctau100")  + f.Get("ggZH_HToPhiPhi_ctau100")
    
                try: 
                    sig0   += f.Get("ttH_HToPhiPhi_ctau0")
                    sig100 += f.Get("ttH_HToPhiPhi_ctau100")
                except NotImplementedError:
                    print("no ttH contribution")
                try:
                    sig0   += f.Get("WH_HToPhiPhi_ctau0")
                    sig100 += f.Get("WH_HToPhiPhi_ctau100")
                except NotImplementedError:
                    print('no WH contribution')
                sigScale = 0.1
            else:
                sig0 = f.Get("sig0")
                sig100 = f.Get("sig100")
                loose_sb    = f.Get(loose_name)
                tight_sb    = f.Get(tight_name)
                sigScale = 1
        fOut.cd()
        
        if v == 'Z':
            loose_sb.RebinX(5)
            tight_sb.RebinX(5)

        scale = tight_sb.Integral()/loose_sb.Integral()
        
        loose_sb.Write("data_loose", ROOT.TObject.kOverwrite)
        tight_sb.Write("data_tight", ROOT.TObject.kOverwrite)

        #Unscale the histogram to the the proper Poisson errors
        loose.Scale(1/scale)       
        
        
        compareSideband(loose_sb, tight_sb, m, "closure_{}H_{}_m{}_sideband_{}".format(v,l,m,options.era), "m(#gamma#gamma) > 65 GeV, "+headers[v][l]+" m_{#Phi} = "+"{} GeV".format(m))
        
        
        #sigScale = 0.37327*0.1
        
    
        #compareSignal(loose, sig0, sig100, m, "closure_{}H_{}_m{}_signal_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m),v)
        if not options.blind:
            #NOTE: loose is unscaled
            compareAll(loose, tight, sig0, sig100, m,l, "{}H_{}_m{}_data_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m),v,scale,sigScale=sigScale)

fOut.Close()
