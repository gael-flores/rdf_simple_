# Draw closure plots for a given mass hypothesis
import sys
sys.path.append("/uscms/home/gfavila/nobackup/rdf_simple_/")
import common.tdrstyle as tdrstyle
import common.CMS_lumi as CMS_lumi
sty = tdrstyle.setTDRStyle()
import ROOT, math, os
#from python.plotEditor import plotter as plotEditor
#pe = plotEditor()
ROOT.gROOT.SetBatch(True)

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
dir_to_save = f"/uscms/home/gfavila/nobackup/rdf_simple_/datacards_{date}/plots/"

intLumi = {'2018': 59830,
           '2017': 41480,
           '2016': 36310,
           'Run2': 137620
}   

lumi = {'2018': '59830',
        '2017': '41480',
        '2016': '36310',
        'Run2': '137620'
}   

lumifb = {'2018': "59.83",
          '2017': '41.48',
          '2016': '36.31'}


def fetchError(q, n):
    l = 0
    if n!=0:
        l = ROOT.Math.chisquared_quantile_c(1.0-q, 2.0*n)/2.0
    u = ROOT.Math.chisquared_quantile_c(q, 2.0*n+2)/2.0
    return [l,u]

# Poisson points for histogram
def getPoisson(h):
    q = (1-0.6827)/2.0
    gRate = ROOT.TGraphAsymmErrors()
    n = 0
    for i in range(1, h.GetNbinsX()+1):
        thresh = h.GetBinCenter(i)
        N = h.GetBinContent(i)
        #if N == 0:
        #    continue
        gRate.SetPoint(n, thresh, N)
        error = fetchError(q, N)
        gRate.SetPointError(n, 0, 0, (N-error[0]), (error[1]-N))
        n+=1
    return gRate

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

    errors = getPoisson2(loose, tight.Integral()/loose.Integral())
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


# Compare background and signal distributions for 0 and 100 mm ctau
def compareSignalW(background, WHsignal0, WHsignal100, ttHsignal0,ttHsignal100, m, name, header,v, blind= False):
    
    background.SetLineColor(ROOT.kAzure+5)
    background.SetLineWidth(3)
    background.SetFillColor(ROOT.kAzure+5)
    background.SetFillStyle(1001)

    WHsignal0.SetLineColor(ROOT.kRed)
    WHsignal0.SetLineWidth(3)
    WHsignal0.SetFillStyle(0)
    WHsignal0.Add(background)

    WHsignal100.SetLineColor(ROOT.kRed)
    WHsignal100.SetLineWidth(3)
    WHsignal100.SetFillStyle(0)
    WHsignal100.Add(background)
    WHsignal100.SetLineStyle(2)
    
    ttHsignal0.SetLineColor(ROOT.kBlue)
    ttHsignal0.SetLineWidth(3)
    ttHsignal0.SetFillStyle(0)

    ttHsignal100.SetLineColor(ROOT.kBlue)
    ttHsignal100.SetLineWidth(3)
    ttHsignal100.SetFillStyle(0)
    ttHsignal100.Add(background)
    ttHsignal100.SetLineStyle(2)

    #Stack both signals on top of background
    stk = ROOT.THStack()
    stk.Add(background)
    stk.Add(ttHsignal0)
    stk.SetMaximum(1.2*max(stk.GetMaximum(), WHsignal0.GetMaximum()))
    #will be on same canvas
    c = ROOT.TCanvas("c", '', 800, 600)
    c.SetRightMargin(0.05)
    stk.Draw("hist")
    WHsignal0.Draw("hist,same")
    ttHsignal100.Draw("hist,same")
    WHsignal100.Draw("hist,same")

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
    if v == 'W':
        #separate WH and ttH - don't do Z
        leg.AddEntry(WHsignal0, "WH Signal, c#tau = 0 mm", "lf")
        leg.AddEntry(WHsignal100, "WH Signal, c#tau = 100 mm", "lf")
        leg.AddEntry(ttHsignal0, "ttH Signal, c#tau = 0 mm", "lf")
        leg.AddEntry(ttHsignal100, "ttH Signal, c#tau = 100 mm", "lf")
    elif v == 'Z':
        leg.AddEntry(signal0, "Signal, c#tau = 0 mm", "lf")
        leg.AddEntry(signal100, "Signal, c#tau = 100 mm", "lf")
    br = ROOT.TH1D()
    br.SetLineWidth(0)
    br.SetMarkerSize(0)
    br.SetFillStyle(0)
    leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)=0.05", "l")
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
    c.SaveAs(name+".png")
    c.SaveAs(name+".pdf")

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
        CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = lumifb[options.era], extraText = "Work in progress")
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
    leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)=0.05", "l")
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
def compareAll(background, data, signal0, signal100, m, name, header,v, blind= False):
    
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

    data_pois = getPoisson(data)
    data_pois.SetLineColor(ROOT.kBlack)
    data_pois.SetMarkerColor(ROOT.kBlack)
    data_pois.SetMarkerStyle(20)

    stk = ROOT.THStack()
    stk.Add(background)
    stk.Add(signal0)
    stk.SetMaximum(1.75*max(stk.GetMaximum(), signal100.GetMaximum(), (data.GetBinContent(data.GetMaximumBin())+data.GetBinError(data.GetMaximumBin()))))
    c = ROOT.TCanvas("c", '', 800, 600)
    c.SetRightMargin(0.05)
    stk.Draw("hist")
    signal100.Draw("hist,same")
    data_pois.Draw("p,same")
    if drawprelim:
        CMS_lumi.CMS_lumi(c, 4, 0, relPosX=0.077, lumi_13TeV = lumifb[options.era], extraText = "Preliminary")
    errors = ROOT.TGraphErrors(background.GetNbinsX())
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
    leg.AddEntry(data_pois, "Observed", "lp")
    if v == 'W':
        leg.AddEntry(signal0, "Signal, c#tau = 0 mm", "lf")
        leg.AddEntry(signal100, "Signal, c#tau = 100 mm", "lf")
    elif v == 'Z':
        leg.AddEntry(signal0, "Signal, c#tau = 0 mm", "lf")
        leg.AddEntry(signal100, "Signal, c#tau = 100 mm", "lf")
    br = ROOT.TH1D()
    br.SetLineWidth(0)
    br.SetMarkerSize(0)
    br.SetFillStyle(0)
    leg.AddEntry(br, "BR(H#rightarrow#Phi#Phi)#timesBR(#Phi#rightarrow#gamma#gamma)=0.05", "l")
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
    c.SaveAs(name+".png")
    c.SaveAs(name+".pdf")


headers = {}
headers['W'] = {}
headers['W']['ELE'] = "W#rightarrowe#nu"
headers['W']['MU'] = "W#rightarrow#mu#nu"
headers['Z'] = {}
headers['Z']['ELE'] = "Z#rightarrowee"
headers['Z']['MU'] = "Z#rightarrow#mu#mu"

fOut = ROOT.TFile("validation_{}.root".format(options.era), "UPDATE")
binsLxy = 110

#for v in ['W', 'Z']:
for v in ['W','Z']:
    for l in ['ELE','MU']:
        m = int(options.mass)
        if v == 'Z' and m == 15 and l == 'ELE':
            binsM = 2
        elif v == 'Z' and m ==20 and l == 'ELE':
            binsM = 2
        else:
            binsM =4
        f = ROOT.TFile("datacardInputs_{}H_{}_m{}_{}.root".format(v,l,m,options.era))
        fOut.cd()
        loose_sb = f.Get("loose_sideband")
        tight_sb = f.Get("tight_sideband")
        if v == 'Z':
            loose_sb.RebinX(5)
            tight_sb.RebinX(5)
        loose = f.Get("background")
        tight = f.Get("data_obs")
        compareSideband(loose_sb, tight_sb, m, "closure_{}H_{}_m{}_sideband_{}".format(v,l,m,options.era), "m(#gamma#gamma) > m_{H}/2, "+headers[v][l]+" m_{#Phi} = "+"{} GeV".format(m))
        
        if v =='Z':
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
        else:
            sig0   = f.Get("ZH_HToPhiPhi_ctau0")    + f.Get("ggZH_HToPhiPhi_ctau0")   + f.Get("ttH_HToPhiPhi_ctau0")   + f.Get("WH_HToPhiPhi_ctau0")
            sig100 = f.Get("ZH_HToPhiPhi_ctau100")  + f.Get("ggZH_HToPhiPhi_ctau100") + f.Get("ttH_HToPhiPhi_ctau100")  + f.Get("WH_HToPhiPhi_ctau100")
        sig0.Scale(0.1)
        sig100.Scale(0.1)

        compareSignal(loose, sig0, sig100, m, "closure_{}H_{}_m{}_signal_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m),v)
        if not options.blind:
            compareAll(loose, tight, sig0, sig100, m, "{}H_{}_m{}_data_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m))
        """

        if v == 'Z':
            sig0 = f.Get("ZH_HToPhiPhi_ctau0")  + f.Get("ggZH_HToPhiPhi_ctau0")
            sig0.Scale(0.1)
            sig100 = f.Get("ZH_HToPhiPhi_ctau100")  + f.Get("ggZH_HToPhiPhi_ctau100")
            sig100.Scale(0.1)
            compareSignalZ(loose, sig0, sig100, m, "closure_{}H_{}_m{}_signal_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m),v)
            if not options.blind:
                compareAll(loose, tight, sig0, sig100, m, "{}H_{}_m{}_data_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m))
        
        elif v == 'W':
            WHsig0  = f.Get("WH_HToPhiPhi_ctau0")
            ttHsig0 = f.Get("ttH_HToPhiPhi_ctau0")

            WHsig0.Scale(0.1)
            ttHsig0.Scale(0.1)
            
            WHsig100  = f.Get("WH_HToPhiPhi_ctau100")
            ttHsig100 = f.Get("ttH_HToPhiPhi_ctau100")
            WHsig100.Scale(0.1)
            ttHsig100.Scale(0.1)
            compareSignal(loose, WHsig0, WHsig100,ttHsig0,ttHsig100, m, "closure_{}H_{}_m{}_signal_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m),v)

            if not options.blind:
                compareAll(loose, tight, sig0, sig100, m, "{}H_{}_m{}_data_{}".format(v,l,m,options.era), "Signal Region, "+headers[v][l]+" m_{#Phi}="+"{} GeV".format(m))
        """
fOut.Close()
