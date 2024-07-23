import ROOT, sys, os, array
sys.path.append("/home/tyler/DDP/rdf_simple/")
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.Declare('#include "common/datacardHelpers.h"')
import optparse
parser = optparse.OptionParser()
parser.add_option("-e", "--era", dest = "era", default = "2018", help = "2016, 2017, 2018")
parser.add_option("-d", "--date", dest = "date", default = "03_18_24", help = "current date (to name directory)")
parser.add_option("-p", "--prod", dest = "prod", default = "03_26_24", help = "sample production date")
parser.add_option("-b", "--blind", dest = "blind", action="store_true", default = False, help = "Blind data, sets data_obs = expected")
parser.add_option("-v", "--v", dest = "v", default = 'Z', help = "V category (Z or W)")
parser.add_option("-l", "--lep", dest = "l", default = 'MU', help = "lepton category (ELE or MU)")
parser.add_option("-m", "--model", dest = "model", default = "SM", help = "SM or model_indep")
(options,args) = parser.parse_args()

from python.VHTools.VHcuts import * 
from startVH import *

model_sfx = ""
modelIndependent = False
if options.model != "SM":
    model_sfx = "_modelIndependent"
    modelIndependent = True

era = options.era
date = options.date
prod = options.prod
v = options.v
l = options.l

binsM = 4
binsLxy = 110
lxyMin = -10
lxyMax = 100

if not os.path.isdir("datacards_{}".format(date)):
    os.mkdir("datacards_{}".format(date))
os.chdir("datacards_{}".format(date))

def fillZeroBins(hist, thresh):
    for i in range(1, hist.GetNbinsX()+1):
        if hist.GetBinContent(i) <= 0:
            hist.SetBinContent(i,thresh)
    return hist

xsecUnc = {'W': [-0.005, 0.007],
           'Z': [-0.031, 0.038]}

pdfUnc = {'W': 1.018,
          'Z': 1.016}

# Which lifetime to use for ctau reweighing
lifetimes = {0: 0, 3: 10, 10: 10, 14: 10, 20: 20, 32: 20, 50: 50, 70: 50, 100: 100, 316: 1000, 1000: 1000}

# Binning quantiles
quantiles = {'W': array.array('d', [.2, .4, .6, .8]),
             'Z': array.array('d', [.25, .5, .75])}

# Getting the required plotters
plotters = getPlotters(era, prod, "/home/tyler/samples/DDP/rdf_samples", modelIndependent = modelIndependent)
dataEMU = plotters['dataEMU']
signal = plotters['signal']

for m in masses:
    print("Making card for {}->{}, m={}".format(v,l,m))         
    f = ROOT.TFile("datacardInputs_{}H_{}_m{}_{}{}.root".format(v,l,m,era,model_sfx), "UPDATE")

    # negative lxy sidebands for loose (anti-ID) and tight (ID)
    loose_sb = dataEMU[ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))
    tight_sb = dataEMU[ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['sr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))
    
    # normalization scale factor
    scale = tight_sb.Integral() / loose_sb.Integral() if loose_sb.Integral() > 0 else 1.0

    loose_sb.Write("loose_sideband", ROOT.TObject.kOverwrite)
    tight_sb.Write("tight_sideband", ROOT.TObject.kOverwrite)

    # Get predicted background from antiID region
    background = dataEMU[ana[v][l]].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m], "1", ("data_loose_{}_{}_m{}".format(v,l,m), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
    background.Scale(scale)
    background.Write("data_loose", ROOT.TObject.kOverwrite)
    background = unfoldTH2(background)
    binsX = rebinTH1(background, quantiles[v], binsM, binsLxy, lxyMin, lxyMax)
    background.Rebin(len(binsX) - 1, "data_loose_rebinned", binsX)
    background = ROOT.gDirectory.Get("data_loose_rebinned")

    background = fillZeroBins(background, 0.001)
    background.Write("background", ROOT.TObject.kOverwrite)

    # Get observed data
    data_obs = dataEMU[ana[v][l]].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), cuts['sr'][v][l][m], "1", ("data_tight_{}_{}_m{}".format(v,l,m), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
    data_obs = unfoldTH2(data_obs)
    data_obs.Rebin(len(binsX) - 1, "data_tight_rebinned", binsX)
    data_obs = ROOT.gDirectory.Get("data_tight_rebinned")
    # Only write if not blind, else use expected background
    if not options.blind:
        data_obs.Write("data_obs", ROOT.TObject.kOverwrite)
    else:
        data_obs = ROOT.TH1D("data_obs", "", len(binsX) - 1, binsX)
        for i in range(1, data_obs.GetNbinsX()+1):
            data_obs.SetBinContent(i, int(round(background.GetBinContent(background.FindBin(background.GetBinCenter(i))))))
        data_obs.Write("data_obs", ROOT.TObject.kOverwrite)

    # Scale factors
    sfs = getSF(scaleFactors[v][l] + scaleFactors['g'][m])

    # Defining photon scale corrections
    for ct in ctaus:

        # energy scale/smearing corrections
        # First get the scaled/smeared energies
        signal[ana[v][l]][m][ct]['2G2Q'].define("Photon_ptSmearUp", "ptSmearUp(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaUp)")
        signal[ana[v][l]][m][ct]['2G2Q'].define("Photon_ptSmearDown", "ptSmearDown(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaDown)")
        signal[ana[v][l]][m][ct]['2G2Q'].define("Photon_energyScaleUp", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
        signal[ana[v][l]][m][ct]['2G2Q'].define("Photon_energyScaleDown", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
        signal[ana[v][l]][m][ct]['2G2Q'].define("Photon_ptScaleUp", "Photon_pt*Photon_energyScaleUp")
        signal[ana[v][l]][m][ct]['2G2Q'].define("Photon_ptScaleDown", "Photon_pt*Photon_energyScaleDown")

        # Then calculate the vertex/raw mass with scaled/smeared values
        for corr in ['Scale', 'Smear']:
            for d in ['Up', 'Down']:
                signal[ana[v][l]][m][ct]['2G2Q'].define("best_2g_{}{}_m{}_info".format(corr,d,m), "getDxy(Photon_pt{}{}, Photon_eta, Photon_phi, Photon_isScEtaEB, Photon_isScEtaEE, best_2g_idx1_m{}, best_2g_idx2_m{}, {})".format(corr, d, m, m, m))
                signal[ana[v][l]][m][ct]['2G2Q'].define("best_2g_dxy_{}{}_m{}".format(corr,d,m), "best_2g_{}{}_m{}_info[0]".format(corr,d,m))
                signal[ana[v][l]][m][ct]['2G2Q'].define("best_2g_raw_mass_{}{}_m{}".format(corr,d,m), "best_2g_{}{}_m{}_info[1]".format(corr,d,m))


    # Now get expected signal for each lifetime
    for ctau in lifetimes:
        # Category name
        cat = "{}H_{}_m{}_ctau{}_{}{}".format(v,l,m,ctau,era,model_sfx)
        
        # Lifetime reweighting
        ct = ctau
        reweight = "(1)"
        if ct not in [0, 10, 20, 50, 100, 1000]:
            ct = min([x for x in [0,10,20,50,100,1000] if x >= ctau], key = lambda x: abs(x-ctau))
            reweight = "(reweight_{}to{})".format(ct, ctau)
            signal[ana[v][l]][m][ct]['2G2Q'].define("reweight_{}to{}".format(ct,ctau), "reweightCT(GenPart_ctau[abs(GenPart_pdgId)==9000006], {}/10., {}/10.)".format(ct, ctau))

        # Nominal signal
        sig = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", lumi[era], ("sig_{}_{}_m{}_ct{}".format(v,l,m,ct), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
        sig = unfoldTH2(sig)
        sig.Rebin(len(binsX) - 1, "HToPhiPhi_ctau{}".format(ctau), binsX)
        sig = ROOT.gDirectory.Get("HToPhiPhi_ctau{}".format(ctau))
        sig.Write("HToPhiPhi_ctau{}".format(ctau), ROOT.TObject.kOverwrite)

        # Energy scale/smearing corrections, scale factor corrections
        sigScaleUp = {}
        sigScaleDown = {}

        sigScaleUp[v] = {}
        sigScaleDown[v] = {}
        if l == 'ELE':
            for sf in ['reco', 'id', 'trig']:
                sfUp, sfDown = getSF(scaleFactors[v][l] + scaleFactors['g'][m], 'Electron_{}'.format(sf))
                sigUp = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfUp+")*("+reweight+")", lumi[era], ("sig_ele{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ct), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigUp = unfoldTH2(sigUp)
                sigUp.Rebin(len(binsX) - 1, "sig_ele{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), binsX)
                sigUp = ROOT.gDirectory.Get("sig_ele{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau))
                sigUp.Write("HToPhiPhi_ctau{}_ele_{}Up".format(ctau, sf), ROOT.TObject.kOverwrite)

                sigScaleUp[v][sf] = fillZeroBins(sigUp, 0)

                sigDown = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfDown+")*("+reweight+")", lumi[era], ("sig_ele{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigDown = unfoldTH2(sigDown)
                sigDown.Rebin(len(binsX) - 1, "sig_ele{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), binsX)
                sigDown = ROOT.gDirectory.Get("sig_ele{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau))
                sigDown.Write("HToPhiPhi_ctau{}_ele_{}Down".format(ctau, sf), ROOT.TObject.kOverwrite)

                sigScaleDown[v][sf] = fillZeroBins(sigDown, 0)

        if l == 'MU':
            for sf in ['reco', 'id', 'trig', 'iso']:
                sfUp, sfDown = getSF(scaleFactors[v][l] + scaleFactors['g'][m], 'Muon_{}'.format(sf))
                sigUp = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfUp+")*("+reweight+")", lumi[era], ("sig_mu{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigUp = unfoldTH2(sigUp)
                sigUp.Rebin(len(binsX) - 1, "sig_mu{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), binsX)
                sigUp = ROOT.gDirectory.Get("sig_mu{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau))
                sigUp.Write("HToPhiPhi_ctau{}_mu_{}Up".format(ctau, sf), ROOT.TObject.kOverwrite)

                sigScaleUp[v][sf] = fillZeroBins(sigUp, 0)

                sigDown = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfDown+")*("+reweight+")", lumi[era], ("sig_mu{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigDown = unfoldTH2(sigDown)
                sigDown.Rebin(len(binsX) - 1, "sig_mu{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), binsX)
                sigDown = ROOT.gDirectory.Get("sig_mu{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau))
                sigDown.Write("HToPhiPhi_ctau{}_mu_{}Down".format(ctau, sf), ROOT.TObject.kOverwrite)

                sigScaleDown[v][sf] = fillZeroBins(sigDown, 0)

        sigScaleUp['g'] = {}
        sigScaleDown['g'] = {}
        for sf in ['id', 'pix']:
            sfUp, sfDown = getSF(scaleFactors[v][l] + scaleFactors['g'][m], "Photon_{}".format(sf))
            sigUp = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfUp+")*("+reweight+")", lumi[era], ("sig_pho{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
            sigUp = unfoldTH2(sigUp)
            sigUp.Rebin(len(binsX) - 1, "sig_pho{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), binsX)
            sigUp = ROOT.gDirectory.Get("sig_pho{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau))
            sigUp.Write("HToPhiPhi_ctau{}_pho_{}Up".format(ctau, sf), ROOT.TObject.kOverwrite)

            sigScaleUp['g'][sf] = fillZeroBins(sigUp, 0)

            sigDown = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfDown+")*("+reweight+")", lumi[era], ("sig_pho{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
            sigDown = unfoldTH2(sigDown)
            sigDown.Rebin(len(binsX) - 1, "sig_pho{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), binsX)
            sigDown = ROOT.gDirectory.Get("sig_pho{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau))
            sigDown.Write("HToPhiPhi_ctau{}_pho_{}Down".format(ctau, sf), ROOT.TObject.kOverwrite)

            sigScaleDown['g'][sf] = fillZeroBins(sigDown, 0)

        for corr in ['Scale', 'Smear']:

            sigUp = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_{}Up_m{}".format(corr,m), "best_2g_dxy_{}Up_m{}".format(corr,m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", lumi[era], ("sig_pho{}Up_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
            sigUp = unfoldTH2(sigUp)
            sigUp.Rebin(len(binsX) - 1, "sig_pho{}Up_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau), binsX)
            sigUp = ROOT.gDirectory.Get("sig_pho{}Up_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau))
            sigUp.Write("HToPhiPhi_ctau{}_pho_{}Up".format(ctau, corr), ROOT.TObject.kOverwrite)

            sigScaleUp['g'][corr] = fillZeroBins(sigUp, 0)

            sigDown = signal[ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_{}Down_m{}".format(corr,m), "best_2g_dxy_{}Down_m{}".format(corr,m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", lumi[era], ("sig_pho{}Down_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
            sigDown = unfoldTH2(sigDown)
            sigDown.Rebin(len(binsX) - 1, "sig_pho{}Down_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau), binsX)
            sigDown = ROOT.gDirectory.Get("sig_pho{}Down_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau))
            sigDown.Write("HToPhiPhi_ctau{}_pho_{}Down".format(ctau, corr), ROOT.TObject.kOverwrite)

            sigScaleDown['g'][corr] = fillZeroBins(sigDown, 0)

        # Now make the cards
        for i in range(data_obs.GetNbinsX()):
            if sig.GetBinContent(i+1) <= 0:
                continue
            card = open("datacard_{}_bin{}.txt".format(cat, i), "w")
            card.write("imax 1\n")
            card.write("jmax 1\n")
            card.write("kmax *\n")
            card.write('-------------------------\n')
            card.write("bin {}_bin{}\n".format(cat, i))
            card.write("observation {}\n".format(data_obs.GetBinContent(i+1)))
            card.write('-------------------------\n')
            card.write("bin\t{}_bin{}\t{}_bin{}\n".format(cat, i, cat, i))
            card.write("process\tHToPhiPhi_ctau{}\tbackground\n".format(ctau))
            card.write("process\t0\t1\n")
            card.write("rate\t{}\t{}\n".format(0.1*sig.GetBinContent(i+1), background.GetBinContent(i+1)))
            card.write('-------------------------\n')
            card.write("CMS_lumi_{}\tlnN\t{}\t-\n".format(era, lumiUnc[era]))
            card.write("CMS_pileup_{}\tlnN\t1.05\t-\n".format(era))
            if options.model == "SM":
                card.write("CMS_xsec_{}\tlnN\t{}/{}\t-\n".format(v, (1.0+xsecUnc[v][0])/1.0, (1.0+xsecUnc[v][1])/1.0))
                card.write("CMS_pdf_{}\tlnN\t{}\t-\n".format(v, pdfUnc[v]))
            card.write("CMS_bkg_{}_norm\tgmN\t{}\t-\t{}\n".format(cat, int(tight_sb.Integral()), background.GetBinContent(i+1)/tight_sb.Integral()))
            card.write("CMS_bkg_{}_bin{}_shape\tlnN\t-\t1.4\n".format(cat, i))

            for sf in sigScaleUp[v]:
                errUp = sigScaleUp[v][sf].GetBinContent(i+1) / sig.GetBinContent(i+1)
                errDown = sigScaleDown[v][sf].GetBinContent(i+1) / sig.GetBinContent(i+1)
                if errUp == 0:
                    errUp = 1.0
                if errDown == 0:
                    errDown = 1.0
                card.write("CMS_{}_{}\tlnN\t{}/{}\t-\n".format(l,sf,errDown, errUp))
            for sf in sigScaleUp['g']:
                errUp = sigScaleUp['g'][sf].GetBinContent(i+1) / sig.GetBinContent(i+1)
                errDown = sigScaleDown['g'][sf].GetBinContent(i+1) / sig.GetBinContent(i+1)
                if errUp == 0:
                    errUp = 1.0
                if errDown == 0:
                    errDown = 1.0
                card.write("CMS_gamma_{}\tlnN\t{}/{}\t-\n".format(sf,errDown,errUp))

            card.close()
