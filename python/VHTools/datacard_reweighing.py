import ROOT, sys, os, array
sys.path.append("/uscms/home/gfavila/nobackup/rdf_simple_/")
ROOT.gROOT.SetBatch(True)

ROOT.gInterpreter.Declare('#include "common/datacardHelpers.h"')
import optparse
parser = optparse.OptionParser()
parser.add_option("-e", "--era", dest = "era", default = "2018", help = "2016, 2017, 2018")
parser.add_option("-a", "--mass", dest = "mass", default = 15,type=int ,help = "15,20,30,40,50,55")
parser.add_option("-d", "--date", dest = "date", default = "03_18_24", help = "current date (to name directory)")
parser.add_option("-p", "--prod", dest = "prod", default = "03_26_24", help = "sample production date")
parser.add_option("-b", "--blind", dest = "blind", action="store_true", default = False, help = "Blind data, sets data_obs = expected")
parser.add_option("-v", "--v", dest = "v", default = 'Z', help = "V category (Z or W)")
parser.add_option("-l", "--lep", dest = "l", default = 'MU', help = "lepton category (ELE or MU)")
parser.add_option("-m", "--model", dest = "model", default = "SM", help = "SM or model_indep")
(options,args) = parser.parse_args()

from python.VHTools.VHcuts import * 
import startVH

model_sfx = ""
modelIndependent = False
if options.model != "SM":
    model_sfx = "_modelIndependent"
    modelIndependent = True

era = options.era
mass = options.mass
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

masses = [mass]#[20,30,40,55]

channels = ['ttH','WH','ZH','ggZH']

# Which lifetime to use for ctau reweighing
lifetimes = {0: 0, 3: 10, 10: 10, 14: 10, 20: 20, 32: 20, 50: 50, 70: 50, 100: 100, 316: 1000, 1000: 1000}

# Binning quantiles
quantiles = {'W': array.array('d', [.2, .4, .6, .8]),
             'Z': array.array('d', [.25, .5, .75])}




for m in masses:
    print("Making card for {}->{}, m={}".format(v,l,m))
    # Getting the required plotters
    plotters = startVH.getPlotters(era, prod, "/uscms/home/gfavila/nobackup/rdf_simple_/DDP/",[m], modelIndependent = modelIndependent)
    dataEMU = plotters['dataEMU']
    signal = plotters['signal']
    
    f = ROOT.TFile("datacardInputs_{}H_{}_m{}_{}{}.root".format(v,l,m,era,model_sfx), "UPDATE")

    # negative lxy sidebands for loose (anti-ID) and tight (ID)
    loose_sb = dataEMU[startVH.ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))
    tight_sb = dataEMU[startVH.ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['sr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))
    
    # normalization scale factor
    scale = tight_sb.Integral() / loose_sb.Integral() if loose_sb.Integral() > 0 else 1.0

    loose_sb.Write("loose_sideband", ROOT.TObject.kOverwrite)
    tight_sb.Write("tight_sideband", ROOT.TObject.kOverwrite)

    # Get predicted background from antiID region
    background = dataEMU[startVH.ana[v][l]].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m], "1", ("data_loose_{}_{}_m{}".format(v,l,m), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
    background.Scale(scale)
    background.Write("data_loose", ROOT.TObject.kOverwrite)
    background = startVH.unfoldTH2(background)
    binsX = startVH.rebinTH1(background, quantiles[v], binsM, binsLxy, lxyMin, lxyMax)
    background.Rebin(len(binsX) - 1, "data_loose_rebinned", binsX)
    background = ROOT.gDirectory.Get("data_loose_rebinned")

    background = fillZeroBins(background, 0.001)
    background.Write("background", ROOT.TObject.kOverwrite)

    # Get observed data
    data_obs = dataEMU[startVH.ana[v][l]].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), cuts['sr'][v][l][m], "1", ("data_tight_{}_{}_m{}".format(v,l,m), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
    data_obs = startVH.unfoldTH2(data_obs)
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
    sfs = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m])
    #TODO import signals from config file, which doesn't exist yet
    # Defining photon scale corrections
    for ct in startVH.ctaus:
        for channel in channels:
            # energy scale/smearing corrections
            # First get the scaled/smeared energies
            signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptSmearUp", "ptSmearUp(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaUp)")
            signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptSmearDown", "ptSmearDown(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaDown)")
            if 'tt' not in channel:
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_energyScaleUp", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_energyScaleDown", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
            signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptScaleUp", "Photon_pt*Photon_energyScaleUp")
            signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptScaleDown", "Photon_pt*Photon_energyScaleDown")

            # Then calculate the vertex/raw mass with scaled/smeared values
            for corr in ['Scale', 'Smear']:
                for d in ['Up', 'Down']:
                    signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("best_2g_{}{}_m{}_info".format(corr,d,m), "getDxy(Photon_pt{}{}, Photon_eta, Photon_phi, Photon_isScEtaEB, Photon_isScEtaEE, best_2g_idx1_m{}, best_2g_idx2_m{}, {})".format(corr, d, m, m, m))
                    signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("best_2g_dxy_{}{}_m{}".format(corr,d,m), "best_2g_{}{}_m{}_info[0]".format(corr,d,m))
                    signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("best_2g_raw_mass_{}{}_m{}".format(corr,d,m), "best_2g_{}{}_m{}_info[1]".format(corr,d,m))

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
            for channel in channels:
                print(channel)
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("reweight_{}to{}".format(ct,ctau), "reweightCT(GenPart_ctau[abs(GenPart_pdgId)==9000006], {}/10., {}/10.)".format(ct, ctau))
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].display("reweight_{}to{}".format(ct,ctau))
            exit()
        # Nominal signal
        sig ={}
        for channel in channels:
            sig[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", startVH.lumi[era], ("sig_{}_{}_m{}_ct{}".format(v,l,m,ct), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
            sig[channel] = startVH.unfoldTH2(sig[channel])
            sig[channel].Rebin(len(binsX) - 1, f"{channel}_HToPhiPhi_ctau{ctau}", binsX)
            sig[channel] = ROOT.gDirectory.Get(f"{channel}_HToPhiPhi_ctau{ctau}")
            sig[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}", ROOT.TObject.kOverwrite)

        # Energy scale/smearing corrections, scale factor corrections
        sigScaleUp = {}
        sigScaleDown = {}

        sigScaleUp[v] = {}
        sigScaleDown[v] = {}
    exit()