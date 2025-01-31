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
#skips ttH and ggZH for model indep
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


#cant be too fine 
if v == 'Z' and mass == 15 and l == 'ELE':
    binsM = 2
elif v == 'Z' and mass ==20 and l == 'ELE':
    binsM = 2
else:
    binsM =4
    
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

#source https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
#asymetric
xsecUnc = {'WH'  : [0.005, -0.007],
           'ZH'  : [0.038, -0.031],
           'ggZH': [0.251, -0.189],
           'ttH' : [0.058, -0.092]}

#symmetric
pdfUnc = {'WH'  : 1.019,
          'ZH'  : 1.016,
          'ggZH': 1.024,
          'ttH' : 1.036}

masses = [mass]
channels = ['ZH','ttH','WH','ggZH']

# Which lifetime to use for ctau reweighing
#lifetimes = {0: 0, 3: 10, 10: 10, 14: 10, 20: 20, 32: 20, 50: 50, 70: 50, 100: 100, 316: 1000, 1000: 1000}
lifetimes = {0:0, 10:10,20:20,50:50, 100:100, 1000:1000}

customBins = {
            'W': {
                'ELE': {
                    15:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    20:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    30:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    50:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                },
                'MU': {
                    15:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    20:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    30:array.array('d' , [0.0,25.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    50:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                }
            },
            'Z': {
                'ELE': {
                    15:array.array('d' , [0.0,36.0,72.0,110.0,220.0]),
                    20:array.array('d' , [0.0,110.0,138.0,165.0,192.0,220.0]),
                    30:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,350.0,370.0,440.0]),
                    50:array.array('d' , [0.0,50.0,110.0,138.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,110.0,138.0,192.0,220.0,248.0,275.0,330.0,360.0,440.0]),
                },
                'MU': {
                    15:array.array('d' , [0.0,50.0,80.0,110.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    20:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    30:array.array('d' , [0.0,50.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,350.0,370.0,440.0]),
                    50:array.array('d' , [0.0,50.0,80.0,110.0,138.0,192.0,220.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                }
            }
        }


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
    
    z_tight_sb = dataEMU[startVH.ana[v][l]].hist1d("Z_mass", cuts['sr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("z_loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))
    z_loose_sb = dataEMU[startVH.ana[v][l]].hist1d("Z_mass", cuts['sr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))
    
    # normalization scale factor
    scale = tight_sb.Integral() / loose_sb.Integral() if loose_sb.Integral() > 0 else 1.0
    loose_sb.Write("loose_sideband", ROOT.TObject.kOverwrite)
    tight_sb.Write("tight_sideband", ROOT.TObject.kOverwrite)
    

    # Get predicted background from antiID region
    background = dataEMU[startVH.ana[v][l]].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m], "1", ("data_loose_{}_{}_m{}".format(v,l,m), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
    background.Scale(scale)
    background.Write("data_loose", ROOT.TObject.kOverwrite)
    background = startVH.unfoldTH2(background)
    binsX = customBins[v][l][m]
    background.Rebin(len(binsX) - 1, "data_loose_rebinned", binsX)
    background = ROOT.gDirectory.Get("data_loose_rebinned")
    background = fillZeroBins(background, 0.001)
    background.Write("background", ROOT.TObject.kOverwrite)

    l = 'ELE'
    z = startVH.ana['Z'][l].hist1d("Z_mass", cuts['Z'][l]+"&&"+cuts['photons'][options.mass]+"&&best_2g_sumID_m30=={}".format(i), startVH.lumi[year], ("z_mass_{}_preselection{}".format(l,sfx), "", 20, 70, 110), SFs = sfs['Z'][l][options.mass], titlex = "m({}{}) [GeV]".format(lep,lep))

    Z_mass_loose_sb = loose_sb = dataEMU[startVH.ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", ("loose_sb_{}_{}_m{}".format(v,l,m), "", 50, -1400, 0))

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
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("reweight_{}to{}".format(ct,ctau), "reweightCT(GenPart_ctau[abs(GenPart_pdgId)==9000006], {}/10., {}/10.)".format(ct, ctau))

        # Nominal signal
        sig ={}
        remove = []
        for channel in channels:
            sig[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", startVH.lumi[era], ("sig_{}_{}_m{}_ct{}".format(v,l,m,ct), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
            if sig[channel] == None:
                print(f'\nRemoving channel {channel}\n')
                remove.append(channel)
            else:
                sig[channel] = startVH.unfoldTH2(sig[channel])
                sig[channel].Rebin(len(binsX) - 1, f"{channel}_HToPhiPhi_ctau{ctau}", binsX)
                sig[channel] = ROOT.gDirectory.Get(f"{channel}_HToPhiPhi_ctau{ctau}")
                sig[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}", ROOT.TObject.kOverwrite)
        
        for channel in remove:
            channels.remove(channel)

        # Energy scale/smearing corrections, scale factor corrections
        sigScaleUp = {}
        sigScaleDown = {}

        sigScaleUp[v] = {}
        sigScaleDown[v] = {}
        if l == 'ELE':
            for sf in ['reco', 'id', 'trig']:
                sfUp, sfDown = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m], 'Electron_{}'.format(sf))
                sigUp = {}
                sigDown = {}
                sigScaleUp[v][sf] = {}
                sigScaleDown[v][sf] = {}
                for channel in channels:
                    sigUp[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfUp+")*("+reweight+")", startVH.lumi[era], ("sig_ele{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ct), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                    sigUp[channel] = startVH.unfoldTH2(sigUp[channel])
                    sigUp[channel].Rebin(len(binsX) - 1, f"sig_ele{sf}Up_{v}_{l}_m{m}_ct{ctau}", binsX)
                    sigUp[channel] = ROOT.gDirectory.Get(f"sig_ele{sf}Up_{v}_{l}_m{m}_ct{ctau}")
                    sigUp[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_ele_{sf}Up", ROOT.TObject.kOverwrite)
                    sigScaleUp[v][sf][channel] = fillZeroBins(sigUp[channel], 0)

                    sigDown[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfDown+")*("+reweight+")", startVH.lumi[era], ("sig_ele{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                    sigDown[channel] = startVH.unfoldTH2(sigDown[channel])
                    sigDown[channel].Rebin(len(binsX) - 1, f"sig_ele{sf}Down_{v}_{l}_m{m}_ct{ctau}", binsX)
                    sigDown[channel] = ROOT.gDirectory.Get(f"sig_ele{sf}Down_{v}_{l}_m{m}_ct{ctau}")
                    sigDown[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_ele_{sf}Down", ROOT.TObject.kOverwrite)

                    sigScaleDown[v][sf][channel] = fillZeroBins(sigDown[channel], 0)

        if l == 'MU':
            for sf in ['reco', 'id', 'trig', 'iso']:
                sfUp, sfDown = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m], 'Muon_{}'.format(sf))
                sigUp={}
                sigDown={}
                sigScaleUp[v][sf] = {}
                sigScaleDown[v][sf] = {}
                for channel in channels:
                    sigUp[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfUp+")*("+reweight+")", startVH.lumi[era], ("sig_mu{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                    sigUp[channel] = startVH.unfoldTH2(sigUp[channel])
                    sigUp[channel].Rebin(len(binsX) - 1, f"sig_mu{sf}Up_{v}_{l}_m{m}_ct{ctau}", binsX)
                    sigUp[channel] = ROOT.gDirectory.Get(f"sig_mu{sf}Up_{v}_{l}_m{m}_ct{ctau}")
                    sigUp[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_mu_{sf}Up", ROOT.TObject.kOverwrite)

                    sigScaleUp[v][sf][channel] = fillZeroBins(sigUp[channel], 0)

                    sigDown[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfDown+")*("+reweight+")", startVH.lumi[era], ("sig_mu{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                    sigDown[channel] = startVH.unfoldTH2(sigDown[channel])
                    sigDown[channel].Rebin(len(binsX) - 1, f"sig_mu{sf}Down_{v}_{l}_m{m}_ct{ctau}", binsX)
                    sigDown[channel] = ROOT.gDirectory.Get(f"sig_mu{sf}Down_{v}_{l}_m{m}_ct{ctau}")
                    sigDown[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_mu_{sf}Down", ROOT.TObject.kOverwrite)

                    sigScaleDown[v][sf][channel] = fillZeroBins(sigDown[channel], 0)

        sigScaleUp['g'] = {}
        sigScaleDown['g'] = {}
        for sf in ['id', 'pix']:
            sfUp, sfDown = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m], "Photon_{}".format(sf))
            sigScaleUp['g'][sf] = {}
            sigScaleDown['g'][sf] = {}
            for channel in channels:
                sigUp[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfUp+")*("+reweight+")", startVH.lumi[era], ("sig_pho{}Up_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigUp[channel] = startVH.unfoldTH2(sigUp[channel])
                sigUp[channel].Rebin(len(binsX) - 1, f"sig_pho{sf}Up_{v}_{l}_m{m}_ct{ctau}", binsX)
                sigUp[channel] = ROOT.gDirectory.Get(f"sig_pho{sf}Up_{v}_{l}_m{m}_ct{ctau}")
                sigUp[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_pho_{sf}Up", ROOT.TObject.kOverwrite)

                sigScaleUp['g'][sf][channel] = fillZeroBins(sigUp[channel], 0)

                sigDown[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_m{}".format(m), "best_2g_dxy_m{}".format(m), "("+cuts['sr'][v][l][m]+")*("+sfDown+")*("+reweight+")", startVH.lumi[era], ("sig_pho{}Down_{}_{}_m{}_ct{}".format(sf,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigDown[channel] = startVH.unfoldTH2(sigDown[channel])
                sigDown[channel].Rebin(len(binsX) - 1, f"sig_pho{sf}Down_{v}_{l}_m{m}_ct{ctau}", binsX)
                sigDown[channel] = ROOT.gDirectory.Get(f"sig_pho{sf}Down_{v}_{l}_m{m}_ct{ctau}")
                sigDown[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_pho_{sf}Down", ROOT.TObject.kOverwrite)

                sigScaleDown['g'][sf][channel] = fillZeroBins(sigDown[channel], 0)
        
        for corr in ['Scale', 'Smear']:
            sigUp = {}
            sigDown = {}
            sigScaleUp['g'][corr] = {}
            sigScaleDown['g'][corr]={}
            for channel in channels:
                sigUp[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_{}Up_m{}".format(corr,m), "best_2g_dxy_{}Up_m{}".format(corr,m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", startVH.lumi[era], ("sig_pho{}Up_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigUp[channel] = startVH.unfoldTH2(sigUp[channel])
                sigUp[channel].Rebin(len(binsX) - 1, f"sig_pho{corr}Up_{v}_{l}_m{m}_ct{ctau}", binsX)
                sigUp[channel] = ROOT.gDirectory.Get(f"sig_pho{corr}Up_{v}_{l}_m{m}_ct{ctau}")
                sigUp[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_pho_{corr}Up", ROOT.TObject.kOverwrite)

                sigScaleUp['g'][corr][channel] = fillZeroBins(sigUp[channel], 0)

                sigDown[channel] = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d("best_2g_raw_mass_{}Down_m{}".format(corr,m), "best_2g_dxy_{}Down_m{}".format(corr,m), "("+cuts['sr'][v][l][m]+")*("+sfs+")*("+reweight+")", startVH.lumi[era], ("sig_pho{}Down_{}_{}_m{}_ct{}".format(corr,v,l,m,ctau), "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
                sigDown[channel] = startVH.unfoldTH2(sigDown[channel])
                sigDown[channel].Rebin(len(binsX) - 1, f"sig_pho{corr}Down_{v}_{l}_m{m}_ct{ctau}", binsX)
                sigDown[channel] = ROOT.gDirectory.Get(f"sig_pho{corr}Down_{v}_{l}_m{m}_ct{ctau}")
                sigDown[channel].Write(f"{channel}_HToPhiPhi_ctau{ctau}_pho_{corr}Down", ROOT.TObject.kOverwrite)

                sigScaleDown['g'][corr][channel] = fillZeroBins(sigDown[channel], 0)

        # Now make the cards
        
        for i in range(data_obs.GetNbinsX()):
            skip_channel = {channel: False for channel in channels}
            for channel in channels:        
                # Check bin content for the current channel
                if sig[channel].GetBinContent(i + 1) <= 0:
                    skip_channel[channel] = True  # Mark channel to be skipped
            kept_channels = []
            for channel in channels:
                if not skip_channel[channel]:
                    kept_channels.append(channel)

            # Process channels that are not marked to skip
            for channel in channels:
                if skip_channel[channel]:
                    continue  # Skip processing for this channel
            
                # Process the channel here if it's not skipped
                print(f"Processing {channel} at bin {i} for lifetime {ctau}")

                card = open("datacard_{}_bin{}.txt".format(cat, i), "w")
                card.write("imax 1\n")
                card.write("jmax *\n")
                card.write("kmax *\n")
                card.write('--------------------------------------------------\n')
                card.write("bin {}_bin{}\n".format(cat, i))
                card.write("observation {}\n".format(data_obs.GetBinContent(i+1)))
                card.write('--------------------------------------------------\n')
                bin = "bin\t"
                process = "process\t"
                iprocess= "process\t"
                rate="rate\t"
                for ichannel,channel in enumerate(kept_channels):
                    bin+=f"{cat}_bin{i}\t"
                    process+=f"{channel}_HToPhiPhi_ctau{ctau}\t"
                    iprocess+=f"{-ichannel}\t"
                    rate+=f"{0.1*sig[channel].GetBinContent(i+1)}\t"
                bin+=f"{cat}_bin{i}\n"
                process+="\tbackground\n"
                iprocess+="1\n"
                rate+=f"{background.GetBinContent(i+1)}\n"
                card.write(bin)
                card.write(process)
                card.write(iprocess)
                card.write(rate)
                card.write('--------------------------------------------------\n')
                card.write(f"CMS.lumi_{era}\tlnN\t"+f"{startVH.lumiUnc[era]}\t"*len(kept_channels)+"-\n")
                card.write(f"CMS_pileup_{era}\tlnN\t"+"1.05\t"*len(kept_channels)+"-\n")
                if options.model == "SM":
                        for ichannel,channel in enumerate(kept_channels):
                            card.write(f"CMS_xsec_{channel}\tlnN\t" + 
                                    "-\t"*ichannel +
                                    f"{(1.0+xsecUnc[channel][0])/1.0}/{(1.0+xsecUnc[channel][1])/1.0}\t" +
                                    "-\t"*(len(kept_channels) - ichannel - 1) +
                                    "-\n")
                            card.write(f"CMS_pdf_{channel}\tlnN\t"+
                                    "-\t"*ichannel +
                                    f"{pdfUnc[channel]}\t" +
                                    "-\t"*(len(kept_channels) - ichannel - 1) +
                                    "-\n")
                        #card.write(f"CMS_xsec_{v}\tlnN\t"+f"{(1.0+xsecUnc[v][0])/1.0}/{(1.0+xsecUnc[v][1])/1.0}\t"*len(channels)+"-\n")
                        #card.write(f"CMS_pdf_{v}\tlnN\t"+f"{pdfUnc[v]}\t"*len(channels)+"-\n")
                card.write(f"CMS_bkg_{cat}_norm\tgmN\t{int(tight_sb.Integral())}\t"+f"-\t"*len(kept_channels)+f"{background.GetBinContent(i+1)/tight_sb.Integral()}\n")
                card.write(f"CMS_bkg_{cat}_bin{i}_shape\tlnN\t"+f"-\t"*len(kept_channels)+"1.4\n")
                
                for sf in sigScaleUp[v]:
                    s = f"CMS_{l}_{sf}\tlnN\t"
                    for ichannel,channel in enumerate(kept_channels):
                        errUp = sigScaleUp[v][sf][channel].GetBinContent(i+1) / sig[channel].GetBinContent(i+1)
                        errDown = sigScaleDown[v][sf][channel].GetBinContent(i+1) / sig[channel].GetBinContent(i+1)
                        if errUp == 0:
                            errUp = 1.0
                        if errDown == 0:
                            errDown = 1.0
                        s+=f"{errDown}/{errUp}\t"
                    s+="-\n"
                    card.write(s)
                for sf in sigScaleUp['g']:
                    s = f"CMS_gamma_{sf}\tlnN\t"
                    for ichannel, channel in enumerate(kept_channels):
                        errUp = sigScaleUp['g'][sf][channel].GetBinContent(i+1) / sig[channel].GetBinContent(i+1)
                        errDown = sigScaleDown['g'][sf][channel].GetBinContent(i+1) / sig[channel].GetBinContent(i+1)
                        if errUp == 0:
                            errUp = 1.0
                        if errDown == 0:
                            errDown = 1.0
                        s+=f"{errDown}/{errUp}\t"
                    s+="-\n"
                    card.write(s)
                card.close()
            print(f'done with lifetime {ctau}')