import ROOT, sys, os, array
sys.path.append("/uscms/home/gfavila/nobackup/rdf_simple_/")
ROOT.gROOT.SetBatch(True)
import numpy as np
from python.VHTools.py_helpers import *
from python.VHTools.config import *
from python.VHTools.VHcuts import * 
from typing import Callable


ROOT.gInterpreter.Declare('#include "/uscms/home/gfavila/nobackup/rdf_simple_/common/datacardHelpers.h"')
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
if v == 'Z' and mass == 7 and l == 'ELE':
    binsM = 2
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

if not os.path.isdir("plots"):
    os.mkdir("plots")


masses = [mass]


lifetimes = {0:0,10:10,20:20,50:50,100:100,1000:1000}
#lifetimes = {0}


#for m in masses:
m = mass

#Rebin according to custom bins
binsX = customBins[v][l][m]

# Scale factors
sfs = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m])
final_selection_cuts = f"{cuts['sr'][v][l][m]}"

def unfoldRebin(histogram,name,binsX = binsX):
    unfolded = startVH.unfoldTH2(histogram)
    unfolded.Rebin(len(binsX) - 1, name, binsX)
    rebinned = ROOT.gDirectory.Get(name)
    return(rebinned)

def data_histogram_setup(binsM:int = binsM, binsLxy:int = binsLxy, lxyMin:float = lxyMin, lxyMax:float = lxyMax, v:str=v,l:str=l,m:int=m) -> Callable:
    def make_histogram(name:str,appliedCuts:str):
        histogram = dataEMU[startVH.ana[v][l]].hist2d(
            f"best_2g_raw_mass_m{m}",
            f"best_2g_dxy_m{m}",
            f"({appliedCuts})", 
            "(1)", 
            (name+f"_{v}_{l}_{m}", "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
        return histogram
    return make_histogram

def monteCarlo_histogram_setup(appliedCuts:str =final_selection_cuts, binsM:int = binsM, binsLxy:int = binsLxy, lxyMin:float = lxyMin, lxyMax:float = lxyMax, v:str=v,l:str=l,m:int=m) -> Callable:
    def make_channel_histogram(ct:int, channel:str, name:str, scaleFactors:list[str], correction:str=""):
        histogram = signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].hist2d(
            f"best_2g_raw_mass{correction}_m{m}",
            f"best_2g_dxy{correction}_m{m}",
            f"({appliedCuts})*({scaleFactors})", 
            startVH.lumi[era], 
            (name, "", binsM, 4, m+5, binsLxy, lxyMin, lxyMax))
        return histogram
    return make_channel_histogram

def processChannels(name,scaleFactors,ct,appliedCuts=final_selection_cuts,histogram_suffix = '',correction="") -> dict[str]:
    sig = {}
    channels_to_remove = []
    global channels
    for channel in channels:
        if (histogram := make_histogram(ct,channel,name,scaleFactors,correction)):
            histogram  = unfoldRebin(histogram,name)
            histogram.Write(f"{channel}_HToPhiPhi_ctau{ct}" + histogram_suffix, ROOT.TObject.kOverwrite)
            sig[channel] = histogram
        else:
            print("Removing channel ",channel)
            channels_to_remove.append(channel)
            continue

    for channel in channels_to_remove:
        channels.remove(channel)
        
    return sig

def processLeptonScaleFactors(l:str, ct:int, v:str=v) -> tuple:
    # Energy scale/smearing corrections, scale factor corrections
    sigScaleUp = {}
    sigScaleDown = {}

    sigScaleUp[v] = {}
    sigScaleDown[v] = {}

    if l == 'ELE':
        lepton = "Electron"
        lepton_lowerCase   = 'ele'
        list_of_scaleFactors = ['reco', 'id', 'trig']
    else:
        lepton = "Muon"
        lepton_lowerCase   = 'mu'
        list_of_scaleFactors = ['reco', 'id', 'trig','iso']
    
    for scaleFactor in list_of_scaleFactors:
        
        scaleFactor_name = f'{lepton}_{scaleFactor}'
        
        sfUp,sfDown = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m], scaleFactor_name)
        
        sigScaleUp[v][scaleFactor] = {}
        sigScaleDown[v][scaleFactor] = {}

        sigUp   = processChannels(f"sig_{lepton_lowerCase}{scaleFactor}Up_{v}_{l}_m{m}_ct{ct}",  sfUp,   ct, histogram_suffix = f"_{lepton_lowerCase}_{scaleFactor}Up")
        sigDown = processChannels(f"sig_{lepton_lowerCase}{scaleFactor}Down_{v}_{l}_m{m}_ct{ct}",sfDown, ct, histogram_suffix = f"_{lepton_lowerCase}_{scaleFactor}Down")
        
        for channel in channels:
            sigScaleUp[v][scaleFactor][channel] = fillZeroBins(sigUp[channel], 0)
            sigScaleDown[v][scaleFactor][channel] = fillZeroBins(sigDown[channel], 0)

    return sigScaleUp, sigScaleDown

def processPhotonScaleFactors(factors_toProcess,ct):
    # Energy scale/smearing corrections, scale factor corrections
    sigScaleUp   = {}
    sigScaleDown = {}

    for scaleFactor in factors_toProcess:
        scaleFactor_name = f"{scaleFactor}"

        sfUp, sfDown = startVH.getSF(startVH.scaleFactors[v][l] + startVH.scaleFactors['g'][m], f"Photon_{scaleFactor}")
        sigScaleUp[scaleFactor] = {}
        sigScaleDown[scaleFactor] = {}

        sigUp   = processChannels(f"sig_pho{scaleFactor}Up_{v}_{l}_m{m}_ct{ct}",  sfUp,   ct, histogram_suffix = f"_pho_{scaleFactor}Up")
        sigDown = processChannels(f"sig_pho{scaleFactor}Down_{v}_{l}_m{m}_ct{ct}",sfDown, ct, histogram_suffix = f"_pho_{scaleFactor}Down")

        for channel in channels:
            sigScaleUp[scaleFactor][channel] = fillZeroBins(sigUp[channel], 0)
            sigScaleDown[scaleFactor][channel] = fillZeroBins(sigDown[channel], 0)

    return sigScaleUp, sigScaleDown

def processPhotonCorrections(corrections:list[str], ct:int) -> tuple:
    # Energy scale/smearing corrections, scale factor corrections
    sigScaleUp   = {}
    sigScaleDown = {}
    for correction in corrections:
        sigScaleUp[correction]  = {}    
        sigScaleDown[correction] = {}

        sigUp   = processChannels(f"sig_pho{correction}Up_{v}_{l}_m{m}_ct{ctau}",  sfs,  ct,  histogram_suffix = f"_pho_{correction}Up",correction=f"_{correction}Up")
        sigDown = processChannels(f"sig_pho{correction}Down_{v}_{l}_m{m}_ct{ctau}",sfs, ct,  histogram_suffix = f"_pho_{correction}Down",correction=f"_{correction}Down")

        for channel in channels:
            sigScaleUp[correction][channel] = fillZeroBins(sigUp[channel], 0)
            sigScaleDown[correction][channel] = fillZeroBins(sigDown[channel], 0)
    
    return sigScaleUp, sigScaleDown


print(f"Making card for {era} {v}->{l}, m={m}")
# Getting the required plotters
dataEMU = startVH.getDataPlotters(era,prod,"/store/user/gfavila/DDP",startVH.ana[v][l],modelIndependent=modelIndependent)
signal  = startVH.getSigPlotters(era,prod,"/store/user/gfavila/DDP",[m],startVH.ana[v][l],modelIndependent=modelIndependent)['signal']


f = ROOT.TFile(f"datacardInputs_{v}H_{l}_m{m}_{era}{model_sfx}.root", "UPDATE")

# negative lxy sidebands for loose (anti-ID) and tight (ID)
loose_sb = dataEMU[startVH.ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['cr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", (f"loose_sb_{v}_{l}_m{m}", "", 50, -1400, 0))
tight_sb = dataEMU[startVH.ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['sr'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1", (f"loose_sb_{v}_{l}_m{m}", "", 50, -1400, 0))
closure_sb = dataEMU[startVH.ana[v][l]].hist1d("best_2g_dxy_m{}".format(m), cuts['cl'][v][l][m]+"&&best_2g_raw_mass_m{}>65".format(m), "1",(f"closure_sb_{v}_{l}_m{m}", "", 50, -1400, 0))

#arrays to store historam integral errors
loose_error = np.array([0.0],dtype=float)
tight_error = np.array([0.0],dtype=float)
closure_error = np.array([0.0],dtype=float)

#integrals
loose_integral = loose_sb.IntegralAndError(1,loose_sb.GetNbinsX(),loose_error)
tight_integral = tight_sb.IntegralAndError(1,tight_sb.GetNbinsX(),tight_error)
closure_integral = closure_sb.IntegralAndError(1,closure_sb.GetNbinsX(),closure_error)

# normalization scale factor
scale = tight_integral / loose_integral if loose_integral > 0 else 1.0
closure_scale = loose_integral / closure_integral if closure_integral > 0 else 1.0


#Error propagation, histograms are uncorrelated since there are no overlapping events
scale_error = scale * np.sqrt((loose_error[0]/loose_integral)**2 + (tight_error[0]/tight_integral)**2)
closure_scale_error = closure_scale * np.sqrt((loose_error[0]/loose_integral)**2 + (closure_error[0]/closure_integral)**2)

loose_sb.Write("loose_sideband", ROOT.TObject.kOverwrite)
tight_sb.Write("tight_sideband", ROOT.TObject.kOverwrite)
closure_sb.Write("closure_sideband", ROOT.TObject.kOverwrite)


make_data_histogram = data_histogram_setup()

# Get predicted background from antiID region
background = make_data_histogram("data_loose", cuts['cr'][v][l][m])
background.Write("data_loose", ROOT.TObject.kOverwrite)
background = unfoldRebin(background, "data_loose_rebinned")
background.Scale(scale)
background.Write("background", ROOT.TObject.kOverwrite)

#closure
closure_2bad = make_data_histogram("closure_2bad", cuts['cl'][v][l][m])
closure_2bad.Write("closure_loose", ROOT.TObject.kOverwrite)
closure_2bad = unfoldRebin(closure_2bad, "closure_loose_rebinned")
#before scaling, get poisson errors then return scaled python arrays
closure_loose_bin_centers, closure_loose_array, closure_loose_err_low, closure_loose_err_up, closure_loose_bin_widths = getPoisson_arrays(closure_2bad , closure_scale)
closure_2bad.Scale(closure_scale)
closure_2bad.Write("closure_loose_scaled", ROOT.TObject.kOverwrite)

#this is the same as the background histogram. Keep for now, will do clean up later
closure_1bad = make_data_histogram("closure_1bad", cuts['cr'][v][l][m])
closure_1bad.Write("closure_tight_2d", ROOT.TObject.kOverwrite)
closure_1bad = unfoldRebin(closure_1bad,"closure_tight_rebinned")
closure_tight_bin_centers, closure_tight_array, closure_tight_err_low, closure_tight_err_up, closure_tight_bin_widths = getPoisson_arrays(closure_1bad)
closure_1bad.Write("closure_tight", ROOT.TObject.kOverwrite)


#sigmas 
sigmas = []
deltas_2bad = []
deltas_1bad = []
for i in range(len(binsX)-1):
    N_2bad = closure_loose_array[i]
    N_1bad = closure_tight_array[i]
    
    if N_2bad == 0 and N_1bad == 0:
        sigma_x = 0
    
    if N_2bad > N_1bad:
        delta_2bad = closure_loose_err_low[i]
        delta_1bad = closure_tight_err_up[i]
    else:
        delta_2bad = closure_loose_err_up[i]
        delta_1bad = closure_tight_err_low[i]
    
    deltas_1bad.append(delta_1bad)
    deltas_2bad.append(delta_2bad)
    
    if (N_2bad - N_1bad)**2 - delta_2bad**2 - delta_1bad**2 < 0:
        sigma_x = 0
    else:
        sigma_x = np.sqrt((N_2bad - N_1bad)**2 - delta_2bad**2 - delta_1bad**2)
    
    sigmas.append(sigma_x)

uncertainties = []
for i in range(len(closure_loose_array)):
    if closure_loose_array[i] != 0 and sigmas[i] != 0.0:
        uncertainties.append(1 + sigmas[i]/closure_loose_array[i])
    else:
        uncertainties.append("nan")


# Get observed data
data_obs = make_data_histogram("data_tight",cuts['sr'][v][l][m])
data_obs = unfoldRebin(data_obs, "data_tight_rebinned")
# Only write if not blind, else use expected background
if not options.blind:
    data_obs.Write("data_obs", ROOT.TObject.kOverwrite)
else:
    data_obs = ROOT.TH1D("data_obs", "", len(binsX) - 1, binsX)
    for i in range(1, data_obs.GetNbinsX()+1):
        data_obs.SetBinContent(i, int(round(background.GetBinContent(background.FindBin(background.GetBinCenter(i))))))
    data_obs.Write("data_obs", ROOT.TObject.kOverwrite)

#TODO import signals from config file, which doesn't exist yet
# Defining photon scale corrections
for ct in startVH.ctaus:
    for channel in channels:
        # energy scale/smearing corrections
        # First get the scaled/smeared energies
        signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptSmearUp", "ptSmearUp(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaUp)")
        signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptSmearDown", "ptSmearDown(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaDown)")
        if 'tt' not in channel:
            try:
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_energyScaleUp", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
            except:
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].redefine("Photon_energyScaleUp", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
            try:
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_energyScaleDown", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
            except:
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].redefine("Photon_energyScaleDown", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, sample_isMC)".format(era = "2016preVFP" if options.era=="2016" else options.era))
        signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptScaleUp", "Photon_pt*Photon_energyScaleUp")
        signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("Photon_ptScaleDown", "Photon_pt*Photon_energyScaleDown")

        # Then calculate the vertex/raw mass with scaled/smeared values
        for corr in ['Scale', 'Smear']:
            for d in ['Up', 'Down']:
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("best_2g_{}{}_m{}_info".format(corr,d,m), "getDxy(Photon_pt{}{}, Photon_eta, Photon_phi, Photon_isScEtaEB, Photon_isScEtaEE, best_2g_idx1_m{}, best_2g_idx2_m{}, {})".format(corr, d, m, m, m))
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("best_2g_dxy_{}{}_m{}".format(corr,d,m), "best_2g_{}{}_m{}_info[0]".format(corr,d,m))
                signal[channel][startVH.ana[v][l]][m][ct]['2G2Q'].define("best_2g_raw_mass_{}{}_m{}".format(corr,d,m), "best_2g_{}{}_m{}_info[1]".format(corr,d,m))

make_histogram = monteCarlo_histogram_setup()

# Now get expected signal for each lifetime
for ctau in lifetimes:
    # Category name
    cat = "{}H_{}_m{}_ctau{}_{}{}".format(v,l,m,ctau,era,model_sfx)
    
    # Nominal signal
    sig = processChannels(f"sig_{v}_{l}_m{m}_ct{ct}",sfs,ctau)
    
    sigScaleUp = {}
    sigScaleDown = {}

    Lepton_sigScaleUp, Lepton_sigScaleDown = processLeptonScaleFactors(l,ctau)
    print(Lepton_sigScaleUp)
    sigScaleUp.update(Lepton_sigScaleUp)
    sigScaleDown.update(Lepton_sigScaleDown)

    Photon_corrScaleUp, Photon_corrScaleDown = processPhotonCorrections(['Scale', 'Smear'],ctau)
    Photon_sigScaleUp, Photon_sigScaleDown  = processPhotonScaleFactors(['id', 'pix'],ctau)

    
    Photon_corrScaleUp.update(Photon_sigScaleUp)
    Photon_corrScaleDown.update(Photon_sigScaleDown)

    
    sigScaleUp['g'] = Photon_corrScaleUp
    sigScaleDown['g'] = Photon_corrScaleDown

    # Now make the cards
    print(channels)
    for ibin in range(data_obs.GetNbinsX()):
        skip_channel  = {channel:True if sig[channel].GetBinContent(ibin+1) <= 0 else False for channel in channels}
        kept_channels = [channel for channel in channels if not skip_channel[channel]]
        
        for channel in channels:
            if skip_channel[channel]:
                continue
    
            print(f"Processing {channel} at bin {ibin} for lifetime {ctau}, era {era}")
            card = open(f"datacard_{cat}_bin{ibin}.txt", "w")
            card.write("#MC BR(Higgs->Phi Phi)(Phi->gamma gamma) = 0.5\n")
            card.write("#Signal processes scaled by 0.1\n")
            card.write('#BR_limit = r_limit * 0.05\n')
            card.write('\n')
            card.write("imax 1\n")
            card.write("jmax *\n")
            card.write("kmax *\n")
            card.write('--------------------------------------------------\n')
            card.write("bin {}_bin{}\n".format(cat, ibin))
            card.write("observation {}\n".format(data_obs.GetBinContent(ibin+1)))
            card.write('--------------------------------------------------\n')
            bin = "bin\t"
            process = "process\t"
            iprocess= "process\t"
            rate="rate\t"
            for ichannel,channel in enumerate(kept_channels):
                bin+=f"{cat}_bin{ibin}\t"
                process+=f"{channel}_HToPhiPhi_ctau{ctau}\t"
                iprocess+=f"{-ichannel}\t"
                rate+=f"{0.1*sig[channel].GetBinContent(ibin+1)}\t"
            bin+=f"{cat}_bin{ibin}\n"
            process+="\tbackground\n"
            iprocess+="1\n"
            rate+=f"{background.GetBinContent(ibin+1)}\n" 
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
                            f"{(1.0+xsecUnc[channel][0])/1.0}/{(1.0+xsecUnc[channel][1])/1.0}\t" + #uncertainties are in percentages
                            "-\t"*(len(kept_channels) - ichannel - 1) +
                            "-\n")
                    card.write(f"CMS_pdf_{channel}\tlnN\t"+
                            "-\t"*ichannel +
                            f"{pdfUnc[channel]}\t" +
                            "-\t"*(len(kept_channels) - ichannel - 1) +
                            "-\n")
            card.write(f"CMS_bkg_{cat}_bin{ibin}_norm\tgmN\t{round(background.GetBinContent(ibin+1)/scale)}\t"+f"-\t"*len(kept_channels)+f"{scale}\n")
            card.write(f"CMS_bkg_{cat}_shape\tlnN\t"+f"-\t"*len(kept_channels)+f"{1 + scale_error/scale}"+"\n")

            if uncertainties[ibin] == 'nan':
                continue
            else:
                card.write(f"CMS_bkg_{cat}_bin{ibin}_closure\tlnN\t"+f"-\t"*len(kept_channels)+f"{uncertainties[ibin]}"+"\n")


            for sf in sigScaleUp[v]:
                s = f"CMS_{l}_{sf}\tlnN\t"
                for ichannel,channel in enumerate(kept_channels):
                    errUp = sigScaleUp[v][sf][channel].GetBinContent(ibin+1) / sig[channel].GetBinContent(ibin+1)
                    errDown = sigScaleDown[v][sf][channel].GetBinContent(ibin+1) / sig[channel].GetBinContent(ibin+1)
                    if errUp == 0:
                        errUp = 1.0
                    if errDown == 0:
                        errDown = 1.0
                    if errDown == 1.0 and errUp == 1.0:
                        s+= "-\t"
                    else:
                        s+=f"{errDown}/{errUp}\t"
                s+="-\n"
                card.write(s)
            for sf in sigScaleUp['g']:
                s = f"CMS_gamma_{sf}\tlnN\t"
                for ichannel, channel in enumerate(kept_channels):
                    errUp = sigScaleUp['g'][sf][channel].GetBinContent(ibin+1) / sig[channel].GetBinContent(ibin+1)
                    errDown = sigScaleDown['g'][sf][channel].GetBinContent(ibin+1) / sig[channel].GetBinContent(ibin+1)
                    if errUp == 0:
                        errUp = 1.0
                    if errDown == 0:
                        errDown = 1.0
                    if errDown == 1.0 and errUp == 1.0:
                        s+= "-\t"
                    else:
                        s+=f"{errDown}/{errUp}\t"
                s+="-\n"
                card.write(s)
            card.close()
    print(f'done with lifetime {ctau} {era} {v} {l} category')
    ROOT.gSystem.Exit(0)