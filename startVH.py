from common.plotter import *
from python.VHTools.VHcuts import *
import ROOT
import glob
import array
import pandas as pd
import subprocess

ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')

#TODO make config file to carry all dictionaries
# Run 2 luminosity
lumi = {'2018': "59830",
        '2017': "41480",
        '2016': "36310",
        'Run2': '137620',
        'HEM': "38750",
        'preHEM': "21080"}

# Run 2 luminosity (integers) (too lazy to manually convert strings to ints)
intLumi = {'2018': 59830,
           '2017': 41480,
           '2016': 36310,
           'Run2': 137620,
           'HEM': 38750,
           'preHEM': 21080}

# Run 2 luminosity uncertainty
lumiUnc = {'2018': 1.025,
           '2017': 1.023,
           '2016': 1.012}

xsecs = {'ggZ': "0.1057", #units in pb
        'Z': "0.8696",
        'Wplus': "0.84",
        'Wminus': "0.5328",
        'tt':"0.5071"}

BRs = {"Z": "(.03363+.03366+.033696)", #e nue, mu numu, tau nutau
       "W": "(.1063+.1071+.1138)", #mu numu, e nue, tau nutau
       "ttSemiLeptonic": "2*(1 - (.1110 + .1140 + .107))*(.1110 + .1140 + .107)" #e nue b,mu numu b, tau nutau b. semileptonic decay
       }

# Phi masses and lifetimes
#masses= [15,20, 30,40, 50, 55]
ctaus = [0, 10, 20, 50, 100, 1000]
#ctaus = [1000]

# Cuts for W->e+nu 2018 HEM
cutsHEM = "(Electron_eta[W_l1_idx]>-1.3 || Electron_eta[W_l1_idx]<-3.0) && (Electron_phi[W_l1_idx]>-0.85 || Electron_phi[W_l1_idx]<-1.57)"

def isValidFile(fname, tree):
    f = ROOT.TFile(fname)
    valid = True
    try:
        t = f.Get(tree)
        t.GetEntry(0)
        t.Photon_pt
    except:
        f.Close()
        return False
    f.Close()
    return True

# Temp function to use official cut based ID instead of custom ID
# Remove once ntuples get reprocessed with correct ID (not needed after 03_26_24)
def redoPhotonID(plotter,era,isMC = False):
    plotter.define("PassPhIso" , "passPhIso(Photon_vidNestedWPBitmap)")
    for m in masses:
        plotter.define("best_2g_looseID1_m{}".format(m), "Photon_cutBased[best_2g_idx1_m{}]>0".format(m))
        plotter.define("best_2g_looseID2_m{}".format(m), "Photon_cutBased[best_2g_idx2_m{}]>0".format(m))
        plotter.redefine("best_2g_sumID_m{}".format(m), "best_2g_looseID1_m{m}+best_2g_looseID2_m{m}".format(m=m))
        plotter.filter("PassPhIso[best_2g_idx1_m{m}]==1 && PassPhIso[best_2g_idx2_m{m}]==1".format(m=m))
    if isMC:
        plotter.define("pho_SFs_id", "scaleFactors_2d(Photon_eta, Photon_pt, PHO_ID_{era}_sf, PHO_ID_{era}_binsX, PHO_ID_{era}_binsY, sample_isMC, Photon_cutBased>0)".format(era=era))
        plotter.redefine("Photon_idSF_val", "pho_SFs_id[0]")
        plotter.redefine("Photon_idSF_unc", "pho_SFs_id[1]")
        plotter.define("pho_SFs_pix", "getPixelSeedSF(Photon_isScEtaEB, Photon_isScEtaEE, hasPix_UL{era}_sf, sample_isMC, !Photon_pixelSeed)".format(era=era))
        plotter.redefine("Photon_pixSF_val", "pho_SFs_pix[0]")
        plotter.redefine("Photon_pixSF_unc", "pho_SFs_pix[1]")

#TODO: this function doesn't properly open the files yet
def getFiles(query,sampleDir,sampleType,era,prod):
    ser = pd.Series(subprocess.check_output(['xrdfs', 'root://cmseos.fnal.gov', 'ls', f"{sampleDir}/{sampleType}{era}_{prod}/"], text=True).split("\n"))
    return(list('root://cmseos.fnal.gov/' + ser[ser.str.contains(query)]))

def getPlotters(era,prod,sampleDir,mass_list,modelIndependent=False):
    dataEMU = {}     # merged_plotter of data
    dyPlotters = {}  # merged_plotter of dy+jets
    wjPlotters = {}  # merged_plotter of w+jets
    signal = {}      # merged_plotter of signal
    global masses
    global scaleFactors
    masses = mass_list
    for ana in ['zmm2g', 'zee2g', 'wen2g', 'wmn2g']:
        singleMuSamples = glob.glob("{d}/DATA{era}_{prod}/SingleMuon_*.root".format(d=sampleDir,era=era, prod=prod))
        #singleMuSamples = getFiles("SingleMuon_",sampleDir,sampleType="DATA",era=era,prod="03_26_24")
        EGSamples = glob.glob("{d}/DATA{era}_{prod}/EGamma_*.root".format(d=sampleDir,era=era, prod=prod))
        #EGSamples = getFiles("EGamma_",sampleDir,sampleType="DATA",era=era,prod="03_26_24")
        if era != '2018':
            EGSamples = glob.glob("{d}/DATA{era}_{prod}/SingleElectron_*.root".format(d=sampleDir,era=era, prod=prod))
            #EGSamples = getFiles("SingleElectron_",sampleDir,sampleType="DATA",era=era,prod="03_26_24")
        singleMuPlotters = []
        egPlotters = []

        for sample in singleMuSamples:
            if not isValidFile(sample, ana):
                continue
            singleMuPlotters.append(rdf_plotter(sample, tree=ana))
            redoPhotonID(singleMuPlotters[-1], era)

        for sample in EGSamples:
            if not isValidFile(sample, ana):
                continue
            egPlotters.append(rdf_plotter(sample, tree=ana))
            redoPhotonID(egPlotters[-1], era)
            if era == '2018' and ana == 'wen2g':
                egPlotters[-1].defaultCuts = "((run>=319077&&(Electron_eta[W_l1_idx]>-1.3||Electron_eta[W_l1_idx]<-3.0)&&(Electron_phi[W_l1_idx]>-0.87||Electron_phi[W_l1_idx]<-1.57))||(run<319077))"
        if 'e' in ana:
            dataEMU[ana] = merged_plotter(egPlotters)
        else:
            dataEMU[ana] = merged_plotter(singleMuPlotters)
        dySamples_nJ = glob.glob("{d}/MC{era}_{prod}/DYJetsToLL_?J_NLO_*.root".format(d=sampleDir,era=era, prod=prod))
        #dySamples_nJ = getFiles("DYJetsToLL_\dJ_NLO",sampleDir,sampleType="MC",era=era,prod="03_26_24")
        dyPlotters_nJ = []
        for sample in dySamples_nJ:
            if not isValidFile(sample, ana):
                continue
            dyPlotters_nJ.append(rdf_plotter(sample, True, tree = ana))
            if era == '2016':
                if 'APV' in sample or 'preVFP' in sample:
                    dyPlotters_nJ[-1].addCorrectionFactor("0.54", "flat")
                    redoPhotonID(dyPlotters_nJ[-1], '2016preVFP', isMC=True)
                else:
                    dyPlotters_nJ[-1].addCorrectionFactor("0.46", "flat")
                    redoPhotonID(dyPlotters_nJ[-1], '2016postVFP', isMC=True)
            else:
                redoPhotonID(dyPlotters_nJ[-1], era)
            if era == '2018' and ana=='wen2g':
                dyPlotters_nJ[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                dyPlotters_nJ.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                redoPhotonID(dyPlotters_nJ[-1], era, isMC=True)
                dyPlotters_nJ[-1].addCorrectionFactor(str(38750./59830), "flat")

        dyPlotters[ana] = merged_plotter(dyPlotters_nJ)
        dyPlotters[ana].setFillProperties(1001, ROOT.kAzure+5)
        dyPlotters[ana].setLineProperties(1, ROOT.kAzure+5, 3)

        wjSamples_nJ = glob.glob("{d}/MC{era}_{prod}/WJetsToLNu_?J_NLO_*.root".format(d=sampleDir,era=era, prod=prod))
        #wjSamples_nJ = getFiles("WJetsToLNu_\dJ_NLO",sampleDir,sampleType="MC",era=era,prod="03_26_24")
        wjPlotters_nJ = []
        for sample in wjSamples_nJ:
            if not isValidFile(sample, ana):
                continue
            wjPlotters_nJ.append(rdf_plotter(sample, True, tree = ana))
            if era == '2016':
                if 'APV' in sample or 'preVFP' in sample:
                    wjPlotters_nJ[-1].addCorrectionFactor("0.54", "flat")
                    redoPhotonID(wjPlotters_nJ[-1], "2016preVFP", isMC=True)
                else:
                    wjPlotters_nJ[-1].addCorrectionFactor("0.46", "flat")
                    redoPhotonID(wjPlotters_nJ[-1], "2016postVFP", isMC=True)
            else:
                redoPhotonID(wjPlotters_nJ[-1], era, isMC=True)
            if era == '2018' and ana=='wen2g':
                wjPlotters_nJ[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                wjPlotters_nJ.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                redoPhotonID(wjPlotters_nJ[-1], era, isMC=True)
                wjPlotters_nJ[-1].addCorrectionFactor(str(38750./59830), "flat")

        wjPlotters[ana] = merged_plotter(wjPlotters_nJ)
        wjPlotters[ana].setFillProperties(1001, ROOT.kAzure-9)
        wjPlotters[ana].setLineProperties(1, ROOT.kAzure-9, 3)

    # Make signal plotters, one merge_plotter for each category
    #for sig in ['total']:
    for sig in ['ZH','ggZH','WH','ttH']:
    
        signal[sig] = {}
        for ana in ['zmm2g', 'zee2g', 'wen2g', 'wmn2g']:
            signal[sig][ana]  = {}        
            for m in masses:
                signal[sig][ana][m] = {}
                for ct in ctaus:
                    signal[sig][ana][m][ct] = {}
                    for br in ['2G2Q']:
                        sigPlotters = []
                        if sig == 'ttH':
                            V = ['tt']
                        elif sig == 'WH':
                            V = ['Wplus','Wminus']
                        elif sig == 'ZH':
                            V = ['Z']
                        elif sig == 'ggZH':
                            V = ['ggZ']
                        elif sig == 'total':
                            V = ['Z','ggZ','Wplus','Wminus','tt']
                        for v in V:
                            samples = glob.glob("{d}/MC{era}_{prod}/{v}H{br}_M{m}_ctau{ct}_UL{yr}*.root".format(d=sampleDir,era=era,prod=prod,v=v,br=br,m=m,ct=ct,yr=era[-2:]))
                            #samples = getFiles(f"{v}H{br}_M{m}_ctau{ct}_UL{era[-2:]}",sampleDir,sampleType="MC",era=era,prod="03_26_24")
                            for sample in samples:
                                if not isValidFile(sample, ana):
                                    continue
                                sigPlotters.append(rdf_plotter(sample, True, tree = ana))
                                # xsec*BR weight
                                weight = "(1)"
                                if not modelIndependent:
                                    weight +="*"+xsecs[v]
                                if "Z" in v:
                                    weight+="*"+BRs['Z']
                                elif 'W' in v:
                                    weight+="*"+BRs['W']
                                elif 'tt' in v:
                                    weight+="*"+BRs['ttSemiLeptonic']
                                sigPlotters[-1].addCorrectionFactor(weight, "flat")

                                #lumi/era normalizations
                                if era == '2016':
                                    if 'APV' in sample or 'preVFP' in sample:
                                        sigPlotters[-1].addCorrectionFactor("0.54", "flat")
                                        redoPhotonID(sigPlotters[-1], "2016preVFP", isMC=True)
                                    else:
                                        sigPlotters[-1].addCorrectionFactor("0.46", "flat")
                                        redoPhotonID(sigPlotters[-1], "2016postVFP", isMC=True)
                                else:
                                    redoPhotonID(sigPlotters[-1], era, isMC=True)
                                if era=='2018' and ana=='wen2g':
                                    #Monte Carlo normalized to unaffected luminosity of 21.08fb-1
                                    sigPlotters[-1].addCorrectionFactor(str(21080./59830), "flat")
                                    #Monte carlo normalized to affected luminosity of 38.75fb-1 veto events is affected region defined by HEM cuts
                                    sigPlotters.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                                    sigPlotters[-1].addCorrectionFactor(str(38750./59830), "flat")
                                    redoPhotonID(sigPlotters[-1], era, isMC=True)

                        signal[sig][ana][m][ct][br] = merged_plotter(sigPlotters)
                        signal[sig][ana][m][ct][br].setFillProperties(0, ROOT.kWhite)
                        signal[sig][ana][m][ct][br].setLineProperties(1, ROOT.kRed, 3)
                    

    # Scale factor stuff
    scaleFactors = {}
    scaleFactors['W'] = {'ELE': ['Electron_recoSF_val[W_l1_idx]', 'Electron_idSF_val[W_l1_idx]', 'Electron_trigSF_val[W_l1_idx]'],
                     'MU': ['Muon_recoSF_val[W_l1_idx]', 'Muon_idSF_val[W_l1_idx]', 'Muon_isoSF_val[W_l1_idx]', 'Muon_trigSF_val[W_l1_idx]']
    }

    scaleFactors['Z'] = {'ELE': ['Electron_recoSF_val[Z_idx[0]]', 'Electron_recoSF_val[Z_idx[1]]', 'Electron_idSF_val[Z_idx[0]]', 'Electron_idSF_val[Z_idx[1]]', 'Electron_trigSF_val[Z_idx[0]]'],
                     'MU': ['Muon_recoSF_val[Z_idx[0]]', 'Muon_recoSF_val[Z_idx[1]]', 'Muon_idSF_val[Z_idx[0]]', 'Muon_idSF_val[Z_idx[1]]', 'Muon_isoSF_val[Z_idx[0]]', 'Muon_isoSF_val[Z_idx[1]]', 'Muon_trigSF_val[Z_idx[0]]']
    }

    scaleFactors['g'] = {}
    for m in masses:
        scaleFactors['g'][m] = ['Photon_idSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_idSF_val[best_2g_idx2_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx2_m{}]'.format(m)]

    plots = {}
    plots['dataEMU'] = dataEMU
    plots['signal'] = signal
    plots['DYJets'] = dyPlotters
    plots['WJets'] = wjPlotters
    return plots

def getSF(scaleFactors, sfToVary = ''):
    if sfToVary == '':
        return "(" + ")*(".join(scaleFactors) + ")"
    toVary = list(filter(lambda x: sfToVary in x, scaleFactors))
    scaleFactors = list(filter(lambda x: sfToVary not in x, scaleFactors))
    sfUp = ["{}+{}".format(sf, sf.replace("val", "unc")) for sf in toVary]
    sfDown = ["{}-{}".format(sf, sf.replace("val", "unc")) for sf in toVary]
    return ["(" + ")*(".join(scaleFactors + sfUp) + ")", "(" + ")*(".join(scaleFactors + sfDown) + ")"]


# Useful dictionaries
ana = {}
ana['W'] = {'MU': "wmn2g",
            'ELE': 'wen2g'}
ana['Z'] = {'MU': 'zmm2g',
            'ELE': 'zee2g'}

headers = {}
headers['W'] = {'ELE': "W#rightarrowe#nu",
                'MU': "W#rightarrow#mu#nu"}
headers['Z'] = {'ELE': "Z#rightarrowee",
                'MU': "Z#rightarrow#mu#mu"}
                

# Histogram methods for data cards
def unfoldTH2(hist):
    binsX = hist.GetNbinsX()
    binsY = hist.GetNbinsY()
    hOut = ROOT.TH1D("hOut", "", binsX*binsY, 0, binsX*binsY)
    for x in range(binsX):
        for y in range(binsY):
            hOut.SetBinContent(1+y+x*binsY, hist.GetBinContent(x+1, y+1))
            hOut.SetBinError(1+y+x*binsY, hist.GetBinError(x+1, y+1))
    return hOut

def rebinTH1(hist, quantiles, binsM, binsLxy, lxyMin, lxyMax,m): #if highest mass and lowest statistics - use coarser quantiles
    q = {}
    bins = array.array('d')
    # Rebin each row of unfolded histogram using getquantiles
    for binM in range(binsM):
        quants = quantiles
        if binM == 0:
            quants = array.array('d', [0.5])
        elif binM == 3 and m == 15:
            quants = array.array('d', [0.5])
        bins.append(binsLxy * binM)
        q[binM] = array.array('d', [0]*len(quants))
        tmpTH1 = ROOT.TH1D("tmpTH1_m{}".format(binM), "", binsLxy, lxyMin, lxyMax)
        for i in range(1, binsLxy + 1):
            tmpTH1.SetBinContent(i, hist.GetBinContent(i + binsLxy*binM))
        tmpTH1.GetQuantiles(len(quants), q[binM], quants)
        for i in range(len(q[binM])):
            bins.append(binsLxy*binM + tmpTH1.GetXaxis().FindBin(q[binM][i]))
    bins.append(binsM*binsLxy)
    out = array.array('d')
    # Remove repeated bin entries
    for b in bins:
        if b in out:
            continue
        out.append(b)
    return outgit 