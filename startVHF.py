from common.plotter import *
from python.VHTools.VHFcuts import *
from python.VHTools.config import *
from scipy.stats import beta
import ROOT,sys,os
import glob
import array
import pandas as pd
import subprocess
import matplotlib
import datetime

ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')

cuts = {}
# Tight + veto lepton cuts
cuts['W'] = {'MU': "(Muon_isTrigger[W_l1_idx] && Muon_nloose==1 && Electron_nloose==0)",
             'ELE': "(Electron_isTrigger[W_l1_idx] && Muon_nloose==0 && Electron_nloose==1)"}

ctaus = [0,10,20,50,100,1000]

def isValidFile(fname, tree):
    f = ROOT.TFile.Open(fname)  # ✅ works for root:// and local paths
    if not f or f.IsZombie():
        return False
    try:
        t = f.Get(tree)
        if not t:
            f.Close()
            return False
        t.GetEntry(0)  # test if it’s readable
        _ = t.Photon_pt  # test access
    except Exception:
        f.Close()
        return False
    f.Close()
    return True


# Temp function to use official cut based ID instead of custom ID
# Remove once ntuples get reprocessed with correct ID (not needed after 03_26_24)
def redoPhotonID(plotter,era,ana,isMC = False):
    return
    branches = plotter.rdf.GetColumnNames()
        
    plotter.define("loose_electron", "Electron_pt>15&&abs(Electron_scEta)<2.4 && ((abs(Electron_scEta)>1.479 && abs(Electron_dxy)<0.1 && abs(Electron_dz)<0.2)||(abs(Electron_scEta)<=1.479 && abs(Electron_dxy)<0.05 && abs(Electron_dz)<0.1)) &&Electron_cutBased>1")
    
    plotter.define("tight_electron", "loose_electron&&Electron_cutBased==4")
    
    plotter.define("veto_electron", "Electron_pt>10&&abs(Electron_scEta)<2.4 && ((abs(Electron_scEta)>1.479 && abs(Electron_dxy)<0.1 && abs(Electron_dz)<0.2)||(abs(Electron_scEta)<=1.479 && abs(Electron_dxy)<0.05 && abs(Electron_dz)<0.1)) && (tight_electron==0) && (loose_electron==0)")
    
    plotter.define("Electron_nloose", "Sum(loose_electron)")
    plotter.redefine("Electron_nveto", "Sum(veto_electron)")
    
    if isMC:
        plotter.define("pho_SFs_id", "scaleFactors_2d(Photon_eta, Photon_pt, PHO_ID_{era}_sf, PHO_ID_{era}_binsX, PHO_ID_{era}_binsY, sample_isMC, Photon_cutBased>0)".format(era=era))
        plotter.redefine("Photon_idSF_val", "pho_SFs_id[0]")
        plotter.redefine("Photon_idSF_unc", "pho_SFs_id[1]")
        plotter.define("pho_SFs_pix", "getPixelSeedSF(Photon_isScEtaEB, Photon_isScEtaEE, hasPix_UL{era}_sf, sample_isMC, !Photon_pixelSeed)".format(era=era))
        plotter.redefine("Photon_pixSF_val", "pho_SFs_pix[0]")
        plotter.redefine("Photon_pixSF_unc", "pho_SFs_pix[1]")

def getFiles(query,sampleDir,sampleType,era,prod):
    ser = pd.Series(subprocess.check_output(['xrdfs', 'root://cmseos.fnal.gov', 'ls', f"{sampleDir}/{sampleType}{era}_{prod}/"], text=True).split("\n"))
    matched = ser[ser.str.contains(f'/{query}', regex=True)]
    return list('root://cmseos.fnal.gov/' + matched)

def getPlotters(era,prod,sampleDir,mass_list,modelIndependent=False,which_ana='all'):
    dataEMU = {}     # merged_plotter of data
    dyPlotters = {}  # merged_plotter of dy+jets
    wjPlotters = {}  # merged_plotter of w+jets
    signal = {}      # merged_plotter of signal
    global masses
    masses = mass_list
    #TODO: make option for only running on single ana instead of making everything always
    if which_ana == 'all':
        analyses = ['zmm2g', 'zee2g', 'wen2g', 'wmn2g']
    else:
        analyses = [which_ana]
    for ana in analyses:
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
            redoPhotonID(singleMuPlotters[-1], era,ana)

        for sample in EGSamples:
            if not isValidFile(sample, ana):
                continue
            egPlotters.append(rdf_plotter(sample, tree=ana))
            redoPhotonID(egPlotters[-1], era,ana)
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
                    redoPhotonID(dyPlotters_nJ[-1], '2016preVFP',ana, isMC=True)
                else:
                    dyPlotters_nJ[-1].addCorrectionFactor("0.46", "flat")
                    redoPhotonID(dyPlotters_nJ[-1], '2016postVFP', ana,isMC=True)
            else:
                redoPhotonID(dyPlotters_nJ[-1], era,ana)
            if era == '2018' and ana=='wen2g':
                dyPlotters_nJ[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                dyPlotters_nJ.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                redoPhotonID(dyPlotters_nJ[-1], era, ana,isMC=True)
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
                    redoPhotonID(wjPlotters_nJ[-1], "2016preVFP", ana,isMC=True)
                else:
                    wjPlotters_nJ[-1].addCorrectionFactor("0.46", "flat")
                    redoPhotonID(wjPlotters_nJ[-1], "2016postVFP",ana, isMC=True)
            else:
                redoPhotonID(wjPlotters_nJ[-1], era,ana, isMC=True)
            if era == '2018' and ana=='wen2g':
                wjPlotters_nJ[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                wjPlotters_nJ.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                redoPhotonID(wjPlotters_nJ[-1], era,ana, isMC=True)
                wjPlotters_nJ[-1].addCorrectionFactor(str(38750./59830), "flat")

        wjPlotters[ana] = merged_plotter(wjPlotters_nJ)
        wjPlotters[ana].setFillProperties(1001, ROOT.kAzure-9)
        wjPlotters[ana].setLineProperties(1, ROOT.kAzure-9, 3)

    # Make signal plotters, one merge_plotter for each category
    #TODO: make option for running on all signals or total
    #for sig in ['total']:
    for sig in ['ZH','ggZH','WH','ttH']: #NOTE: if running on 7 GeV samples, there's no ttH contribution
    
        signal[sig] = {}
        if which_ana == 'all':
            analyses = ['zmm2g', 'zee2g', 'wen2g', 'wmn2g']
        else:
            analyses = [which_ana]
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
                                        redoPhotonID(sigPlotters[-1], "2016preVFP",ana, isMC=True)
                                    else:
                                        sigPlotters[-1].addCorrectionFactor("0.46", "flat")
                                        redoPhotonID(sigPlotters[-1], "2016postVFP",ana, isMC=True)
                                else:
                                    redoPhotonID(sigPlotters[-1], era,ana, isMC=True)
                                if era=='2018' and ana=='wen2g':
                                    #Monte Carlo normalized to unaffected luminosity of 21.08fb-1
                                    sigPlotters[-1].addCorrectionFactor(str(21080./59830), "flat")
                                    #Monte carlo normalized to affected luminosity of 38.75fb-1 veto events is affected region defined by HEM cuts
                                    sigPlotters.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                                    #Michalis: THIS PLOTTER NEEDS the weight
                                    sigPlotters[-1].addCorrectionFactor(weight, "flat")

                                    sigPlotters[-1].addCorrectionFactor(str(38750./59830), "flat")
                                    redoPhotonID(sigPlotters[-1], era,ana, isMC=True)

                        signal[sig][ana][m][ct][br] = merged_plotter(sigPlotters)
                        signal[sig][ana][m][ct][br].setFillProperties(0, ROOT.kWhite)
                        signal[sig][ana][m][ct][br].setLineProperties(1, ROOT.kRed, 3)
                    

    plots = {}
    plots['dataEMU'] = dataEMU
    plots['signal'] = signal
    plots['DYJets'] = dyPlotters
    plots['WJets'] = wjPlotters
    return plots


def getSigPlotters(era,prod,sampleDir,mass_list,analysis,modelIndependent=False,combine_signals=False, which_lifetime = 'all'):
    """
    Inputs
    ======
    analysis (STR): which analysis to make plotter for. Options: wmn2g, wen2g, zee2g, zmm2g
    combine_signals (BOOL): whether or not to combine signals (WH,ttH,ZH,ggZH) into 'total' signal
    which lifetime (INT): Which lifetime to make plotter for. If no argument is added, will create plotter with all lifetimes.

    Returns
    =======
    Dict corresponding to plotters. 
    Ex:plotter['signal']['total']['wmn2g'][m][ct]]
    """
    signal = {}      # merged_plotter of signal
    global masses
    masses = mass_list
    # Make signal plotters, one merge_plotter for each category
    if combine_signals:
        signals = ['total']
    else:
        signals = ['ZH','ggZH','WH','ttH']  #NOTE: if running on 7 GeV samples, there's no ttH contribution
    for sig in signals:
        signal[sig] = {}
        for ana in [analysis]:
            signal[sig][ana]  = {}        
            for m in masses:
                print(f"\nSignal samples used")
                signal[sig][ana][m] = {}
                if which_lifetime == 'all':
                    ctaus = [0,10,20,50,100,1000]
                else:
                    ctaus = [which_lifetime]
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
                            #samples = glob.glob(f"{sampleDir}/MC{era}_{prod}/{v}H{br}_M{m}_ctau{ct}_UL{era[-2:]}*.root")
                            samples = getFiles(f"{v}H{br}_M{m}_ctau{ct}_UL{era[-2:]}",sampleDir,sampleType="MC",era=era,prod=prod)
                            for sample in samples:
                                print(sample)
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
                                        redoPhotonID(sigPlotters[-1], "2016preVFP",ana, isMC=True)
                                    else:
                                        redoPhotonID(sigPlotters[-1], "2016postVFP",ana, isMC=True)
                                else:
                                    redoPhotonID(sigPlotters[-1], era, ana, isMC=True)
                                if era=='2018' and ana=='wen2g':
                                    #Monte Carlo normalized to unaffected luminosity of 21.08fb-1
                                    sigPlotters[-1].addCorrectionFactor(str(21080./59830), "flat")
                                    #Monte carlo normalized to affected luminosity of 38.75fb-1 veto events is affected region defined by HEM cuts
                                    sigPlotters.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                                    #Michalis: THIS PLOTTER NEEDS the weight
                                    sigPlotters[-1].addCorrectionFactor(weight, "flat")

                                    sigPlotters[-1].addCorrectionFactor(str(38750./59830), "flat")
                                    redoPhotonID(sigPlotters[-1], era,ana, isMC=True)

                        signal[sig][ana][m][ct][br] = merged_plotter(sigPlotters)
                        signal[sig][ana][m][ct][br].setFillProperties(0, ROOT.kWhite)
                        signal[sig][ana][m][ct][br].setLineProperties(1, ROOT.kRed, 3)
    plots = {}
    plots['signal'] = signal
    return plots

def getDYPlotters(era,prod,sampleDir,mass_list:list[int],analysis,modelIndependent:bool=False,which_mc:list = ['DYJetsToLL']):  
    dyPlotters = {}  # merged_plotter of w+jets
    global masses #TODO clean this mess up
    masses = mass_list
    for ana in [analysis]:
        Samples = []
        for sample_mc in which_mc:
            Samples+= getFiles(sample_mc,sampleDir,sampleType="MC",era=era,prod=prod)
        print('\nDY samples used:')
        for sample in Samples:
            print(sample)
        print("\n")
        dyPlotters_nJ = []
        for sample in Samples:
            if not isValidFile(sample, ana):
                continue
            dyPlotters_nJ.append(rdf_plotter(sample, True, tree = ana))
            if era == '2016':
                if 'APV' in sample or 'preVFP' in sample:
                    redoPhotonID(dyPlotters_nJ[-1], "2016preVFP",ana, isMC=True)
                else:
                    redoPhotonID(dyPlotters_nJ[-1], "2016postVFP", ana,isMC=True)
            else:
                redoPhotonID(dyPlotters_nJ[-1], era,ana, isMC=True)
            if era == '2018' and ana=='wen2g':
                dyPlotters_nJ[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                dyPlotters_nJ.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                redoPhotonID(dyPlotters_nJ[-1], era, ana,isMC=True)
                dyPlotters_nJ[-1].addCorrectionFactor(str(38750./59830), "flat")
        dyPlotters[ana] = merged_plotter(dyPlotters_nJ)
        dyPlotters[ana].setFillProperties(1001, ROOT.kAzure-9)
        dyPlotters[ana].setLineProperties(1, ROOT.kAzure-9, 3)
    plots = {}
    plots['DY'] = dyPlotters
    return plots

def getWJPlotters(era,prod,sampleDir,mass_list:list[int],analysis,modelIndependent:bool=False,which_mc:list = ['WGToLNuG','TTJets','WJetsToLNu','WGG']):
    wjPlotters = {}  # merged_plotter of w+jets
    global masses #TODO clean this mess up
    masses = mass_list
    for ana in [analysis]:
        Samples = []
        for sample_mc in which_mc:
            Samples+= getFiles(sample_mc,sampleDir,sampleType="MC",era=era,prod=prod)
        print('\nWJ samples used:')
        for sample in Samples:
            print(sample)
        print("\n")
        wjPlotters_nJ = []
        for sample in Samples:
            if not isValidFile(sample, ana):
                continue
            wjPlotters_nJ.append(rdf_plotter(sample, True, tree = ana))
            if era == '2016':
                if 'APV' in sample or 'preVFP' in sample:
                    redoPhotonID(wjPlotters_nJ[-1], "2016preVFP",ana, isMC=True)
                else:
                    redoPhotonID(wjPlotters_nJ[-1], "2016postVFP", ana,isMC=True)
            else:
                redoPhotonID(wjPlotters_nJ[-1], era,ana, isMC=True)
            if era == '2018' and ana=='wen2g':
                wjPlotters_nJ[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                wjPlotters_nJ.append(rdf_plotter(sample, True, tree = ana, defaultCuts = cutsHEM))
                redoPhotonID(wjPlotters_nJ[-1], era, ana,isMC=True)
                wjPlotters_nJ[-1].addCorrectionFactor(str(38750./59830), "flat")
        wjPlotters[ana] = merged_plotter(wjPlotters_nJ)
        wjPlotters[ana].setFillProperties(1001, ROOT.kAzure-9)
        wjPlotters[ana].setLineProperties(1, ROOT.kAzure-9, 3)
    plots = {}
    plots['WJets'] = wjPlotters
    return plots

def getDataPlotters(era:str, prod:str, sampleDir: str,mass_list:list[int], modelIndependent:bool = False, analysis:str="ee2g") -> merged_plotter:

    Samples = getFiles("JetHT",sampleDir,sampleType="DATA",era=era,prod=prod)
    
    Plotters = []
    
    print("\nJetHT samples")
    for sample in Samples:
        print(sample)
        if not isValidFile(sample, analysis):
            continue
        Plotters.append(rdf_plotter(sample, tree=analysis, report = "Report_" + analysis))
        
    dataEMU = merged_plotter(Plotters) 

    return dataEMU

def getPlotter(sample,sampleDir,sampleType,eras,prod,analysis)->merged_plotter:
    plotters=[]
    for era in eras:
        #if we are doing data overrtide the search
        files=[]
        if sampleType=='DATA':
            files = getFiles('JetHT',sampleDir,sampleType,era,prod)
        else: #MC means query
            files = getFiles(sample,sampleDir,sampleType,era,prod)
        print("Files Used")
        for f in files:
            print(f)
            plotters.append(rdf_plotter(f,isMC=(sampleType=='MC'), tree=analysis, report = "Report_" + analysis))
            #scale with the luminosity
            if sampleType=='MC':
                plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                plotters[-1].addCorrectionFactor('1000', "flat") #to convert to pb-1
    p = merged_plotter(plotters)
    p.define("Photon_SieiePassLevel","passLevelSieie(Photon_vidNestedWPBitmap)")
    p.define("Photon_PhIsoPassLevel","passLevelPhIso(Photon_vidNestedWPBitmap)")
    p.define("Photon_HOEPassLevel","passLevelHOE(Photon_vidNestedWPBitmap)")
    p.define("Photon_NeuIsoPassLevel","passLevelNeuIso(Photon_vidNestedWPBitmap)")
    p.define("Photon_ChIsoPassLevel","passLevelChIso(Photon_vidNestedWPBitmap)")

    return p


def calculate_fake_rate(sampleDir,prod,eras=['2016','2017','2018'],ptbins=[25.,30.,35.,40.,50.,60.,80.,150.],etabins=[-2.5,-2.0,-1.57,-1.44,-0.8,0.8,1.44,1.57,2.0,2.5],arrayName="fake_rate_VHF",doMCClosure=False,outdir='.'):
    if len(eras) == 3:
        era = "Run2"
    elif len(eras) == 1:
        era = eras[0]
    else:
        print("what the heck")
        exit()
    def clopper_pearson(k, n, confidence=0.68):
        """
        Compute the Clopper-Pearson interval for a binomial proportion.
        
        k: number of successes (array-like)
        n: number of trials (array-like)
        confidence: confidence level (0.95 for 95%)
        """
        alpha = 1 - confidence
    
        # Lower bound: alpha/2 quantile of Beta(k, n-k+1)
        # We use np.maximum to handle the edge case where k=0
        lower = beta.ppf(alpha / 2, k, n - k + 1)
        lower = np.nan_to_num(lower, nan=0.0)
        
        # Upper bound: 1-alpha/2 quantile of Beta(k+1, n-k)
        # We use np.nan_to_num to handle the edge case where k=n
        upper = beta.ppf(1 - alpha / 2, k + 1, n - k)
        upper = np.nan_to_num(upper, nan=1.0)
    
        return lower, upper
    if doMCClosure:
        wjets=getPlotter('WJetsToLNu_',sampleDir,'MC',eras,prod,'gamma')
        vjets=getPlotter('DYJetsToLL_',sampleDir,'MC',eras,prod,'gamma')
        fr = merged_plotter([wjets,vjets])
    else:    
        fr = getPlotter('nothing',sampleDir,'DATA',eras,prod,'gamma')
    
    fr.define("denominator_mask","Photon_preselection==1")
    fr.define("numerator_mask","Photon_preselection ==1 && Photon_cutBased>0")

    fr.define("Photon_pt_denom","Photon_pt[denominator_mask]")
    fr.define("Photon_eta_denom","Photon_eta[denominator_mask]")
    fr.define("Photon_pt_num","Photon_pt[numerator_mask]")
    fr.define("Photon_eta_num","Photon_eta[numerator_mask]")

    cuts_denom = 'Sum(Photon_preselection==1)==1'
    cuts_num = '&&'.join([cuts_denom,'Sum(Photon_cutBased>0)>0'])

    xedges,yedges,denominator,w2_denom = fr.array2d('Photon_pt_denom','Photon_eta_denom',cuts_denom,('denom','denom',len(ptbins)-1,np.array(ptbins),len(etabins)-1,np.array(etabins)))
    xedges,yedges,numerator,w2_num = fr.array2d('Photon_pt_num','Photon_eta_num',cuts_num,('num','num',len(ptbins)-1,np.array(ptbins),len(etabins)-1,np.array(etabins)))


    #fr.display(['Photon_phi','Photon_eta','denominator_mask','numerator_mask'])
    #cuts_denom = 'Photon_preselection>0 && Photon_electronIdx != -1'
    #cuts_num = '&&'.join([cuts_denom,'Photon_cutBased>0']) 

    result = np.zeros_like(numerator)
    fake_rate = np.divide(numerator,denominator,out=result,where=(denominator !=0))
    fake_rate_down,fake_rate_up = clopper_pearson(numerator,denominator)

    #make a plot
    fig, ax = plt.subplots()
    mesh = ax.pcolormesh(xedges, yedges, fake_rate, cmap='viridis', edgecolors='white', linewidth=0.5)
    plt.colorbar(mesh, label='fake rate')
    # 4. Add labels and uncertainties (The "ROOT TEXT" part)
    # Calculate centers for text placement
    x_centers = (xedges[:-1] + xedges[1:]) / 2
    y_centers = (yedges[:-1] + yedges[1:]) / 2
    #put the values on the plot but also make a string
    st=f'std::vector<std::vector<std::vector<float> > > {arrayName}_vals = {{'
    for i, x_val in enumerate(x_centers):
        st=st+'{'
        for j, y_val in enumerate(y_centers):
            st=st+'{'+','.join([str(fake_rate[j, i]),str(fake_rate_up[j, i]),str(fake_rate_down[j, i])])+'}'
            if j!=(len(y_centers)-1):
                st=st+','
            val = fake_rate[j, i] # Note: pcolormesh uses (row, col) which is (y, x)
            errUp = fake_rate_up[j, i]-val
            errDwn = val-fake_rate_down[j, i]                    
#            label = f"{val:.3f}+\n{errUp:.3f}\n-{errDwn:.3f}"
#            ax.text(x_val, y_val, label, color="black", ha="center", va="center", fontsize=9)
            ax.set_xlabel(r"$\gamma p_{T}$ (GeV)")
            ax.set_ylabel(r"$\gamma \eta$")
        st=st+'}'
        if i!=(len(x_centers)-1):
            st=st+','
    st=st+'};'
    if doMCClosure:
        mh.cms.label(data=False, lumi=lumifb[era],ax=ax, loc=0)
    else:
        mh.cms.label("Preliminary", data=True, lumi=lumifb[era], ax=ax, loc=0)
    #print it in C++ format
    #note that we remove the last edge on how the code is defined to work
    st=st+'\n'+f'std::vector<float> {arrayName}_xbins = {{'+','.join([str(x) for x in xedges[:-1]])+'};\n'+f'std::vector<float> {arrayName}_ybins = {{'+','.join([str(y) for y in yedges[:-1]])+'};\n'
    plt.savefig(f'{outdir}/{arrayName}.png', dpi=400, bbox_inches='tight')
    return st


def getSF(scaleFactors, sfToVary = ''):
    if sfToVary == '':
        return "(" + ")*(".join(scaleFactors) + ")"
    toVary = list(filter(lambda x: sfToVary in x, scaleFactors))
    scaleFactors = list(filter(lambda x: sfToVary not in x, scaleFactors))
    sfUp = ["{}+{}".format(sf, sf.replace("val", "unc")) for sf in toVary]
    sfDown = ["{}-{}".format(sf, sf.replace("val", "unc")) for sf in toVary]
    return ["(" + ")*(".join(scaleFactors + sfUp) + ")", "(" + ")*(".join(scaleFactors + sfDown) + ")"]

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


#debugging code
if __name__ == '__main__':
    try:
        era = '2018'
        prod = '01_23_26'
        sampleDir = "/store/user/gfavila/DDP" #EOS directory
        m = 20
        test_plotters = False
        date = datetime.date.today()
        date = str(date).split('-')
        date = date[1]+date[2]+date[0][2:]
        outdir = f"VHF_{date}"
        if not os.path.isdir(outdir):
            os.mkdir(outdir)

        if test_plotters:

            print('testing data plotters')
            data_plotters = getDataPlotters(era,prod,sampleDir,mass_list = [20],analysis="gamma")
            list_of_plotters = data_plotters.plotters
            
            print("\nIndividual RDFs")
            for plotter in list_of_plotters:
                print(plotter)
            
            print("\nDefault cuts")
            for plotter in list_of_plotters:
                print(plotter.defaultCuts)
            
            print("\nReports")
            for plotter in list_of_plotters:
                report:dict[str] = plotter.readReport()
                for key, value in report.items():
                    print(f"{key}: {value}")
                print("\n")

        print("Running Calculation of Fake Rates")
        eras=["2016","2017","2018"]
        #eras=["2018"]
        with open(outdir+'/vhfFakeRates.h', "w") as file:
            file.write("#ifndef FAKERATES\n")
            file.write("#define FAKERATES\n")
            for e in eras:
                print("running for era ",e)
                fr=calculate_fake_rate(sampleDir,prod,[e],arrayName=f"fake_rate_{e}",outdir=outdir)
                file.write(fr)
            #write the average all run fake rate for studies
            print("running for Run2")
            fr=calculate_fake_rate(sampleDir,prod,eras,arrayName="fake_rate",outdir=outdir)
            file.write(fr)
            #print("running for MC")
            ##write the average all run MC fake rate for studies
            #fr=calculate_fake_rate(sampleDir,prod,eras,arrayName="fake_rate_MC",doMCClosure=True,outdir=outdir)
            #file.write(fr)
            file.write("#endif\n")
    
    except Exception as e:
        import traceback
        traceback.print_exc()
    finally:
        ROOT.gSystem.Exit(0)
