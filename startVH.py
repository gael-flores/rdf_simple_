from common.plotter import *
from python.VHTools.VHcuts import *
from python.VHTools.config import *
import ROOT
import glob
import array
import pandas as pd
import subprocess

ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')



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
    branches = plotter.rdf.GetColumnNames()
    if "Photon_passPhIso" not in branches:
        plotter.define("Photon_passPhIso" , "passPhIso(Photon_vidNestedWPBitmap)")
    for m in [15,20,30,40,50,55]:
        plotter.define(f'DeltaR_g1g2_m{m}', f'DeltaR(Photon_eta[best_2g_idx1_m{m}],Photon_eta[best_2g_idx2_m{m}], Photon_phi[best_2g_idx1_m{m}] , Photon_phi[best_2g_idx2_m{m}])')
        plotter.define(f'DeltaPhi_g1g2_m{m}', f'abs(DeltaPhi(Photon_phi[best_2g_idx1_m{m}] , Photon_phi[best_2g_idx2_m{m}]))')
        if "Jet_pt" in branches and "Jet_eta" in branches and "Jet_phi" in branches:
            plotter.define(f'best_2g_g1_closestJet_Idx_m{m}', f"idx_closest_jet(Photon_eta[best_2g_idx1_m{m}],Photon_phi[best_2g_idx1_m{m}],Jet_pt,Jet_eta,Jet_phi)[0]")
            plotter.define(f'best_2g_g2_closestJet_Idx_m{m}', f"idx_closest_jet(Photon_eta[best_2g_idx2_m{m}],Photon_phi[best_2g_idx2_m{m}],Jet_pt,Jet_eta,Jet_phi)[0]")
            plotter.define(f'best_2g_g1_DeltaR_closestJet_m{m}',f"DeltaR(Photon_eta[best_2g_idx1_m{m}],Jet_eta[best_2g_g1_closestJet_Idx_m{m}],Photon_phi[best_2g_idx1_m{m}],Jet_phi[best_2g_g1_closestJet_Idx_m{m}])")
            plotter.define(f'best_2g_g2_DeltaR_closestJet_m{m}',f"DeltaR(Photon_eta[best_2g_idx2_m{m}],Jet_eta[best_2g_g2_closestJet_Idx_m{m}],Photon_phi[best_2g_idx2_m{m}],Jet_phi[best_2g_g2_closestJet_Idx_m{m}])")
            plotter.define(f"Jet1_puId_m{m}",f"Jet_puId[best_2g_g1_closestJet_Idx_m{m}]")
            plotter.define(f"Jet2_puId_m{m}",f"Jet_puId[best_2g_g2_closestJet_Idx_m{m}]")
            plotter.define(f'Jetfrac1_m{m}',f'Jet_chFPV0EF[best_2g_g1_closestJet_Idx_m{m}]/(Jet_chHEF[best_2g_g1_closestJet_Idx_m{m}] + Jet_chEmEF[best_2g_g1_closestJet_Idx_m{m}] + Jet_chFPV0EF[best_2g_g1_closestJet_Idx_m{m}])')
            plotter.define(f'Jetfrac2_m{m}',f'Jet_chFPV0EF[best_2g_g2_closestJet_Idx_m{m}]/(Jet_chHEF[best_2g_g2_closestJet_Idx_m{m}] + Jet_chEmEF[best_2g_g2_closestJet_Idx_m{m}] + Jet_chFPV0EF[best_2g_g2_closestJet_Idx_m{m}])')
        if "Photon_passPhIso" not in branches:
            print('what the heck')
            plotter.define(f"best_2g_looseID1_m{m}", f"Photon_cutBased[best_2g_idx1_m{m}]>0")
            plotter.define(f"best_2g_looseID2_m{m}", f"Photon_cutBased[best_2g_idx2_m{m}]>0")
            plotter.redefine(f"best_2g_sumID_m{m}", f"best_2g_looseID1_m{m}+best_2g_looseID2_m{m}")
            plotter.filter(f"Photon_passPhIso[best_2g_idx1_m{m}]==1 && Photon_passPhIso[best_2g_idx2_m{m}]==1") #all photons required to pass isolation
        if "Tau_pt" in branches:
            plotter.define(f'best_2g_g1_closestTau_Idx_m{m}', f"idx_closest_tau(Photon_eta[best_2g_idx1_m{m}],Photon_phi[best_2g_idx1_m{m}],Tau_pt,Tau_eta,Tau_phi)[0]")
            plotter.define(f'best_2g_g2_closestTau_Idx_m{m}', f"idx_closest_tau(Photon_eta[best_2g_idx2_m{m}],Photon_phi[best_2g_idx2_m{m}],Tau_pt,Tau_eta,Tau_phi)[0]")
            plotter.define(f'best_2g_g1_DeltaR_closestTau_m{m}',f"DeltaR(Photon_eta[best_2g_idx1_m{m}],Tau_eta[best_2g_g1_closestTau_Idx_m{m}],Photon_phi[best_2g_idx1_m{m}],Tau_phi[best_2g_g1_closestTau_Idx_m{m}])")
            plotter.define(f'best_2g_g2_DeltaR_closestTau_m{m}',f"DeltaR(Photon_eta[best_2g_idx2_m{m}],Tau_eta[best_2g_g2_closestTau_Idx_m{m}],Photon_phi[best_2g_idx2_m{m}],Tau_phi[best_2g_g2_closestTau_Idx_m{m}])")
            plotter.define(f"Tau_mass1_m{m}",f"Tau_mass[best_2g_g1_closestTau_Idx_m{m}]")
            plotter.define(f"Tau_mass2_m{m}",f"Tau_mass[best_2g_g2_closestTau_Idx_m{m}]")
            plotter.define(f"Tau_decayMode1_m{m}",f"Tau_decayMode[best_2g_g1_closestTau_Idx_m{m}]")
            plotter.define(f"Tau_decayMode2_m{m}",f"Tau_decayMode[best_2g_g2_closestTau_Idx_m{m}]")
    if ana =='wen2g':
        plotter.define('ept','Electron_pt[W_l1_idx]')
        plotter.define('mInv_egamma','calculate_lgamma_mass(Electron_pt[W_l1_idx], Electron_eta[W_l1_idx], Electron_phi[W_l1_idx], Electron_mass[W_l1_idx],Photon_pt, Photon_eta, Photon_phi)')
        plotter.define('mInv_egammagamma','calculate_lgammagamma_mass(Electron_pt[W_l1_idx], Electron_eta[W_l1_idx], Electron_phi[W_l1_idx], Muon_mass[W_l1_idx],Photon_pt[best_2g_idx1_m20], Photon_eta[best_2g_idx1_m20], Photon_phi[best_2g_idx1_m20],Photon_pt[best_2g_idx2_m20], Photon_eta[best_2g_idx2_m20], Photon_phi[best_2g_idx2_m20])')
        plotter.define('DeltaR_closest_photon', 'deltaR_lgamma(Electron_eta[W_l1_idx], Electron_phi[W_l1_idx],Photon_eta[Photon_preselection], Photon_phi[Photon_preselection])[0]')
        #plotter.define("loose_electron", "Electron_pt>15&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
        #plotter.define("tight_electron", "loose_electron&&Electron_cutBased>3")
        #plotter.define("veto_electron", "Electron_pt>5&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&Electron_lostHits<2&&Electron_convVeto&&(tight_electron==0)&&(loose_electron==0)")
        #plotter.redefine("Electron_nveto", "Sum(veto_electron)")
    elif ana =='wmn2g':
        plotter.define('mpt','Muon_pt[W_l1_idx]')
        plotter.define('mInv_mugamma','calculate_lgamma_mass(Muon_pt[W_l1_idx], Muon_eta[W_l1_idx], Muon_phi[W_l1_idx], Muon_mass[W_l1_idx],Photon_pt, Photon_eta, Photon_phi)')
        plotter.define('mInv_mugammagamma','calculate_lgammagamma_mass(Muon_pt[W_l1_idx], Muon_eta[W_l1_idx], Muon_phi[W_l1_idx], Muon_mass[W_l1_idx],Photon_pt[best_2g_idx1_m20], Photon_eta[best_2g_idx1_m20], Photon_phi[best_2g_idx1_m20],Photon_pt[best_2g_idx2_m20], Photon_eta[best_2g_idx2_m20], Photon_phi[best_2g_idx2_m20])')
        plotter.define('DeltaR_closest_photon', 'deltaR_lgamma(Muon_eta[W_l1_idx], Muon_phi[W_l1_idx],Photon_eta[Photon_preselection], Photon_phi[Photon_preselection])[0]')
        #plotter.redefine("veto_muon", "Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&(loose_muon==0)&&(tight_muon==0)")
        #plotter.redefine("Muon_nveto", "Sum(veto_muon)")
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
                                        sigPlotters[-1].addCorrectionFactor("0.54", "flat")
                                        redoPhotonID(sigPlotters[-1], "2016preVFP",ana, isMC=True)
                                    else:
                                        sigPlotters[-1].addCorrectionFactor("0.46", "flat")
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
                    dyPlotters_nJ[-1].addCorrectionFactor("0.54", "flat")
                    redoPhotonID(dyPlotters_nJ[-1], "2016preVFP",ana, isMC=True)
                else:
                    dyPlotters_nJ[-1].addCorrectionFactor("0.46", "flat")
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
                    wjPlotters_nJ[-1].addCorrectionFactor("0.54", "flat")
                    redoPhotonID(wjPlotters_nJ[-1], "2016preVFP",ana, isMC=True)
                else:
                    wjPlotters_nJ[-1].addCorrectionFactor("0.46", "flat")
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

def getDataPlotters(era:str, prod:str, sampleDir: str,mass_list:list[int],analysis:str ,modelIndependent:bool = False) -> merged_plotter:
    if 'e' in analysis:
        if era != '2018':
            EGSamples = getFiles("SingleElectron_",sampleDir,sampleType="DATA",era=era,prod=prod)
        else:
            EGSamples = getFiles("EGamma_",sampleDir,sampleType="DATA",era=era,prod=prod)
        
        egPlotters = []
        
        print("\nEG samples")
        for sample in EGSamples:
            print(sample)
            if not isValidFile(sample, analysis):
                continue
            egPlotters.append(rdf_plotter(sample, tree=analysis, report = "Report_" + analysis))
            redoPhotonID(egPlotters[-1], era,analysis)
            
            if era == '2018' and analysis == 'wen2g':
                egPlotters[-1].defaultCuts = "((run>=319077&&(Electron_eta[W_l1_idx]>-1.3||Electron_eta[W_l1_idx]<-3.0)&&(Electron_phi[W_l1_idx]>-0.87||Electron_phi[W_l1_idx]<-1.57))||(run<319077))"
        dataEMU = merged_plotter(egPlotters)            
    
    else:
        singleMuSamples = getFiles(f"SingleMuon_",sampleDir,sampleType="DATA",era=era,prod=prod)
        singleMuPlotters = []
    
        print("\nSingle muon samples:")
        for sample in singleMuSamples:
            print(sample)
            if not isValidFile(sample, analysis):
                continue
            singleMuPlotters.append(rdf_plotter(sample, tree=analysis, report = "Report_" + analysis))
            redoPhotonID(singleMuPlotters[-1], era, analysis)    
        
        dataEMU = merged_plotter(singleMuPlotters)

    return dataEMU

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
        era = '2016'
        prod = '11_14_25'
        sampleDir = "/store/user/gfavila/DDP" #EOS directory
        m = 20
        v = 'Z'
        l = 'MU'
        print('testing data plotters')
        data_plotters = getDataPlotters(era,prod,sampleDir,ana[v][l])
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
            

        ROOT.gSystem.Exit(0)
    except Exception as e:
        import traceback
        traceback.print_exc()
        ROOT.gSystem.Exit(0)
