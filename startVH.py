import ROOT
import glob
import array
import pandas as pd
import subprocess
from scipy.stats import beta
import matplotlib
matplotlib.use('TkAgg')
from common.plotter import *
from common.datacard_maker import cnc_datacard_maker
import gc

ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')

ROOT.gInterpreter.Declare('#include "common/vhFakeRates.h"')        
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')        
        
ROOT.ROOT.EnableImplicitMT()

analysis_status='Preliminary'

masses = [15,20,30,40,50,55]
lifetimes = [0,10,20,50,100,1000]
signal_colors={0:'red',
               10:'lightcoral',
               20:'chocolate',
               50:'darkorange',
               100:'yellowgreen',
               1000:'lightseagreen'}
analyses  = ['wmn2g','wen2g','zmm2g','zee2g']
signals=['ZH','ggZH','WH','ttH']

lumifb = {'2018': "59.83",
          '2017': '41.48',
          '2016': '36.31',
          'Run2': '137.62'}

center_of_mass = {'2018': 13,
                  '2017': 13,
                  '2016': 13,
                  'Run2': 13}


# Standard Analysis Cuts
cuts = {}
# Tight + veto lepton cuts
cuts['W'] = {'MU': "(Muon_isTrigger[W_l1_idx] && Muon_nloose==1 && Electron_nloose==0)",
             'ELE': "(Electron_isTrigger[W_l1_idx] && Muon_nloose==0 && Electron_nloose==1)"}
            
cuts['Z'] = {'MU': "((Muon_isTrigger[Z_idx[0]] || Muon_isTrigger[Z_idx[1]]) && Muon_nloose==2 && Electron_nloose==0)",
             'ELE': "((Electron_isTrigger[Z_idx[0]] || Electron_isTrigger[Z_idx[1]]) && Electron_nloose==2 && Muon_nloose==0)"}
# Photon pt and kinematic cuts
cuts['pt'] = {}
cuts['photons'] = {}
for m in masses:
    cuts['pt'][m] = "(Photon_pt[best_2g_idx1_m{m}] > 35 && Photon_pt[best_2g_idx2_m{m}] > 25)".format(m=m)
    cuts['photons'][m] = "(best_2g_pt_m{m} > 20 && best_2g_raw_mass_m{m} > 4)".format(m=m)

# w->e+nu misID cuts
cuts['misID'] = {}
for v in ['W', 'Z']:
    cuts['misID'][v] = {}
    for l in ['ELE', 'MU']:
        cuts['misID'][v][l] = {}
        for m in masses:
            if (v == 'W' and l == 'ELE'):
                cuts['misID'][v][l][m] = "(abs(best_2g_misID1_m{m}-91) > 10 && abs(best_2g_misID2_m{m}-91) > 10 && abs(best_2g_misID3_m{m} - 91) > 15)".format(m=m)
            else:
                cuts['misID'][v][l][m] = "(1)"
# Z->ll FSR cuts
#cuts['fsr'] = {}
#for v in ['W', 'Z']:
#    cuts['fsr'][v] = {}
#    for m in masses:
#        cuts['fsr'][v][m] = "(1)"
#        if (v=='W'):
#            cuts['fsr'][v][m] = "(1)"
#        else:
#            cuts['fsr'][v][m] = "(Photon_isFSR[best_2g_idx1_m{m}]==0 && Photon_isFSR[best_2g_idx2_m{m}]==0)".format(m=m)

            
#redefine the cuts here based on the analysis for easier parsing
for analysis in ['wmn2g','wen2g','zmm2g','zee2g']:
    cuts[analysis]={}
    if analysis=='wmn2g':
        v='W'
        l='MU'
    elif analysis=='zmm2g':
        v='Z'
        l='MU'
    elif analysis=='wen2g':
        v='W'
        l='ELE'
    elif analysis=='zee2g':
        v='Z'
        l='ELE'
    for m in masses:
        cuts[analysis][m] = {} 
        cuts[analysis][m]['presr']='&&'.join([cuts[v][l],
                                           cuts['pt'][m],
                                           cuts['photons'][m],
                                           cuts['misID'][v][l][m],
                                           f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)"])
        cuts[analysis][m]['precr']='&&'.join([cuts[v][l],
                                              cuts['pt'][m],
                                              cuts['photons'][m],
                                              cuts['misID'][v][l][m],
                                              f"(Photon_passCutBasedID[best_2g_idx1_m{m}]>0)"])                                              
        cuts[analysis][m]['precr_abcd']='&&'.join([cuts[v][l],
                                              cuts['pt'][m],
                                              cuts['photons'][m],
                                              cuts['misID'][v][l][m],
                                              f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==1)"])


        cuts[analysis][m]['sr']='&&'.join([cuts[analysis][m]['presr'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['cr']='&&'.join([cuts[analysis][m]['precr'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['cr_abcd']='&&'.join([cuts[analysis][m]['precr_abcd'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['ssb']='&&'.join([cuts[analysis][m]['presr'],f'(best_2g_raw_mass_m{m}>62.5)'])
        cuts[analysis][m]['csb_abcd']='&&'.join([cuts[analysis][m]['precr_abcd'],f'(best_2g_raw_mass_m{m}>62.5)'])                  
# Cuts for signal efficiencies:
effCuts = {}
effCuts['HLT'] = {'ELE': {'2018': "(HLT_Ele32_WPTight_Gsf)",
                          '2017': "(HLT_Ele32_WPTight_Gsf_L1DoubleEG)",
                          '2016': "(HLT_Ele27_WPTight_Gsf)"},
                  'MU': {'2018': "(HLT_IsoMu24)",
                         '2017': "(HLT_IsoMu27)",
                         '2016': "(HLT_IsoMu24)"}
              }

# Cuts for W->e+nu 2018 HEM
cutsHEM = "(Electron_eta[W_l1_idx]>-1.3 || Electron_eta[W_l1_idx]<-3.0) && (Electron_phi[W_l1_idx]>-0.85 || Electron_phi[W_l1_idx]<-1.57)"








binning={'wen2g': {7: [((4.0, 6.0), [-10.0, 40.0, 70.0, 110.0]), ((6.0, 8.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((8.0, 10.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((10.0, 12.0), [-10.0, 20.0, 110.0])], 15: [((4.0, 8.0), [-10.0, 40.0, 70.0, 110.0]), ((8.0, 12.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((12.0, 16.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((16.0, 20.0), [-10.0, 20.0, 110.0])], 20: [((4.0, 9.25), [-10.0, 40.0, 70.0, 110.0]), ((9.25, 14.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((14.5, 19.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.75, 25.0), [-10.0, 20.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 40.0, 70.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 20.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 70.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 70.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}, 'wmn2g': {7: [((4.0, 6.0), [-10.0, 40.0, 70.0, 110.0]), ((6.0, 8.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((8.0, 10.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((10.0, 12.0), [-10.0, 20.0, 110.0])], 15: [((4.0, 8.0), [-10.0, 40.0, 70.0, 110.0]), ((8.0, 12.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((12.0, 16.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((16.0, 20.0), [-10.0, 20.0, 110.0])], 20: [((4.0, 9.25), [-10.0, 40.0, 70.0, 110.0]), ((9.25, 14.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((14.5, 19.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.75, 25.0), [-10.0, 20.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 15.0, 70.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 70.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 20.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 70.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 70.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}, 'zee2g': {7: [((4.0, 8.0), [-10.0, 26.0, 62.0, 110.0]), ((8.0, 12.0), [-10.0, 110.0])], 15: [((4.0, 12.0), [-10.0, 26.0, 62.0, 110.0]), ((12.0, 20.0), [-10.0, 110.0])], 20: [((4.0, 14.5), [-10.0, 110.0]), ((14.5, 25.0), [-10.0, 18.0, 45.0, 72.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 40.0, 70.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 10.0, 30.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}, 'zmm2g': {7: [((4.0, 6.0), [-10.0, 40.0, 70.0, 110.0]), ((6.0, 8.0), [-10.0, 45.0, 72.0, 110.0]), ((8.0, 10.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((10.0, 12.0), [-10.0, 20.0, 110.0])], 15: [((4.0, 8.0), [-10.0, 40.0, 70.0, 110.0]), ((8.0, 12.0), [-10.0, 45.0, 72.0, 110.0]), ((12.0, 16.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((16.0, 20.0), [-10.0, 20.0, 110.0])], 20: [((4.0, 9.25), [-10.0, 40.0, 70.0, 110.0]), ((9.25, 14.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((14.5, 19.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.75, 25.0), [-10.0, 20.0, 110.0])], 30: [((4.0, 11.75), [-10.0, 40.0, 110.0]), ((11.75, 19.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((19.5, 27.25), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((27.25, 35.0), [-10.0, 20.0, 110.0])], 40: [((4.0, 14.25), [-10.0, 40.0, 70.0, 110.0]), ((14.25, 24.5), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((24.5, 34.75), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((34.75, 45.0), [-10.0, 10.0, 30.0, 110.0])], 50: [((4.0, 16.75), [-10.0, 40.0, 70.0, 110.0]), ((16.75, 29.5), [-10.0, 18.0, 72.0, 110.0]), ((29.5, 42.25), [-10.0, 45.0, 72.0, 110.0]), ((42.25, 55.0), [-10.0, 20.0, 110.0])], 55: [((4.0, 18.0), [-10.0, 40.0, 70.0, 110.0]), ((18.0, 32.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((32.0, 46.0), [-10.0, 18.0, 45.0, 72.0, 110.0]), ((46.0, 60.0), [-10.0, 20.0, 110.0])]}}


# Scale factor stuff 
leptonSF = {}
leptonSF['wmn2g'] = "*".join(['Muon_recoSF_val[W_l1_idx]', 'Muon_idSF_val[W_l1_idx]', 'Muon_isoSF_val[W_l1_idx]', 'Muon_trigSF_val[W_l1_idx]'])
leptonSF['wen2g'] = "*".join(['Electron_recoSF_val[W_l1_idx]', 'Electron_idSF_val[W_l1_idx]', 'Electron_trigSF_val[W_l1_idx]'])
leptonSF['zee2g'] = "*".join(['Electron_recoSF_val[Z_idx[0]]', 'Electron_recoSF_val[Z_idx[1]]', 'Electron_idSF_val[Z_idx[0]]', 'Electron_idSF_val[Z_idx[1]]', 'Electron_trigSF_val[Z_idx[0]]'])
leptonSF['zmm2g'] = "*".join(['Muon_recoSF_val[Z_idx[0]]', 'Muon_recoSF_val[Z_idx[1]]', 'Muon_idSF_val[Z_idx[0]]', 'Muon_idSF_val[Z_idx[1]]', 'Muon_isoSF_val[Z_idx[0]]', 'Muon_isoSF_val[Z_idx[1]]', 'Muon_trigSF_val[Z_idx[0]]'])

    
photonSF={}
for m in masses:
    photonSF[m] = "*".join(['Photon_idSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_idSF_val[best_2g_idx2_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx2_m{}]'.format(m)])




def convertOldCustomBins(customBins):
    from ctypes import c_int, byref
    d={}
    for v in ['W','Z']:
        for l in ['ELE','MU']:
            if v =='W':
                if l=='ELE':
                    ana='wen2g'
                else:    
                    ana='wmn2g'
            else: 
                if l=='ELE':
                    ana='zee2g'
                else:    
                    ana='zmm2g'
            d[ana]={}        
            for mass in [7,15,20,30,40,50,55]:
                #create a fake histogram
                if v == 'Z' and mass == 7 and l == 'ELE':
                    binsM = 2
                elif v == 'Z' and mass == 15 and l == 'ELE':
                    binsM = 2
                elif v == 'Z' and mass ==20 and l == 'ELE':
                    binsM = 2
                else:
                    binsM =4
                models=[]    
                h=ROOT.TH2D("test","test",binsM,4,mass+5,110,-10,100)
                hU = ROOT.TH1D("test2","test2",binsM*110,0,binsM*110)
                #loop on the mass bins
                for i in range(1,h.GetNbinsX()+1):
                    xedges=(h.GetXaxis().GetBinLowEdge(i),h.GetXaxis().GetBinUpEdge(i))
                    #now we need to identify the y bins associated with this xbin
                    yedges=[]
                    for b in customBins[v][l][mass][:-1]:
                        binglobal = hU.GetXaxis().FindBin(b)
                        bx = int((binglobal-1)/110)+1
                        if bx!=i:
                            continue
                        by = int((binglobal-1)%110)+1
                        yedges.append(h.GetYaxis().GetBinLowEdge(by))
                    yedges.append(110.)
                    models.append((xedges,yedges))
                d[ana][mass]=models
    print(d)            
                        
                    
                               


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







def getFiles(query,sampleDir,sampleType,era,prod):
    if '/store' in sampleDir:
        ser = pd.Series(subprocess.check_output(['xrdfs', 'root://cmseos.fnal.gov', 'ls', f"{sampleDir}/{sampleType}{era}_{prod}/"], text=True).split("\n"))
        return(list('root://cmseos.fnal.gov/' + ser[ser.str.contains(query)]))
    else:
        ser = pd.Series(subprocess.check_output(['ls', f"{sampleDir}/{sampleType}{era}_{prod}/"], text=True).split("\n"))
        return(list(f"{sampleDir}/{sampleType}{era}_{prod}/" +ser[ser.str.contains(query)]))

def getPlotter(sample,sampleDir,sampleType,eras,prod,analysis):
    plotters=[]
    for era in eras:
        #if we are doing data overrtide the search
        files=[]
        if sampleType=='DATA':
            if analysis in ['wen2g','zee2g']:
                if era=='2018':
                    files = getFiles('EGamma',sampleDir,sampleType,era,prod)
                else:
                    files = getFiles('SingleElectron',sampleDir,sampleType,era,prod)
            else:
                files = getFiles('SingleMuon',sampleDir,sampleType,era,prod)
        else: #MC means query
            files = getFiles(sample,sampleDir,sampleType,era,prod)
        for f in files:
            plotters.append(rdf_plotter(f,isMC=(sampleType=='MC'), tree=analysis, report = "Report_" + analysis))
            #scale with the luminosity
            if sampleType=='MC':
                plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1
            
            #Deal with the HEM cuts
            if era == '2018' and analysis == 'wen2g' and sampleType=='DATA':
                plotters[-1].defaultCuts = "((run>=319077&&(Electron_eta[W_l1_idx]>-1.3||Electron_eta[W_l1_idx]<-3.0)&&(Electron_phi[W_l1_idx]>-0.87||Electron_phi[W_l1_idx]<-1.57))||(run<319077))"
            elif era=='2018' and  analysis == 'wen2g' and sampleType=='MC':   
                plotters[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                plotters.append(rdf_plotter(f, True, tree = analysis, defaultCuts = cutsHEM))
                plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1                
                plotters[-1].addCorrectionFactor(str(38750./59830), "flat")
    p = merged_plotter(plotters)

    return p


def calculate_fake_rate(sampleDir,prod,eras=['2016','2017','2018'],ptbins=[25.,30.,35.,40.,50.,60.,80.,150.],etabins=[-2.5,-2.0,-1.57,-1.44,-0.8,0.8,1.44,1.57,2.0,2.5],arrayName="fake_rate_VH",doMCClosure=False,outdir='.',file_extension='png'):
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
        wjets=getPlotter('WJetsToLNu_HT',sampleDir,'MC',eras,prod,'wmugamma')
        vjets=getPlotter('DYJetsToLL_M50_LO',sampleDir,'MC',eras,prod,'wmugamma')
        tt=getPlotter('TTJets',sampleDir,'MC',eras,prod,'wmugamma')
        fr = merged_plotter([wjets,vjets,tt])
    else:    
        fr = getPlotter('nothing',sampleDir,'DATA',eras,prod,'wmugamma')

        
    cuts_denom = "&&".join([cuts['W']['MU'],'(Photon_preselection>0)'])
#    cuts_denom = '(Photon_preselection)'
    cuts_num = '&&'.join([cuts_denom,'(Photon_cutBased>0)'])
    xedges,yedges,denominator,w2_denom = fr.array2d('Photon_pt','Photon_eta',cuts_denom,('denom','denom',len(ptbins)-1,np.array(ptbins),len(etabins)-1,np.array(etabins)))
    xedges,yedges,numerator,w2_num = fr.array2d('Photon_pt','Photon_eta',cuts_num,('num','num',len(ptbins)-1,np.array(ptbins),len(etabins)-1,np.array(etabins)))
    result = np.zeros_like(numerator)
    fake_rate = np.divide(numerator,denominator,out=result,where=(denominator !=0))

    fake_rate_down,fake_rate_up = clopper_pearson(numerator,denominator)
    #make a plot    
    fig, ax = plt.subplots()

    mesh = ax.pcolormesh(xedges, yedges, fake_rate, cmap='plasma', edgecolors='white', linewidth=0.5)
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
        mh.cms.label(data=False, ax=ax, loc=0)
    else:
        mh.cms.label("Preliminary", data=True, lumi=(lumifb[eras[0]] if len(eras)==1 else lumifb['Run2']), ax=ax, loc=0)
    #print it in C++ format
    #note that we remove the last edge on how the code is defined to work
    st=st+'\n'+f'std::vector<float> {arrayName}_xbins = {{'+','.join([str(x) for x in xedges[:-1]])+'};\n'+f'std::vector<float> {arrayName}_ybins = {{'+','.join([str(y) for y in yedges[:-1]])+'};\n'
    plt.savefig(f'{outdir}/{arrayName}.{file_extension}', dpi=400, bbox_inches='tight')
    return st


    



def getSignalPlotter(sampleDir,prod,eras,analysis,mass,lifetime,signals=['ZH','ggZH','WH','ttH'],modelIndependent=False):
    xsecs = {'ggZ': "0.1057", #units in pb
             'Z': "0.8696",
             'Wplus': "0.84",
             'Wminus': "0.5328",
             'tt':"0.5071"} 

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


    #BRs from PDG
    BRs = {"Z": "(.03363+.03366+.033696)",  #e⁺e⁻, μ⁺μ⁻, τ⁺τ⁻
           "W": "(.1063+.1071+.1138)",      #μν, eν, τν
           "ttSemiLeptonic": "2*(1 - (.1110 + .1140 + .107))*(.1110 + .1140 + .107)" #tt̄->W⁺W⁻bb̄(100%) -> W->(hadronic)W->(leptonic)bb̄->(jets) * 2
           }

    V=[]
    br='2G2Q'
    for sig in signals:
        if sig == 'ttH':
            V.extend(['tt'])
        elif sig == 'WH':
            V.extend(['Wplus','Wminus'])
        elif sig == 'ZH':
            V.extend(['Z'])
        elif sig == 'ggZH':
            V.extend(['ggZ'])
    plotters=[]
    for era in eras:
        for sig in V:
            fs=getFiles(f"{sig}H{br}_M{mass}_ctau{lifetime}_{era}",sampleDir,"MC",era,prod)
            for f in fs:
                plotters.append(rdf_plotter(f, tree=analysis,isMC=True, report = "Report_" + analysis))
                plotters[-1].addCorrectionFactor(lumifb[era], "flat")                
                plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1                
                weight = "(1)"
                if not modelIndependent:
                    weight +="*"+xsecs[sig]
                    if sig in ['Z','ggZ']:
                        weight+="*"+BRs['Z']
                    elif sig in ['Wplus','Wminus']:
                        weight+="*"+BRs['W']
                    elif sig in ['tt']:
                        weight+="*"+BRs['ttSemiLeptonic']
                    plotters[-1].addCorrectionFactor(weight, "flat")
                    
                #Deal with the HEM cuts
                if era=='2018' and  analysis == 'wen2g':   
                    plotters[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                    plotters.append(rdf_plotter(f, isMC=True, tree = analysis, defaultCuts = cutsHEM))
                    plotters[-1].addCorrectionFactor(lumifb[era], "flat")
                    plotters[-1].addCorrectionFactor('1000', "flat") #to conevrt to pb-1                             
                    plotters[-1].addCorrectionFactor(str(38750./59830), "flat")
                    plotters[-1].addCorrectionFactor(weight, "flat")
    p = merged_plotter(plotters)
    return p

def getAnalysis(sampleDir,prod,ana,era='Run2',masses=masses,lifetimes=lifetimes,signals=['ZH','ggZH','WH','ttH'],modelIndependent=False,br=0.01,background_method="fakerate"):
    analysis={}



    
    if era=='Run2':
        eras=['2016','2017','2018']
    else:
        eras=[era]
    analysis['data']=getPlotter('nothing',sampleDir,'DATA',eras,prod,ana)        
    #now create a background plotter with the sideband method we are talking about
    analysis['bkg']={}
    for m in masses:
        if background_method=="abcd":
            analysis['bkg'][m]=abcd_plotter(cuts[ana][m]['sr'],
                                                  cuts[ana][m]['cr_abcd'],
                                                  cuts[ana][m]['ssb'],
                                                  cuts[ana][m]['csb_abcd'],analysis['data'].plotters)
        elif background_method=="fakerate":
            if era!='Run2':
                st = f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0),(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_{era}_vals,fake_rate_{era}_xbins,fake_rate_{era}_ybins)'
            else:#assume run2
                st = f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0),(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_vals,fake_rate_xbins,fake_rate_ybins)'
                
            analysis['bkg'][m]=fakerate_plotter(cuts[ana][m]['presr'],
                                                  cuts[ana][m]['precr'],
                                                  analysis['data'].plotters,
                                                  'fakeRate',
                                                  st
                                                  )                                                  
            #define systematics
            analysis['bkg'][m].define("fakeRate_val","fakeRate[0]")           
            analysis['bkg'][m].define("fakeRate_up","fakeRate[1]")
            analysis['bkg'][m].define("fakeRate_down","fakeRate[2]")
        analysis['bkg'][m].setFillProperties(1001, ROOT.kAzure+5)
        analysis['bkg'][m].setLineProperties(1, ROOT.kAzure+5, 3)        


    #For these MC sampleswe do not apply photon scale factors since they are not used for data MC/comparison
    #We do it for the signals though
    analysis['wjets']=getPlotter('WJetsToLNu_HT',sampleDir,'MC',eras,prod,ana)
    analysis['wjets'].addCorrectionFactor(leptonSF[ana],'flat')
    
    analysis['zjets']=getPlotter('DYJetsToLL_M50_LO',sampleDir,'MC',eras,prod,ana)
    analysis['zjets'].addCorrectionFactor(leptonSF[ana],'flat')

    analysis['tt']=getPlotter('TTJets',sampleDir,'MC',eras,prod,ana)
    analysis['tt'].addCorrectionFactor(leptonSF[ana],'flat')
    
    analysis['tt'].setFillProperties(1001, ROOT.kAzure-2)
    analysis['tt'].setLineProperties(1, ROOT.kAzure-2, 3)

    analysis['wg']=getPlotter('WGToLNuG',sampleDir,'MC',eras,prod,ana)
    analysis['wg'].addCorrectionFactor(leptonSF[ana],'flat')
    
    analysis['wgg']=getPlotter('WGG',sampleDir,'MC',eras,prod,ana)
    analysis['wgg'].addCorrectionFactor(leptonSF[ana],'flat')

    
    analysis['signal']={}
    for m in masses:
        analysis['signal'][m]={}
        i=0
        for ct in lifetimes:
            analysis['signal'][m][ct]={}
            analysis['signal'][m][ct]['sum']=getSignalPlotter(sampleDir,prod,eras,ana,m,ct,signals,modelIndependent)
            analysis['signal'][m][ct]['sum'].addCorrectionFactor(leptonSF[ana],'flat')
            analysis['signal'][m][ct]['sum'].addCorrectionFactor(photonSF[m],'flat')            
            analysis['signal'][m][ct]['sum'].addCorrectionFactor(str(br),'flat')
            for signal in signals:
                analysis['signal'][m][ct][signal]=getSignalPlotter(sampleDir,prod,eras,ana,m,ct,[signal],modelIndependent)
                analysis['signal'][m][ct][signal].addCorrectionFactor(leptonSF[ana],'flat')
                analysis['signal'][m][ct][signal].addCorrectionFactor(photonSF[m],'flat')            
                analysis['signal'][m][ct][signal].addCorrectionFactor(str(br),'flat')
                
    return analysis
            



def runAction(sampleDir,prod,action='fakerate_closure',masses=masses,outputDir='VHresults',era='Run2',analyses=analyses,signals=['ZH','ggZH','WH','ttH'],lifetimes=lifetimes,signal_br=0.01,blinded=False,file_extension='png'):



    if era=='Run2':
        eras=['2016','2017','2018']
    else:
        eras=[era]


    #ACTION: Kinematic Fit plots
    if action=="kinfit_plots":
        for m in masses:
            ana='wmn2g'
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=signals,lifetimes=lifetimes)
            stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
            stack.stack=False
            for ctau in lifetimes:
                stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r'$c\tau=$'+f"{ctau} mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
            stack.hist1d(f"best_2g_raw_mass_m{m}",cuts[ana][m]['sr'],model=('a','a',100,8,m+2),alpha=1,xlabel=r"$m_{\gamma\gamma}$",xunits="GeV",show=False,legend_loc='upper left')
            plt.savefig(f'{outputDir}/kinfit_mass_m{m}.{file_extension}', dpi=400, bbox_inches='tight')
            stack.hist1d(f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],model=('a','a',90,-10,80),alpha=1,xlabel=r"$d_{xy}$",xunits="cm",show=False,logscale=False)
            plt.savefig(f'{outputDir}/kinfit_dxy{m}.{file_extension}', dpi=400, bbox_inches='tight')
            
        
        

        
    #ACTION: Calculate fake rates
    elif action=="fakerate_calc":
        print("Running Calculation of Fake rates")
        with open('common/vhFakeRates.h', "w") as file:
            file.write("#ifndef FAKERATES\n")
            file.write("#define FAKERATES\n")
            for e in eras:
                fr=calculate_fake_rate(sampleDir,prod,[e],arrayName=f"fake_rate_{e}",outdir=outputDir,file_extension=file_extension)
                file.write(fr)
            #write the average all run fake rate for studies
            fr=calculate_fake_rate(sampleDir,prod,eras,arrayName="fake_rate",outdir=outputDir,file_extension=file_extension)
            file.write(fr)
            #write the average all run MC fake rate for studies
            fr=calculate_fake_rate(sampleDir,prod,eras,arrayName="fake_rate_MC",outdir=outputDir,doMCClosure=True,file_extension=file_extension)
            file.write(fr)
            file.write("#endif\n")
    


    #ACTION: fake rate MC Closure
    elif action=="fakerate_closure":
        print("Running MC Closure of Fake rates")        
        for ana in analyses:
            #create a plotter that has all MC as data
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=[],lifetimes=[])
            plotters=[analysis['wjets'],analysis['zjets'],analysis['tt']]
            #create a fake rate MC plotter
            for m in masses:
                print(f"Running {ana} m={m} GeV") 
                bkg=fakerate_plotter(cuts[ana][m]['sr'],
                                            cuts[ana][m]['cr'],
                                            plotters,
                                            'fakeRate_MC',
                                            f'fake_rate(Photon_pt[best_2g_idx1_m{m}],Photon_eta[best_2g_idx1_m{m}],Photon_pt[best_2g_idx2_m{m}],Photon_eta[best_2g_idx2_m{m}],(Photon_cutBased[best_2g_idx1_m{m}]>0),(Photon_cutBased[best_2g_idx2_m{m}]>0),fake_rate_MC_vals,fake_rate_MC_xbins,fake_rate_MC_ybins)')                                                  
                #make a stack plotter and plot the stack
                stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
                stack.add_plotter(bkg,label='Fake rate estimate',typeP='data',error_mode='poisson_bootstrap')               
                stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['zjets'],label='DY+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
                #draw a plot
                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",show=False)
                plt.savefig(f'{outputDir}/fakerate_closure_{ana}_{m}.{file_extension}', dpi=400, bbox_inches='tight')
                stack=None
                fr_plotter=None
            analysis=None
            plotters=None
                


    #ACTION: ABCD MC Closure               
    elif action=="abcd_closure":
        print("Running MC Closure of ABCD Method")
        
        for ana in analyses:
            #create a plotter that has all MC as data
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,signals=[],lifetimes=[])
            if ana in ['wmn2g','wen2g']:
                plotters=[analysis['wjets'],analysis['zjets'],analysis['tt']]
            else:
                plotters=[analysis['zjets'],analysis['tt']]                
            for m in masses:
                print(f"Running {ana} m={m} GeV")
                #create a ABCD MC plotter
                bkg=abcd_plotter(cuts[ana][m]['sr'],
                                 cuts[ana][m]['cr_abcd'],
                                 cuts[ana][m]['ssb'],
                                 cuts[ana][m]['csb_abcd'],plotters)
                #make a stack plotter and plot the stack
                stack=mplhep_plotter(com=center_of_mass[era],data=False,lumi=None)
                stack.add_plotter(bkg,label='ABCD estimate',typeP='data',error_mode='poisson_bootstrap')              #                stack.add_plotter(analysis['data'],label='Data',typeP='data',error_mode='poisson',color='red')               
                stack.add_plotter(analysis['wjets'],label='W+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['zjets'],label='DY+jets',typeP='background',error_mode='w2')
                stack.add_plotter(analysis['tt'],label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
                #draw a plot
                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",show=False)
                plt.savefig(f'{outputDir}/abcd_closure_{ana}_{m}.{file_extension}', dpi=400, bbox_inches='tight')
                stack=None
                fr_plotter=None
            analysis=None
            plotters=None
                
    #ACTION: Control region plots
    elif action=="control_region_plots":
        print("Make Plots of the control region for data and signal for both background estimate methods")
        for ana in analyses:
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,br=signal_br,signals=signals,lifetimes=[100])
            for m in masses:
                print(f"Running {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era])
                stack.add_plotter(analysis['signal'][m][100]['sum'],label=r'$m_{\phi}$='+f"{m} GeV,"+r" $c\tau =$ 100 mm",typeP='signal',error_mode='w2',color='red')               
                stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['cr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",show=False)
                plt.savefig(f'{outputDir}/control_region_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
                stack=None
            analysis=None
            
        
    


            
    #ACTION: Final Plots              
    elif action=="final_plots":
        print("Make final Plots")
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
            
        
        for ana in analyses:
            #create a plotter that has all MC as data
            analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=era,br=signal_br,signals=mySignals[ana],lifetimes=lifetimes)
            for m in masses:
                print(f"Running {ana} m={m} GeV")
                stack=mplhep_plotter(label=analysis_status,data=True,lumi=lumifb[era],com=center_of_mass[era])
                stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='poisson_bootstrap')
                for ctau in lifetimes:
                    stack.add_plotter(analysis['signal'][m][ctau]['sum'],label=r'$m_{\phi}$='+f"{m} GeV,"+r" $c\tau =$ 0 mm",typeP='signal',error_mode='w2',color=signal_colors[ctau])
                if blinded==False:
                    stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
                #draw a plot
                stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr'],binning[ana][m],alpha=1.0,xlabel=r"$d_{xy}$",xunits="cm",show=False,legend_ax=0,legend_loc='upper left')
                if blinded:
                    plt.savefig(f'{outputDir}/blinded_prefit_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
                else:
                    plt.savefig(f'{outputDir}/prefit_{ana}_{m}_{era}.{file_extension}', dpi=400, bbox_inches='tight')
                stack=None
            analysis=None
            
    #ACTION: Make datacards             
    elif action=="make_datacards":
        print("Make Datacards")
        lumiUnc = {'2018': 1.025,
                   '2017': 1.023,
                   '2016': 1.012}
        xsecUnc = {'WH'  : [1+0.005, 1-0.007],
                   'ZH'  : [1+0.038, 1-0.031],
                   'ggZH': [1+0.251, 1-0.189],
                   'ttH' : [1+0.058, 1-0.092]}

        #symmetric
        pdfUnc = {'WH'  : 1.019,
                  'ZH'  : 1.016,
                  'ggZH': 1.024,
                  'ttH' : 1.036}
        mySignals={
            'wmn2g':['WH','ttH'],
            'wen2g':['WH','ttH'],
            'zmm2g':['ZH','ggZH'],
            'zee2g':['ZH','ggZH']}
            
        for ana in analyses:
            for e in eras:
                for m in masses:                    
                    for ctau in lifetimes:
                        analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era=e,br=signal_br,masses=[m],signals=signals,lifetimes=[ctau])                        
                        print(f"Making Datacards for {ana} in era {era} for m={m} GeV and ctau={ctau} mm")
                        for ibinx,bin_setup in enumerate(binning[ana][m]):
                            mass_min = bin_setup[0][0]
                            mass_max=  bin_setup[0][1]
                            for ibiny,biny in enumerate(bin_setup[1][:-1]):
                                dxy_min=bin_setup[1][ibiny]
                                dxy_max=bin_setup[1][ibiny+1]
                                print(f"Bin {mass_min}<=m<{mass_max} {dxy_min}<= dxy <={dxy_max}")
                                
                                cutstring = cuts[ana][m]['sr']+f"*(best_2g_raw_mass_m{m}>={mass_min}&&best_2g_raw_mass_m{m}<{mass_max})*(best_2g_dxy_m{m}>={dxy_min}&&best_2g_dxy_m{m}<{dxy_max})"

                                dcm = cnc_datacard_maker(outDir=outputDir,binname=f"{ana}_m{m}_ctau{ctau}_era{e}_binx{ibinx}_biny{ibiny}",cuts=cutstring)
                                dcm.add('data','data',analysis['data'],{})

                                for signal in mySignals[ana]:                                    
                                    signalUncertainties={
                                        f'CMS_lumi_{e}':{'type':'adhoc','kind':'lnN','value':lumiUnc[e]},                                        
                                        f'CMS_{signal}_xsec':{'type':'adhoc','kind':'lnN','value':f"{xsecUnc[signal][1]}/{xsecUnc[signal][0]}"},
                                        'CMS_pdf':{'type':'adhoc','kind':'lnN','value':pdfUnc[signal]}
                                    }
                                    #add lepton ID SF, some string manipulations in place
                                    leptonSFs = leptonSF[ana]
                                    individuals = list(set([x.split('SF_val')[0] for x in leptonSFs.split("*")]))
                                    for unc in individuals:
                                        l = unc.replace('Muon','mu').replace('Electron','ele')
                                        signalUncertainties[f'CMS_{l}_{e}']={'type':'weightAsymm','weightUp':leptonSFs.replace(unc+'SF_val',unc+'SF_up'),'weightDown':leptonSFs.replace(unc+'SF_val',unc+'SF_down'),'weightOrig':leptonSFs}

                                    #Same for photons
                                    photonSFs = photonSF[m]
                                    individuals = list(set([x.split('SF_val')[0] for x in photonSFs.split("*")]))

                                    for unc in individuals:
                                        l = unc.replace('Photon','gamma')
                                        signalUncertainties[f'CMS_{l}_{e}']={'type':'weightAsymm','weightUp':photonSFs.replace(unc+'SF_val',unc+'SF_up'),'weightDown':photonSFs.replace(unc+'SF_val',unc+'SF_down'),'weightOrig':photonSFs}

                                    #now evaluate the systematics because of photon scale and resolution
                                    #first define the new photon momenta
                                    analysis['signal'][m][ctau][signal].define("Photon_ptSmearUp", "ptSmearUp(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaUp)")
                                    analysis['signal'][m][ctau][signal].define("Photon_ptSmearDown", "ptSmearDown(Photon_pt, Photon_eta, Photon_phi, Photon_dEsigmaDown)")
                                    analysis['signal'][m][ctau][signal].define("Photon_ptScaleUp", "Photon_pt*photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{era}_val, PHO_scaledown_{era}_bins, 1)".format(era = "2016preVFP" if e=="2016" else e))
                                    analysis['signal'][m][ctau][signal].define("Photon_ptScaleDown", "Photon_pt*photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaleup_{era}_val, PHO_scaledown_{era}_bins, 1)".format(era = "2016preVFP" if e=="2016" else e))
                                    for unc in ['Scale','Smear']:
                                        for direction in ['Up','Down']:
                                            #calculate the mass and the dxy with new pt
                                            analysis['signal'][m][ctau][signal].define(f"best_2g_{unc}{direction}_m{m}_info", f"kinfit_systematics(Photon_pt{unc}{direction}, Photon_eta, Photon_phi, Photon_isScEtaEB, Photon_isScEtaEE, best_2g_idx1_m{m}, best_2g_idx2_m{m}, {m})")
                                            analysis['signal'][m][ctau][signal].define(f"best_2g_dxy_{unc}{direction}_m{m}", f"best_2g_{unc}{direction}_m{m}_info[0]")
                                            analysis['signal'][m][ctau][signal].define(f"best_2g_raw_mass_{unc}{direction}_m{m}", f"best_2g_{unc}{direction}_m{m}_info[1]")
                                        systName=unc.replace('Scale','scale').replace('Smear','res')     
                                        signalUncertainties[f'CMS_gamma_{systName}_{e}']={'type':'replication','originals':['Photon_pt',f'best_2g_raw_mass_m{m}',f'best_2g_dxy_m{m}'],
                                                                                      'replacementsUp':[f"Photon_pt{unc}Up",f"best_2g_raw_mass_{unc}Up_m{m}",f"best_2g_dxy_{unc}Up_m{m}"],
                                                                                      'replacementsDown':[f"Photon_pt{unc}Down",f"best_2g_raw_mass_{unc}Down_m{m}",f"best_2g_dxy_{unc}Down_m{m}"]}
                                    dcm.add(signal,'signal',analysis['signal'][m][ctau][signal],uncertainties=signalUncertainties)
                                    
                                    
                                dcm.add('bkg','background',analysis['bkg'][m],uncertainties={
                                    f"CMS_DDP_{ana}_m{m}_{e}_binx{ibinx}_biny{ibiny}_CR_stats":{'type':'statAsym'},
                                    f"CMS_DDP_fakeRateUnc_{e}":{'type':'weightAsymm','weightUp':'fakeRate_up','weightDown':'fakeRate_down','weightOrig':'fakeRate_val'}
                                    
                                },error_mode='poisson_bootstrap')
                                #write only bins that have bkg>0 or data>0
                                signalTot = sum([dcm.rates[s] for s in  mySignals[ana]])
                                print('total signal :',signalTot)
                                if (dcm.rates['bkg']+dcm.data)>0 and signalTot>0:
                                    dcm.write()
                                    
                        analysis=None
                        gc.collect()
                        #This forced garbage collection enables better running in slow machines, eg your laptop
        





#debugging code
#if __name__ == '__main__':
#
#    sampleDir='/tank/ddp/DDP'
#    prod='2025_12_23'
#    ana='wmn2g'
#    m=20
#    analysis=getAnalysis(sampleDir,prod,ana,background_method='fakerate',era='2016',br=0.01,signals=['WH','ttH'],lifetimes=[100])
#    print(f"Running {ana} m={m} GeV")
#    stack=mplhep_plotter()
#    stack.add_plotter(analysis['bkg'][m],label='Background',typeP='background',error_mode='w2')               
#    stack.add_plotter(analysis['signal'][m][100],label=r'$m_{\phi}$='+f"{m} GeV,"+r" $c\tau =$ 100 mm",typeP='signal',error_mode='w2',color='red')               
#    stack.add_plotter(analysis['data'],label="Data",typeP='data',error_mode='poisson')               
#    stack.define("loose_muon", "(Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pfIsoId>1&&Muon_looseId>0)")
#    stack.define("tight_muon", "loose_muon&&Muon_tightId>0&&Muon_pfIsoId>3")
#   
#    stack.define("overlap_muon", "(Muon_pt>3&&abs(Muon_eta)<2.4&&Muon_softId>0)")
#    stack.define("Photon_overlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[overlap_muon], Muon_eta[overlap_muon],0.8,0.8)")
#    stack.define("Muon_nloose", "Sum(loose_muon)")
#    stack.unrolledCustom(f"best_2g_raw_mass_m{m}",f"best_2g_dxy_m{m}",cuts[ana][m]['sr']+"&&(Photon_overlap[best_2g_idx1_m20]==0 && Photon_overlap[best_2g_idx2_m20]==0)",binning[ana][m],alpha=1.0,xlabel=r'$d_{xy}$',xunits='cm',show=True)



    
#    ana=getAnalysis('/tank/ddp/DDP','2025_12_23','wmn2g',background_method='fakerate',era='Run2')
#    ana['wjets'].define("Photon_prompt","prompt_photon(Photon_genPartFlav,Photon_genPartIdx,GenPart_pdgId,GenPart_genPartIdxMother)")
#    ana['wg'].define("Photon_prompt","prompt_photon(Photon_genPartFlav,Photon_genPartIdx,GenPart_pdgId,GenPart_genPartIdxMother)")
#    ana['wgg'].define("Photon_prompt","prompt_photon(Photon_genPartFlav,Photon_genPartIdx,GenPart_pdgId,GenPart_genPartIdxMother)")
#    ana['wjets'].addCorrectionFactor('((Photon_prompt[best_2g_idx1_m20]+Photon_prompt[best_2g_idx2_m20])==0)','flat')
#    ana['wg'].addCorrectionFactor('((Photon_prompt[best_2g_idx1_m20]+Photon_prompt[best_2g_idx2_m20])==1)','flat')
#    ana['wgg'].addCorrectionFactor('((Photon_prompt[best_2g_idx1_m20]+Photon_prompt[best_2g_idx2_m20])==2)','flat')

#    abcd=abcd_plotter(cuts['wmn2g'][20]['sr'],
#                     cuts['wmn2g'][20]['cr'],
#                     cuts['wmn2g'][20]['ssb'],
#                     cuts['wmn2g'][20]['csb'],[ana['wjets'],ana['wg'],ana['wgg']])

#    fr=fakerate_plotter(cuts['wmn2g'][20]['sr'],
#                                cuts['wmn2g'][20]['cr'],
#                                [ana['wjets'],ana['wg'],ana['wgg']],
#                                'fakeRate_MC',
#                                f'fake_rate(Photon_pt[best_2g_idx1_m20],Photon_eta[best_2g_idx1_m20],Photon_pt[best_2g_idx2_m20],Photon_eta[best_2g_idx2_m20],(Photon_cutBased[best_2g_idx1_m20]>0),(Photon_cutBased[best_2g_idx2_m20]>0),fake_rate_MC_vals,fake_rate_MC_xbins,fake_rate_MC_ybins)')                                                  


#    bkg_plotter=merged_plotter([ana['wjets'],ana['wg'],ana['wgg']])
    


#    stack=mplhep_plotter()
#    stack.add_plotter(ana['bkg'][20],name='bkg',label='Background',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['signal'][20][100],name='sig100',label=r'$m_{\phi}=20$ GeV, $c\tau=100$ mm',typeP='signal',error_mode='w2',color='red')
#    stack.add_plotter(ana['signal'][20][1000],name='sig1000',label=r'$m_{\phi}=20$ GeV, $c\tau=1000$ mm',typeP='signal',error_mode='w2',color='orange')    
#    stack.add_plotter(ana['data'],name='data',label='Data',typeP='data',error_mode='w2')
    
#    stack.add_plotter(ana['tt'],name='ttjets',label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['zjets'],name='zjets',label='DY+jets',typeP='background',error_mode='w2')   
#    stack.add_plotter(bkg_plotter,name='wjets',label='W+jets',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['wjets'],name='wjets',label='0 Prompt',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['wg'],name='wg',label='1 Prompt',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['wgg'],name='wgg',label='2 Prompt',typeP='background',error_mode='w2')
#    stack.add_plotter(fr,name='fr',label='Fake Rate estimate',typeP='data',error_mode='w2')
#    stack.add_plotter(abcd,name='fr',label='ABCD estimate',typeP='data',error_mode='w2')



#    stack.define("loose_muon", "Muon_pt>5&&Muon_looseId==1&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pfIsoId>1")
#    stack.define("tight_muon", "loose_muon&&Muon_tightId&&Muon_pfIsoId>3")
#    stack.define("veto_muon", "(Muon_softId>0)")

#    stack.define("drlg1","deltaR_lgamma(Photon_eta[best_2g_idx1_m20], Photon_phi[best_2g_idx1_m20],Muon_eta[veto_muon],Muon_phi[veto_muon])")
#    stack.define("drlg2","deltaR_lgamma(Photon_eta[best_2g_idx2_m20], Photon_phi[best_2g_idx2_m20],Muon_eta[veto_muon],Muon_phi[veto_muon])")



#    stack.unrolledCustom("best_2g_raw_mass_m20","best_2g_dxy_m20",cuts['wmn2g'][20]['sr']+"*(drlg1>1&&drlg2>1)",binning['wmn2g'][20],alpha=0.5,xlabel=r"$d_{xy}$",xunits="cm")



#    ana=getAnalysis('/tank/ddp/DDP','2025_12_23','wmn2g',background_method='fakerate',era='Run2')
#    stack=mplhep_plotter()
#    stack.add_plotter(ana['data'],name='data',label='Data',typeP='data',error_mode='poisson')
#    stack.add_plotter(ana['signal'][20][100],name='signal',label='Signal',typeP='signal',color='red',error_mode='w2')
#    stack.stack=False
    
    
#    stack.add_plotter(ana['bkg'][20],name='MC',label='MC',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['wjets'],name='wjets',label='W+jets',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['zjets'],name='zjets',label='Z+jets',typeP='background',error_mode='w2')
#    stack.add_plotter(ana['tt'],name='ttjets',label=r'$t\bar{t}$+jets',typeP='background',error_mode='w2')



#    mc=merged_plotter([ana['wjets'],ana['zjets'],ana['tt']])

#    stack.define("loose_muon", "Muon_pt>5&&Muon_looseId==1&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pfIsoId>1")
#    stack.define("tight_muon", "loose_muon&&Muon_tightId&&Muon_pfIsoId>3")
#    stack.define("veto_muon", "(Muon_softId>0)")

#    stack.define("drlg1","deltaR_lgamma(Photon_eta[best_2g_idx1_m20], Photon_phi[best_2g_idx1_m20],Muon_eta[veto_muon],Muon_phi[veto_muon])")
#    stack.define("drlg2","deltaR_lgamma(Photon_eta[best_2g_idx2_m20], Photon_phi[best_2g_idx2_m20],Muon_eta[veto_muon],Muon_phi[veto_muon])")


#    stack.define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[veto_muon], Muon_eta[veto_muon])")
#    stack.define("overlap", "(Photon_muOverlap[best_2g_idx1_m20]||Photon_muOverlap[best_2g_idx2_m20])")
#    stack.hist1d('drlg1',cuts['wmn2g'][20]['sr'],('a','a',4,0,4),xlabel=r"$Delta R(\mu,\gamma)$",alpha=0.5)
    
    

                    
