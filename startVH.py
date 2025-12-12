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
            if analysis in ['wenu2g','zee2g']:
                if era==2018:
                    files = getFiles('EGamma',sampleDir,sampleType,era,prod)
                else:
                    files = getFiles('SingleElectron',sampleDir,sampleType,era,prod)
            else:
                files = getFiles('SingleMuon',sampleDir,sampleType,era,prod)
        else: #MC means query
            files = getFiles(sample,sampleDir,sampleType,era,prod)
        for f in files:
            plotters.append(rdf_plotter(f, tree=analysis, report = "Report_" + analysis))
            #Deal with the HEM cuts
            if era == '2018' and analysis == 'wen2g' and sampleType=='DATA':
                plotters[-1].defaultCuts = "((run>=319077&&(Electron_eta[W_l1_idx]>-1.3||Electron_eta[W_l1_idx]<-3.0)&&(Electron_phi[W_l1_idx]>-0.87||Electron_phi[W_l1_idx]<-1.57))||(run<319077))"
            elif era=='2018' and  analysis == 'wen2g' and sampleType=='MC':   
                plotters[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                plotters.append(rdf_plotter(f, True, tree = analysis, defaultCuts = cutsHEM))
                plotters[-1].addCorrectionFactor(str(38750./59830), "flat")
    return merged_plotter(plotters)


def getSignalPlotter(sampleDir,prod,eras,analysis,mass,lifetime,signals=['ZH','ggZH','WH','ttH'],modelIndependent=False):
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
                plotters.append(rdf_plotter(f, tree=analysis, report = "Report_" + analysis))
                weight = "(1)"
                if not modelIndependent:
                    weight +="*"+xsecs[v]
                    if "Z" in v:
                        weight+="*"+BRs['Z']
                    elif 'W' in v:
                        weight+="*"+BRs['W']
                    elif 'tt' in v:
                        weight+="*"+BRs['ttSemiLeptonic']
                    plotters[-1].addCorrectionFactor(weight, "flat")

                #Deal with the HEM cuts
                elif era=='2018' and  analysis == 'wen2g':   
                    plotters[-1].addCorrectionFactor(str(21080.0/59830), "flat")
                    plotters.append(rdf_plotter(f, True, tree = analysis, defaultCuts = cutsHEM))
                    plotters[-1].addCorrectionFactor(str(38750./59830), "flat")
                    plotters[-1].addCorrectionFactor(weight, "flat")
    return merged_plotter(plotters)


def getAnalysis(sampleDir,prod,ana,eras=['2016','2017','2018'],masses=[15, 20, 30, 40, 50, 55],lifetimes=[0, 10, 20, 50, 100, 1000],signals=['ZH','ggZH','WH','ttH'],modelIndependent=False):
    analysis={}
    analysis['data']=getPlotter('nothing',sampleDir,'DATA',eras,prod,ana)
    analysis['wjets']=getPlotter('WJetsToLNu_HT',sampleDir,'MC',eras,prod,ana)
    analysis['zjets']=getPlotter('DYJetsToLL_M50_LO',sampleDir,'MC',eras,prod,ana)
    analysis['tt']=getPlotter('TTJets',sampleDir,'MC',eras,prod,ana)
    analysis['signal']={}
    for m in masses:
        analysis['signal'][m]={}
        for ct in lifetimes:
            analysis['signal'][m][ct]=getSignalPlotter(sampleDir,prod,eras,ana,m,ct,signals,modelIndependent)
    return analysis       
            
        

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
        analysis = getAnalysis("/tank/ddp/DDP","2025_12_12",'wmn2g')

    except Exception as e:
        import traceback
        traceback.print_exc()
        ROOT.gSystem.Exit(0)
