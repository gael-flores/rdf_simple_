import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')


def zmumuH(data):
    #Declare dataframe
    dataframe =ROOT.RDataFrame(data)

    #pass HLT
    zmm = dataframe.Filter('HLT_IsoMu27==1','passed HLT')
    #Apply Muon ID (no ISO)
    zmm = zmm.Define("good_muons",'Muon_looseId==1').Filter("Sum(good_muons)>1","At least two good muons")
    #Look for + and - muons
    zmm = zmm.Filter("Sum(Muon_charge[good_muons]==1)>0 && Sum(Muon_charge[good_muons]==-1)>0","At least one opposite sign muon pair")
    #create the best Zmumu candidate 
    zmm = zmm.Define("Zmumu_idx", "best_z(Muon_pt[good_muons], Muon_eta[good_muons], Muon_phi[good_muons], Muon_mass[good_muons], Muon_charge[good_muons])")
    #Apply Thresholds to the muon pts
    zmm = zmm.Filter("(Muon_pt[Zmumu_idx[0]]>8&&Muon_pt[Zmumu_idx[1]]>8)&&(Muon_pt[Zmumu_idx[0]]>27||Muon_pt[Zmumu_idx[1]]>27)","Muons above threshold")
    #require just a super cluster
    zmm = zmm.Filter("nPhoton>0","At least one super cluster")
    
    #Apply photon ID (no ISO)
#    zmm =zmm.Define("Photon_IDNoIso","((Photon_isScEtaEB&&Photon_hoe<0.04596&&Photon_sieie<0.0106)||(Photon_isScEtaEE&&Photon_hoe<0.0590&&Photon_sieie<0.0272))&&(Photon_electronVeto==1)&&(Photon_eta>-2.5&&Photon_eta<2.5)").Filter('Sum(Photon_IDNoIso)>0','At least one ID No Iso Photon')

    zmm =zmm.Define("Photon_IDNoIso","((Photon_isScEtaEB&&Photon_hoe<0.2&&Photon_sieie<0.05)||(Photon_isScEtaEE&&Photon_hoe<0.2&&Photon_sieie<0.05))&&(Photon_electronVeto==1)&&(Photon_eta>-2.5&&Photon_eta<2.5)").Filter('Sum(Photon_IDNoIso)>0','At least one ID No Iso Photon')

    #Correct Photon Isolation for both muons and photons
    zmm=zmm.Define("Photon_corrIso","correct_gammaIso_for_muons_and_photons(Zmumu_idx,Muon_pt,Muon_eta,Muon_phi,Photon_pt,Photon_eta,Photon_phi,Photon_pfRelIso03_all,Photon_IDNoIso)")
    #Define ID and isolation
    zmm=zmm.Define("Photon_ID","Photon_IDNoIso==1 &&Photon_corrIso<0.15")
    #At least one Photon
    zmm=zmm.Filter('Sum(Photon_ID==1)>0','At least one ID photon')
    #FSR recovery
    zmm=zmm.Define("fsr","fsr_recovery(Zmumu_idx,Muon_pt[good_muons], Muon_eta[good_muons], Muon_phi[good_muons], Muon_mass[good_muons],Photon_pt,Photon_eta,Photon_phi,Photon_ID)");
    #Correct Muon isolation for photons
    zmm=zmm.Define("corr_muon_iso","correct_muoniso_for_photons(Muon_pt, Muon_eta,Muon_phi,Muon_pfRelIso03_all,Photon_pt,Photon_eta,Photon_phi,fsr)");
    zmm=zmm.Filter("(corr_muon_iso[Zmumu_idx[0]]<0.15) &&(corr_muon_iso[Zmumu_idx[1]]<0.15)","Muon Isolation with FSR recovery")
    #Now cut on the ll+gamma mass around the Z 
    zmm=zmm.Define("llgamma_mass","calculate_llgamma_mass(Zmumu_idx,Muon_pt[good_muons], Muon_eta[good_muons], Muon_phi[good_muons], Muon_mass[good_muons],Photon_pt,Photon_eta,Photon_phi, fsr)").Filter('llgamma_mass>75&&llgamma_mass<115','Z+gamma mass near the Z')

    #Exacty four good photons
    zmm4g=zmm.Define('good_photons','fsr==0&&Photon_ID==1').Filter('Sum(good_photons)==4','Exactly 4 good no FSR photons')

#    report=zmm2g.Report()
#    report.Print()
    zmm4g.Snapshot("zmm4g","zmm4g.root")


def analysis(data):
    zmumuH(data)
