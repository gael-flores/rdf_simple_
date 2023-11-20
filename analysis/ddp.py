import ROOT
from common.pyhelpers import load_meta_data
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')


def zmumuH(data,phi_mass=[5,10,20,30]):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)

    ####################
    #ANALYSIS CODE HERE#        
    ####################

    #pass HLT
    zmm = dataframe.Filter('HLT_passed','passed HLT')
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
    zmm =zmm.Define("Photon_IDNoIso","((Photon_isScEtaEB&&Photon_hoe<0.04596&&Photon_sieie<0.0106)||(Photon_isScEtaEE&&Photon_hoe<0.0590&&Photon_sieie<0.0272))&&(Photon_electronVeto==1)&&(Photon_eta>-2.5&&Photon_eta<2.5)").Filter('Sum(Photon_IDNoIso)>0','At least one ID No Iso Photon')

#    zmm =zmm.Define("Photon_IDNoIso","((Photon_isScEtaEB&&Photon_hoe<0.2&&Photon_sieie<0.05)||(Photon_isScEtaEE&&Photon_hoe<0.2&&Photon_sieie<0.05))&&(Photon_electronVeto==1)&&(Photon_eta>-2.5&&Photon_eta<2.5)").Filter('Sum(Photon_IDNoIso)>0','At least one ID No Iso Photon')

    #Correct Photon Isolation for both muons and photons
    zmm=zmm.Define("Photon_corrIso","correct_gammaIso_for_muons_and_photons(Zmumu_idx,Muon_pt,Muon_eta,Muon_phi,Photon_pt,Photon_eta,Photon_phi,Photon_pfRelIso03_all,Photon_IDNoIso)")
    #Define ID and isolation
    zmm=zmm.Define("Photon_ID","Photon_IDNoIso==1 &&Photon_corrIso<0.1")
    #At least one Photon
    zmm=zmm.Filter('Sum(Photon_ID==1)>0','At least one ID photon')
    #FSR recovery
    zmm=zmm.Define("fsr","fsr_recovery(Zmumu_idx,Muon_pt[good_muons], Muon_eta[good_muons], Muon_phi[good_muons], Muon_mass[good_muons],Photon_pt,Photon_eta,Photon_phi,Photon_ID)");
    #Correct Muon isolation for photons
    zmm=zmm.Define("corr_muon_iso","correct_muoniso_for_photons(Muon_pt, Muon_eta,Muon_phi,Muon_pfRelIso03_all,Photon_pt,Photon_eta,Photon_phi,fsr)");
    zmm=zmm.Filter("(corr_muon_iso[Zmumu_idx[0]]<0.15) &&(corr_muon_iso[Zmumu_idx[1]]<0.15)","Muon Isolation with FSR recovery")
    #Now cut on the ll+gamma mass around the Z 
    zmm=zmm.Define("llgamma_mass","calculate_llgamma_mass(Zmumu_idx,Muon_pt[good_muons], Muon_eta[good_muons], Muon_phi[good_muons], Muon_mass[good_muons],Photon_pt,Photon_eta,Photon_phi, fsr)").Filter('llgamma_mass>75&&llgamma_mass<115','Z+gamma mass near the Z')
    #Define good photons
    zmm=zmm.Define('good_photons','fsr==0&&Photon_ID==1')

    ##### gamma gamma +X analysis ######
    ####################################

    #Al least Two good photons
    zmm2g=zmm.Filter('Sum(good_photons)>1','At least 2 good no FSR photons')
    for mass in phi_mass:
        zmm2g=zmm2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEE[good_photons], Photon_isScEtaEB[good_photons],{})'.format(float(mass)))
        zmm2g=zmm2g.Define('best_2g_gamma1_pt_m{}'.format(mass),'raw_best_2g_m{}[0]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_eta_m{}'.format(mass),'raw_best_2g_m{}[1]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_phi_m{}'.format(mass),'raw_best_2g_m{}[2]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_pt_m{}'.format(mass),'raw_best_2g_m{}[3]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_eta_m{}'.format(mass),'raw_best_2g_m{}[4]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_phi_m{}'.format(mass),'raw_best_2g_m{}[5]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))

    
    ##### 3 gamma analysis ######
    #Targets events where one two photons are merged
    ####################################
    zmm3g=zmm.Filter('Sum(good_photons)==3','Exactly three no FSR photons')
    for mass in phi_mass:
        zmm3g=zmm3g.Define('raw_best_3g_m{}'.format(mass),'best_3gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEE[good_photons], Photon_isScEtaEB[good_photons],{})'.format(float(mass)))
        zmm3g=zmm3g.Define('best_3g_phi_gamma1_pt_m{}'.format(mass),'raw_best_3g_m{}[0]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_gamma1_eta_m{}'.format(mass),'raw_best_3g_m{}[1]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_gamma1_phi_m{}'.format(mass),'raw_best_3g_m{}[2]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_gamma2_pt_m{}'.format(mass),'raw_best_3g_m{}[3]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_gamma2_eta_m{}'.format(mass),'raw_best_3g_m{}[4]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_gamma2_phi_m{}'.format(mass),'raw_best_3g_m{}[5]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_gamma3_pt_m{}'.format(mass),'raw_best_3g_m{}[9]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_gamma3_eta_m{}'.format(mass),'raw_best_3g_m{}[10]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_gamma3_phi_m{}'.format(mass),'raw_best_3g_m{}[11]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_dxy_m{}'.format(mass),'raw_best_3g_m{}[6]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_valid_m{}'.format(mass),'raw_best_3g_m{}[7]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_phi_mass_m{}'.format(mass),'raw_best_3g_m{}[8]'.format(mass))
        zmm3g=zmm3g.Define('best_3g_raw_mass_m{}'.format(mass),'raw_best_3g_m{}[12]'.format(mass))




    ##### 4 gamma  analysis       ######
    ####################################
    #Four good photons
    zmm4g=zmm.Filter('Sum(good_photons)>3','At least 4 good no FSR photons')
    #Pick the best combination dependening on the mass
    for mass in phi_mass:
        zmm4g=zmm4g.Define('raw_best_4g_m{}'.format(mass),'best_4gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEE[good_photons], Photon_isScEtaEB[good_photons],{})'.format(float(mass)))
        zmm4g=zmm4g.Define('best_4g_phi1_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[0]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[1]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[2]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[3]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[4]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[5]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_dxy_m{}'.format(mass),'raw_best_4g_m{}[6]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[8]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[9]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[10]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[11]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[12]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[13]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_dxy_m{}'.format(mass),'raw_best_4g_m{}[14]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi1_mass_m{}'.format(mass),'raw_best_4g_m{}[16]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_phi2_mass_m{}'.format(mass),'raw_best_4g_m{}[17]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_uncorr_mass_m{}'.format(mass),'raw_best_4g_m{}[18]'.format(mass))
        zmm4g=zmm4g.Define('best_4g_corr_mass_m{}'.format(mass),'raw_best_4g_m{}[19]'.format(mass))
     
        
    actions.append(zmm2g.Snapshot('zmm2g','zmm2g.root',"best_2g.*|sample_.*"))
    actions.append(zmm3g.Snapshot('zmm3g','zmm3g.root',"best_3g.*|sample_.*"))
    actions.append(zmm4g.Snapshot('zmm4g','zmm4g.root',"best_4g.*|sample_.*"))
    r=zmm2g.Report()
    r.Print()
    return actions
    

def analysis(data):
    phi_mass=[5,10,20,30]

    actions = []
    actions.append(zmumuH(data,phi_mass))

    

    
