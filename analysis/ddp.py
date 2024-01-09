import ROOT
from common.pyhelpers import load_meta_data
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
def zeeH(data,phi_mass=[5,10,20,30]):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################

    return actions

def wenuH(data,phi_mass=[5,10,20,30]):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################

    return actions

def wmunuH(data,phi_mass=[5,10,20,30]):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################

    return actions

def ggH(data,phi_mass=[5,10,20,30]):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    #pass HLT
    ggH = dataframe.Filter('HLT_passed','passed HLT')
    
    #your code here 
    ggH=ggH.Filter('nPhoton>2','At least three photons')
 
# Implement matching only for signal efficiency studies    
#    ggH=ggH.Define('genPhotonDR','minMatchDR(Photon_eta,Photon_phi,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)')
#    ggH=ggH.Filter('Sum(genPhotonDR<0.2)>2','At least three matched photons')

    ggH=ggH.Define("Photon_IDNoIso","((Photon_isScEtaEB&&Photon_hoe<0.04596&&Photon_sieie<0.0106)||(Photon_isScEtaEE&&Photon_hoe<0.0590&&Photon_sieie<0.0272))&&(Photon_electronVeto==1)&&(Photon_eta>-2.5&&Photon_eta<2.5)")
    
    #Correct Photon Isolation for photons
    ggH=ggH.Define("Photon_corrIso","correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,Photon_pfRelIso03_all,Photon_IDNoIso)")
    #Define ID and isolation
    ggH=ggH.Define("Photon_ID","Photon_IDNoIso==1 &&Photon_corrIso<0.1")
        
    #Anti-photon RDF
    ggH3g_antiID = ggH.Filter("Sum(Photon_ID==1)>1 && Sum(Photon_ID==0)>0","at least 2 good and 1 bad photon")

    #At least three Photons
    ggH=ggH.Filter('Sum(Photon_IDNoIso)>2','At least 3 ID No Iso Photon')
    ggH=ggH.Filter('Sum(Photon_ID==1)>2','At least 3 ID photon')

    #Filter out muons above 10Gev and electrons above 15GeV
    ggH=ggH.Define('electron_veto','Electron_cutBased==2 && Electron_pt>15').Filter('Sum(electron_veto==1)==0','veto loose id electrons above 15GeV')
    ggH=ggH.Define('muon_veto','Muon_looseId==1 && Muon_pt>10').Filter('Sum(muon_veto==1)==0','veto loose id muons above 10GeV')
    
    ggH3g_antiID=ggH3g_antiID.Define('electron_veto','Electron_cutBased==2 && Electron_pt>15').Filter('Sum(electron_veto==1)==0','veto loose id electrons above 15GeV')
    ggH3g_antiID=ggH3g_antiID.Define('muon_veto','Muon_looseId==1 && Muon_pt>10').Filter('Sum(muon_veto==1)==0','veto loose id muons above 10GeV')
    
    #exactly 3 photons
    ggH3g=ggH.Filter('Sum(Photon_ID==1)==3','exactly 3 ID photon')
    ggH3g=ggH3g.Define('good_photons','Photon_ID==1')
    ggH3g=ggH3g.Define("m_3g", "InvariantMass(Photon_pt[good_photons], Photon_eta[good_photons], Photon_phi[good_photons], Photon_mass[good_photons])")
    for mass in phi_mass:
        ggH3g=ggH3g.Define('raw_best_3g_m{}'.format(mass),"best_3gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})".format(float(mass)))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_pt_m{}'.format(mass),'raw_best_3g_m{}[0]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_eta_m{}'.format(mass),'raw_best_3g_m{}[1]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_phi_m{}'.format(mass),'raw_best_3g_m{}[2]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_pt_m{}'.format(mass),'raw_best_3g_m{}[3]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_eta_m{}'.format(mass),'raw_best_3g_m{}[4]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_phi_m{}'.format(mass),'raw_best_3g_m{}[5]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_pt_m{}'.format(mass),'raw_best_3g_m{}[9]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_eta_m{}'.format(mass),'raw_best_3g_m{}[10]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_phi_m{}'.format(mass),'raw_best_3g_m{}[11]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_dxy_m{}'.format(mass),'raw_best_3g_m{}[6]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_valid_m{}'.format(mass),'raw_best_3g_m{}[7]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_mass_m{}'.format(mass),'raw_best_3g_m{}[8]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_raw_mass_m{}'.format(mass),'raw_best_3g_m{}[12]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_corr_mass_m{}'.format(mass),'raw_best_3g_m{}[13]'.format(mass))
        
    #blinding region for data samples only
    ggH3g=ggH3g.Define('non_MC_cut','sample_isMC==0 && best_3g_raw_mass_m30<30|best_3g_raw_mass_m30>140')
    ggH3g=ggH3g.Filter('sample_isMC==1 | non_MC_cut==1','blinding data samples')
   
    #exactly 3 photons: two good and one failing ID
    ggH3g_antiID=ggH3g_antiID.Filter('Sum(Photon_ID==1)==2 && Sum(Photon_ID==0)==1','exactly 2 ID photon and 1 fail ID photon')
    ggH3g_antiID=ggH3g_antiID.Define('good_photons','Photon_ID==1|Photon_ID==0')
    ggH3g_antiID=ggH3g_antiID.Define("m_3g", "InvariantMass(Photon_pt[good_photons], Photon_eta[good_photons], Photon_phi[good_photons], Photon_mass[good_photons])")
    for mass in phi_mass:
        ggH3g_antiID=ggH3g_antiID.Define('raw_best_3g_m{}'.format(mass),"best_3gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})".format(float(mass)))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_gamma1_pt_m{}'.format(mass),'raw_best_3g_m{}[0]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_gamma1_eta_m{}'.format(mass),'raw_best_3g_m{}[1]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_gamma1_phi_m{}'.format(mass),'raw_best_3g_m{}[2]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_gamma2_pt_m{}'.format(mass),'raw_best_3g_m{}[3]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_gamma2_eta_m{}'.format(mass),'raw_best_3g_m{}[4]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_gamma2_phi_m{}'.format(mass),'raw_best_3g_m{}[5]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_gamma3_pt_m{}'.format(mass),'raw_best_3g_m{}[9]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_gamma3_eta_m{}'.format(mass),'raw_best_3g_m{}[10]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_gamma3_phi_m{}'.format(mass),'raw_best_3g_m{}[11]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_dxy_m{}'.format(mass),'raw_best_3g_m{}[6]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_valid_m{}'.format(mass),'raw_best_3g_m{}[7]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_phi_mass_m{}'.format(mass),'raw_best_3g_m{}[8]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_raw_mass_m{}'.format(mass),'raw_best_3g_m{}[12]'.format(mass))
        ggH3g_antiID=ggH3g_antiID.Define('best_3g_corr_mass_m{}'.format(mass),'raw_best_3g_m{}[13]'.format(mass))
        
    #blinding region for data samples only
    ggH3g_antiID=ggH3g_antiID.Define('non_MC_cut','sample_isMC==0 && best_3g_raw_mass_m30<30|best_3g_raw_mass_m30>140')
    ggH3g_antiID=ggH3g_antiID.Filter('sample_isMC==1 | non_MC_cut==1','blinding data samples')

    #at least 4 photons
    ggH4g=ggH.Filter('Sum(Photon_ID==1)>3','at least 4 ID photon')
    ggH4g=ggH4g.Define('good_photons','Photon_ID==1')
    ggH4g=ggH4g.Define("m_4g", "InvariantMass(Photon_pt[good_photons], Photon_eta[good_photons], Photon_phi[good_photons], Photon_mass[good_photons])")
    for mass in phi_mass:
        ggH4g=ggH4g.Define('raw_best_4g_m{}'.format(mass),"best_4gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})".format(float(mass)))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[0]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[1]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[2]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[3]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[4]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[5]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_dxy_m{}'.format(mass),'raw_best_4g_m{}[6]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[8]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[9]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[10]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[11]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[12]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[13]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_dxy_m{}'.format(mass),'raw_best_4g_m{}[14]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_mass_m{}'.format(mass),'raw_best_4g_m{}[16]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_mass_m{}'.format(mass),'raw_best_4g_m{}[17]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_uncorr_mass_m{}'.format(mass),'raw_best_4g_m{}[18]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_corr_mass_m{}'.format(mass),'raw_best_4g_m{}[19]'.format(mass))
        
    #blinding region for data samples only
    ggH4g=ggH4g.Define('non_MC_cut','sample_isMC==0 && best_4g_uncorr_mass_m30<90|best_4g_uncorr_mass_m30>150')
    ggH4g=ggH4g.Filter('sample_isMC==1 | non_MC_cut==1','blinding data samples')
    
    
    actions.append(ggH4g.Snapshot('ggH4g','ggH4g.root',"best_4g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH3g.Snapshot('ggH3g','ggH3g.root',"best_3g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH3g_antiID.Snapshot('ggH3g_antiID','ggH3g_antiID.root',"best_3g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))

    #r=ggH.Report()
    #r.Print()

    return actions

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
        zmm2g=zmm2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})'.format(float(mass)))
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
        zmm3g=zmm3g.Define('raw_best_3g_m{}'.format(mass),'best_3gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})'.format(float(mass)))
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
        zmm4g=zmm4g.Define('raw_best_4g_m{}'.format(mass),'best_4gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})'.format(float(mass)))
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
#    r=zmm2g.Report()
#    r.Print()
    return actions
    

def analysis(data):
    phi_mass=[5,10,20,30]

    actions = []
    actions.append(zmumuH(data,phi_mass))
    actions.append(zeeH(data,phi_mass))
    actions.append(wmunuH(data,phi_mass))
    actions.append(wenuH(data,phi_mass))
    actions.append(ggH(data,phi_mass))

    

    
