import ROOT
from common.pyhelpers import load_meta_data
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

# Common Object ID:
def muonAna(dataframe):

    # Common Muon ID definitions (No isolation)
    muons = dataframe.Define("loose_muon", "Muon_looseId==1&&abs(Muon_eta)<2.5&&Muon_dxy<0.2&&abs(Muon_dz)<0.5&&Muon_pt>10")
    muons = muons.Define("tight_muon", "loose_muon&&Muon_tightId")

    return muons

def electronAna(dataframe):

    # Common Electron ID definitions
    electrons = dataframe.Define("loose_electrons", "Electron_pt>15&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&Electron_dxy<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
    return electrons

# Must run muon + electron analyzer first to do overlap with loose leptons
def photonAna(dataframe):

    # Overlap with loose leptons
    photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi, Muon_eta, loose_muon)")
    photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi, Electron_eta, loose_electrons)")
    photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")

    # Photon Preselection criteria
    photons = photons.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&!Photon_overlap&&(Photon_isScEtaEE||Photon_isScEtaEB)")

    # Common Photon ID definitions (No isolation)
    photons = photons.Define("Photon_IDNoIso","((Photon_isScEtaEB&&Photon_hoe<0.04596&&Photon_sieie<0.0106)||(Photon_isScEtaEE&&Photon_hoe<0.0590&&Photon_sieie<0.0272))").Filter('Sum(Photon_IDNoIso)>0','At least one ID No Iso Photon')

    return photons    

def makeZ(dataframe, lepton):
    # pT, eta, phi, mass, deltaR, deltaPhi
    Zs = dataframe.Define("Z_idx", "best_z({L}_pt[loose_{l}], {L}_eta[loose_{l}], {L}_phi[loose_{l}], {L}_mass[loose_{l}], {L}_charge[loose_{l}])".format(L=lepton, l=lepton.lower()))
    Zs = Zs.Define("bestZ_info", "best_Z_info({l}_pt, {l}_eta, {l}_phi, {l}_mass, Z_idx)".format(l=lepton))
    Zs = Zs.Define("Z_pt", "bestZ_info[0]")
    Zs = Zs.Define("Z_eta", "bestZ_info[1]")
    Zs = Zs.Define("Z_phi", "bestZ_info[2]")
    Zs = Zs.Define("Z_mass", "bestZ_info[3]")
    Zs = Zs.Define("Z_deltaR", "bestZ_info[4]")
    Zs = Zs.Define("Z_deltaPhi", "bestZ_info[5]")
    return Zs

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
    ggH = dataframe['Events'].Filter('HLT_passed','passed HLT')
    
    #your code here 
    ggH=ggH.Filter('nPhoton>2','At least three photons')

    #Filter out muons above 10Gev and electrons above 15GeV
    ggH=electronAna(ggH)
    ggH=muonAna(ggH)
    ggH=photonAna(ggH)
    
    ggH=ggH.Filter("Sum(loose_muon==1)==0",'muon veto')
    ggH=ggH.Filter("Sum(loose_electrons==1)==0",'electron veto')
# Implement matching only for signal efficiency studies    
#    ggH=ggH.Define('genPhotonDR','minMatchDR(Photon_eta,Photon_phi,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)')
#    ggH=ggH.Filter('Sum(genPhotonDR<0.2)>2','At least three matched photons')

    #Correct Photon Isolation for photons
    ggH=ggH.Define("Photon_corrIso","correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,Photon_pfRelIso03_all,Photon_IDNoIso)")
    
    #Define ID and isolation
    ggH=ggH.Define("Photon_ID","Photon_IDNoIso==1 &&Photon_corrIso<0.1")
    ggH=ggH.Define("Photon_antiID","Photon_ID==0 && (Photon_isScEtaEB|Photon_isScEtaEE)")
        
    #Anti-photon RDF
    ggH_antiID = ggH.Filter("Sum(Photon_ID==1)>1 && Sum(Photon_antiID==1)>0 ","at least 2 good and 1 bad photon")

    #At least three Photons
    ggH=ggH.Filter('Sum(Photon_IDNoIso)>2','At least 3 ID No Iso Photon')
    ggH=ggH.Filter('Sum(Photon_ID==1)>2','At least 3 ID photon')
    
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
    ggH3g_antiID=ggH_antiID.Filter('Sum(Photon_ID==1)==2 && Sum(Photon_antiID==1)==1','exactly 2 ID photon and 1 fail ID photon')
    ggH3g_antiID=ggH3g_antiID.Define('good_photons','Photon_ID==1|Photon_antiID==1')
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

    #at least 4 photons with exactly one anti ID photon
    ggH4g_antiID=ggH_antiID.Filter('Sum(Photon_ID==1)>2 && Sum(Photon_antiID==1)==1','at least 3 ID photon and exactly one fail ID photon')
    ggH4g_antiID=ggH4g_antiID.Define('good_photons','Photon_ID==1|Photon_antiID==1')
    ggH4g_antiID=ggH4g_antiID.Define("m_4g", "InvariantMass(Photon_pt[good_photons], Photon_eta[good_photons], Photon_phi[good_photons], Photon_mass[good_photons])")
    for mass in phi_mass:
        ggH4g_antiID=ggH4g_antiID.Define('raw_best_4g_m{}'.format(mass),"best_4gamma(Photon_pt[good_photons],Photon_eta[good_photons],Photon_phi[good_photons],Photon_isScEtaEB[good_photons], Photon_isScEtaEE[good_photons],{})".format(float(mass)))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[0]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[1]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[2]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[3]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[4]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[5]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_dxy_m{}'.format(mass),'raw_best_4g_m{}[6]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[8]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[9]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[10]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[11]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[12]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[13]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_dxy_m{}'.format(mass),'raw_best_4g_m{}[14]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi1_mass_m{}'.format(mass),'raw_best_4g_m{}[16]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_phi2_mass_m{}'.format(mass),'raw_best_4g_m{}[17]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_uncorr_mass_m{}'.format(mass),'raw_best_4g_m{}[18]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_corr_mass_m{}'.format(mass),'raw_best_4g_m{}[19]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_idx1_m{}'.format(mass),'raw_best_4g_m{}[20]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_idx2_m{}'.format(mass),'raw_best_4g_m{}[21]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_idx3_m{}'.format(mass),'raw_best_4g_m{}[22]'.format(mass))
        ggH4g_antiID=ggH4g_antiID.Define('best_4g_idx4_m{}'.format(mass),'raw_best_4g_m{}[23]'.format(mass))
    
    ggH4g_antiID=ggH4g_antiID.Define('best_4g_includes_failID_m30','Photon_antiID[best_4g_idx1_m30]==1|Photon_antiID[best_4g_idx2_m30]==1|Photon_antiID[best_4g_idx3_m30]==1|Photon_antiID[best_4g_idx4_m30]==1')
    ggH4g_antiID=ggH4g_antiID.Filter('best_4g_includes_failID_m30==1','fit keeps one fail ID photon')

    actions.append(ggH4g.Snapshot('Events','ggH4g.root',"gen.*|best_4g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH3g.Snapshot('Events','ggH3g.root',"gen.*|best_3g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH3g_antiID.Snapshot('Events','ggH3g_antiID.root',"gen.*|best_3g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH4g_antiID.Snapshot('Events','ggH4g_antiID.root',"gen.*|best_4g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))

    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, "ggH4g.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, "ggH3g.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, "ggH3g_antiID.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, "ggH4g_antiID.root", "", opts))

    #r=ggH.Report()
    #r=ggH4g_antiID.Report()
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
    zmm = dataframe['Events'].Filter('HLT_passed','passed HLT')
    #Apply Lepton ID (no ISO)
    zmm = muonAna(zmm)
    zmm = electronAna(zmm)
    #Look for + and - muons
    zmm = zmm.Filter("Sum(Muon_charge[loose_muon]==1)>0 && Sum(Muon_charge[loose_muon]==-1)>0","At least one opposite sign muon pair")
    #create the best Zmumu candidate and filter
    zmm = makeZ(zmm, "Muon")  
    zmm = zmm.Filter("Z_mass>70&&Z_mass<110")
    #Apply Thresholds to the muon pts
    ptThresh = 28 if data['era']=='2017' else 25
    zmm = zmm.Filter("(tight_muon[Z_idx[0]]||tight_muon[Z_idx[1]])&&(Muon_pt[Z_idx[0]]>{p}||Muon_pt[Z_idx[1]]>{p})".format(p=ptThresh))

    #require just a super cluster
    zmm = zmm.Filter("nPhoton>0","At least one super cluster")
    
    #Apply photon ID (no ISO)
    zmm = photonAna(zmm)

    # FSR Recovery with loose muons and preselection photons
    zmm=zmm.Define("Photon_isFSR","fsr_recovery(Z_idx,Muon_pt[loose_muon], Muon_eta[loose_muon], Muon_phi[loose_muon], Muon_mass[loose_muon],Photon_pt,Photon_eta,Photon_phi,Photon_preselection)");
    # Correct isolation for FSR
    zmm = zmm.Define("Muon_pfRelIso04_fsrCorr", "correct_muoniso_for_photons(Muon_pt, Muon_eta,Muon_phi,Muon_pfRelIso04_all,Photon_pt,Photon_eta,Photon_phi,Photon_isFSR)")
    zmm = zmm.Define("Photon_pfRelIso03_fsrCorr", "correct_gammaIso_for_muons(Z_idx, Muon_pt, Muon_eta, Muon_phi, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_isFSR)")
    # Defining muon ID after FSR
    zmm = zmm.Define("Muon_isLoose_corr", "loose_muon&&Muon_pfRelIso04_fsrCorr<.25")
    zmm = zmm.Define("Muon_isTight_corr", "tight_muon&&Muon_pfRelIso04_fsrCorr<.15")
    #Now cut on the ll+gamma mass around the Z 
    zmm=zmm.Define("llgamma_mass","calculate_llgamma_mass(Z_idx,Muon_pt[loose_muon], Muon_eta[loose_muon], Muon_phi[loose_muon], Muon_mass[loose_muon],Photon_pt,Photon_eta,Photon_phi, Photon_isFSR)")
    
    # Define loose photons: passing preselection and not FSR tagged
    zmm = zmm.Define("Photon_isLoose", "Photon_preselection==1&&Photon_isFSR==0")

    #At least one Photon
    zmm = zmm.Filter('Sum(Photon_isLoose==1)>0','At least one ID photon')

    ##### gamma gamma +X analysis ######
    ####################################

    #Al least Two good photons
    zmm2g=zmm.Filter('Sum(Photon_isLoose==1)>1','At least 2 good no FSR photons')
    for mass in phi_mass:
        zmm2g=zmm2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt[Photon_isLoose],Photon_eta[Photon_isLoose],Photon_phi[Photon_isLoose],Photon_isScEtaEB[Photon_isLoose], Photon_isScEtaEE[Photon_isLoose],{})'.format(float(mass)))
        zmm2g=zmm2g.Define('best_2g_gamma1_pt_m{}'.format(mass),'raw_best_2g_m{}[0]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_eta_m{}'.format(mass),'raw_best_2g_m{}[1]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_phi_m{}'.format(mass),'raw_best_2g_m{}[2]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_pt_m{}'.format(mass),'raw_best_2g_m{}[3]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_eta_m{}'.format(mass),'raw_best_2g_m{}[4]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_phi_m{}'.format(mass),'raw_best_2g_m{}[5]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        zmm2g=zmm2g.Define("Photon_corrIso_m{}".format(mass), "correct_gammaIso_for_photons(best_2g_idx1_m{m}, best_2g_idx2_m{m}, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all)".format(m=mass))
        zmm2g = zmm2g.Define("Photon_ID_m{}".format(mass), "Photon_isLoose&&Photon_corrIso_m{}<0.1&&Photon_IDNoIso".format(mass))
        zmm2g=zmm2g.Define('best_2g_sumID_m{}'.format(mass), 'Photon_ID_m{m}[best_2g_idx1_m{m}]+Photon_ID_m{m}[best_2g_idx2_m{m}]'.format(m=mass))
    
    ##### 3 gamma analysis ######
    #Targets events where one two photons are merged
    ####################################
    zmm3g=zmm.Filter('Sum(Photon_isLoose)==3','Exactly three no FSR photons')
    for mass in phi_mass:
        zmm3g=zmm3g.Define('raw_best_3g_m{}'.format(mass),'best_3gamma(Photon_pt[Photon_isLoose],Photon_eta[Photon_isLoose],Photon_phi[Photon_isLoose],Photon_isScEtaEB[Photon_isLoose], Photon_isScEtaEE[Photon_isLoose],{})'.format(float(mass)))
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
    zmm4g=zmm.Filter('Sum(Photon_isLoose)>3','At least 4 good no FSR photons')
    #Pick the best combination dependening on the mass
    for mass in phi_mass:
        zmm4g=zmm4g.Define('raw_best_4g_m{}'.format(mass),'best_4gamma(Photon_pt[Photon_isLoose],Photon_eta[Photon_isLoose],Photon_phi[Photon_isLoose],Photon_isScEtaEB[Photon_isLoose], Photon_isScEtaEE[Photon_isLoose],{})'.format(float(mass)))
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
     
        
    actions.append(zmm2g.Snapshot('Events','zmm2g.root',"best_2g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*"))
    actions.append(zmm3g.Snapshot('Events','zmm3g.root',"best_3g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*"))
    actions.append(zmm4g.Snapshot('Events','zmm4g.root',"best_4g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*"))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, 'zmm2g.root', "", opts))
        actions.append(dataframe[tree].Snapshot(tree, 'zmm3g.root', "", opts))
        actions.append(dataframe[tree].Snapshot(tree, 'zmm4g.root', "", opts))
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

    

    
