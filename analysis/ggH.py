import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data


cols = "best_2g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Electron_.*|^Pileup_.*"

def photonAna(dataframe):

    # Overlap with loose leptons
    photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[loose_muon], Muon_eta[loose_muon])")
    photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron])")
    photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")

    # Photon Preselection criteria
    photons = photons.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&!Photon_overlap&&(Photon_isScEtaEE||Photon_isScEtaEB)")

    # Common Photon ID definitions (No isolation)
    photons = photons.Define("Photon_IdNoIso","((Photon_isScEtaEB&&Photon_hoe<0.04596&&Photon_sieie<0.0106)||(Photon_isScEtaEE&&Photon_hoe<0.0590&&Photon_sieie<0.0272))")

    return photons    

def ggH(data,phi_mass,sample):
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
    ggH=ggH.Filter("Sum(loose_electron==1)==0",'electron veto')
# Implement matching only for signal efficiency studies    
#    ggH=ggH.Define('genPhotonDR','minMatchDR(Photon_eta,Photon_phi,GenIsolatedPhoton_eta,GenIsolatedPhoton_phi)')
#    ggH=ggH.Filter('Sum(genPhotonDR<0.2)>2','At least three matched photons')

    #Correct Photon Isolation for photons
    ggH=ggH.Define("Photon_corrIso","correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,Photon_pfRelIso03_all,Photon_IdNoIso)")
    
    #Define ID and isolation
    ggH=ggH.Define("Photon_ID","Photon_IdNoIso==1 &&Photon_corrIso<0.1")
    ggH=ggH.Define("Photon_antiID","Photon_ID==0 && (Photon_isScEtaEB|Photon_isScEtaEE)")
        
    #Anti-photon RDF
    ggH_antiID = ggH.Filter("Sum(Photon_ID==1)>1 && Sum(Photon_antiID==1)>0 ","at least 2 good and 1 bad photon")

    #At least three Photons
    ggH=ggH.Filter('Sum(Photon_IdNoIso)>2','At least 3 ID No Iso Photon')
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

    actions.append(ggH4g.Snapshot('Events',sample+'_ggH4g.root',"gen.*|best_4g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH3g.Snapshot('Events',sample+'_ggH3g.root',"gen.*|best_3g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH3g_antiID.Snapshot('Events',sample+'_ggH3g_antiID.root',"gen.*|best_3g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))
    actions.append(ggH4g_antiID.Snapshot('Events',sample+'_ggH4g_antiID.root',"gen.*|best_4g.*|sample_.*|.*LHE.*|Pileup.*|^PV.*|run|event|luminosity|Block|genWeight"))

    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+"_ggH4g.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, sample+"_ggH3g.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, sample+"_ggH3g_antiID.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, sample+"_ggH4g_antiID.root", "", opts))

    #r=ggH.Report()
    #r=ggH4g_antiID.Report()
    #r.Print()

    return actions

def analysis(data,sample):
    phi_mass=[7,15,20,30,40,50,55]
    actions.extend(ggH(data,phi_mass,sample))
    return actions

