import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')
opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data

cols = "best_2g.*|sample_.*|^Photon_.*|^Muon_.*|^Z.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Electron_.*|^Pileup_.*"

# Muon trigger[era][par], par = ['name', 'bits', 'pt']
# Name = branch name in tree
# Bits = trigger bits for HLT path
# pt = pt threshold for trigger
muTrig = {'2018': [{'name': 'HLT_IsoMu24', 'bits': 2+8, 'pt': 24}],
          '2017': [{'name': 'HLT_IsoMu27', 'bits': 2+8, 'pt': 27}],
          '2016postVFP': [{'name': 'HLT_IsoMu24', 'bits': 2+8, 'pt': 24},
                          {'name': 'HLT_IsoTkMu24', 'bits': 1+8, 'pt': 24}],
          '2016preVFP': [{'name': 'HLT_IsoMu24', 'bits': 2+8, 'pt': 24},
                          {'name': 'HLT_IsoTkMu24', 'bits': 1+8, 'pt': 24}]}

eleTrig = {'2018': [{'name': 'HLT_Ele32_WPTight_Gsf', 'bits': 2, 'pt': 32}],
           '2017': [{'name': 'HLT_Ele32_WPTight_Gsf', 'bits': 2+1024, 'pt': 32}],
           '2016': [{'name': 'HLT_Ele27_WPTight_Gsf', 'bits': 2, 'pt': 27}]}
          
# Common Object ID:
def muonAna(dataframe, era = '2018'):

    # Common Muon ID definitions (No isolation)
    muons = dataframe.Define("loose_muon", "Muon_looseId==1&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pt>10&&Muon_pfIsoId>1")
    muons = muons.Define("tight_muon", "loose_muon&&Muon_tightId&&Muon_pfIsoId>3")
    muons = muons.Define("veto_muon", "Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&(loose_muon==0)&&(tight_muon==0)")
    muons = muons.Define("Muon_nloose", "Sum(loose_muon)")
    muons = muons.Define("Muon_ntight", "Sum(tight_muon)")
    muons = muons.Define("Muon_nveto", "Sum(veto_muon)")

    for trigger in muTrig[era]:
        muons = muons.Define("Muon_pass{}".format(trigger['name']), "matchTrigger(Muon_eta, Muon_phi, Muon_pdgId, TrigObj_eta, TrigObj_phi, TrigObj_pt, TrigObj_id, TrigObj_filterBits, {}, {})".format(trigger['bits'], trigger['pt']))
    muons = muons.Define("Muon_isTrigger", "||".join(["Muon_pass{}".format(trig['name']) for trig in muTrig[era]]))

    muons = muons.Define("mu_SFs_reco", 'scaleFactors_2d(abs(Muon_eta), Muon_pt, MU_RECO_{era}_sf, MU_RECO_{era}_binsX, MU_RECO_{era}_binsY, sample_isMC, Muon_pt>10)'.format(era=era))
    muons = muons.Define("Muon_recoSF_val", "mu_SFs_reco[0]")
    muons = muons.Define("Muon_recoSF_unc", "mu_SFs_reco[1]")

    muons = muons.Define("mu_SFs_id", 'scaleFactors_2d(abs(Muon_eta), Muon_pt, MU_ID_{era}_sf, MU_ID_{era}_binsX, MU_ID_{era}_binsY, sample_isMC, tight_muon)'.format(era=era))
    muons = muons.Define("Muon_idSF_val", "mu_SFs_id[0]")
    muons = muons.Define("Muon_idSF_unc", "mu_SFs_id[1]")

    muons = muons.Define("mu_SFs_iso", 'scaleFactors_2d(abs(Muon_eta), Muon_pt, MU_ISO_{era}_sf, MU_ISO_{era}_binsX, MU_ISO_{era}_binsY, sample_isMC, tight_muon)'.format(era=era))
    muons = muons.Define("Muon_isoSF_val", "mu_SFs_iso[0]")
    muons = muons.Define("Muon_isoSF_unc", "mu_SFs_iso[1]")

    muons = muons.Define("mu_SFs_trig", 'scaleFactors_3d(Muon_charge, Muon_eta, Muon_pt, MU_TRIG_{era}_sf, MU_TRIG_{era}_binsX, MU_TRIG_{era}_binsY, MU_TRIG_{era}_binsZ, sample_isMC, Muon_isTrigger)'.format(era=era))
    muons = muons.Define("Muon_trigSF_val", "mu_SFs_trig[0]")
    muons = muons.Define("Muon_trigSF_unc", "mu_SFs_trig[1]")
    

    return muons

def electronAna(dataframe, era = '2018'):

    # Common Electron ID definitions
    electrons = dataframe.Define("loose_electron", "Electron_pt>15&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
    electrons = electrons.Define("tight_electron", "loose_electron&&Electron_cutBased>3")
    electrons = electrons.Define("veto_electron", "Electron_pt>55&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&(tight_electron==0)&&(loose_electron==0)")
    electrons = electrons.Define("Electron_nloose", "Sum(loose_electron)")
    electrons = electrons.Define("Electron_ntight", "Sum(tight_electron)")
    electrons = electrons.Define("Electron_nveto", "Sum(veto_electron)")

    for trigger in eleTrig[era]:
        electrons = electrons.Define("Electron_pass{}".format(trigger['name']), "matchTrigger(Electron_eta, Electron_phi, Electron_pdgId, TrigObj_eta, TrigObj_phi, TrigObj_pt, TrigObj_id, TrigObj_filterBits, {}, {})".format(trigger['bits'], trigger['pt']))
    electrons = electrons.Define("Electron_isTrigger", "||".join(["Electron_pass{}".format(trig['name']) for trig in eleTrig[era]]))

    electrons = electrons.Define("ele_SFs_id", 'scaleFactors_2d(Electron_eta, Electron_pt, ELE_ID_{era}_sf, ELE_ID_{era}_binsX, ELE_ID_{era}_binsY, sample_isMC, tight_electron)'.format(era=era))
    electrons = electrons.Define("Electron_idSF_val", "ele_SFs_id[0]")
    electrons = electrons.Define("Electron_idSF_unc", "ele_SFs_id[1]")

    electrons = electrons.Define("ele_SFs_reco", 'scaleFactors_eleReco(Electron_eta, Electron_eta, ELE_RECO_ptBelow20_{era}_sf, ELE_RECO_ptBelow20_{era}_binsX, ELE_RECO_ptBelow20_{era}_binsY,  ELE_RECO_ptAbove20_{era}_sf, ELE_RECO_ptAbove20_{era}_binsX, ELE_RECO_ptAbove20_{era}_binsY, sample_isMC, tight_electron)'.format(era=era))
    electrons = electrons.Define("Electron_recoSF_val", "ele_SFs_id[0]")
    electrons = electrons.Define("Electron_recoSF_unc", "ele_SFs_id[1]")
    
    electrons = electrons.Define("ele_SFs_trig", "scaleFactors_2d(Electron_eta, Electron_pt, ELE_TRIG_{era}_sf, ELE_TRIG_{era}_binsX, ELE_TRIG_{era}_binsY, sample_isMC, Electron_pt>35)".format(era=era))
    electrons = electrons.Define("Electron_trigSF_val", "ele_SFs_trig[0]")
    electrons = electrons.Define("Electron_trigSF_unc", "ele_SFs_trig[1]")

    return electrons

# Must run muon + electron analyzer first to do overlap with loose leptons
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

def makeZ(dataframe, lepton):
    # pT, eta, phi, mass, deltaR, deltaPhi
    Zs = dataframe.Define("Z_idx", "best_z({L}_pt[tight_{l}], {L}_eta[tight_{l}], {L}_phi[tight_{l}], {L}_mass[tight_{l}], {L}_charge[tight_{l}])".format(L=lepton, l=lepton.lower()))
    Zs = Zs.Define("bestZ_info", "best_Z_info({L}_pt[tight_{l}], {L}_eta[tight_{l}], {L}_phi[tight_{l}], {L}_mass[tight_{l}], Z_idx)".format(L=lepton, l=lepton.lower()))
    Zs = Zs.Define("Z_pt", "bestZ_info[0]")
    Zs = Zs.Define("Z_eta", "bestZ_info[1]")
    Zs = Zs.Define("Z_phi", "bestZ_info[2]")
    Zs = Zs.Define("Z_mass", "bestZ_info[3]")
    Zs = Zs.Define("Z_deltaR", "bestZ_info[4]")
    Zs = Zs.Define("Z_deltaPhi", "bestZ_info[5]")

    return Zs

def makeZ_fsr(dataframe, lepton):
    Zgs = dataframe.Define("Zg_idx", "best_zg({L}_pt[tight_{l}], {L}_eta[tight_{l}], {L}_phi[tight_{l}], {L}_mass[tight_{l}], {L}_charge[tight_{l}], Photon_pt, Photon_eta, Photon_phi)".format(L=lepton, l=lepton.lower()))
    Zgs = Zgs.Define("bestZg_info", "best_Zg_info({L}_pt[tight_{l}], {L}_eta[tight_{l}], {L}_phi[tight_{l}], {L}_mass[tight_{l}], Photon_pt, Photon_eta, Photon_phi, Zg_idx)".format(L=lepton, l=lepton.lower()))
    Zgs = Zgs.Define("Zg_pt", "bestZg_info[0]")
    Zgs = Zgs.Define("Zg_eta", "bestZg_info[1]")
    Zgs = Zgs.Define("Zg_phi", "bestZg_info[2]")
    Zgs = Zgs.Define("Zg_mass", "bestZg_info[3]")
    Zgs = Zgs.Define("Zg_deltaR_l1g", "bestZg_info[4]")
    Zgs = Zgs.Define("Zg_deltaR_l2g", "bestZg_info[5]")
    Zgs = Zgs.Define("Zg_deltaPhi_l1g", "bestZg_info[6]")
    Zgs = Zgs.Define("Zg_deltaPhi_l2g", "bestZg_info[7]")
    Zgs = Zgs.Define("Zg_l1_idx", "Zg_idx[0]")
    Zgs = Zgs.Define("Zg_l2_idx", "Zg_idx[1]")
    Zgs = Zgs.Define("Zg_g_idx", "Zg_idx[2]")
    
    return Zgs

def zeeH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################

    return actions

def wenuH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################

    return actions

def wmunuH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################

    return actions



def zmumuH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)

    ####################
    #ANALYSIS CODE HERE#        
    ####################

    #pass HLT
    zmm = dataframe['Events'].Filter('HLT_passed','passed HLT')
    zmm = zmm.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{},sample_isMC)".format(data['era']))
    #Apply Lepton ID (no ISO)
    zmm = muonAna(zmm, data['era'])
    zmm = electronAna(zmm, data['era'])
    #Look for + and - muons
    zmm = zmm.Filter("Sum(Muon_charge[tight_muon]==1)>0 && Sum(Muon_charge[tight_muon]==-1)>0","At least one opposite sign muon pair, both passing tight ID and preselection")

    # Separate dataframe for Z FSR tagging for photon ID
    zmmg = zmm.Filter("nPhoton==1", "Exactly 1 photon")
    zmmg = photonAna(zmmg) # Run photon analyzer but no cuts on ID
    zmmg = makeZ_fsr(zmmg, "Muon")
    zmmg = zmmg.Filter("Zg_mass > 70 && Zg_mass < 110", "Dimuon + photon mass between 70-110 GeV")
    actions.append(zmmg.Snapshot("Events", sample+"_zmmg.root", cols))
    actions.append(dataframe['Runs'].Snapshot("Runs", sample+"_zmmg.root", "", opts))

    #create the best Zmumu candidate and filter
    zmm = makeZ(zmm, "Muon")  
    zmm = zmm.Filter("Z_mass>70&&Z_mass<110", "Dimuon mass between 70-110 GeV")  

    #Apply Thresholds to the muon pts and cut on muon pf iso
    ptThresh = 28 if data['era']=='2017' else 25
    zmm = zmm.Filter("(Muon_pt[Z_idx[0]]>{p}||Muon_pt[Z_idx[1]]>{p})".format(p=ptThresh), "At least one muon in Z with pt>{}".format(ptThresh))

    #require just a super cluster
    zmm = zmm.Filter("nPhoton>0","At least one photon")
    
    #Apply photon ID (no ISO)
    zmm = photonAna(zmm)
    #actions.append(zmm.Snapshot("Events", "zmmg.root", cols, opts))

    # FSR Recovery with loose muons and preselection photons
    zmm=zmm.Define("Photon_isFSR","fsr_recovery(Z_idx,Muon_pt, Muon_eta, Muon_phi, Muon_mass,Photon_pt,Photon_eta,Photon_phi,Photon_preselection)");
    # Correct isolation for FSR
    zmm = zmm.Define("Photon_pfRelIso03_fsrCorr", "correct_gammaIso_for_muons(Z_idx, Muon_pt, Muon_eta, Muon_phi, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_isFSR)")
    
    
    #Now cut on the ll+gamma mass around the Z 
    zmm=zmm.Define("llgamma_mass","calculate_llgamma_mass(Z_idx,Muon_pt[loose_muon], Muon_eta[loose_muon], Muon_phi[loose_muon], Muon_mass[loose_muon],Photon_pt,Photon_eta,Photon_phi, Photon_isFSR)")

    #zmm = zmm.Filter("Sum(Photon_preselection==1)>1", "At least 2 photons passing preselection")

    # Define loose photons: passing preselection and not FSR tagged
    zmm = zmm.Define("Photon_isLoose", "Photon_preselection==1&&Photon_isFSR==0")

    ##### gamma gamma +X analysis ######
    ####################################

    #Al least Two good photons
    zmm2g=zmm.Filter('Sum(Photon_isLoose==1)>1','At least 2 non FSR photons passing preselection')
    for mass in phi_mass:
        zmm2g=zmm2g.Define('raw_best_2g_m{}'.format(mass),'best_2gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB, Photon_isScEtaEE, Photon_isLoose, Photon_IdNoIso, Photon_pfRelIso03_fsrCorr, {})'.format(float(mass)))
        zmm2g=zmm2g.Define('best_2g_gamma1_pt_m{}'.format(mass),'raw_best_2g_m{}[0]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_eta_m{}'.format(mass),'raw_best_2g_m{}[1]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_phi_m{}'.format(mass),'raw_best_2g_m{}[2]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma1_id_m{}'.format(mass), 'raw_best_2g_m{}[13]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_pt_m{}'.format(mass),'raw_best_2g_m{}[3]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_eta_m{}'.format(mass),'raw_best_2g_m{}[4]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_phi_m{}'.format(mass),'raw_best_2g_m{}[5]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_gamma2_id_m{}'.format(mass), 'raw_best_2g_m{}[14]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_deltaPhi_m{}'.format(mass), 'raw_best_2g_m{}[11]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_deltaR_m{}'.format(mass), 'raw_best_2g_m{}[12]'.format(mass))
        zmm2g=zmm2g.Define('best_2g_pt_m{}'.format(mass), 'raw_best_2g_m{}[15]'.format(mass))
        zmm2g=zmm2g.Define("Photon_corrIso_m{}".format(mass), "correct_gammaIso_for_photons(best_2g_idx1_m{m}, best_2g_idx2_m{m}, Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all)".format(m=mass))
        zmm2g = zmm2g.Define("Photon_ID_m{}".format(mass), "Photon_isLoose&&Photon_corrIso_m{}<0.1&&Photon_IdNoIso".format(mass))
        zmm2g = zmm2g.Define("best_2g_sumID_m{}".format(mass), "raw_best_2g_m{m}[13]+raw_best_2g_m{m}[14]".format(m=mass))
        #zmm2g=zmm2g.Define('best_2g_sumID_m{}'.format(mass), 'Photon_ID_m{m}[best_2g_idx1_m{m}]+Photon_ID_m{m}[best_2g_idx2_m{m}]'.format(m=mass))
    
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
     
        

    actions.append(zmm2g.Snapshot('Events',sample+'_zmm2g.root',cols))
    actions.append(zmm3g.Snapshot('Events',sample+'_zmm3g.root',cols))
    actions.append(zmm4g.Snapshot('Events',sample+'_zmm4g.root',cols))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+'_zmm2g.root', "", opts))
        actions.append(dataframe[tree].Snapshot(tree, sample+'_zmm3g.root', "", opts))
        actions.append(dataframe[tree].Snapshot(tree, sample+'_zmm4g.root', "", opts))
#    r=zmm2g.Report()
#    r.Print()

    return actions
    


def analysis(data,sample):
    #phi_mass=[5,10,20,30]
    phi_mass=[7,15,20,30,40,50,55]
    actions = []
    actions.extend(zmumuH(data,phi_mass,sample))
    actions.extend(zeeH(data,phi_mass,sample))
    actions.extend(wmunuH(data,phi_mass,sample))
    actions.extend(wenuH(data,phi_mass,sample))
    return actions

    

    
