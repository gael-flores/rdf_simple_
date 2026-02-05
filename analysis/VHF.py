import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')
ROOT.gInterpreter.Declare('#include "common/signalEfficiency.h"')
opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data


cols = ".*best_2g.*|sample_.*|^Fsr.*|^PV_.*|^Photon_.*|^Electron_.*|nElectron|^Muon_.*|nMuon|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*|^MET_.*|^run.*|^luminosityBlock.*"


def muonAna(dataframe, era = '2018'):

    # Common Muon ID definitions (No isolation)
    muons = dataframe.Define("loose_muon", "(Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pfIsoId>1&&Muon_looseId>0)")
    muons = muons.Define("tight_muon", "(loose_muon&&Muon_tightId>0&&Muon_pfIsoId>3)")
    muons = muons.Define("overlap_muon", "(Muon_pt>3&&abs(Muon_eta)<2.4&&Muon_softId>0)")        

    
    muons = muons.Define("Muon_ntight", "Sum(tight_muon)")
    muons = muons.Define("Muon_nloose", "Sum(loose_muon)") 

    return muons

def electronAna(dataframe, era = '2018'):

    # Common Electron ID definitions
    electrons = dataframe.Define("Electron_scEta","Electron_eta + Electron_deltaEtaSC")
    electrons = electrons.Define("loose_electron", "Electron_pt>5&&abs(Electron_eta)<2.5&&(abs(Electron_scEta)>1.57||abs(Electron_scEta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
    electrons = electrons.Define("tight_electron", "loose_electron&&Electron_cutBased>3")
    electrons = electrons.Define("Electron_ntight", "Sum(tight_electron)")
    electrons = electrons.Define("Electron_nloose", "Sum(loose_electron)")
    
    return electrons

def photonAna(dataframe, era = '2018'):
    
    photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[overlap_muon], Muon_eta[overlap_muon],0.8,0.8)")
    photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron],0.5,0.4)")
    photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")
    
    # Photon Preselection criteria
    photons = photons.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&(!Photon_overlap)&&(Photon_isScEtaEE||Photon_isScEtaEB)&&passPhIso(Photon_vidNestedWPBitmap)")

    # Common Photon ID definitions (No isolation)
    photons = photons.Define("Photon_passCutBasedID","Photon_preselection>0 && Photon_cutBased>0")

    #For scale factors redefine them so the uncertainty is actually a variation up or down instead of just uncertainty
    photons = photons.Define("pho_SFs_id", "scaleFactors_2d(Photon_eta, Photon_pt, PHO_ID_{era}_sf, PHO_ID_{era}_binsX, PHO_ID_{era}_binsY, sample_isMC, Photon_passCutBasedID)".format(era=era))
    photons = photons.Define("Photon_idSF_val", "pho_SFs_id[0]")
    photons = photons.Define("Photon_idSF_up", "pho_SFs_id[1]+pho_SFs_id[0]")
    photons = photons.Define("Photon_idSF_down", "pho_SFs_id[0]-pho_SFs_id[1]")

    photons = photons.Define("pho_SFs_pix", "getPixelSeedSF(Photon_isScEtaEB, Photon_isScEtaEE, hasPix_UL{}_sf, sample_isMC, !Photon_pixelSeed)".format(era))
    photons = photons.Define("Photon_pixSF_val", "pho_SFs_pix[0]")
    photons = photons.Define("Photon_pixSF_up", "pho_SFs_pix[0]+pho_SFs_pix[1]")
    photons = photons.Define("Photon_pixSF_down", "pho_SFs_pix[0]-pho_SFs_pix[1]")

    photons = photons.Define("Photon_energyScaleUp", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{}_val, PHO_scaledown_{}_bins, sample_isMC)".format(era, era))
    photons = photons.Define("Photon_energyScaleDown", "photonEnergyScale(Photon_eta, Photon_seedGain, PHO_scaledown_{}_val, PHO_scaleup_{}_bins, sample_isMC)".format(era, era))
    return photons   

def genAna(dataframe):
    gen = dataframe.Define("GenPart_ctau", "getctau(GenPart_dx, GenPart_dy, GenPart_dz, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    gen = gen.Define("GenPart_ecalEta", "getGenScEta(GenPart_vx, GenPart_vy, GenPart_vz, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    gen = gen.Define("GenPart_ecalPhi", "getGenScPhi(GenPart_vx, GenPart_vy, GenPart_vz, GenPart_pt, GenPart_eta, GenPart_phi, GenPart_mass)")
    gen = gen.Define("GenPart_isSignal", "isGenSignal(GenPart_pdgId, GenPart_genPartIdxMother)")
    gen = gen.Define("GenPart_nSignal", "Sum(GenPart_isSignal)")
    return gen

def jet2g(data,phi_mass,sample):
    print("Here!")
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    #pass HLT
    JetHT = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    JetHT = JetHT.Filter('HLT_passed','passed_HLT')
    if data['isMC']:
        JetHT = JetHT.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{},sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            JetHT = genAna(JetHT)

    JetHT = muonAna(JetHT, data['era'])
    JetHT = electronAna(JetHT, data['era'])
    JetHT = photonAna(JetHT,data['era'])

    #At least Two good photons
    jet2g = JetHT.Filter("Sum(Photon_preselection==1)>1", "at_least_2_preselection_photons")
    
    for mass in phi_mass:
        #electrons assumed to be massless
        print("running kinematic fit")
        jet2g=jet2g.Define(f'raw_best_2g_m{mass}',f'best_2gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB, Photon_isScEtaEE, Photon_preselection, Photon_cutBased, {float(mass)})')
        jet2g=jet2g.Define('best_2g_gamma1_pt_m{}'.format(mass),'raw_best_2g_m{}[0]'.format(mass))
        jet2g=jet2g.Define('best_2g_gamma1_eta_m{}'.format(mass),'raw_best_2g_m{}[1]'.format(mass))
        jet2g=jet2g.Define('best_2g_gamma1_phi_m{}'.format(mass),'raw_best_2g_m{}[2]'.format(mass))
        jet2g=jet2g.Define('best_2g_gamma1_id_m{}'.format(mass), 'raw_best_2g_m{}[13]'.format(mass)) #Photon cut based ID > 0
        jet2g=jet2g.Define('best_2g_gamma2_pt_m{}'.format(mass),'raw_best_2g_m{}[3]'.format(mass))
        jet2g=jet2g.Define('best_2g_gamma2_eta_m{}'.format(mass),'raw_best_2g_m{}[4]'.format(mass))
        jet2g=jet2g.Define('best_2g_gamma2_phi_m{}'.format(mass),'raw_best_2g_m{}[5]'.format(mass))
        jet2g=jet2g.Define('best_2g_gamma2_id_m{}'.format(mass), 'raw_best_2g_m{}[14]'.format(mass))  #photon cut based ID > 0
        jet2g=jet2g.Define('best_2g_dxy_m{}'.format(mass),'raw_best_2g_m{}[6]'.format(mass))
        jet2g=jet2g.Define('best_2g_valid_m{}'.format(mass),'raw_best_2g_m{}[7]'.format(mass))
        jet2g=jet2g.Define('best_2g_raw_mass_m{}'.format(mass),'raw_best_2g_m{}[8]'.format(mass))
        jet2g=jet2g.Define('best_2g_idx1_m{}'.format(mass),'raw_best_2g_m{}[9]'.format(mass))
        jet2g=jet2g.Define('best_2g_idx2_m{}'.format(mass),'raw_best_2g_m{}[10]'.format(mass))
        jet2g=jet2g.Define('best_2g_deltaPhi_m{}'.format(mass), 'raw_best_2g_m{}[11]'.format(mass))
        jet2g=jet2g.Define('best_2g_deltaR_m{}'.format(mass), 'raw_best_2g_m{}[12]'.format(mass))
        jet2g=jet2g.Define('best_2g_pt_m{}'.format(mass), 'raw_best_2g_m{}[15]'.format(mass))
        jet2g=jet2g.Define('best_2g_eta_m{}'.format(mass), 'raw_best_2g_m{}[16]'.format(mass))
        jet2g=jet2g.Define('best_2g_phi_m{}'.format(mass), 'raw_best_2g_m{}[17]'.format(mass))
        
        jet2g = jet2g.Define("best_2g_sumID_m{}".format(mass), f"raw_best_2g_m{mass}[13]+raw_best_2g_m{mass}[14]")

    actions.append(jet2g.Snapshot('jet2g',sample+'.root',cols,opts))
    report = ROOT.RDataFrame(1)
    r = jet2g.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_jet2g", sample+'.root', "", opts))
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+'.root', "", opts))

    return actions

def gamma(data,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    ####################
    #ANALYSIS CODE HERE#        
    ####################
    gamma = dataframe['Events'].Filter("isGoodLumi", "passed_lumiFilter")
    gamma = gamma.Filter('HLT_passed', 'passed_HLT')
    
    if data['isMC']:
        gamma = gamma.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{}, sample_isMC)".format(data['era']))
        if data['customNanoAOD']:
            gamma = genAna(gamma)
    gamma = muonAna(gamma, data['era'])
    gamma = electronAna(gamma, data['era'])

    gamma = photonAna(gamma, data['era'])
    gamma = gamma.Filter("Sum(Photon_preselection==1)>0", "at_least_1_preselected_photon")

    actions.append(gamma.Snapshot('gamma', sample+".root", cols,opts))
    report = ROOT.RDataFrame(1)
    r = gamma.Report()
    for cut in r:
        report = report.Define("report_{}_all".format(cut.GetName()), "{}".format(cut.GetAll()))
        report = report.Define("report_{}_pass".format(cut.GetName()), "{}".format(cut.GetPass()))
    actions.append(report.Snapshot("Report_gamma", sample+'.root', "", opts))

    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, sample+".root", "", opts))        
    return actions

def analysis(data,sample):
    phi_mass=[20,50]
    actions = []
    actions.extend(jet2g(data,phi_mass,sample))
    actions.extend(gamma(data,sample))
    return actions

