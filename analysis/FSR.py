import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')
opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data

cols = "sample_.*|^Photon_.*|^Muon_.*|^Z.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*"

muTrig = {'2018': [{'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_subleading", 'bits': 16, 'pt': 8},
                   {'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_leading", 'bits': 16, 'pt': 17}],
          '2017': [{'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_subleading", 'bits': 16, 'pt': 8},
                   {'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_leading", 'bits': 16, 'pt': 17}],
          '2016postVFP': [{'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_subleading", 'bits': 16, 'pt': 8},
                          {'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_leading", 'bits': 16, 'pt': 17}],
          '2016preVFP': [{'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_subleading", 'bits': 16, 'pt': 8},
                         {'name' : "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_leading", 'bits': 16, 'pt': 17}],          
          }

def muonAna(dataframe, era = '2018'):

    muons = dataframe.Define("loose_muon", "Muon_looseId==1&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pt>10&&Muon_pfIsoId>1")
    muons = muons.Define("med_muon", "loose_muon&&Muon_mediumId&&Muon_pfIsoId>2")
    muons = muons.Define("veto_muon", "Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&(loose_muon==0)")
    muons = muons.Define("Muon_nloose", "Sum(loose_muon)")
    muons = muons.Define("Muon_nmedium", "Sum(med_muon)")
    muons = muons.Define("Muon_nveto", "Sum(veto_muon)")
    for trigger in muTrig[era]:
        muons = muons.Define("Muon_pass{}".format(trigger['name']), "matchTrigger(Muon_eta, Muon_phi, Muon_pdgId, TrigObj_eta, TrigObj_phi, TrigObj_pt, TrigObj_id, TrigObj_filterBits, {}, {})".format(trigger['bits'], trigger['pt']))

    return muons


def photonAna(dataframe):

    # Photon Preselection criteria
    photons = dataframe.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&(Photon_isScEtaEE||Photon_isScEtaEB)")

    # Common Photon ID definitions (No isolation)
    photons = photons.Define("Photon_IdNoIso","((Photon_isScEtaEB&&Photon_hoe<0.04596&&Photon_sieie<0.0106)||(Photon_isScEtaEE&&Photon_hoe<0.0590&&Photon_sieie<0.0272))")

    return photons


def makeZ_fsr(dataframe, lepton):
    Zgs = dataframe.Define("Zg_idx", "best_zg({L}_pt, {L}_eta, {L}_phi, {L}_mass, {L}_charge, med_{l}, Photon_pt, Photon_eta, Photon_phi)".format(L=lepton, l=lepton.lower()))
    Zgs = Zgs.Define("bestZg_info", "best_Zg_info({L}_pt, {L}_eta, {L}_phi, {L}_mass, Photon_pt, Photon_eta, Photon_phi, Zg_idx)".format(L=lepton, l=lepton.lower()))
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
    Zgs = Zgs.Define("Zg_Z_mass", "bestZg_info[8]")
    return Zgs

def Zmmg(data, sample):
    actions = []
    dataframe = load_meta_data(data)
    zmm = dataframe['Events'].Filter('HLT_passed', "passed HLT")
    if data['isMC']:
        zmm = zmm.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{}, sample_isMC)".format(data['era']))
    zmm = muonAna(zmm, data['era'])
    zmm = photonAna(zmm)

    zmmg = zmm.Filter("nPhoton>0", "At least 1 photon")
    zmmg = zmmg.Filter("Sum(Muon_charge[med_muon] == 1) > 0 && Sum(Muon_charge[med_muon] == -1) > 0", "At least one opposite sign muon pair passing medium ID")
    
    zmmg = makeZ_fsr(zmmg, "Muon")

    # Applying same cuts as done for existing pixel veto here:
    # https://indico.cern.ch/event/890335/contributions/3937761/attachments/2070580/3475894/Eleveto_UL2018.pdf
    zmmg = zmmg.Filter("Zg_mass > 70 && Zg_mass < 110", "Dimuon + photon mass between 70-110 GeV")
    zmmg = zmmg.Filter("Zg_mass + Zg_Z_mass < 180", "Reject ISR photon")
    zmmg = zmmg.Filter("Zg_deltaR_l1g < 0.8 || Zg_deltaR_l2g < 0.8", "Photon close to muon")
    zmmg = zmmg.Filter("Zg_deltaR_l1g > 0.1 && Zg_deltaR_l2g > 0.1", "Photon not too close to muon")
    zmmg = zmmg.Filter("Muon_pt[Zg_l1_idx] > 20 && Muon_pt[Zg_l2_idx] > 10", "Leading/subleading muon pt > 20/10")
    actions.append(zmmg.Snapshot("Events", sample+"_zmmg.root", cols))
    actions.append(dataframe['Runs'].Snapshot("Runs", sample+"_zmmg.root", "", opts))

    return actions

def analysis(data, sample):
    actions = []
    actions.extend(Zmmg(data, sample))
    return actions
