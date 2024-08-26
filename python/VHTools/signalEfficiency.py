import ROOT, sys, glob
sys.path.append("/home/tyler/DDP/rdf_simple/")
ROOT.gInterpreter.Declare('#include "interface/signalEfficiency.h"')
ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
from analysis.VH import *
ROOT.ROOT.EnableImplicitMT()

actions = []
histProxies = []

for v in ['W', 'Z']:
    for m in [15, 20, 30, 40, 50, 55]:
        rdf = ROOT.RDataFrame("Events", "/home/tyler/samples/DDP/rdf_samples/signalSamples/2018/*{}*H_HTo2LongLivedTo*_MH-125_MFF-{}*.root".format(v,m)) #Input uses raw signal nanoAOD samples (not processed using analyzer)
        rdf = rdf.Define("sample_isMC", "1")
        rdf = genAna(rdf)
        rdf = muonAna(rdf, era="2018")
        rdf = electronAna(rdf, era="2018")
        rdf = photonAna(rdf, era='2018')
        rdf = rdf.Redefine("Photon_preselection", "abs(Photon_eta)<2.5")
        rdf = rdf.Filter("GenPart_nSignal>0")
        rdf = rdf.Define("GenPart_dxy", "sqrt(GenPart_vx*GenPart_vx+GenPart_vy*GenPart_vy)") # Vertex Dxy
        rdf = rdf.Define("GenPart_decayDxy", "sqrt(GenPart_dx*GenPart_dx+GenPart_dy*GenPart_dy)") # Decay length dxy
        rdf = rdf.Define("Photon_pfRelIso03_corr", "correct_gammaIso(Photon_pt, Photon_eta, Photon_phi, Photon_pfRelIso03_all, Photon_preselection)")
        rdf = rdf.Define("Photon_id_corr", "Photon_IdNoIso&&Photon_pfRelIso03_corr<0.1")
        rdf = rdf.Define("Photon_looseId_noPhIso", "passIDNoIso(Photon_vidNestedWPBitmap)")
        # Define good gen signal photons:
        #   Decay inside ecal (dxy, dz)
        #   Decays to photons
        #   pt > 20
        #   |eta_sc| < 2.5
        for var in ['pt', 'eta', 'phi', 'mass', 'genPartIdxMother', 'ecalEta', 'ecalPhi']:
            rdf = rdf.Define("GenSignal_{}".format(var), "GenPart_{}[GenPart_isSignal&&GenPart_dxy<140&&abs(GenPart_vx)<330&&abs(GenPart_ecalEta)<2.5&&GenPart_pt>20]".format(var))

            
        rdf = rdf.Define("nGenSignal", "Sum(GenSignal_pt>0)")
        rdf = rdf.Define("GenSignal_energy", "getEnergy(GenSignal_pt, GenSignal_eta, GenSignal_phi, GenSignal_mass)")
        rdf = rdf.Define("GenPart_decayToSignal", "decayToSignal(GenPart_pdgId, GenSignal_genPartIdxMother)")

        rdf = rdf.Define("Photon_energy", "getEnergy(Photon_pt, Photon_eta, Photon_phi, Photon_mass)")

        rdf = rdf.Define("GenSignal_recoIdx_all", "matchSignal(GenSignal_energy, GenSignal_ecalEta, GenSignal_ecalPhi, Photon_energy, Photon_eta, Photon_phi, 0.3)")
        rdf = rdf.Define("GenPart_nMatchedDaughters_all", "pairMatchedToReco(GenPart_pdgId, GenSignal_recoIdx_all, GenSignal_genPartIdxMother)")

        rdf = rdf.Define("GenSignal_recoIdx_preselection", "matchSignal(GenSignal_energy, GenSignal_ecalEta, GenSignal_ecalPhi, Photon_energy[Photon_preselection], Photon_eta[Photon_preselection], Photon_phi[Photon_preselection], 0.3)")
        rdf = rdf.Define("GenPart_nMatchedDaughters_preselection", "pairMatchedToReco(GenPart_pdgId, GenSignal_recoIdx_preselection, GenSignal_genPartIdxMother)")


        rdf = rdf.Define("GenSignal_recoIdx_looseId", "matchSignal(GenSignal_energy, GenSignal_ecalEta, GenSignal_ecalPhi, Photon_energy[Photon_cutBased>0&&Photon_preselection], Photon_eta[Photon_cutBased>0&&Photon_preselection], Photon_phi[Photon_cutBased>0&&Photon_preselection], 0.3)")
        rdf = rdf.Define("GenPart_nMatchedDaughters_looseId", "pairMatchedToReco(GenPart_pdgId, GenSignal_recoIdx_looseId, GenSignal_genPartIdxMother)")

        rdf = rdf.Define("GenSignal_recoIdx_idNoIso", "matchSignal(GenSignal_energy, GenSignal_ecalEta, GenSignal_ecalPhi, Photon_energy[Photon_looseId_noPhIso&&Photon_preselection], Photon_eta[Photon_looseId_noPhIso&&Photon_preselection], Photon_phi[Photon_looseId_noPhIso&&Photon_preselection], 0.3)")
        rdf = rdf.Define("GenPart_nMatchedDaughters_idNoIso", "pairMatchedToReco(GenPart_pdgId, GenSignal_recoIdx_idNoIso, GenSignal_genPartIdxMother)")

        rdf = rdf.Define("GenSignal_recoIdx_id_corr", "matchSignal(GenSignal_energy, GenSignal_ecalEta, GenSignal_ecalPhi, Photon_energy[Photon_id_corr&&Photon_preselection], Photon_eta[Photon_id_corr&&Photon_preselection], Photon_phi[Photon_id_corr&&Photon_preselection], 0.3)")
        rdf = rdf.Define("GenPart_nMatchedDaughters_id_corr", "pairMatchedToReco(GenPart_pdgId, GenSignal_recoIdx_id_corr, GenSignal_genPartIdxMother)")


        rdf = rdf.Define("GenPhi_decayDxy", "GenPart_decayDxy[GenPart_decayToSignal==2]")

        for var in ['HOE', 'Sieie', 'ChIso', 'NeuIso', 'PhIso']:
            rdf = rdf.Define("Photon_pass{}".format(var), "pass{}(Photon_vidNestedWPBitmap)".format(var))
            rdf = rdf.Define("GenSignal_recoIdx_pass{}".format(var), "matchSignal(GenSignal_energy, GenSignal_ecalEta, GenSignal_ecalPhi, Photon_energy[Photon_pass{var}&&Photon_preselection], Photon_eta[Photon_pass{var}&&Photon_preselection], Photon_phi[Photon_pass{var}&&Photon_preselection], 0.3)".format(var=var))
            rdf = rdf.Define("GenPart_nMatchedDaughters_pass{}".format(var), "pairMatchedToReco(GenPart_pdgId, GenSignal_recoIdx_pass{}, GenSignal_genPartIdxMother)".format(var))
            rdf = rdf.Define("GenPhi_nMatchedDaughters_pass{}".format(var), "GenPart_nMatchedDaughters_pass{}[GenPart_decayToSignal==2]".format(var))

        for var in ['all', 'preselection', 'looseId', 'idNoIso', 'id_corr']:
            rdf = rdf.Define("GenPhi_nMatchedDaughters_{}".format(var), "GenPart_nMatchedDaughters_{}[GenPart_decayToSignal==2]".format(var))

        histProxies.append(rdf.Histo1D(("denom_{}_m{}".format(v, m), "", 20, 0, 100), "GenPhi_decayDxy"))
        for var in ['all', 'preselection', 'looseId', 'idNoIso', 'id_corr', 'passHOE', 'passSieie', 'passChIso', 'passNeuIso', 'passPhIso']:
            rdf = rdf.Define("GenPhi_dxy_{}".format(var), "GenPhi_decayDxy[GenPhi_nMatchedDaughters_{}>1]".format(var))
            histProxies.append(rdf.Histo1D(("num_{}_m{}_{}".format(v,m,var), "", 20, 0, 100), "GenPhi_dxy_{}".format(var)))

ROOT.RDF.RunGraphs(histProxies)

histos = [histoproxy.GetValue() for histoproxy in histProxies]

f = ROOT.TFile("photonID.root", "RECREATE")

for h in histos:
    h.Write()

f.Close()
