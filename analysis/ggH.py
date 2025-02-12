import ROOT
ROOT.gInterpreter.Declare('#include "analysis/ddp_vertex.h"')
ROOT.gInterpreter.Declare('#include "common/scaleFactors.h"')

opts = ROOT.RDF.RSnapshotOptions()
opts.fMode = "UPDATE"
opts.fOverwriteIfExists = True

from common.pyhelpers import load_meta_data


#cols = "best_3g.*|best_4g.*|sample_.*|^Photon_.*|^Muon_.*|^Z_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Electron_.*|^Pileup_.*|^run.*"
cols = "best_.*|sample_.*|^Photon_.*|Weight.*|^Gen.*|^weight.*|^TrigObj_.*|^event.*|^Pileup_.*|^run.*|gen.*|.*LHE.*|^PV.*|luminosity|Block|genWeight|HLT_passed"

iso = {'2017': 'Photon_pfRelIso03_all',
       '2018': 'Photon_pfRelIso03_all',
       '2022': 'Photon_pfRelIso03_all',
       '2023': 'Photon_pfRelIso03_all_Fall17V2'}

#indices of raw_best_4g_m{} determined by "result" output of best_4gamma functions
#in ddp_vertex. we saved the indices of best 4 in ddp_vertex.h as 24,25,26,27
ID_dict1_ggH4g={
 'Photon_preselection':'==1'}
ID_dict2_ggH4g={
 'Photon_IdNoIso':'==1'}

ID1_ggH4g = " && ".join(f"{branch}[raw_best_4g_m{{m}}[{i}]]{b}" for branch, b in ID_dict1_ggH4g.items() for i in range(24,28))
ID2_ggH4g = " && ".join(f"{branch}[raw_best_4g_m{{m}}[{i}]]{b}" for branch, b in ID_dict2_ggH4g.items() for i in range(24,28))
ID3_ggH4g = " && ".join(f"((Photon_isScEtaEB[raw_best_4g_m{{m}}[{i}]]==1 && Photon_corrIso_m{{m}}[raw_best_4g_m{{m}}[{i}]]<0.25)||(Photon_isScEtaEE[raw_best_4g_m{{m}}[{i}]]==1 && Photon_corrIso_m{{m}}[raw_best_4g_m{{m}}[{i}]]<0.3))" for i in range(24,28))

#only preselection and ID rules change for ggH4g_1bad. EE and EB logic identical
ID1_1bad = "(((Photon_preselection[raw_best_4g_m{m}[24]])+(Photon_preselection[raw_best_4g_m{m}[25]])+(Photon_preselection[raw_best_4g_m{m}[26]])+(Photon_preselection[raw_best_4g_m{m}[27]]))==3)"
ID2_1bad = "(((Photon_IdNoIso[raw_best_4g_m{m}[24]])+(Photon_IdNoIso[raw_best_4g_m{m}[25]])+(Photon_IdNoIso[raw_best_4g_m{m}[26]])+(Photon_IdNoIso[raw_best_4g_m{m}[27]]))==3)"

#strings to pass into Define in phi mass loop
best_4g_ID_str_ggH4g_1bad = "{} && {} && {}".format(ID1_1bad, ID2_1bad, ID3_ggH4g)
best_4g_ID_str_ggH4g="{} && {} && {}".format(ID1_ggH4g, ID2_ggH4g, ID3_ggH4g)
    
# Common Object ID:
def muonAna(dataframe):

    # Common Muon ID definitions (No isolation)
    muons = dataframe.Define("loose_muon", "Muon_looseId==1&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&Muon_pt>10&&Muon_pfIsoId>1")
    muons = muons.Define("tight_muon", "loose_muon&&Muon_tightId&&Muon_pfIsoId>3")
    muons = muons.Define("veto_muon", "Muon_pt>5&&abs(Muon_eta)<2.4&&abs(Muon_dxy)<0.2&&abs(Muon_dz)<0.5&&(loose_muon==0)&&(tight_muon==0)")
    muons = muons.Define("Muon_nloose", "Sum(loose_muon)")
    muons = muons.Define("Muon_ntight", "Sum(tight_muon)")
    muons = muons.Define("Muon_nveto", "Sum(veto_muon)")
    return muons

def electronAna(dataframe):

    # Common Electron ID definitions
    electrons = dataframe.Define("loose_electron", "Electron_pt>15&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&Electron_cutBased>0")
    electrons = electrons.Define("tight_electron", "loose_electron&&Electron_cutBased>3")
    electrons = electrons.Define("veto_electron", "Electron_pt>55&&abs(Electron_eta)<2.5&&(abs(Electron_eta)>1.57||abs(Electron_eta)<1.44)&&abs(Electron_dxy)<0.2&&abs(Electron_dz)<0.2&&Electron_lostHits<2&&Electron_convVeto&&(tight_electron==0)&&(loose_electron==0)")
    electrons = electrons.Define("Electron_nloose", "Sum(loose_electron)")
    electrons = electrons.Define("Electron_ntight", "Sum(tight_electron)")
    electrons = electrons.Define("Electron_nveto", "Sum(veto_electron)")
    return electrons
    
def photonAna(dataframe):

    # Overlap with loose leptons
    photons = dataframe.Define("Photon_muOverlap", "overlapClean(Photon_phi, Photon_eta, Muon_phi[loose_muon], Muon_eta[loose_muon])")
    photons = photons.Define("Photon_eleOverlap", "overlapClean(Photon_phi, Photon_eta, Electron_phi[loose_electron], Electron_eta[loose_electron])")
    photons = photons.Define("Photon_overlap", "Photon_muOverlap||Photon_eleOverlap")

    # Photon Preselection criteria
    photons = photons.Define("Photon_preselection", "Photon_pt>20&&!Photon_pixelSeed&&abs(Photon_eta)<2.5&&(abs(Photon_eta)>1.57||abs(Photon_eta)<1.44)&&!Photon_overlap&&(Photon_isScEtaEE||Photon_isScEtaEB)")
    #photons = photons.Define("Photon_preselection", "(Photon_isScEtaEE||Photon_isScEtaEB)")

    # Common Photon ID definitions (No isolation)
    photons = photons.Define("Photon_IdNoIso","((Photon_isScEtaEB&&Photon_hoe<0.3&&Photon_sieie<0.035)||(Photon_isScEtaEE&&Photon_hoe<0.2&&Photon_sieie<0.045))")

    return photons    

def save_report(df, report_name, sample, opts, actions):
        report = ROOT.RDataFrame(1)  # Create a dummy dataframe with one entry
        r = df.Report()
        for cut in r:
            # Define new columns for total and passing entries for each cut
            report = report.Define(f"report_{cut.GetName()}_all", f"{cut.GetAll()}")
            report = report.Define(f"report_{cut.GetName()}_pass", f"{cut.GetPass()}")
        # Append the snapshot action to the actions list
        actions.append(report.Snapshot(report_name, f"{sample}.root", "", opts))

def ggH(data,phi_mass,sample):
    actions=[]

    #Declare dataframe and load all meta data 
    dataframe =load_meta_data(data)
    #pass HLT
    #ggH = dataframe['Events'].Filter('HLT_passed','passed_HLT')
    ggH=dataframe["Events"]


    if data["isMC"]:
        ggH = ggH.Define("Pileup_weight", "getPUweight(Pileup_nPU, puWeight_UL{}, sample_isMC)".format(data["era"]))

    

    #Filter out muons above 10Gev and electrons above 15GeV
    ggH=electronAna(ggH)
    ggH=muonAna(ggH)
    
    ggH=ggH.Filter("Sum(loose_muon==1)==0",'muon_veto')
    ggH=ggH.Filter("Sum(loose_electron==1)==0",'electron_veto')

    #Add photon preselection and common photon ID definitions 
    ggH=photonAna(ggH)
    
    #filtering for all snapshots done here:
    ggH3g=ggH.Filter('nPhoton>2','at_least_3_photons')
    ggH3g=ggH3g.Filter('Sum(Photon_preselection==1)==3','exactly_3_preselected_photons')
    ggH4g=ggH.Filter('nPhoton>3','at_least_4_photons')
    ggH4g=ggH4g.Filter('Sum(Photon_preselection==1)>3','at_least_3_preselected_photons')
    ggH4g_1bad=ggH.Filter('nPhoton>3','at_least_4_photons')
    ggH4g_1bad=ggH4g_1bad.Filter('Sum(Photon_preselection==1)==3','exactly_3_preselected_photons')

    for mass in phi_mass:
        ggH3g=ggH3g.Define('Photon_corrIso_m{}'.format(mass),'correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,{},Photon_preselection)'.format(iso[data['era']]))
        ggH3g=ggH3g.Define('raw_best_3g_m{}'.format(mass),"best_3gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso,Photon_corrIso_m{},{})".format(mass,float(mass)))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_pt_m{}'.format(mass),'raw_best_3g_m{}[0]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_eta_m{}'.format(mass),'raw_best_3g_m{}[1]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_phi_m{}'.format(mass),'raw_best_3g_m{}[2]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma1_id_m{}'.format(mass),'raw_best_3g_m{}[9]'.format(mass)) #id criteria: pass loose nonIsoID, corrected isolation <0.1
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_pt_m{}'.format(mass),'raw_best_3g_m{}[3]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_eta_m{}'.format(mass),'raw_best_3g_m{}[4]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_phi_m{}'.format(mass),'raw_best_3g_m{}[5]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_gamma2_id_m{}'.format(mass),'raw_best_3g_m{}[10]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_pt_m{}'.format(mass),'raw_best_3g_m{}[12]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_eta_m{}'.format(mass),'raw_best_3g_m{}[13]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_phi_m{}'.format(mass),'raw_best_3g_m{}[14]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_gamma3_id_m{}'.format(mass),'raw_best_3g_m{}[11]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_dxy_m{}'.format(mass),'raw_best_3g_m{}[6]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_valid_m{}'.format(mass),'raw_best_3g_m{}[7]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_phi_mass_m{}'.format(mass),'raw_best_3g_m{}[8]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_raw_mass_m{}'.format(mass),'raw_best_3g_m{}[15]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_corr_mass_m{}'.format(mass),'raw_best_3g_m{}[16]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_deltaPhi_m{}'.format(mass),'raw_best_3g_m{}[17]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_deltaR_m{}'.format(mass),'raw_best_3g_m{}[18]'.format(mass))
        ggH3g=ggH3g.Define('best_3g_sumID_m{}'.format(mass),'raw_best_3g_m{m}[9]+raw_best_3g_m{m}[10]+raw_best_3g_m{m}[11]'.format(m=mass))
        ggH3g=ggH3g.Define('Photon_ID_m{}'.format(mass),'Photon_preselection && Photon_IdNoIso && ((Photon_isScEtaEE==1 && Photon_corrIso_m{m}<0.3)||(Photon_isScEtaEB==1 && Photon_corrIso_m{m}<0.25))'.format(m=mass))
        ggH3g=ggH3g.Define('non_MC_cut_m{}'.format(mass),'sample_isMC==0 && best_3g_raw_mass_m{m}< 30 | best_3g_raw_mass_m{m}>140'.format(m=mass))
        ggH3g=ggH3g.Define('best_3g_idx1_m{}'.format(mass),'raw_best_3g_m{}[19]'.format(mass)) #index of photon1
        ggH3g=ggH3g.Define('best_3g_idx2_m{}'.format(mass),'raw_best_3g_m{}[20]'.format(mass)) #index of photon 2
        ggH3g=ggH3g.Define('best_3g_idx3_m{}'.format(mass),'raw_best_3g_m{}[21]'.format(mass)) #index of photon 3
        

    for mass in phi_mass:
        ggH4g=ggH4g.Define('Photon_corrIso_m{}'.format(mass),'correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,{},Photon_preselection)'.format(iso[data['era']]))
        ggH4g=ggH4g.Define('raw_best_4g_m{}'.format(mass),"best_4gamma(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso,Photon_corrIso_m{},{})".format(mass,float(mass)))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[0]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[1]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[2]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma1_id_m{}'.format(mass),'raw_best_4g_m{}[20]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[3]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[4]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[5]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_gamma2_id_m{}'.format(mass),'raw_best_4g_m{}[21]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_dxy_m{}'.format(mass),'raw_best_4g_m{}[6]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[8]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[9]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[10]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma1_id_m{}'.format(mass),'raw_best_4g_m{}[22]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[11]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[12]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[13]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_gamma2_id_m{}'.format(mass),'raw_best_4g_m{}[23]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_dxy_m{}'.format(mass),'raw_best_4g_m{}[14]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi1_mass_m{}'.format(mass),'raw_best_4g_m{}[16]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_phi2_mass_m{}'.format(mass),'raw_best_4g_m{}[17]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_uncorr_mass_m{}'.format(mass),'raw_best_4g_m{}[18]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_corr_mass_m{}'.format(mass),'raw_best_4g_m{}[19]'.format(mass))
        ggH4g=ggH4g.Define('best_4g_sumID_m{}'.format(mass),'raw_best_4g_m{m}[20]+raw_best_4g_m{m}[21]+raw_best_4g_m{m}[22]+raw_best_4g_m{m}[23]'.format(m=mass))
        ggH4g=ggH4g.Define('Photon_ID_m{}'.format(mass),'Photon_preselection && Photon_IdNoIso && ((Photon_isScEtaEE==1 && Photon_corrIso_m{m}<0.3)||(Photon_isScEtaEB==1 && Photon_corrIso_m{m}<0.25))'.format(m=mass))
        ggH4g=ggH4g.Define('non_MC_cut_m{}'.format(mass),'sample_isMC==0 && best_4g_uncorr_mass_m{m}<90|best_4g_uncorr_mass_m{m}>150'.format(m=mass))
        ggH4g=ggH4g.Define('best_4g_idx1_m{}'.format(mass),'raw_best_4g_m{}[24]'.format(mass)) #index of photon 1
        ggH4g=ggH4g.Define('best_4g_idx2_m{}'.format(mass),'raw_best_4g_m{}[25]'.format(mass)) #index of photon 2
        ggH4g=ggH4g.Define('best_4g_idx3_m{}'.format(mass),'raw_best_4g_m{}[26]'.format(mass)) #index of photon 3
        ggH4g=ggH4g.Define('best_4g_idx4_m{}'.format(mass),'raw_best_4g_m{}[27]'.format(mass)) #index of photon 4
        ggH4g=ggH4g.Define('best_4g_ID_m{}'.format(mass),best_4g_ID_str_ggH4g.format(m=mass))
        
        

    for mass in phi_mass:
        ggH4g_1bad=ggH4g_1bad.Define('Photon_corrIso_m{}'.format(mass),'correct_gammaIso(Photon_pt,Photon_eta,Photon_phi,{},Photon_preselection)'.format(iso[data['era']]))
        ggH4g_1bad=ggH4g_1bad.Define('raw_best_4g_m{}'.format(mass),"best_4gamma_1bad(Photon_pt,Photon_eta,Photon_phi,Photon_isScEtaEB,Photon_isScEtaEE,Photon_preselection,Photon_IdNoIso,Photon_corrIso_m{},{})".format(mass,float(mass)))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[0]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[1]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[2]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma1_id_m{}'.format(mass),'raw_best_4g_m{}[20]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_idx1_m{}'.format(mass),'raw_best_4g_m{}[24]'.format(mass)) #index of photon 1
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[3]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[4]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[5]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_gamma2_id_m{}'.format(mass),'raw_best_4g_m{}[21]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_idx2_m{}'.format(mass),'raw_best_4g_m{}[25]'.format(mass)) #index of photon 2
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_dxy_m{}'.format(mass),'raw_best_4g_m{}[6]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_valid_m{}'.format(mass),'raw_best_4g_m{}[7]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma1_pt_m{}'.format(mass),'raw_best_4g_m{}[8]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma1_eta_m{}'.format(mass),'raw_best_4g_m{}[9]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma1_phi_m{}'.format(mass),'raw_best_4g_m{}[10]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma1_id_m{}'.format(mass),'raw_best_4g_m{}[22]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_idx3_m{}'.format(mass),'raw_best_4g_m{}[26]'.format(mass)) #index of photon 3
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma2_pt_m{}'.format(mass),'raw_best_4g_m{}[11]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma2_eta_m{}'.format(mass),'raw_best_4g_m{}[12]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma2_phi_m{}'.format(mass),'raw_best_4g_m{}[13]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_gamma2_id_m{}'.format(mass),'raw_best_4g_m{}[23]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_idx4_m{}'.format(mass),'raw_best_4g_m{}[27]'.format(mass)) #index of photon 4
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_dxy_m{}'.format(mass),'raw_best_4g_m{}[14]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_valid_m{}'.format(mass),'raw_best_4g_m{}[15]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi1_mass_m{}'.format(mass),'raw_best_4g_m{}[16]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_phi2_mass_m{}'.format(mass),'raw_best_4g_m{}[17]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_uncorr_mass_m{}'.format(mass),'raw_best_4g_m{}[18]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_corr_mass_m{}'.format(mass),'raw_best_4g_m{}[19]'.format(mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_sumID_m{}'.format(mass),'raw_best_4g_m{m}[20]+raw_best_4g_m{m}[21]+raw_best_4g_m{m}[22]+raw_best_4g_m{m}[23]'.format(m=mass))
        ggH4g_1bad=ggH4g_1bad.Define('Photon_ID_m{}'.format(mass),'Photon_preselection && Photon_IdNoIso && ((Photon_isScEtaEE==1 && Photon_corrIso_m{m}<0.3)||(Photon_isScEtaEB==1 && Photon_corrIso_m{m}<0.25))'.format(m=mass))
        ggH4g_1bad=ggH4g_1bad.Define('non_MC_cut_m{}'.format(mass),'sample_isMC==0 && best_4g_uncorr_mass_m{m}<90|best_4g_uncorr_mass_m{m}>150'.format(m=mass))
        ggH4g_1bad=ggH4g_1bad.Define('best_4g_ID_m{}'.format(mass),best_4g_ID_str_ggH4g_1bad.format(m=mass))


        
    ggH4g=ggH4g.Filter('sample_isMC==1 | non_MC_cut_m30==1','blinding_data_samples')
    ggH3g=ggH3g.Filter('sample_isMC==1 | non_MC_cut_m30==1','blinding_data_samples')
    ggH4g_1bad=ggH4g_1bad.Filter('sample_isMC==1 | non_MC_cut_m30==1','blinding_data_samples')


    # Create snapshots for ggH3g and ggH4g
    actions.append(ggH4g.Snapshot('ggH4g', f"{sample}_ggH4g.root", cols, opts))
    actions.append(ggH3g.Snapshot('ggH3g', f"{sample}_ggH3g.root", cols, opts))
    actions.append(ggH4g_1bad.Snapshot('ggH4g_1bad', f"{sample}_ggH4g_1bad.root", cols, opts))

    # Generate and save reports for ggH3g and ggH4g
    save_report(ggH3g, "Report_ggH3g", f"{sample}_ggH3g", opts, actions)
    save_report(ggH4g, "Report_ggH4g", f"{sample}_ggH4g", opts, actions)
    save_report(ggH4g_1bad, "Report_ggH4g_1bad", f"{sample}_ggH4g_1bad", opts, actions)

    # Create snapshots for the 'Runs' tree in both files
    for tree in ['Runs']:
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_ggH4g.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_ggH3g.root", "", opts))
        actions.append(dataframe[tree].Snapshot(tree, f"{sample}_ggH4g_1bad.root", "", opts))

    return actions

def analysis(data,sample):
    actions=[]
    phi_mass=[7,15,20,30,40,50,55]
    actions.extend(ggH(data,phi_mass,sample))
    return actions

