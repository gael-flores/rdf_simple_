mass_list = [15,20,30,40,50,55]


# Standard Analysis Cuts
cuts = {}

# Tight + veto lepton cuts
cuts['W'] = {'MU': "(Muon_isTrigger[W_l1_idx] && Muon_nveto == 0 && Electron_nveto==0 && Muon_nloose==1 && Electron_nloose==0)",
             'ELE': "(Electron_isTrigger[W_l1_idx] && Electron_nloose==1 &&Electron_nveto==0 && Muon_nveto == 0 && Muon_nloose==0)"}
             
cuts['Z'] = {'MU': "((Muon_isTrigger[Z_idx[0]] || Muon_isTrigger[Z_idx[1]]) && Muon_nloose==2 && Muon_nveto == 0 && Electron_nveto==0 && Electron_nloose==0)",
             'ELE': "((Electron_isTrigger[Z_idx[0]] || Electron_isTrigger[Z_idx[1]]) && Electron_nloose==2 && Electron_nveto == 0 && Muon_nveto == 0 && Muon_nloose==0)"}


# Photon pt and kinematic cuts
cuts['pt'] = {}
cuts['photons'] = {}
for m in mass_list:
    cuts['pt'][m] = "(Photon_pt[best_2g_idx1_m{m}] > 35 && Photon_pt[best_2g_idx2_m{m}] > 25)".format(m=m)
    cuts['photons'][m] = "(best_2g_pt_m{m} > 20 && best_2g_raw_mass_m{m} > 4)".format(m=m)

# w->e+nu misID cuts
cuts['misID'] = {}
for v in ['W', 'Z']:
    cuts['misID'][v] = {}
    for l in ['ELE', 'MU']:
        cuts['misID'][v][l] = {}
        for m in mass_list:
            if (v == 'W' and l == 'ELE'):
                cuts['misID'][v][l][m] = "(abs(best_2g_misID1_m{m}-91) > 10 && abs(best_2g_misID2_m{m}-91) > 10 && abs(best_2g_misID3_m{m} - 91) > 15)".format(m=m)
            else:
                cuts['misID'][v][l][m] = "(1)"

# Z->ll FSR cuts
cuts['fsr'] = {}
for v in ['W', 'Z']:
    cuts['fsr'][v] = {}
    for m in mass_list:
        cuts['fsr'][v][m] = "(1)"
        if (v=='W'):
            cuts['fsr'][v][m] = "(1)"
        else:
            cuts['fsr'][v][m] = "(Photon_isFSR[best_2g_idx1_m{m}]==0 && Photon_isFSR[best_2g_idx2_m{m}]==0)".format(m=m)


#redefine the cuts here based on the analysis for easier parsing
for analysis in ['wmn2g','wen2g','zmm2g','zee2g']:
    cuts[analysis]={}
    if analysis=='wmn2g':
        v='W'
        l='MU'
    elif analysis=='zmm2g':
        v='Z'
        l='MU'
    elif analysis=='wen2g':
        v='W'
        l='ELE'
    elif analysis=='zee2g':
        v='Z'
        l='ELE'
    for m in mass_list:
        cuts[analysis][m] = {} 
        cuts[analysis][m]['presr']='&&'.join([cuts[v][l],
                                           cuts['pt'][m],
                                           cuts['photons'][m],
                                           cuts['misID'][v][l][m],
                                           cuts['fsr'][v][m],
                                           f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==2)"])
        cuts[analysis][m]['precr']='&&'.join([cuts[v][l],
                                              cuts['pt'][m],
                                              cuts['photons'][m],
                                              cuts['misID'][v][l][m],
                                              cuts['fsr'][v][m],
                                              f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==1)"])
                                              
#                                              f"(Photon_passCutBasedID[best_2g_idx1_m{m}]>0)"])                                              
#                                           f"(Photon_passCutBasedID[best_2g_idx1_m{m}]>0)"])

        cuts[analysis][m]['precl']='&&'.join([cuts[v][l],
                                           cuts['pt'][m],
                                           cuts['photons'][m],
                                           cuts['misID'][v][l][m],
                                           cuts['fsr'][v][m],
                                           f"((Photon_passCutBasedID[best_2g_idx1_m{m}]+Photon_passCutBasedID[best_2g_idx2_m{m}])==0)"])
        cuts[analysis][m]['sr']='&&'.join([cuts[analysis][m]['presr'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['cr']='&&'.join([cuts[analysis][m]['precr'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['cl']='&&'.join([cuts[analysis][m]['precl'],f'(best_2g_dxy_m{m}>-10)'])
        cuts[analysis][m]['ssb']='&&'.join([cuts[analysis][m]['presr'],f'(best_2g_raw_mass_m{m}>62.5)'])
        cuts[analysis][m]['csb']='&&'.join([cuts[analysis][m]['precr'],f'(best_2g_raw_mass_m{m}>62.5)'])
        cuts[analysis][m]['lsb']='&&'.join([cuts[analysis][m]['precl'],f'(best_2g_raw_mass_m{m}>62.5)'])
        

            

# Cuts for signal efficiencies:
effCuts = {}
effCuts['HLT'] = {'ELE': {'2018': "(HLT_Ele32_WPTight_Gsf)",
                          '2017': "(HLT_Ele32_WPTight_Gsf_L1DoubleEG)",
                          '2016': "(HLT_Ele27_WPTight_Gsf)"},
                  'MU': {'2018': "(HLT_IsoMu24)",
                         '2017': "(HLT_IsoMu27)",
                         '2016': "(HLT_IsoMu24)"}
              }

# Cuts for W->e+nu 2018 HEM
cutsHEM = "(Electron_eta[W_l1_idx]>-1.3 || Electron_eta[W_l1_idx]<-3.0) && (Electron_phi[W_l1_idx]>-0.85 || Electron_phi[W_l1_idx]<-1.57)"


#analysis definitions
def applyDefinitions(plotter):
    pass
#    for m in mass_list:
#        plotter.define(f"best_2g_gamma1_tauMass_m{m}",f"Photon_tauMass[best_2g_idx1_m{m}]")
#        plotter.define(f"best_2g_gamma2_tauMass_m{m}",f"Photon_tauMass[best_2g_idx2_m{m}]")
#        plotter.define(f"best_2g_gamma1_ID_m{m}",f"Photon_passCutBasedID[best_2g_idx1_m{m}]")
#        plotter.define(f"best_2g_gamma2_ID_m{m}",f"Photon_passCutBasedID[best_2g_idx2_m{m}]")
#        plotter.define(f"best_2g_gamma1_Iso_m{m}",f"Photon_passPhIso[best_2g_idx1_m{m}]")
#        plotter.define(f"best_2g_gamma2_Iso_m{m}",f"Photon_passPhIso[best_2g_idx2_m{m}]")
#        plotter.define(f"best_2g_gamma1_IDandIso_m{m}",f"Photon_IDandIso[best_2g_idx1_m{m}]")
#        plotter.define(f"best_2g_gamma2_IDandIso_m{m}",f"Photon_IDandIso[best_2g_idx2_m{m}]")
