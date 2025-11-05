mass_list = [7,15,20,30,40,50,55]


# Standard Analysis Cuts
cuts = {}

# Tight + veto lepton cuts
cuts['W'] = {'MU': "Muon_isTrigger[W_l1_idx] && Muon_nloose == 1 && Muon_nveto == 0 && Electron_nloose==0",
             'ELE': "Electron_isTrigger[W_l1_idx] && Electron_nloose==1 && Muon_nloose == 0"}
             
cuts['Z'] = {'MU': "(Muon_isTrigger[Z_idx[0]] || Muon_isTrigger[Z_idx[1]]) && Muon_nloose == 2 && Electron_nloose==0",
             'ELE': "(Electron_isTrigger[Z_idx[0]] || Electron_isTrigger[Z_idx[1]]) && Electron_nloose == 2 && Muon_nloose == 0"}

# Photon pt and kinematic cuts
cuts['pt'] = {}
cuts['photons'] = {}
for m in mass_list:
    cuts['pt'][m] = "Photon_pt[best_2g_idx1_m{m}] > 35 && Photon_pt[best_2g_idx2_m{m}] > 25".format(m=m)
    cuts['photons'][m] = "best_2g_pt_m{m} > 20 && best_2g_raw_mass_m{m} > 4".format(m=m)

# w->e+nu misID cuts
cuts['misID'] = {}
for v in ['W', 'Z']:
    cuts['misID'][v] = {}
    for l in ['ELE', 'MU']:
        cuts['misID'][v][l] = {}
        for m in mass_list:
            if (v == 'W' and l == 'ELE'):
                cuts['misID'][v][l][m] = "abs(best_2g_misID1_m{m}-91) > 10 && abs(best_2g_misID2_m{m}-91) > 10 && abs(best_2g_misID3_m{m} - 91) > 15".format(m=m)
            else:
                cuts['misID'][v][l][m] = "1"

# Z->ll FSR cuts
cuts['fsr'] = {}
for v in ['W', 'Z']:
    cuts['fsr'][v] = {}
    for m in mass_list:
        cuts['fsr'][v][m] = "1"
        if (v=='W'):
            cuts['fsr'][v][m] = "1"
        else:
            cuts['fsr'][v][m] = "Photon_isFSR[best_2g_idx1_m{m}]==0 && Photon_isFSR[best_2g_idx2_m{m}]==0".format(m=m)


cuts['sr'] = {} # Signal region
cuts['cr'] = {} # anti-ID control region
cuts['cl'] = {} # anti-ID anti-ID control region
for v in ['W', 'Z']:
    cuts['sr'][v] = {}
    cuts['cr'][v] = {}
    cuts['cl'][v] = {}
    for l in ['ELE', 'MU']:
        cuts['sr'][v][l] = {}
        cuts['cr'][v][l] = {}
        cuts['cl'][v][l] = {}
        for m in mass_list:
            cuts['sr'][v][l][m] = "(" + cuts[v][l] + ")&&(" + cuts['pt'][m] + ")&&(" + cuts['photons'][m] + ")&&(" + \
                                        cuts['misID'][v][l][m] + ")&&(" + cuts['fsr'][v][m] + ")&&(" +"best_2g_sumID_m{m}==2 && best_2g_valid_m{m})".format(m=m)
            cuts['cr'][v][l][m] = "(" + cuts[v][l] + ")&&(" + cuts['pt'][m] + ")&&(" + cuts['photons'][m] + ")&&(" + \
                                        cuts['misID'][v][l][m] + ")&&(" + cuts['fsr'][v][m] + ")&&(" +"best_2g_sumID_m{m}==1 && best_2g_valid_m{m})".format(m=m)            
            cuts['cl'][v][l][m] = "(" + cuts[v][l] + ")&&(" + cuts['pt'][m] + ")&&(" + cuts['photons'][m] + ")&&(" + \
                                        cuts['misID'][v][l][m] + ")&&(" + cuts['fsr'][v][m] + ")&&("  +"best_2g_sumID_m{m}==0 && best_2g_valid_m{m})".format(m=m)

# Cuts for signal efficiencies:
effCuts = {}
effCuts['HLT'] = {'ELE': {'2018': "HLT_Ele32_WPTight_Gsf",
                          '2017': "HLT_Ele32_WPTight_Gsf_L1DoubleEG",
                          '2016': "HLT_Ele27_WPTight_Gsf"},
                  'MU': {'2018': "HLT_IsoMu24",
                         '2017': "HLT_IsoMu27",
                         '2016': "HLT_IsoMu24"}
              }

# Cuts for W->e+nu 2018 HEM
cutsHEM = "(Electron_eta[W_l1_idx]>-1.3 || Electron_eta[W_l1_idx]<-3.0) && (Electron_phi[W_l1_idx]>-0.85 || Electron_phi[W_l1_idx]<-1.57)"
