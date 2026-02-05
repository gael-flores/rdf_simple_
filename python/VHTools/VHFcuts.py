mass_list = [20,50]

# Standard Analysis Cuts
cuts = {}

# Photon pt and kinematic cuts
cuts['pt'] = {}
cuts['photons'] = {}
for m in mass_list:
    cuts['pt'][m] = "Photon_pt[best_2g_idx1_m{m}] > 35 && Photon_pt[best_2g_idx2_m{m}] > 25".format(m=m)
    cuts['photons'][m] = f"best_2g_pt_m{m} > 20 && best_2g_raw_mass_m{m} > 4"

cuts['sr'] = {} # Signal region
cuts['cr'] = {} # anti-ID control region
cuts['cl'] = {} # anti-ID anti-ID control region
cuts['fr'] = {} #fake rate control region
for m in mass_list:
    cuts['sr'][m] = f"({cuts['pt'][m]})&&({cuts['photons'][m]})&&(best_2g_sumID_m{m}==2 && best_2g_valid_m{m})"
    cuts['cr'][m] = f"({cuts['pt'][m]})&&({cuts['photons'][m]})&&(best_2g_sumID_m{m}==1 && best_2g_valid_m{m})"
    cuts['cl'][m] = f"({cuts['pt'][m]})&&({cuts['photons'][m]})&&(best_2g_sumID_m{m}==0 && best_2g_valid_m{m})"
    cuts['fr'][m] = f"({cuts['pt'][m]})&&({cuts['photons'][m]})&&(Photon_passCutBasedID[best_2g_idx1_m{m}]>0)"