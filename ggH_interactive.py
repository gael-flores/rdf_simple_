from common.plotter import *
from math import sqrt
luminosity=87.09

#make a plotter for the antiID data that will be used to estimate the background
#background = rdf_plotter('data_ggH4g_antiID.root',False,'0.0009943081222885356','ggH4g_antiID')
#make a plotter for the signal
data_sideband_cut="((best_4g_corr_mass_m30<110)||(best_4g_corr_mass_m30 > 140))"
endcap_branches={
"Photon_hoe":"<0.2",
"Photon_sieie":"<0.045",
"Photon_corrIso_m30":"<0.3"
}

barrel_branches={
"Photon_hoe":"<0.3",
"Photon_sieie":"<0.035",
"Photon_corrIso_m30":"<0.25"
}



barrel_loose_cut = "&&".join(f"{branch}[best_4g_idx{i}_m30]{cut}" for branch, cut in barrel_branches.items() for i in range(1,5))
endcap_loose_cut = "&&".join(f"{branch}[best_4g_idx{i}_m30]{cut}" for branch, cut in endcap_branches.items() for i in range(1,5))


signal  = rdf_plotter('ggH_M30_ctau0_ggH4g.root',True,'(56*(3e-4))','ggH4g','1','Report_ggH4g')
sideband_data = rdf_plotter('EGamma_2018_all_ggH4g.root',False,'1','ggH4g',data_sideband_cut, "Report_ggH4g")
background    = rdf_plotter('EGamma_2018_all_ggH4g_1bad.root',False,'0.9307568438003','ggH4g_1bad','1', "Report_ggH4g_1bad")

#set the colors
background.setFillProperties(1001,ROOT.kAzure-9)
signal.setFillProperties(0,ROOT.kWhite)
signal.setLineProperties(1,ROOT.kRed,3)


#create the stack and add the plotters above
stack = combined_plotter()
#stack.add_plotter(sideband_data,"bkg","sideband data",typeP = "background")
stack.add_plotter(background,"bkg","background",typeP = "background")
stack.add_plotter(signal,"signal","signal",typeP = "signal")

#stack.add_plotter(data,"data","data",typeP = "data")














 
