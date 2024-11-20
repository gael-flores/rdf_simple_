from common.plotter import *

luminosity=87.09

#make a plotter for the antiID data that will be used to estimate the background
#background = rdf_plotter('data_ggH4g_antiID.root',False,'0.0009943081222885356','ggH4g_antiID')
#make a plotter for the signal
signal = rdf_plotter('ggH_M30_ctau0_ggH4g.root',True,'1','ggH4g','1','Report_ggH4g')
#make a plotter for the blinded data 
#data = rdf_plotter('data_ggH4g.root',False,'1','ggH4g')

#set the colors
#background.setFillProperties(1001,ROOT.kAzure-9)
signal.setFillProperties(0,ROOT.kWhite)
signal.setLineProperties(1,ROOT.kRed,3)


#create the stack and add the plotters above
stack = combined_plotter()
#stack.add_plotter(background,"bkg","background",typeP = "background")
stack.add_plotter(signal,"signal","ggH #rightarrow #phi#phi#rightarrow 4#gamma ",typeP = "signal")
#stack.add_plotter(data,"data","data",typeP = "data")





 
