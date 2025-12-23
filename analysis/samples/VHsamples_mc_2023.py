
#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

samples_mc_2023 = {
# 2023 DY+Jets MC
'DYJetsToLL_M50_LO_2023postBPIX': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer23BPixNanoAODv12-Pilot_130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2023',
                         'jobs': 8,
                         'sigma': 5558.0}, 

'DYJetsToLL_M50_2023preBPIX': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer23NanoAODv12-Pilot_130X_mcRun3_2023_realistic_v14-v3/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2023',
                         'jobs': 8,
                         'sigma': 5645.0},


#2023 tt+Jets MC
'TTto2L2Nu_2Jets_2023preBPIX': {'dataset': '/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 94.28}, 

'TTtoLminusNu2Q_2Jets_2023postBPIX': {'dataset': '/TTtoLminusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 189.8}, 

'TTtoLplusNu2Q_2Jets_2023postBPIX': {'dataset': '/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 190.5}, 

'TTto2L2Nu_2Jets_2023preBPIX': {'dataset': '/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 95.87}, 

'TTtoLminusNu2Q_2Jets_2023preBPIX': {'dataset': '/TTtoLminusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 192.2}, 

'TTtoLplusNu2Q_2Jets_2023preBPIX': {'dataset': '/TTtoLplusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 192.4}, 

#2023 W+Jets MC
'WToLNu_2Jets_0J_2023postBPIX': {'dataset': '/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 55810.}, 

'WToLNu_2Jets_1J_2023postBPIX': {'dataset': '/WtoLNu-2Jets_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 9542.}, 

'WToLNu_2Jets_2J_2023postBPIX': {'dataset': '/WtoLNu-2Jets_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 3632.}, 

'WToLNu_2Jets_0J_2023preBPIX': {'dataset': '/WtoLNu-2Jets_0J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 55740.},

'WToLNu_2Jets_1J_2023preBPIX': {'dataset': '/WtoLNu-2Jets_1J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 9484.}, 

'WToLNu_2Jets_2J_2023preBPIX': {'dataset': '/WtoLNu-2Jets_2J_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v14-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 3593.}, 

#2023 W+Gamma MC
'WGToLNuG_1Jets_PTG10to100_2023postBPIX': {'dataset': '/WGtoLNuG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 667.6},

'WGToLNuG_1Jets_PTG100to200_2023postBPIX': {'dataset': '/WGtoLNuG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 2.213},

'WGToLNuG_1Jets_PTG10to100_2023preBPIX': {'dataset': '/WGtoLNuG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 669.2}, 

'WGToLNuG_1Jets_PTG100to200_2023preBPIX': {'dataset': '/WGtoLNuG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 2.204},

#2018 tt+Gamma
'TTG_1Jets_PTG10to100_2023postBPIX': {'dataset': '/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 4.200},

'TTG_1Jets_PTG100to200_2023postBPIX': {'dataset': '/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 0.4083},

'TTG_1Jets_PTG200_2023postBPIX': {'dataset': '/TTG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 0.1270},

'TTG_1Jets_PTG10to100_2023preBPIX': {'dataset': '/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 4.197},

'TTG_1Jets_PTG100to200_2023preBPIX': {'dataset': '/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 0.4054},

'TTG_1Jets_PTG200_2023preBPIX': {'dataset': '/TTG-1Jets_PTG-200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2023',
                       'jobs': 8,
                       'sigma': 0.1277},

#2023 tt+GammaGamma
'TTGG_LO_2023preBPIX': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2023',
                           'jobs':1,
                           'sigma':0.02391},

'TTGG_LO_2023postBPIX': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2023',
                           'jobs':1,
                           'sigma':0.02391}

}
