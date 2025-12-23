
#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

samples_mc_2022 = {
# 2022 DY+Jets MC
'DYJetsToLL_M50_LO_2022postEE': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22EENanoAODv11-forPOG_126X_mcRun3_2022_realistic_postEE_v1-v1/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2022',
                         'jobs': 8,
                         'sigma': 5731.0},

'DYJetsToLL_M50_LO_2022preEE': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22NanoAODv11-forPOG_126X_mcRun3_2022_realistic_v2-v1/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2022',
                         'jobs': 8,
                         'sigma': 5558.0},

# 2022 W+Jets MC
'WJetsToLNu_LO_2022postEE': {'dataset': '/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter22NanoAOD-122X_mcRun3_2021_realistic_v9-v1/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 56120.},

#2022 tt+Jets MC
'TTto2L2Nu_2Jets_2022postEE': {'dataset': '/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 95.49}, 

'TTtoLminusNu2Q_2Jets_2022postEE': {'dataset': '/TTtoLminusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 190.4}, 

'TTtoLplusNu2Q_2Jets_2022postEE': {'dataset': '/TTtoLplusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 192.6}, 


'TTto2L2Nu_2Jets_2022preEE': {'dataset': '/TTto2L2Nu-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 96.71}, 

'TTtoLminusNu2Q_2Jets_2022preEE': {'dataset': '/TTtoLminusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 190.9}, 

'TTtoLplusNu2Q_2Jets_2022postEE': {'dataset': '/TTtoLplusNu2Q-2Jets_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 191.5}, 


#2022 W+Gamma MC
'WGtoLNuG_1Jets_PTG100to200_2022postEE': {'dataset':'/WGtoLNuG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':2.219}, 

'WGtoLNuG_1Jets_PTG100to200_2022preEE': {'dataset':'/WGtoLNuG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':2.221}, 

'WGtoLNuG_1Jets_PTG10to100_2022postEE': {'dataset':'/WGtoLNuG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':662.2}, 

'WGtoLNuG_1Jets_PTG10to100_2022preEE': {'dataset':'/WGtoLNuG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':662.2},

'WGtoLNuG_1Jets_PTG200to400_2022postEE': {'dataset':'/WGtoLNuG-1Jets_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.2908}, 

'WGtoLNuG_1Jets_PTG200to400_2022preEE': {'dataset':'/WGtoLNuG-1Jets_PTG-200to400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.2908}, 

'WGtoLNuG_1Jets_PTG400to600_2022postEE': {'dataset':'/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.02231},

'WGtoLNuG_1Jets_PTG400to600_2022preEE': {'dataset':'/WGtoLNuG-1Jets_PTG-400to600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.02231},

'WGtoLNuG_1Jets_PTG600_2022postEE': {'dataset':'/WGtoLNuG-1Jets_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.004907}, 

'WGtoLNuG_1Jets_PTG600_2022preEE': {'dataset':'/WGtoLNuG-1Jets_PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.004907},


#2022 tt+Gamma

'TTGJets_1Jets_PTG-100to200_2022preEE': {'dataset':'/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.4114},

'TTGJets_1Jets_PTG-100to200_2022postEE': {'dataset':'/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.4072},

'TTGJets_1Jets_PTG-10to100_2022preEE': {'dataset':'/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':4.216},

'TTGJets_1Jets_PTG-10to100_2022postEE': {'dataset':'/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v4/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':4.237},


#2022 tt+GammaGamma
'TTGG_LO_2022preEE': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.02391},

'TTGG_LO_2022postEE': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.02391}

}