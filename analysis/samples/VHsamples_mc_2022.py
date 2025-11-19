
#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

samples_mc_2022 = {
# 2022 DY+Jets MC
'DYJetsToLL_M50_LO_EE_2022': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22EENanoAODv11-forPOG_126X_mcRun3_2022_realistic_postEE_v1-v1/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2022',
                         'jobs': 8,
                         'sigma': 5558.0},


'DYJetsToLL_M50_LO_2022': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer22NanoAODv11-forPOG_126X_mcRun3_2022_realistic_v2-v1/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2022',
                         'jobs': 8,
                         'sigma': 5558.0},

# 2022 W+Jets MC
'WJetsToLNu_LO_2022': {'dataset': '/WJetsToLNu_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter22NanoAOD-122X_mcRun3_2021_realistic_v9-v1/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2022',
                       'jobs': 8,
                       'sigma': 56120.},

#2022 W+Gamma MC - no datasets appear on DAS for 2022


#2022 W+GammaGamma - no datasets appear on DAS for 2022


#2022 tt+Gamma

'TTGJets_1Jets_PTG-100to200_2022': {'dataset':'/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.4114},

'TTGJets_1Jets_PTG-100to200_EE_2022': {'dataset':'/TTG-1Jets_PTG-100to200_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v3/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.4114},

'TTGJets_1Jets_PTG-10to100_2022': {'dataset':'/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22NanoAODv12-130X_mcRun3_2022_realistic_v5-v1/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':4.216},

'TTGJets_1Jets_PTG-10to100_EE_2022': {'dataset':'/TTG-1Jets_PTG-10to100_TuneCP5_13p6TeV_amcatnloFXFXold-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v4/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':4.216},


#2022 tt+GammaGamma
'TTGG_LO_2022': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer22NanoAODv11-126X_mcRun3_2022_realistic_v2-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.02391},

'TTGG_LO_EE_2022': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer22EENanoAODv12-130X_mcRun3_2022_realistic_postEE_v6-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2022',
                           'jobs':1,
                           'sigma':0.02391},


#2022 tt+Jets MC - no datasets appear on DAS for 2022

}
