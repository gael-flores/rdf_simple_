
#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

samples_mc_2023 = {
# 2023 DY+Jets MC
'DYJetsToLL_M50_LO_BPix_2023': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer23BPixNanoAODv12-Pilot_130X_mcRun3_2023_realistic_postBPix_v2-v3/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2023',
                         'jobs': 8,
                         'sigma': 5558.0},

'DYJetsToLL_M50_2023': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Summer23NanoAODv12-Pilot_130X_mcRun3_2023_realistic_v14-v3/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2023',
                         'jobs': 8,
                         'sigma': 5558.0},

#2023 tt+GammaGamma
'TTGG_LO_2023': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer23NanoAODv12-130X_mcRun3_2023_realistic_v15-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2023',
                           'jobs':1,
                           'sigma':0.02391},

'TTGG_LO_BPix_2023': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/Run3Summer23BPixNanoAODv12-130X_mcRun3_2023_realistic_postBPix_v6-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2023',
                           'jobs':1,
                           'sigma':0.02391},

}
