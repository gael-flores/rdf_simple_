
#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

samples_mc_2024 = {
#2024 DY+Jets MC
'DYJetsToLL_M50_NoPU_v6_v2_2024': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/Run3Winter24NanoAOD-NoPU_Pilot_133X_mcRun3_2024_realistic_v6-v2/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2024',
                         'jobs': 8,
                         'sigma': 5669.0}, 

'DYJetsToLL_M50_v26-v2_2024': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13p6TeV-madgraphMLM-pythia8/RunIII2024Summer24NanoAOD-Pilot2024wmLHEGS_140X_mcRun3_2024_realistic_v26-v2/NANOAODSIM',
                         'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                         'veto_triggers': [],
                         'era': '2024',
                         'jobs': 8,
                         'sigma': 5634.0}, 

#2024 W+Gamma MC
'WGtoLNuG_1Jets_PTG-100_2024': {'dataset': '/WGtoLNuG-1Jets_Bin-PTG-100_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2024',
                       'jobs': 8,
                       'sigma': 2.551},

'WGtoLNuG_1Jets_PTG-200_2024': {'dataset': '/WGtoLNuG-1Jets_Bin-PTG-200_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2024',
                       'jobs': 8,
                       'sigma': 0.3176}, 

'WGtoLNuG_1Jets_PTG-400_2024': {'dataset': '/WGtoLNuG-1Jets_Bin-PTG-400_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2024',
                       'jobs': 8,
                       'sigma': 0.02642}, 

'WGtoLNuG_1Jets_PTG-600_2024': {'dataset': '/WGtoLNuG-1Jets_Bin-PTG-600_TuneCP5_13p6TeV_amcatnloFXFX-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM',
                       'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                       'veto_triggers': [],
                       'era': '2024',
                       'jobs': 8,
                       'sigma': 0.004754}, 

#2024 tt+GammaGamma
'TTGG_LO_2024': {'dataset': '/TTGG_TuneCP5_13p6TeV_madgraph-madspin-pythia8/RunIII2024Summer24NanoAODv15-150X_mcRun3_2024_realistic_v2-v2/NANOAODSIM',
                           'triggers':['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                           'veto_triggers':[],
                           'era':'2024',
                           'jobs':1,
                           'sigma':0.02391}

}
