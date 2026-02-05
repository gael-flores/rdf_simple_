
#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

samples_mc_2018 = {
# 2018 DY+Jets MC
'DYJetsToLL_M50_LO_2018': {'dataset':'/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
                           'triggers':["(1)"],
                           'veto_triggers':[],
                           'era': '2018',
                           'jobs':8,
                           'sigma':1921.8*3},
'DYJetsToLL_M50_LO_ext1_2018': {'dataset': '/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1_ext1-v1/NANOAODSIM',
                           'triggers':["(1)"],
                           'veto_triggers':[],
                           'era': '2018',
                           'jobs':8,
                           'sigma': 5321.0},
'DYJetsToLL_M10to50_LO_2018': {'dataset': '/DYJetsToLL_M-10to50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
                           'triggers':["(1)"],
                           'veto_triggers':[],
                           'era': '2018',
                           'jobs':8,
                           'sigma': 15810.0},
'WJetsToLNu_LO_2018': {'dataset': '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
                       'triggers': ["(1)"],
                       'veto_triggers': [],
                       'era': '2018',
                       'jobs': 8,
                       'sigma': 3*20508.9},
}
