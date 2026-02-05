#source for listed cross sections: https://xsecdb-xsdb-official.app.cern.ch/xsdb
#source for calculating cross sections using GenXSecAnalyzer: https://twiki.cern.ch/twiki/bin/view/CMSPublic/GeneratorsHATSatLPC2019InDevelopment#Using_GenXSecAnalyzer

#multiply preVFP with 0.54,postVFP by 0.46
samples_mc_2016 = {
'DYJetsToLL_M50_LO_2016preVFP': {'dataset':'/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
                           'triggers':["(1)"],
                           'veto_triggers':[],
                           'era': '2016preVFP',
                           'jobs':8,
                                 'sigma':1921.8*3*0.54},
'WJetsToLNu_LO_2016preVFP': {'dataset': '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODAPVv9-106X_mcRun2_asymptotic_preVFP_v11-v1/NANOAODSIM',
                       'triggers': ["(1)"],
                       'veto_triggers': [],
                       'era': '2016preVFP',
                       'jobs': 8,
                       'sigma': 3*20508.9*0.54},
# postVFP samples
'DYJetsToLL_M50_LO_2016postVFP': {'dataset':'/DYJetsToLL_M-50_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
                           'triggers':["(1)"],
                           'veto_triggers':[],
                           'era': '2016postVFP',
                           'jobs':8,
                                  'sigma':1921.8*3*0.46},
'WJetsToLNu_LO_2016postVFP': {'dataset': '/WJetsToLNu_TuneCP5_13TeV-madgraphMLM-pythia8/RunIISummer20UL16NanoAODv9-106X_mcRun2_asymptotic_v17-v1/NANOAODSIM',
                       'triggers': ["(1)"],
                       'veto_triggers': [],
                       'era': '2016postVFP',
                       'jobs': 8,
                       'sigma': 3*20508.9*0.46},
}
