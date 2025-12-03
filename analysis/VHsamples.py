samples = {}

# Import samples for all eras
for era in ['2016','2017','2018', '2022']:
    for s in ['mc','signal','data']:
        cmd = """
from analysis.samples.VHsamples_{s}_{era} import *
for s in samples_{s}_{era}:
   if s not in samples:
        samples[s] = samples_{s}_{era}[s]
""".format(s=s, era=era)

        exec(cmd)


#samples['TTLep_NLO_2018'] = {'dataset': '/TTTo2L2Nu_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
#                             'triggers': ['(HLT_IsoMu24||HLT_Ele32_WPTight_Gsf)'],
#                             'veto_triggers': [],
#                             'jobs': 16,
#                             'sigma': 831.76*((3*0.108)**2)}

#samples['TTSemi_NLO_2018'] = {'dataset': '/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIISummer20UL18NanoAODv9-106X_upgrade2018_realistic_v16_L1v1-v1/NANOAODSIM',
#                              'triggers': ['(HLT_IsoMu24||HLT_Ele32_WPTight_Gsf)'],
#                              'veto_triggers': [],
#                              'jobs': 16,
#                              'sigma': 831.76*2*(3*0.108)*(1-3*0.108)}

