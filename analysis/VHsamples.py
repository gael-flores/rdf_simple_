samples = {}

# Import samples for all eras
for era in ['2018', '2017', '2016']:
    for s in ['signal', 'mc', 'data']:

        cmd = """
from analysis.samples.VHsamples_{s}_{era} import *
for s in samples_{s}_{era}:
   if s not in samples:
        samples[s] = samples_{s}_{era}[s]
""".format(s=s, era=era)

        #cmd = "import analysis.samples.VHsamples_{s}_{era} as samples_{s}_{era}".format(s=s, era=era)
        exec(cmd)

#import pdb; pdb.set_trace()
'''
samples={}

for s in samples_data_2018:
    if s not in samples:
        samples[s] = samples_data_2018[s]

for s in samples_signal_2018:
    if s not in samples:
        samples[s] = samples_signal_2018[s]

for s in samples_mc_2018:
    if s not in samples:
        samples[s] = samples_mc_2018[s]

for s in samples_data_2017:
    if s not in samples:
        samples[s] = samples_data_2017[s]

for s in samples_mc_2017:
    if s not in samples:
        samples[s] = samples_mc_2017[s]

for s in samples_sig_2017:
    if s not in samples:
        samples[s] = samples_signal_2017[s]

for s in samples_data_2016:
    if s not in samples:
        samples[s] = samples_data_2016[s]

for s in samples_mc_2016:
    if s not in samples:
        samples[s] = samples_mc_2016[s]

for s in samples_mc_2016:
    if s not in samples:
        samples[s] = samples_mc_2016[s]
'''
