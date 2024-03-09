samples = {}

# Import samples for all eras
for era in ['2018', '2017', '2016']:
    for s in ['signal', 'mc', 'data']:
#    for s in ['data']:
        cmd = """
from analysis.samples.VHsamples_{s}_{era} import *
for s in samples_{s}_{era}:
   if s not in samples:
        samples[s] = samples_{s}_{era}[s]
""".format(s=s, era=era)

        exec(cmd)
