from optparse import OptionParser
import ROOT
import os
import time
import importlib
import subprocess
from datetime import date

parser = OptionParser()

parser.add_option("-a", "--analysis", dest="analysis",
                  help="analysis to run",type='str',default='VH')

parser.add_option("-o", "--eosdir", dest="eos",default='root://cmseos.fnal.gov//store/user/bachtis/analysis',
                  help="EOS output Directory")
parser.add_option("-O", "--localdir", dest="localdir",default='/uscmst1b_scratch/lpc1/3DayLifetime/bachtis/',
                  help="Local directory to store the output")

parser.add_option("-y", "--years", dest="years",default='2016,2017,2018',
                  help="Years used in the analysis")
parser.add_option("-d", "--datasets", dest="datasets",default='EGamma,SingleMuon,SingleElectron,JetHT,MET',
                  help="Primary datasets used in the analysis")

(options, args) = parser.parse_args()

samp = importlib.import_module('analysis.{}samples'.format(options.analysis))


root = options.eos.split('/store')[0]
path = options.eos.split(root)[1]

result = subprocess.run(["xrdfs",root, "ls", path], capture_output=True, text=True, check=True)

all_files = result.stdout.split('\n')

years=options.years.split(',')
pds=options.datasets.split(',')


today = (date.today()).strftime("%Y_%m_%d")

#loop on the samples
for d,info in samp.samples.items():    
    #find all the files on eos, for this sample
    matches = [s for s in all_files if ('/'+d) in s]
    #ifg a single file just copy it.
    if len(matches)==0:
        print(f"Files not found for sample {d}")
    elif len(matches)==1:
        os.system(f"xrdcp {root}/{matches[0]} {options.localdir}{d}.root")
    else: #many files
        new_files = [root+item for item in matches]
        os.system(f"hadd {options.localdir}{d}.root "+" ".join(new_files))

#Now deal with extenstions
result = subprocess.run(["ls",options.localdir], capture_output=True, text=True, check=True)
files=result.stdout.split('\n')
#find files with extentions
filesWithExt= [item for item in files if 'ext' in item]

#replace the extention with *
fs4 = [s.replace('_ext4_','*') for s in filesWithExt]
fs3 = [s.replace('_ext3_','*') for s in fs4]
fs2 = [s.replace('_ext2_','*') for s in fs3]
fs1 = [s.replace('_ext1_','*') for s in fs2]
fs0 = [s.replace('_ext_','*') for s in fs1]
wildcard_str = list(set(fs0))
for s in wildcard_str:
    os.system(f"hadd  {options.localdir}tmp.root {options.localdir}{s}")
    os.system(f"rm {options.localdir}{s}")
    fname=s.replace("*",'_')
    os.system(f"mv {options.localdir}tmp.root {options.localdir}{fname}")
#finally categorize files
for year in years:
    os.system(f"mkdir {options.localdir}MC{year}_{today}")
    os.system(f"mkdir {options.localdir}DATA{year}_{today}")    
    for pd in pds:
        os.system(f'mv {options.localdir}{pd}*{year}*.root {options.localdir}DATA{year}_{today}/')
    os.system(f'mv {options.localdir}*{year}*.root {options.localdir}MC{year}_{today}/')
