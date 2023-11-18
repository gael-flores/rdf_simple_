import ROOT
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
from common.pyhelpers import *

from optparse import OptionParser

parser = OptionParser()
#parser.add_option("-s", "--samples", dest="samples",
#                  help="List of sample nick names")
(options, args) = parser.parse_args()

from analysis.ddp import *
from analysis.ddpSamples import analysis_samples

toProcess=[]
for s in args:
    toProcess.append(analysis_samples[s])
data=createDataSet(toProcess)

analysis(data)
 
