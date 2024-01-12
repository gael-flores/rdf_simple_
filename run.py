import ROOT
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
#ROOT.ROOT.DisableImplicitMT()
ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
from common.pyhelpers import *

from optparse import OptionParser

parser = OptionParser()
parser.add_option("-s", "--splitFactor", dest="splitFactor",
                  help="Split the sample in N sub samples",type='int',default=0)
parser.add_option("-p", "--processPart", dest="processPart",
                  help="process the i-th sample if you decide to split ",type='int',default=0)

(options, args) = parser.parse_args()

from analysis.ddp import *
from analysis.ddpSamples import analysis_samples


#root://cmsxrootd.fnal.gov//

data=createDataSet(analysis_samples[args[0]],options.splitFactor,options.processPart,'root://cmsxrootd.fnal.gov//')
analysis(data)
 
