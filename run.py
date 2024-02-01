import ROOT
import importlib
# Enable multi-threading
ROOT.ROOT.EnableImplicitMT()
#ROOT.ROOT.DisableImplicitMT()
ROOT.gInterpreter.Declare('#include "common/chelpers.h"')
from common.pyhelpers import *
from optparse import OptionParser



def load_meta_data(data):
    dataframe = {}
    #Declare dataframe
    dataframe['Events'] =ROOT.RDataFrame('Events',data['files'])
    #read the HLT string from the sample
    #dataframe=dataframe.DefinePerSample('HLTstring','rdfsampleinfo_.GetS("trigger")')
    dataframe['Events']=dataframe['Events'].Define('HLT_passed',data['trigger']) #need to fix to suppport trigger/amples!!!!
    dataframe['Events']=dataframe['Events'].Define('sample_isMC',str(data['isMC']))
    #get the list of meta keys
    meta_keys =data.keys()
    for key in meta_keys:
        if key=='files' or key=='isMC' or key=='trigger' or key=='jobs':
            continue;
        dataframe['Events']=dataframe['Events'].Define('sample_'+key,str(data[key]))
    dataframe['Runs'] = ROOT.RDataFrame("Runs", data['files'])
    return dataframe




parser = OptionParser()
parser.add_option("-a", "--analysis", dest="analysis",
                  help="Split the sample in N sub samples",type='str',default='VH')
parser.add_option("-s", "--splitFactor", dest="splitFactor",
                  help="Split the sample in N sub samples",type='int',default=0)
parser.add_option("-p", "--processPart", dest="processPart",
                  help="process the i-th sample if you decide to split ",type='int',default=0)
(options, args) = parser.parse_args()



an = importlib.import_module('analysis.{}'.format(options.analysis))
samp = importlib.import_module('analysis.{}samples'.format(options.analysis))

actions=[]
for sample in args:
    data=createDataSet(samp.samples[sample],options.splitFactor,options.processPart,'root://cmsxrootd.fnal.gov/')
    actions.extend(an.analysis(data,sample))
ROOT.RDF.RunGraphs(actions)

