import subprocess
import ROOT

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
        if key=='files' or key=='isMC' or key=='trigger':
            continue;
        dataframe['Events']=dataframe['Events'].Define('sample_'+key,str(data[key]))
    dataframe['Runs'] = ROOT.RDataFrame("Runs", data['files'])
    return dataframe


def loadSample(info,locator='root://cms-xrd-global.cern.ch//'):
    meta={}
    files=[]
    print(info['dataset'])
    if 'local:' in info['dataset']:
        files=info['dataset'].split('local:')[1].split(',')
    else:
        p = subprocess.Popen('/cvmfs/cms.cern.ch/common/dasgoclient -dasmaps=./ -query="file dataset={}"'.format(info['dataset']), shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            files.append(locator+(line.decode('ASCII').split('\n')[0]))
        retval = p.wait()
    triggerStr=[]
    for t in info['triggers']:
        triggerStr.append('{}==1'.format(t))
    for t in info['veto_triggers']:
        triggerStr.append('{}==0'.format(t))
    triggerDecision = '&&'.join(triggerStr)
    meta['trigger']=triggerDecision
    if not ('SIM' in info['dataset'] or 'local' in info['dataset']):
        meta['isMC']=0
    else:
        meta['isMC']=1
    for tag, data in info.items():
        if tag not in ['dataset','triggers','veto_triggers']:
            meta[tag]=data
    meta['files'] = files        
    return meta


def split(a, n):
    k, m = divmod(len(a), n)
    return (a[i*k+min(i, m):(i+1)*k+min(i+1, m)] for i in range(n))


def createDataSet(sample,splitFactor=0,processPart=0,locator='root://cms-xrd-global.cern.ch//'):
    spec = loadSample(sample,locator)
    if splitFactor==0:
        return spec
    else:
        #split the sample
        fileNames = spec['files']
        splitFileNames = list(split(fileNames,splitFactor))[processPart]
        print('split file names')
        for s in splitFileNames:
            print(s)
        spec['files']=splitFileNames
        return spec
        
