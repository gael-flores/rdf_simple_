import subprocess
import ROOT

def loadSample(info,locator='root://cms-xrd-global.cern.ch//'):
    files=[]
    print(info['dataset'])
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
    metaData = ROOT.RDF.Experimental.RMetaData()
    if 'SIM' in info['dataset']:
        metaData.Add('isMC',1)
    else:
        metaData.Add('isMC',0)
    for tag, data in info.items():
        if tag not in ['dataset','triggers','veto_triggers']:
            metaData.Add(tag,data)
    sample = ROOT.RDF.Experimental.RSample(info['dataset'],'Events',files,metaData)
    return sample


def createDataSet(samples,locator='root://cms-xrd-global.cern.ch//'):
    spec = ROOT.RDF.Experimental.RDatasetSpec()
    for s in samples:
        spec.AddSample(loadSample(s,locator))
    return spec
        
