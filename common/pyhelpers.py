import subprocess
import ROOT



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
    if 'era' in info.keys():
        meta['era'] = info['era']
    else:
        meta['era'] = '2018' # Assume 2018 if not included
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
        
