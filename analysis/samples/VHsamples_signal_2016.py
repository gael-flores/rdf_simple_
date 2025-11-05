samples_signal_2016 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [15, 20, 30, 40, 50, 55]:
        for ct in [0, 10, 20, 50, 100, 1000]:
            for br in ['2G2Q']:
                for run in ['preVFP', 'postVFP']:
                    samples_signal_2016['{}H{}_M{}_ctau{}_UL16{}'.format(v,br,m,ct,run)] = {
                        'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/tlam/DDP/UL2016{}/{}H_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}cm.root'.format(run,v,br, m, ct),
                        'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                        'veto_triggers': [],
                        'era': '2016{}'.format(run),
                        'sigma': 1.0,
                        'customNanoAOD': True}

#ttH samples
for m in [15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100]:
        for br in ['2G2Q']:
            for run in ['preVFP', 'postVFP']:
                samples_signal_2016['ttH{}_M{}_ctau{}_UL16{}'.format(br,m,ct,run)] = {
                    'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/gfavila/DDP/UL2016{}/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm_v2.root'.format(run, br, m, ct),
                    'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                    'veto_triggers': [],
                    'era': '2016{}'.format(run),
                    'sigma': 1.0,
                    'customNanoAOD': True}

for m in [15, 20, 30, 40, 50, 55]:
    for ct in [1000]:
        for br in ['2G2Q']:
            for run in ['preVFP', 'postVFP']:
                samples_signal_2016['ttH{}_M{}_ctau{}_UL16{}'.format(br,m,ct,run)] = {
                    'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/tlam/DDP/UL2016{}/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm.root'.format(run, br, m, ct),
                    'triggers': ['(HLT_IsoMu24||HLT_IsoTkMu24||HLT_Ele27_WPTight_Gsf)'],
                    'veto_triggers': [],
                    'era': '2016{}'.format(run),
                    'sigma': 1.0,
                    'customNanoAOD': True}


