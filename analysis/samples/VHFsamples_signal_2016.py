samples_signal_2016 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [20,50]:
        for ct in [0, 100, 1000]:
            for br in ['2G2Q']:
                for run in ['preVFP', 'postVFP']:
                    samples_signal_2016[f'{v}H{br}_M{m}_ctau{ct}_UL16{run}'] = {
                        'dataset': f'local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2016{run}/{v}H_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}cm.root',
                        'triggers': ["(1)"],
                        'veto_triggers': [],
                        'era': f'2016{run}',
                        'sigma': 1.0,
                        'customNanoAOD': True}

#ttH samples
#0-100 mm samples
for m in [20,50]:
    for ct in [0, 100]:
        for br in ['2G2Q']:
            for run in ['preVFP', 'postVFP']:
                samples_signal_2016[f'ttH{br}_M{m}_ctau{ct}_UL16{run}'] = {
                    'dataset': f'local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2016{run}/ttH_ttSemiLep_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}mm_v2.root',
                    'triggers': ["(1)"],
                    'veto_triggers': [],
                    'era': f'2016{run}',
                    'sigma': 1.0,
                    'customNanoAOD': True}

#1000 mm samples
for m in [20,50]:
    for ct in [1000]:
        for br in ['2G2Q']:
            for run in ['preVFP', 'postVFP']:
                samples_signal_2016[f'ttH{br}_M{m}_ctau{ct}_UL16{run}'] = {
                    'dataset': f'local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2016{run}/ttH_ttSemiLep_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}mm.root',
                    'triggers': ["(1)"],
                    'veto_triggers': [],
                    'era': f'2016{run}',
                    'sigma': 1.0,
                    'customNanoAOD': True}


