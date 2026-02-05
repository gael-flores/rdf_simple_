samples_signal_2017 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [20,50]:
        for ct in [0,100, 1000]:
            for br in ['2G2Q']:
                samples_signal_2017[f'{v}H{br}_M{m}_ctau{ct}_UL17'] = {
                    'dataset': f'local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2017/{v}H_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}cm.root',
                    'triggers': ["(1)"],
                    'veto_triggers': [],
                    'era': '2017',
                    'sigma': 1.0,
                    'customNanoAOD': True}

# ttH samples:
#0-100mm samples
for m in [20,50]:
    for ct in [0,100]:
        for br in ['2G2Q']:
            samples_signal_2017[f"ttH{br}_M{m}_ctau{ct}_UL17"] = {
                'dataset': f"local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2017/ttH_ttSemiLep_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}mm_v2.root",
                'triggers': ["(1)"],
                'veto_triggers': [],
                'era': '2017',
                'sigma': 1.0,
                'customNanoAOD': True}

# ttH samples:
for m in [20,50]:
    for ct in [1000]:
        for br in ['2G2Q']:
            samples_signal_2017[f"ttH{br}_M{m}_ctau{ct}_UL17"] = {
                'dataset': f"local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2017/ttH_ttSemiLep_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}mm.root",
                'triggers': ["(1)"],
                'veto_triggers': [],
                'era': '2017',
                'sigma': 1.0,
                'customNanoAOD': True}
