samples_signal_2018 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [20,50]:
        for ct in [0, 100, 1000]:
            for br in ['2G2Q']:
                samples_signal_2018[f"{v}H{br}_M{m}_ctau{ct}_UL18"] = {
                    'dataset': f"local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2018/{v}H_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}cm.root",
                    'triggers': ["(1)"],
                    'veto_triggers': [],
                    'era': '2018',
                    'sigma': 1.0,
                    'customNanoAOD': True}

# ttH samples:
for m in [20,50]:
    for ct in [0, 100]:
        for br in ['2G2Q']:
            samples_signal_2018[f"ttH{br}_M{m}_ctau{ct}_UL18"] = {
                'dataset': f"local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2018/ttH_ttSemiLep_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}mm_v2.root",
                'triggers': ["(1)"],
                'veto_triggers': [],
                'era': '2018',
                'sigma': 1.0,
                'customNanoAOD': True}

# ttH samples:
for m in [20,50]:
    for ct in [1000]:
        for br in ['2G2Q']:
            samples_signal_2018[f"ttH{br}_M{m}_ctau{ct}_UL18"] = {
                'dataset': f"local:root://cmsxrootd.fnal.gov//store/user/gfloresa/DDP/UL2018/ttH_ttSemiLep_HTo2LongLivedTo{br}_MH-125_MFF-{m}_ctau-{ct}mm.root",
                'triggers': ["(1)"],
                'veto_triggers': [],
                'era': '2018',
                'sigma': 1.0,
                'customNanoAOD': True}
