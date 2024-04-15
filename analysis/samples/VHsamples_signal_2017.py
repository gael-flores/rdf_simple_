samples_signal_2017 = {}
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [15, 20, 30, 40, 50, 55]:
        for ct in [0, 10, 20, 50, 100, 1000]:
            for br in ['2G2Q', '4G']:
                samples_signal_2017['{}H{}_M{}_ctau{}_UL17'.format(v, br, m, ct)] = {
                    'dataset': 'local:root://cmsxrootd.fnal.gov//store/user/tlam/DDP/UL2017/{}H_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}cm.root'.format(v, br, m, ct),
                    'triggers': ["(HLT_IsoMu27||HLT_Ele32_WPTight_Gsf_L1DoubleEG)"],
                    'veto_triggers': [],
                    'era': '2017',
                    'sigma': 1.0,
                    'customNanoAOD': True}
