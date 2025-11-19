samples_signal_2022 = {}

# TODO: if everything is working as expected, move the samples to a public directory

# ggZH, WplusH, WminusH, ZH samples: 
for v in ['ggZ', 'Z', 'Wplus', 'Wminus']:
    for m in [7, 15, 20, 30, 40, 50, 55]:
        for ct in [0, 10, 20, 50, 100, 1000]:
            for br in ['2G2Q']:
                samples_signal_2022["{}H{}_M{}_ctau{}".format(v, br, m, ct)] = {
                    'dataset': "local:root://eosuser.cern.ch//eos/user/e/ebenjami/b-sandbox/CMSSW_12_4_11/src/php-plots/LLP-RUN3-SAMPLES-og/{}H/2g2d/Run3-Summer22EE-{}HToHTo2LongLivedTo{}_mass{}_ctau{}_GEN.root".format(v, v, '2G2D', ms, ct),
                    'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                    'veto_triggers': [],
                    'era': '2022',
                    'sigma': 1.0,
                    'customNanoAOD': True}


# ttH samples:
for m in [7, 15, 20, 30, 40, 50, 55]:
    for ct in [0, 10, 20, 50, 100, 1000]:
        for br in ['2G2Q']:
            samples_signal_2022["ttH{}_M{}_ctau{}".format(br, m, ct)] = {
                'dataset': "local:root://eosuser.cern.ch//eos/user/e/ebenjami/b-sandbox/CMSSW_12_4_11/src/php-plots/LLP-RUN3-SAMPLES-og/ttH/2g2d/Run3-Summer22EE-ttToHTo2LongLivedTo{}_mass{}_ctau{}_GEN.root".format('2G2D', ms, ct),
                'triggers': ['(HLT_IsoMu24||HLT_Ele30_WPTight_Gsf)'],
                'veto_triggers': [],
                'era': '2022',
                'sigma': 1.0,
                'customNanoAOD': True}

## ttH samples:
#or m in [15, 20, 30, 40, 50, 55]:
 #   for ct in [1000]:
 #       for br in ['2G2Q']:
 #           samples_signal_2018["ttH{}_M{}_ctau{}_UL18".format(br, m, ct)] = {
 #               'dataset': "local:root://cmsxrootd.fnal.gov//store/user/tlam/DDP/UL2018/ttH_ttSemiLep_HTo2LongLivedTo{}_MH-125_MFF-{}_ctau-{}mm.root".format(br, m, ct),
 #               'triggers': ['(HLT_IsoMu24||HLT_Ele32_WPTight_Gsf)'],
 #               'veto_triggers': [],
 #               'era': '2018',
 #               'sigma': 1.0,
 #               'customNanoAOD': True}
