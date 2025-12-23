import array

mass_list = [7,15,20,30,40,50,55]
ctaus = [0,10,20,50,100,1000]
channels = ['WH','ZH','ggZH','ttH']
#channels = ['total']

# Run 2 luminosity
lumi = {'2018': "59830",
        '2017': "41480",
        '2016': "36310",
        'Run2': '137620',
        'HEM': "38750",
        'preHEM': "21080"}

lumifb = {'2018': "59.83",
          '2017': '41.48',
          '2016': '36.31',
          'Run2': '137.62'}

# Run 2 luminosity (integers) (too lazy to manually convert strings to ints)
intLumi = {'2018': 59830,
           '2017': 41480,
           '2016': 36310,
           'Run2': 137620,
           'HEM': 38750,
           'preHEM': 21080}

# Run 2 luminosity uncertainty
lumiUnc = {'2018': 1.025,
           '2017': 1.023,
           '2016': 1.012}

xsecs = {'ggZ': "0.1057", #units in pb
        'Z': "0.8696",
        'Wplus': "0.84",
        'Wminus': "0.5328",
        'tt':"0.5071"} 

#source https://twiki.cern.ch/twiki/bin/view/LHCPhysics/CERNYellowReportPageAt13TeV
#asymetric
xsecUnc = {'WH'  : [0.005, -0.007],
           'ZH'  : [0.038, -0.031],
           'ggZH': [0.251, -0.189],
           'ttH' : [0.058, -0.092]}

#symmetric
pdfUnc = {'WH'  : 1.019,
          'ZH'  : 1.016,
          'ggZH': 1.024,
          'ttH' : 1.036}


#BRs from PDG
BRs = {"Z": "(.03363+.03366+.033696)",  #e⁺e⁻, μ⁺μ⁻, τ⁺τ⁻
       "W": "(.1063+.1071+.1138)",      #μν, eν, τν
       "ttSemiLeptonic": "2*(1 - (.1110 + .1140 + .107))*(.1110 + .1140 + .107)" #tt̄->W⁺W⁻bb̄(100%) -> W->(hadronic)W->(leptonic)bb̄->(jets) * 2
       }

# Scale factor stuff 
scaleFactors = {}
scaleFactors['W'] = {'ELE': ['Electron_recoSF_val[W_l1_idx]', 'Electron_idSF_val[W_l1_idx]', 'Electron_trigSF_val[W_l1_idx]'],
                     'MU': ['Muon_recoSF_val[W_l1_idx]', 'Muon_idSF_val[W_l1_idx]', 'Muon_isoSF_val[W_l1_idx]', 'Muon_trigSF_val[W_l1_idx]']
}

scaleFactors['Z'] = {'ELE': ['Electron_recoSF_val[Z_idx[0]]', 'Electron_recoSF_val[Z_idx[1]]', 'Electron_idSF_val[Z_idx[0]]', 'Electron_idSF_val[Z_idx[1]]', 'Electron_trigSF_val[Z_idx[0]]'],
                     'MU': ['Muon_recoSF_val[Z_idx[0]]', 'Muon_recoSF_val[Z_idx[1]]', 'Muon_idSF_val[Z_idx[0]]', 'Muon_idSF_val[Z_idx[1]]', 'Muon_isoSF_val[Z_idx[0]]', 'Muon_isoSF_val[Z_idx[1]]', 'Muon_trigSF_val[Z_idx[0]]']
}

scaleFactors['g'] = {}
scaleFactors['FakeRate'] = {}
for m in mass_list:
    scaleFactors['g'][m] = ['Photon_idSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_idSF_val[best_2g_idx2_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx1_m{}]'.format(m), 'Photon_pixSF_val[best_2g_idx2_m{}]'.format(m)]
    scaleFactors['FakeRate'][m] = [f"Photon_fakeRateSF_val[antiID_photon_idx_m{m}]"]




# Cuts for W->e+nu 2018 HEM
#cutsHEM = "(Electron_eta[W_l1_idx]>-1.3 || Electron_eta[W_l1_idx]<-3.0) && (Electron_phi[W_l1_idx]>-0.85 || Electron_phi[W_l1_idx]<-1.57)"

# Useful dictionaries
ana = {}
ana['W'] = {'MU': "wmn2g",
            'ELE': 'wen2g'}
ana['Z'] = {'MU': 'zmm2g',
            'ELE': 'zee2g'}
            
#Latex Headers
headers = {'W':{'ELE': "W#rightarrowe#nu",
                'MU' : "W#rightarrow#mu#nu",
                'l'  : "W#rightarrowl#nu"
                },
           'Z': {'ELE': "Z#rightarrowee",
                'MU' : "Z#rightarrow#mu#mu",
                'l'  : "Z#rightarrowll"
                }}


python_headers = {}
python_headers['W'] = {'ELE': r"$W\rightarrow e\nu$",
                'MU': r"$W\rightarrow\mu\nu$"}
python_headers['Z'] = {'ELE': r"$Z\rightarrow ee$",
                'MU': r"$Z\rightarrow\mu\mu$"}


#custom bins for unrolled 2d histograms
customBins = {
            'W': {
                'ELE': {
                     7:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    15:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    20:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    30:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    50:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                },
                'MU': {
                     7:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    15:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    20:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    30:array.array('d' , [0.0,25.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    50:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                }
            },
            'Z': {
                'ELE': {
                     7:array.array('d' , [0.0,36.0,72.0,110.0,220.0]),
                    15:array.array('d' , [0.0,36.0,72.0,110.0,220.0]),
                    20:array.array('d' , [0.0,110.0,138.0,165.0,192.0,220.0]),
                    30:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,350.0,370.0,440.0]),
                    50:array.array('d' , [0.0,50.0,110.0,138.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,110.0,138.0,192.0,220.0,248.0,275.0,330.0,360.0,440.0]),
                },
                'MU': {
                     7:array.array('d' , [0.0,50.0,80.0,110.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    15:array.array('d' , [0.0,50.0,80.0,110.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    20:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    30:array.array('d' , [0.0,50.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                    40:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,350.0,370.0,440.0]),
                    50:array.array('d' , [0.0,50.0,80.0,110.0,138.0,192.0,220.0,275.0,302.0,330.0,360.0,440.0]),
                    55:array.array('d' , [0.0,50.0,80.0,110.0,138.0,165.0,192.0,220.0,248.0,275.0,302.0,330.0,360.0,440.0]),
                }
            }
        }


 
