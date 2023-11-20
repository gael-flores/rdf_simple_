# rdf_simple
RDataFrame analysis framework for CMS data.
Examples using Displaced Di-photon analysis 
## Quick start
- Set up the environment
  ```
  voms_proxy_init --voms cms --valid 100:00
  source setupo.sh
  ```
- Add samples in analysis/ddpSamples.py
- Modify the analysis code in the analysis/ddp.py. The running script is the `analysis()` function 
- Run interactive job ion a sample
  ```
  python run.py ZH4G_M30_ctau1000

  ```
- Submit datasets to condor
```
python runCondor.py -o root://cmseos.fnal.gov//store/user/<your username>/analysis  ZH4G_M30_ctau1000 SingleMuon_Run2018A ... etc
```
Make sure that you have write access to the EOS area you define with -o . This is where the output files will be stored 


