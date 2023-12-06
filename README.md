# rdf_simple
RDataFrame analysis framework for CMS data.
Examples using Displaced Di-photon analysis 
## Quick start
- Set up the environment
  ```
  voms_proxy_init --voms cms --valid 100:00
  source setup.sh
  ```
- Add samples in analysis/ddpSamples.py
- Modify the analysis code in the analysis/ddp.py. The running script is the `analysis()` function 
- Run interactive job in a sample
  ```
  python3 run.py ZH4G_M30_ctau1000

  ```
- Submit datasets to condor
```
python3 runCondor.py -o -s 4 root://cmseos.fnal.gov//store/user/<your username>/analysis  ZH4G_M30_ctau1000 SingleMuon_Run2018A ... etc
```
Make sure that you have write access to the EOS area you define with -o . This is where the output files will be stored 
-s 4 implies that you want to split the sample into 4 jobs. If you do not want to split, ommit it.



