from optparse import OptionParser
import os
import time
parser = OptionParser()
parser.add_option("-o", "--eosdir", dest="eos",default='root://cmseos.fnal.gov//store/user/bachtis/analysis',
                  help="EOS output Directory")
(options, args) = parser.parse_args()

datasets= args
os.system('tar -czvf /tmp/sandbox.tar.gz .')
os.system('mv /tmp/sandbox.tar.gz .')


for d in datasets:     
    shell="""#!/bin/sh
    echo starting script
    tar -xzvf sandbox.tar.gz
    source ./setup.sh
    python run.py {dataset}
    echo python done

    OUTDIR={eos}
    echo "xrdcp output for condor to "
    echo $OUTDIR
    for FILE in *.root
    do
    xrdcp -f ${{FILE}} ${{OUTDIR}}/{dataset}_${{FILE}} 2>&1
    XRDEXIT=$?
    if [[ $XRDEXIT -ne 0 ]]; then
    rm *.root
    echo "exit code $XRDEXIT, failure in xrdcp"
    exit $XRDEXIT
    fi
    rm ${{FILE}}
    done
    """.format(dataset=d,eos=options.eos)
    print(shell)
    f=open("{dataset}_condor.sh".format(dataset=d),"w")
    f.write(shell)
    f.close()
    os.system('chmod +x {dataset}_condor.sh'.format(dataset=d))

    condor="""
    universe = vanilla
    Executable = {dataset}_condor.sh
    should_transfer_files = YES
    Transfer_Input_Files = sandbox.tar.gz
    when_to_transfer_output = ON_EXIT
    Output = {dataset}_$(Cluster)_$(Process).stdout
    Error = {dataset}_$(Cluster)_$(Process).stderr
    Log = {dataset}_$(Cluster)_$(Process).log
    Queue 1
    """.format(dataset=d)
    f=open("{dataset}_condor.jdl".format(dataset=d),"w")
    f.write(condor)
    f.close()
    
    os.system('condor_submit {dataset}_condor.jdl'.format(dataset=d))
    time.sleep(0.1)
