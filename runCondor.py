from optparse import OptionParser
import os
import time
import importlib
parser = OptionParser()

parser.add_option("-a", "--analysis", dest="analysis",
                  help="analysis to run",type='str',default='VH')

parser.add_option("-o", "--eosdir", dest="eos",default='root://cmseos.fnal.gov//store/user/bachtis/analysis',
                  help="EOS output Directory")

(options, args) = parser.parse_args()

samp = importlib.import_module('analysis.{}samples'.format(options.analysis))


datasets= args



os.system('tar --overwrite --exclude="common/__pycache__" --exclude="analysis/__pycache__"  -czvf /tmp/sandbox.tar.gz .')
os.system('mv /tmp/sandbox.tar.gz .')


for d,info in samp.samples.items():    
    if (not ('jobs' in info)):
        shell="""#!/bin/sh
        echo starting script
        tar -xzvf sandbox.tar.gz
        mkdir .dasmaps
        source ./setup.sh
        python3 run.py {dataset}
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
        Output = condor_{dataset}_$(Cluster)_$(Process).stdout
        Error = condor_{dataset}_$(Cluster)_$(Process).stderr
        Log = condor_{dataset}_$(Cluster)_$(Process).log
        Queue 1
        """.format(dataset=d)
        f=open("{dataset}_condor.jdl".format(dataset=d),"w")
        f.write(condor)
        f.close()
        
        os.system('condor_submit {dataset}_condor.jdl'.format(dataset=d))
        time.sleep(0.1)
    else:
        for i in range(0,info['jobs']):
            shell="""#!/bin/sh
            echo starting script
            tar -xzvf sandbox.tar.gz
            mkdir .dasmaps
            source ./setup.sh
            python3 run.py {dataset} -s {splitFactor} -p {part}
            echo python done
        
            OUTDIR={eos}
            echo "xrdcp output for condor to "
            echo $OUTDIR
            for FILE in *.root
            do
            xrdcp -f ${{FILE}} ${{OUTDIR}}/{dataset}_{part}_${{FILE}} 2>&1
            XRDEXIT=$?
            if [[ $XRDEXIT -ne 0 ]]; then
            rm *.root
            echo "exit code $XRDEXIT, failure in xrdcp"
            exit $XRDEXIT
            fi
            rm ${{FILE}}
            done
            """.format(dataset=d,eos=options.eos,splitFactor=info['jobs'],part=i)
            print(shell)
            f=open("{dataset}_{part}_condor.sh".format(dataset=d,part=i),"w")
            f.write(shell)
            f.close()
            os.system('chmod +x {dataset}_{part}_condor.sh'.format(dataset=d,part=i))
            
            condor="""
            universe = vanilla
            Executable = {dataset}_{part}_condor.sh
            should_transfer_files = YES
            Transfer_Input_Files = sandbox.tar.gz
            when_to_transfer_output = ON_EXIT
            Output = condor_{dataset}_{part}_$(Cluster)_$(Process).stdout
            Error = condor_{dataset}_{part}_$(Cluster)_$(Process).stderr
            Log = condor_{dataset}_{part}_$(Cluster)_$(Process).log
            Queue 1
            """.format(dataset=d,part=i)
            f=open("{dataset}_{part}_condor.jdl".format(dataset=d,part=i),"w")
            f.write(condor)
            f.close()
            
            os.system('condor_submit {dataset}_{part}_condor.jdl'.format(dataset=d,part=i))
            time.sleep(0.1)
            
