#!/bin/bash
#$ -cwd
#$ -l h_vmem=1G
# -l nice_cpu=1
# -m ea


# -l h_cpu=2:0:0 
# -l h_rt=3:0:0 



export PATH=/home/miladim/miniconda2/bin:$PATH
export PYTHONNOUSERSITE=True
source activate mutarna

# all parameters for the shell script
script=$1
shift

## after shifting perscript is invoked for each parameter set


if [ -z "$SGE_TASK_ID" ]; then
	echo "No SGE environment found! Set SGE_TASK_ID to 0!"
	let SGE_TASK_ID=0
fi

if [ -z "$JOB_ID" ]; then
        echo "No SGE environment found! Set JOB_ID to 1!"
	        let JOB_ID=1
fi

if ! [ -z "$PBS_JOBID" ]; then
        echo "PBS environment found instead of SGE! Set JOB_ID to PBS_JOBID!"
	        let JOB_ID=$PBS_JOBID
fi


params=$@;
$script $params $(expr $SGE_TASK_ID - 1)

