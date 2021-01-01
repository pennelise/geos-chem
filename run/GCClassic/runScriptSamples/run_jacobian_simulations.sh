#!/bin/bash

#SBATCH -J {RunName}
#SBATCH -c 1
#SBATCH -N 1
#SBATCH -p huce_intel
#SBATCH --mem 2000
#SBATCH -t 0-2:00

### Run directory
RUNDIR=$(pwd -P)

### Get current task ID
x=${SLURM_ARRAY_TASK_ID}

### The task is going to run a series of 4 jobs. For each job:
for iter in {0..3}; do

  ### Select a cluster Id
  xi=$((4*x+iter))

  ### Add zeros to the cluster Id
  if [ $xi -lt 10 ]; then
      xstr="000${xi}"
  elif [ $xi -lt 100 ]; then
      xstr="00${xi}"
  elif [ $xi -lt 1000 ]; then
      xstr="0${xi}"
  else
      xstr="${xi}"
  fi

   ### Run GEOS-Chem in the directory corresponding to the cluster Id
  cd  ${RUNDIR}/{RunName}_${xstr}
  ./{RunName}_${xstr}.run

done

exit 0
