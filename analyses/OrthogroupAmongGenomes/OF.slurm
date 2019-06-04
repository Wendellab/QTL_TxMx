#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 148:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 100   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "FILENAME"   # job name
#SBATCH --mem=100G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-OF.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
ml orthofinder/2.2.0-py2-tf3rsoz

orthofinder -f proteins -t 100
