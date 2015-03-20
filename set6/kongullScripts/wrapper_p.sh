#!/bin/bash

# Name job 'poisson'
#PBS -N wrapper_jc

# Allocate two nodes with 12 processors from the default resources
#PBS -lnodes=1:ppn=1:default

# Expect to run up to 5 minutes
#PBS -lwalltime=00:1:00

# Memory per process
#PBS -lpmem=2000MB

# Run on the freecycle account
#PBS -A freecycle

# Run in the optimist queue by default
#PBS -q optimist

# Join stdout and stderr output to one file
#PBS -j oe

# Change directory to dir with the job script
#cd ${PBS_O_WORKDIR}

# prepare output file
echo "n, np, nt, maxError, used time" > out.txt



# Run with differing amounts of p and t st p*t = 36
qsub -t 1-36 job_varyp.sh

