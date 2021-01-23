#!/bin/csh
#PBS -N wardakgong
##PBS -o output.txt
#PBS -j oe
#PBS -q yossarian
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=32:mem=180GB
#PBS -m n
#PBS -V

#PBS -J 105-200:5

# Runs the simulation on the USYD Physics cluster for a range of `alpha`s

cd "$PBS_O_WORKDIR"
module load Anaconda2-4.0.0
source /usr/physics/nest/bin/nest_vars.csh

python wardakgong2021_sim.py 32 $PBS_ARRAY_INDEX

module unload Anaconda2-4.0.0
exit
