#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=40
#SBATCH --job-name=asffast_nextflow_sbatch
#SBATCH --mem=90gb
#SBATCH --time=02:00:00
#SBATCH --output=asffast_nextflow_sbatch_%j.log

module load python37
export NXF_CLUSTER_SEED=$( shuf -i 0-16777216 -n 1 )

samplename=$( basename $1 )

echo -e "Enqueuing job for data directory: ${samplename}, output directory: ${2}, reference: ${3}"
# Remember to add in singularity container path as -s argument once we get cachedir setup
srun slss-asffast/bin/asffast.py -i "$1" -o "$2" -r "$3" --slurm

echo -e "Job finished for data directory ${samplename}"
