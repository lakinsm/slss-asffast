#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --cpus-per-task=40
#SBATCH --job-name=asffast_nextflow_sbatch
#SBATCH --mem=120gb
#SBATCH --time=04:00:00
#SBATCH --output=asffast_nextflow_sbatch_%j.log

module load python37
module load singularity
export NXF_CLUSTER_SEED=$( shuf -i 0-16777216 -n 1 )

samplename=$( basename $1 )

echo -e "\n$(date): Enqueuing job for data directory: ${samplename}, output directory: ${2}, reference: ${3}"
srun slss-asffast/bin/asffast.py -i "$1" -o "$2" -r "$3" -s singularity_containers/slss-asffast_alpha1.sif -w /lustrefs/scratch_dir/work --slurm
echo -e "$(date): Job finished for data directory ${samplename}\n"
