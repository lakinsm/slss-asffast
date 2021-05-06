#!/usr/bin/env bash
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --partition=cpu_compute
#SBATCH --cpus-per-task=48
#SBATCH --job-name=asffast_nextflow
#SBATCH --mem=160gb
#SBATCH --time=10:00:00
#SBATCH --output=asffast_nextflow_%j.log

module load python37
export NXF_CLUSTER_SEED=$( shuf -i 0-16777216 -n 1 )

echo -e "\n$(date): Enqueuing job"
echo -e "Executable: ${1}"
echo -e "Data directory: ${2}"
echo -e "Output directory: ${3}"
echo -e "Reference: ${4}"
echo -e "Container: ${5}"
echo -e "Working directory: ${6}"
srun "$1" -i "$2" -o "$3" -r "$4" -s "$5" -w "$6" --slurm
echo -e "\n$(date): Job finished for data directory ${2}\n"
