#!/bin/bash
#SBATCH --job-name=polyR30_293
#SBATCH --dependency=singleton
#SBATCH --nodes=1
#SBATCH --partition=qgpu
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8888
#SBATCH -t 24:00:00
#SBATCH -o polyR30_293.out
#SBATCH -e polyR30_293.err

source /home/people/sorbul/.bashrc
module purge
conda activate calvados

echo $SLURM_CPUS_PER_TASK
echo $SLURM_JOB_NODELIST

python run.py --config config.yaml