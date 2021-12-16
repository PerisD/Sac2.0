#!/bin/bash
#SBATCH --job-name="[STRAINNAME]-BUSCO"
#SBATCH --chdir=[WORKINGDIR]
#SBATCH --ntasks=[THREADS]
#SBATCH --time=[TIME]
#SBATCH --mail-type=ALL
#SBATCH --mail-user=david.perisnavarro@gmail.com
#SBATCH --output=SLURM-outputs/Output_%j_%x.out
#SBATCH --error=SLURM-outputs/Output_%j_%x.err


module load miniconda/3
source /storage/home/vlc81/vlc81510/.bashrc
source activate /storage/scratch/vlc81/vlc81510/conda_env/biopython3
which python

PYTHONPATH=/storage/scratch/vlc81/vlc81510/conda_env/biopython3/lib/python3.6/site-packages/:$PYTHONPATH
export PYTHONPATH

pwd; hostname; date

cd [WORKINGDIR]

busco -i [ASSEMBLYPATH] -o [STRAINNAME] -l [BUSCODB] -m geno -c [THREADS] --augustus --augustus_species [SPPTAG] --offline

date



