#! /bin/bash
#SBATCH -J scheduler
#SBATCH -t 5-00:00:00
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -p all
#SBATCH --mem=2gb

cd /cluster/home/selghamr/workflows/ExomeSeq
condaprefix='/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda'

snakemake \
--jobs 6 \
--profile slurm \
--cluster-config slurm/cluster.json \
--conda-frontend conda \
--use-conda \
--use-singularity \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:///cluster/home/selghamr/snakemake_wrappers/' \
--rerun-incomplete
