#!/bin/bash


#SBATCH --nodes=2
#SBATCH --mem=200G 
#SBATCH --cpus-per-task=64
#SBATCH --time=24:00:00
#SBATCH --partition=amd
#SBATCH --export=ALL




module load apptainer

apptainer exec -H $PWD \
--unsquash \
--bind /iridisfs/ddnb/Ahmed/NN/bionemo:/iridisfs/ddnb/Ahmed/NN/bionemo \
--bind /iridisfs/ddnb/Ahmed/data/scalebio_600K_pbmc:/iridisfs/ddnb/Ahmed/data/scalebio_600K_pbmc \
--bind /iridisfs/ddnb/Ahmed/new_cells:/iridisfs/ddnb/Ahmed/new_cells \
--bind /iridisfs/ddnb/Ahmed/BioCreator:/iridisfs/ddnb/Ahmed/BioCreator \
/iridisfs/ddnb/images/bionemo-framework_nightly--bionemo1.sif \
python $1


