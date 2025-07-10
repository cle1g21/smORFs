#!/bin/bash

#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --mem=150G
#SBATCH --time=4:00:00




module load apptainer
module load cuda
apptainer run -H $PWD \
--nv \
--bind "$PWD:/workspace" \
/iridisfs/ddnb/images/bionemo-framework_nightly \
pip install $1

