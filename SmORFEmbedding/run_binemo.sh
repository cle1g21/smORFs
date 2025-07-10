#!/bin/bash

#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --time=6:00:00
#SBATCH --partition=a100
#SBATCH --export=ALL
#SBATCH --gpus=1



module load apptainer
module load cuda
apptainer run -H $PWD \
--nv \
--bind "$PWD:/workspace" \
/iridisfs/ddnb/images/bionemo-framework_nightly \
python $1


