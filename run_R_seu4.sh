#!/bin/bash
#SBATCH --job-name=apptainer_R
#SBATCH --output=job_output_%j.txt
#SBATCH --error=job_error_%j.txt
#SBATCH --nodes=1
#SBATCH --mem=150G
#SBATCH --time=60:00:00

# Load Apptainer
module load apptainer

# Run your R script inside the container
apptainer exec --no-mount fuse \
  -H "$PWD" \
  --bind /lyceum/cle1g21:/lyceum/cle1g21 \
  /lyceum/cle1g21/ddnb_r_full.sif \
  Rscript /lyceum/cle1g21/smORFs/TCGA/BRCA/cellxgene/Gareth_job_submission.r
