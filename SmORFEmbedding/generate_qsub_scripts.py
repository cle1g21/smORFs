"""
generate_qsub_scripts.py
-----------------------
Generates qsub job scripts for each main analysis step in the pipeline.
Each script activates the conda environment and runs the corresponding Python script.
"""
import os

# List of main analysis scripts and their default arguments (if any)
JOBS = [
    ("generate_canonical_fasta.py", ""),
    ("generate_all_transcripts_fasta.py", ""),
    ("summarize_and_classify_fasta.py", "-i data/canonical_genes.fasta"),
    ("analyze_smorfs.py", ""),
    ("embed_fasta_dnabert.py", "-i data/all_transcripts.fasta"),
    # One job per model type for embed_fasta_transformer.py
]

EMBED_TRANSFORMER_MODELS = [
    "dnabert",
    "nucleotide-transformer",
    "esm2",
    "hyenadna",
]

# Ensure output directories exist
os.makedirs("qsub_scripts", exist_ok=True)
os.makedirs("logs", exist_ok=True)

# Template for qsub job script
QSUB_TEMPLATE = """#!/bin/bash
#$ -N {job_name}
#$ -cwd
#$ -V
#$ -l h_rt=24:00:00
#$ -l h_vmem=16G
#$ -o logs/{job_name}.out
#$ -e logs/{job_name}.err

source ~/.bashrc
conda activate canonical-fasta

python {script} {args}
"""

# Generate a qsub script for each job
for script, args in JOBS:
    job_name = os.path.splitext(script)[0]
    qsub_script = QSUB_TEMPLATE.format(job_name=job_name, script=script, args=args)
    script_path = os.path.join("qsub_scripts", f"{job_name}.sh")
    with open(script_path, "w") as f:
        f.write(qsub_script)
    os.chmod(script_path, 0o755)

# Generate a qsub script for each model type for embed_fasta_transformer.py
for model_type in EMBED_TRANSFORMER_MODELS:
    job_name = f"embed_fasta_transformer_{model_type}"
    args = f"--type {model_type} --input data/all_transcripts.fasta"
    qsub_script = QSUB_TEMPLATE.format(job_name=job_name, script="embed_fasta_transformer.py", args=args)
    script_path = os.path.join("qsub_scripts", f"{job_name}.sh")
    with open(script_path, "w") as f:
        f.write(qsub_script)
    os.chmod(script_path, 0o755)

# Generate a qsub script for each protein model for embed_protein_fasta_transformer.py
PROTEIN_MODELS = [
    ("protbert", "--model protbert"),
    ("prott5", "--model prott5"),
]
for model_name, model_arg in PROTEIN_MODELS:
    job_name = f"embed_protein_fasta_transformer_{model_name}"
    args = f"{model_arg}"
    qsub_script = QSUB_TEMPLATE.format(job_name=job_name, script="embed_protein_fasta_transformer.py", args=args)
    script_path = os.path.join("qsub_scripts", f"{job_name}.sh")
    with open(script_path, "w") as f:
        f.write(qsub_script)
    os.chmod(script_path, 0o755)

print("Qsub scripts generated in qsub_scripts/.")
print("Submit a job with: qsub qsub_scripts/<job>.sh") 