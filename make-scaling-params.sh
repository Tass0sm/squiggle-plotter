#!/bin/bash
#SBATCH --job-name=make_vir_scaling_file
#SBATCH --time=05:00:00
#SBATCH --nodes=1 --ntasks-per-node=20
#SBATCH --account=pas1405
#SBATCH --mail-type=ALL

set -x

SCRIPT="/users/PAS1405/tassosm/Desktop/tombo-plot/fit-scaling-params.py"

KMER_MODEL_PATH="/users/PAS1405/tassosm/Desktop/kmer_models/r9.4_180mv_70bps_5mer_RNA/template_median69pA.model"
GUPPY_DIR="/users/PAS1405/tassosm/Desktop/common-data/guppy_vir"
SEQ_SUM_PATH="$GUPPY_DIR/sequencing_summary.txt"
WORKSPACE_DIR="$GUPPY_DIR/workspace"

module load python
source activate tombo-env

mkdir $SLURM_SUBMIT_DIR/scaling-job-results

cd $WORKSPACE_DIR

/usr/bin/time python $SCRIPT $SEQ_SUM_PATH $KMER_MODEL_PATH > $TMPDIR/scaling.txt

cp $TMPDIR/scaling.txt $SLURM_SUBMIT_DIR/scaling-job-results
