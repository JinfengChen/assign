export R_LIBS=/home/jichen/software/BETSY/install/envs/ASSIGN/lib/R/library
ln -s ../Expression_Tables/expression.gene.tpm ./
sbatch --array 1-10 Assign.sh

#combat
export R_LIBS=/home/jichen/software/BETSY/install/envs/BulkRNAseq/lib/R/library/

