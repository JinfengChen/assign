#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=20G
#SBATCH --time=4:00:00
#SBATCH --output=Assign.sh.%A_%a.stdout
#SBATCH -p fast,all
#SBATCH --workdir=./


#export R_LIBS=/home/jichen/software/BETSY/install/envs/BulkRNAseq/lib/R/library/
export R_LIBS=/home/jichen/software/BETSY/install/envs/ASSIGN/lib/R/library

start=`date +%s`

CPU=$SLURM_NTASKS
if [ ! $CPU ]; then
   CPU=2
fi

N=$SLURM_ARRAY_TASK_ID
if [ ! $N ]; then
    N=1
fi

echo "CPU: $CPU"
echo "N: $N"

#FILE=`ls *_1.fastq.gz | grep _1\.fastq\.gz | head -n $N | tail -n 1`
#R1=$FILE
#R2=`echo $R1 | perl -p -e 's/_1\.fastq/_2.fastq/'`
#SAMPLE=${FILE%_1.fastq.gz}
#echo "File: $FILE"
#echo "Read: $R1 and $R2"
#echo "Sample: $SAMPLE"

#if [ ! -e $SAMPLE\.bam ]; then
#echo "mapping $SAMPLE ..."
#fi

#cat Assign.R | /home/jichen/software/BETSY/install/envs/ASSIGN/bin/R --slave
/home/jichen/software/BETSY/install/envs/ASSIGN/bin/Rscript Assign_steps_ABCB1_qsub.R $N

end=`date +%s`
runtime=$((end-start))

echo "Start: $start"
echo "End: $end"
echo "Run time: $runtime"

echo "Done"

