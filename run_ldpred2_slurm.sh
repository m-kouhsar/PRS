#!/bin/bash

#SBATCH --export=ALL
#SBATCH --job-name=LDpred2  # job name
#SBATCH --output=LDpred2.log  # output
#SBATCH --error=LDpred2.err  # errors
#SBATCH -D . # set working directory to .
#SBATCH -p mrcq # submit to the parallel test queue
#SBATCH --time=5:00:00  # walltime
#SBATCH --nodes=1 # specify number of nodes.
#SBATCH --ntasks-per-node=16 # specify number of processors.
#SBATCH --mail-type=END # send email at job completion
#SBATCH --mail-user=mail@sample.co # email address


fileGeno=./EUR.bed
fileGenoRDS=./EUR.rds
fileSumstats=./Height.gwas.txt
fileOut=./EUR

LDPRED2_REF=./ldpred2_ref
SCRIPTDIR=.

module load R

if [ ! -f $fileGenoRDS ]
then
  # Convert plink binary format to rds,bk
  Rscript ${SCRIPTDIR}/createBackingFile.R --file-input $fileGeno --file-output $fileGenoRDS
fi
# Generate PRS using LDPRED2
Rscript ${SCRIPTDIR}/ldpred2.R \
 --ldpred-mode inf \
 --col-stat OR \
 --col-stat-se SE \
 --stat-type OR \
 --geno-file-rds $fileGenoRDS \
 --sumstats $fileSumstats \
 --out $fileOut \
 --ld-file ${LDPRED2_REF}/ldref_hm3_plus/LD_with_blocks_chr@.rds \
 --ld-meta-file ${LDPRED2_REF}/map_hm3_plus.rds \
 --col-A2 A2 \
 --col-bp BP \
 --trait-type q \
 --geno-impute-zero
