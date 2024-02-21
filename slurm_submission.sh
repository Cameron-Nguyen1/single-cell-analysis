#!/bin/bash
#SBATCH -N 1
#SBATCH -n 2
#SBATCH --mem=80g
#SBATCH -o Mylog.out
#SBATCH --mail-type=END,FAIL # notifications for job done & fail
#SBATCH --mail-user=myemail@host.ext
#SBATCH -t 48:00:00

module load r/4.2.2
wd=$(pwd)
samples=/path/to/10x/counts

Rscript ${wd}/sc_mergeqc.r --hs_or_mm mm --samples $samples
Rscript ${wd}/sc_normal.r
Rscript ${wd}/sc_sct.r --hs_or_mm mm
Rscript ${wd}/sc_harm.r --hs_or_mm mm
Rscript ${wd}/sc_azimuth.r --wd $wd --samples $samples