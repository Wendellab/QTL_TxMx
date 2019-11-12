#!/bin/bash

#Submit this script with: sbatch thefilename
# the pronto scheduler is pronto.las.iastate.edu
# note: change the memory, threads, wall, etc

#SBATCH -t 18:00:00   # walltime
#SBATCH -N 1   # number of nodes in this job
#SBATCH -n 15   # total number of processor cores in this job; each node has 272 cores, Nova has 36
#SBATCH -J "mummer"   # job name
#SBATCH --mem=20G # how much memory you need; each box has ~340G (legion) Nova has 190G or 380G
#SBATCH --output=slurm-mummer.out
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=corrinne@iastate.edu   # email address


# LOAD MODULES, INSERT CODE, AND RUN YOUR PROGRAMS HERE
# This file was run several times, commenting out the steps already completed

# compiled mummer4.05b
module load mummer #for any outside calls
module load parallel

parallel -j 15 --tmpdir tmp --colsep '\t' "./nucmer -p align {1} {2}" :::: input.files
ls *.fa.delta | parallel -j 15 './show-coords -b -r -L 15000 {} > {.}.filtered.15k.coords'

cat *.filtered.15k.coords > all.15k.coords
# manually sort columns into bed format

# manually construct Fang.Wang.bed from WangNatureGenetics2017, Table S7 and FangGenomeBiologyGB2017 Table S9
bedtools merge -i Fang.Wang.bed > Fang.Wang.merge.bed
bedtools merge -i all.15k.coords.bed -c 4 -o collapse > all.15k.coords.merge.bed
bedtools intersect -a Fang.Wang.merge.bed -b all.15k.coords.merge.bed -wb > selected.sweep.regions

# manually construct the selected sweep bed for the TM1 used here by taking the coordinates in the name of the "b" region and the matching "a" regions
# to select only the b regions matched by a to create `selected.sweep.regions.ZJC`

curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/Tx-JGI_G.hirsutum_AD1genome/genes/Tx-JGI_G.hirsutum_v1.1.gene.gff3.gz
gunzip Tx-JGI_G.hirsutum_v1.1.gene.gff3.gz

grep gene Tx-JGI_G.hirsutum_v1.1.gene.gff3 | sed '/Gohir[.]1Z/d' > Tx.gene.gff3

bedtools intersect -a selected.sweep.regions.ZJC -b Tx.gene.gff3 -wb | cut -f3 -d '=' > ZJC.genes.selected
