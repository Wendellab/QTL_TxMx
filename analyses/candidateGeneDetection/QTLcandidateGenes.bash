### Candidate gene analysis, QTLs
# Genome used : https://www.cottongen.org/species/Gossypium_hirsutum/jgi-AD1_genome_v1.1
# targets: QTLtargetseqsUpdated.fasta

module load gmap-gsnap/2018-07-04-gtu46xu
module load (bedops)
module load (bedtools)

### 0. fix stupid QTL file
sed 's/\[A\/T\]/W/g' QTLtargetseqsUpdated.fasta | sed 's/\[A\/C\]/M/g' | sed 's/\[A\/G\]/R/g' | sed 's/\[C\/G\]/S/g' | sed 's/\[C\/T\]/Y/g' | sed 's/\[G\/T\]/K/g' | sed 's/\[A/M/g' | sed 's/\/C\]//g' > QTLtargetseqsFixed.fasta # this file is renamed simply QTLtargetseqs.fasta, for simplicity

#### 1. map targets to genome with gmap
# downloaded Jan 4 2019
curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/Tx-JGI_G.hirsutum_AD1genome/assembly/Tx-JGI_G.hirsutum_v1.1.fa.gz
curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/Tx-JGI_G.hirsutum_AD1genome/genes/Tx-JGI_G.hirsutum_v1.1.gene.gff3.gz
gunzip *

gmap_build -D . -d TM1 Tx-JGI_G.hirsutum_v1.1.fa
gmap -D . -d TM1 -t 150 -n 2 -f 4 QTLtargetseqsFixed.fasta > QTLtargetseqsUpdated.gmap.2path.out # too stringent; only finds 63 loci in both subgenomes

makeblastdb -in Tx-JGI_G.hirsutum_v1.1.fa -dbtype nucl -out TM1 

### map targets to D5 genome
mkdir D5_location
cd D5_location
cp ../QTLtargetseqsFixed.fasta .

gmap_build -D . -d D5 Dgenome2_13.fasta
gmap -D . -d D5 -t 150 -n 2 -f 4 QTLtargetseqsFixed.fasta > QTLtargetseqs.D5.gmap.out 

### map targets to A2 genome
mkdir A2_location
cd A2_location
cp ../QTLtargetseqsFixed.fasta .

gmap_build -D . -d A2 A2Du_26.fasta
gmap -D . -d A2 -t 150 -n 2 -f 4 QTLtargetseqsFixed.fasta > QTLtargetseqs.A2.gmap.out 

### recover arabidopsis homologs

























### code for blast, in case I want it later
blastn -db TM1 -query QTLtargetseqsFixed.fasta -out QTLtargetseqsUpdated.blastn.out -evalue 1e-50 -outfmt "6 qseqid sseqid sstart pident evalue bitscore length" -perc_identity 90 -num_threads 150

cp QTLtargetseqsUpdated.blastn.out QTLtargetseqsUpdated.blastn.multiple.out

cut -f1 QTLtargetseqsUpdated.blastn.out | sort | uniq | cut -f2 -d '/' | while read line; do 
num=`grep -c $line QTLtargetseqsUpdated.blastn.out`
if [ "$num" -lt "2" ]
then
    grep $line QTLtargetseqsUpdated.blastn.out >> QTLtargetseqsUpdated.blastn.one_hit.out
	sed -i "/$line/d" QTLtargetseqsUpdated.blastn.multiple.out
fi
done


















#### 2a. convert mapped sam file to bed file

sam2bed < foo.sam2bed > foo.bed

#### 3. extract marker intervals from table 2
grep "first" or "last" from foo.bed
bedtools intersect to figure out what is between them 



