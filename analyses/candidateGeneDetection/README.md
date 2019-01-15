Genomes used:

https://www.cottongen.org/species/Gossypium_hirsutum/jgi-AD1_genome_v1.1
https://www.cottongen.org/species/Gossypium_raimondii/jgi_genome_221

Accessed 01/10/2019

The arboreum genome is from https://www.nature.com/articles/s41588-018-0116-x
We acquired our copy directly from the authors. We have also asked CottonGen to host it.

QTLcandidateGenes.bash is used to run gmap and find marker locations. The download link for the genome is listed within this file. The marker sequences are found in QTLtargetseqs.fasta. gmap results are in QTLtargetseqs.gmap.* 

A preliminary sorted table is found in QTLcandidateGenes.xlsx, which was inspected for each instance of two paths. If the paths were similar in score, each was considered viable and the chromosome which matched the majority of the other markers was used. If the scores/lengths were remarkably different, the longest and/or best score was selected to represent that marker.

Candidate genes are recovered by bedtools v2.27.1, spack module bedtools2/2.27.1-s2mtpsu
```
cat QTL.bed | while read line; do echo $line | tr ' ' '\t' | bedtools intersect -a stdin -b Tx-JGI_G.hirsutum_v1.1.geneONLY.gff3 >> QTL.candidateGenes.recovered

cat QTL.bed | while read line; do echo $line | tr ' ' '\t' | bedtools intersect -a stdin -b Tx-JGI_G.hirsutum_v1.1.geneONLY.gff3 >> QTL.candidateGenes.recovered.named
```
Numbers of Candidate loci recovered per QTL is given by
```
cut -f 4 QTL.candidateGenes.recovered | sort | uniq | while read line; do echo $line `grep -c $line QTL.candidateGenes.recovered` >> QTL.candidateGenes.recovered.number; done
```

Expression information is from Jing's ongoing cis-trans work. Only Tx2094 v Maxxa comparisons are used (@ 10 and 20 dpa)
```
for a in *.txt; do awk '{print FILENAME (NF?"\t":"") $0}' $a | sed 's/"//g' | sed '/NA/d' | sed 's/[.]txt//g' | sed 's/vs/_/g' | sed '/mt_/d' | sed '/Z/d'>> all.expression; done

awk ' { if ( $8 <= 0.005 ) print $0 } ' all.expression > all.expression.padj0.005
cut -f3 -d '=' QTL.candidateGenes.recovered.named | sort | uniq | while read line; do grep $line all.expression.padj0.005 >> candidateGene.expression; done
# of the 13,274 genes under a QTL, 1063 have expression differences between TX2094 and Maxxa
```
