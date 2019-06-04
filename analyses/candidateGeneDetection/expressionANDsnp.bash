for a in *.txt; do awk '{print FILENAME (NF?"\t":"") $0}' $a | sed 's/"//g' | sed '/NA/d' | sed 's/[.]txt//g' | sed 's/vs/_/g' | sed '/mt_/d' | sed '/Z/d'>> all.expression; done

awk ' { if ( $8 <= 0.005 ) print $0 } ' all.expression > all.expression.padj0.005
cut -f3 -d '=' QTL.candidateGenes.recovered.named | sort | uniq | while read line; do grep $line all.expression.padj0.005 >> candidateGene.expression; done
# of the 13,274 genes under a QTL, 1063 have expression differences between TX2094 and Maxxa


cat QTL.bed | while read line; do echo $line | tr ' ' '\t' | bedtools intersect -a stdin -b Tx-JGI_G.hirsutum_v1.1.geneONLY.gff3 >> QTL.candidateGenes.recovered.named


###############################
mkdir Tx26
cp Tx26.snpeff.fa.gz Tx26/sequences.fa.gz
cp Tx26.snpeff.gtf.gz Tx26/genes.gtf.gz
echo "Tx26.genome : Tx26" >> snpEff.config
snpEff build -v Tx26 -gtf22 -dataDir . -c snpEff.config

########## strict SNP annotation
vcftools --vcf maxxa_tx2094_AD4.vcf --min-alleles 2 --max-alleles 2 --max-missing 1 --remove-indv AD4-16 --remove-indv AD4 --remove-indv AD4-W400 --out TxMx.filter.indel --recode --remove-indels

#CHROM POS ID REF ALT QUAL FILTER INFO FORMAT AD4_10 Maxxa SAMN04998813
# AD4 = GEN[0]
# Maxxa = GEN[1]
# TX2094 = GEN[2]
cat ../TxMx.filter.indel.recode.vcf | java -jar SnpSift.jar filter "( isVariant( (GEN[0].GT) ) | isVariant( (GEN[1].GT) ) | isVariant( (GEN[2].GT) ) & ( (GEN[0].GT) = (GEN[2].GT) ) & ( (GEN[0].GT) != (GEN[1].GT)) ) " > TxMx.final.vcf



snpEff ann -t -no-intergenic -no-intron -no-upstream -no-downstream -dataDir . -c snpEff.config Tx26 TxMx.final.vcf > TxMx.snpeff.vcf

cat TxMx.snpeff.vcf | java -jar SnpSift.jar filter "(exists ANN[*].EFFECT)" > TxMx.ann.vcf

grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq > genes.with.variants
grep Gohir TxMx.snpeff.vcf | sed '/synonymous_variant/d' | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq > genes.with.eff.variants
grep Gohir TxMx.snpeff.vcf | grep "missense" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq > genes.with.aa.variants
cut -f3 -d '=' QTL.candidateGenes.recovered.named | sort | uniq > QTL.candidates
comm -12 QTL.candidates genes.with.eff.variants > QTL.genes.with.eff.variants
comm -12 QTL.candidates genes.with.aa.variants > QTL.genes.with.aa.variants


echo subset number >> TxMx.polarized.snps
echo snps >> TxMx.polarized.snps
echo all_snps `grep -c "GT:AD" TxMx.snpeff.vcf` >> TxMx.polarized.snps
echo synonymous `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 2 -d '|' | grep -c "synonymous_variant"` >> TxMx.polarized.snps
echo missense `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 2 -d '|' | grep -c "missense_variant"` >> TxMx.polarized.snps
echo 5_UTR `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 2 -d '|' | grep -c "5_prime_UTR_variant"` >> TxMx.polarized.snps
echo 3_UTR `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 2 -d '|' | grep -c "3_prime_UTR_variant"` >> TxMx.polarized.snps
echo start `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 2 -d '|' | grep -c "start"` >> TxMx.polarized.snps
echo stop `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 2 -d '|' | grep -c "stop"` >> TxMx.polarized.snps
echo >> TxMx.polarized.snps 
echo genes >> TxMx.polarized.snps
echo genes `grep Gohir TxMx.snpeff.vcf | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_syn `grep Gohir TxMx.snpeff.vcf | grep "synonymous" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_aa `grep Gohir TxMx.snpeff.vcf | grep "missense" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_5UTR `grep Gohir TxMx.snpeff.vcf | grep "5_prime_UTR_variant" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_3UTR `grep Gohir TxMx.snpeff.vcf | grep "3_prime_UTR_variant" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_start `grep Gohir TxMx.snpeff.vcf | grep "start" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_stop `grep Gohir TxMx.snpeff.vcf | grep "stop" | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo genes_aff `grep Gohir TxMx.snpeff.vcf | sed '/synonymous_variant/d' | sed 's/^.*ANN=//g' | cut -f 4 -d '|' | sort | uniq | wc -l` >> TxMx.polarized.snps
echo QTL_genes_eff `comm -12 QTL.candidates genes.with.eff.variants | wc -l` >> TxMx.polarized.snps
echo QTL_genes_aa `comm -12 QTL.candidates genes.with.aa.variants | wc -l` >> TxMx.polarized.snps
echo >> TxMx.polarized.snps 
echo compared2expression >> TxMx.polarized.snps
echo expSNP_genes_eff `comm -12 QTL.genes.DGE QTL.genes.with.eff.variants | wc -l` >> TxMx.polarized.snps
echo expSNP_genes_aa `comm -12 QTL.genes.DGE QTL.genes.with.aa.variants | wc -l` >> TxMx.polarized.snps



################## Transcription factors under QTL

curl -O ftp://ftp.bioinfo.wsu.edu/species/Gossypium_hirsutum/Tx-JGI_G.hirsutum_AD1genome/functional/G.hirsutum_Tx-JGI_v1.1_genes2IPR.txt.gz
gunzip G.hirsutum_Tx-JGI_v1.1_genes2IPR.txt.gz
grep "transcription factor" G.hirsutum_Tx-JGI_v1.1_genes2IPR.txt > Gohir.TF
cp SNPaa/snpeff/QTL.candidates .
cut -f1 Gohir.TF | sed 's/[.][0-9]//g' | sed '/Z/d' | sort | uniq > TF.gene.names

comm -12 TF.gene.names QTL.candidates
 
comm -12 TF.gene.names QTL.genes.with.eff.variants
comm -12 TF.gene.names QTL.genes.with.aa.variants

comm -12 TF.gene.names QTL.genes.DGE

########### Candidate genes for phenotypes

cat QTL.fiberColor.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.fiberColor.genes.redundant; done
cat QTL.architecture.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.architecture.genes.redundant; done
cat QTL.fiberLength.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.fiberLength.genes.redundant; done
cat QTL.fiberQuality.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.fiberQuality.genes.redundant; done
cat QTL.flower.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.flower.genes.redundant; done
cat QTL.fruitingHabit.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.fruitingHabit.genes.redundant; done
cat QTL.seed.names | while read line; do grep $line QTL.candidateGenes.recovered.named | cut -f3 -d '=' | sort | uniq >> QTL.seed.genes.redundant; done

for a in *.redundant; do sort $a | uniq > ${a%.redundant}; done

for a in QTL.*.genes; do wc -l $a; done
#2078 QTL.architecture.genes
#3087 QTL.fiberColor.genes
#4329 QTL.fiberLength.genes
#968 QTL.fiberQuality.genes
#1926 QTL.flower.genes
#3266 QTL.fruitingHabit.genes
#3306 QTL.seed.genes

for a in QTL.*.genes; do comm -12 $a QTL.genes.DGE > $a.DGE; wc -l $a.DGE; comm -12 $a TF.gene.names > $a.TF; wc -l $a.TF; comm -12 $a QTL.genes.with.eff.variants > $a.eff; wc -l $a.eff; comm -12 $a QTL.genes.with.aa.variants > $a.aa; wc -l $a.aa; done

# 197 QTL.architecture.genes.DGE
# 13 QTL.architecture.genes.TF
# 957 QTL.architecture.genes.eff
# 594 QTL.architecture.genes.aa
# 236 QTL.fiberColor.genes.DGE
# 16 QTL.fiberColor.genes.TF
# 1298 QTL.fiberColor.genes.eff
# 816 QTL.fiberColor.genes.aa
# 345 QTL.fiberLength.genes.DGE
# 27 QTL.fiberLength.genes.TF
# 1930 QTL.fiberLength.genes.eff
# 1178 QTL.fiberLength.genes.aa
# 106 QTL.fiberQuality.genes.DGE
# 8 QTL.fiberQuality.genes.TF
# 452 QTL.fiberQuality.genes.eff
# 266 QTL.fiberQuality.genes.aa
# 134 QTL.flower.genes.DGE
# 10 QTL.flower.genes.TF
# 889 QTL.flower.genes.eff
# 555 QTL.flower.genes.aa
# 251 QTL.fruitingHabit.genes.DGE
# 12 QTL.fruitingHabit.genes.TF
# 1542 QTL.fruitingHabit.genes.eff
# 982 QTL.fruitingHabit.genes.aa
# 264 QTL.seed.genes.DGE
# 23 QTL.seed.genes.TF
# 1513 QTL.seed.genes.eff
# 921 QTL.seed.genes.aa



### genes from "JJ's list" ###

gmap -D . -d TM1 -t 10 -n 2 -f 4 JJgorai.fasta > cottonGenesOfInterest.gmap.2path.out #Use JJ's Gorai genes to find regions of interest in TM1 genome. 

gff2bed < Tx-JGI_G.hirsutum_v1.1.geneONLY.gff3 >> geneONLY.bed	#converted gff3 to bed file, in an attempt to get rid of an error that was actually induced by bad selection of options

bedtools intersect -a geneONLY.bed -b cottonGenesOfInterest.gmap.2path.out -wo -e -f 0.5 -F 0.5 > Gorai.Gohir.results
cut -f 4 Gorai.Gohir.results | while read line; do grep $line QTL.candidateGenes.recovered.named >> QTL.JJ.candidateGenes.recovered; done

cut -f 13 QTL.JJ.candidateGenes.recovered | sort | uniq | wc -l
# 66 genes are included (out of 402)

### how many are under fiber related QTLs
cat QTL.fiberColor.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.fiberColor.genes.redundant; done
cat QTL.architecture.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.architecture.genes.redundant; done
cat QTL.fiberLength.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.fiberLength.genes.redundant; done
cat QTL.fiberQuality.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.fiberQuality.genes.redundant; done
cat QTL.flower.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.flower.genes.redundant; done
cat QTL.fruitingHabit.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.fruitingHabit.genes.redundant; done
cat QTL.seed.names | while read line; do grep $line QTL.JJ.candidateGenes.recovered | cut -f3 -d '=' | sort | uniq >> QTL.JJ.seed.genes.redundant; done

for a in QTL.JJ.*.redundant; do sort $a | uniq > ${a%.redundant}; done

for a in QTL.JJ*.genes; do wc -l $a; done
# 12 QTL.JJ.architecture.genes
# 11 QTL.JJ.fiberColor.genes
# 34 QTL.JJ.fiberLength.genes
# 13 QTL.JJ.fiberQuality.genes
# 10 QTL.JJ.flower.genes
# 10 QTL.JJ.fruitingHabit.genes
# 9 QTL.JJ.seed.genes


for a in QTL.JJ.*.genes; do comm -12 $a QTL.genes.DGE > $a.DGE; wc -l $a.DGE; comm -12 $a TF.gene.names > $a.TF; wc -l $a.TF; comm -12 $a QTL.genes.with.eff.variants > $a.eff; wc -l $a.eff; comm -12 $a QTL.genes.with.aa.variants > $a.aa; wc -l $a.aa; done

# 2 QTL.JJ.architecture.genes.DGE
# 0 QTL.JJ.architecture.genes.TF
# 0 QTL.JJ.architecture.genes.eff
# 0 QTL.JJ.architecture.genes.aa

# 1 QTL.JJ.fiberColor.genes.DGE
# 0 QTL.JJ.fiberColor.genes.TF
# 3 QTL.JJ.fiberColor.genes.eff
# 1 QTL.JJ.fiberColor.genes.aa

# 2 QTL.JJ.fiberLength.genes.DGE
# 0 QTL.JJ.fiberLength.genes.TF
# 11 QTL.JJ.fiberLength.genes.eff
# 6 QTL.JJ.fiberLength.genes.aa

# 2 QTL.JJ.fiberQuality.genes.DGE
# 0 QTL.JJ.fiberQuality.genes.TF
# 5 QTL.JJ.fiberQuality.genes.eff
# 3 QTL.JJ.fiberQuality.genes.aa

# 2 QTL.JJ.flower.genes.DGE
# 0 QTL.JJ.flower.genes.TF
# 4 QTL.JJ.flower.genes.eff
# 1 QTL.JJ.flower.genes.aa

# 0 QTL.JJ.fruitingHabit.genes.DGE
# 0 QTL.JJ.fruitingHabit.genes.TF
# 2 QTL.JJ.fruitingHabit.genes.eff
# 0 QTL.JJ.fruitingHabit.genes.aa

# 1 QTL.JJ.seed.genes.DGE
# 0 QTL.JJ.seed.genes.TF
# 1 QTL.JJ.seed.genes.eff
# 1 QTL.JJ.seed.genes.aa


# extract from Arabidopsis homologs

grep -f ../QTL.fiberColor.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.fiberColor.Arabidopsis.genes
grep -f ../QTL.architecture.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.architecture.Arabidopsis.genes
grep -f ../QTL.fiberLength.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.fiberLength.Arabidopsis.genes
grep -f ../QTL.flower.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.flower.Arabidopsis.genes
grep -f ../QTL.JJ.architecture.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.architecture.Arabidopsis.genes
grep -f ../QTL.JJ.fiberLength.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.fiberLength.Arabidopsis.genes
grep -f ../QTL.JJ.flower.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.flower.Arabidopsis.genes
grep -f ../QTL.JJ.seed.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.seed.Arabidopsis.genes
grep -f ../QTL.fiberColor.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.fiberColor.Arabidopsis.genes
grep -f ../QTL.fiberQuality.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.fiberQuality.Arabidopsis.genes
grep -f ../QTL.fruitingHabit.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.fruitingHabit.Arabidopsis.genes
grep -f ../QTL.JJ.fiberColor.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.fiberColor.Arabidopsis.genes
grep -f ../QTL.JJ.fiberQuality.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.fiberQuality.Arabidopsis.genes
grep -f ../QTL.JJ.fruitingHabit.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.JJ.fruitingHabit.Arabidopsis.genes
grep -f ../QTL.seed.genes G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.seed.Arabidopsis.genes
grep -f ../QTL.genes.DGE G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.DGE.Arabidopsis.genes
grep -f ../genes.with.aa.variants G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.SNPaa.Arabidopsis.genes
grep -f ../genes.with.eff.variants G.hirsutum_Tx-JGI_v1.1_vs_TAIR10.txt > QTL.SNPeff.Arabidopsis.genes








cat QTL.indiv.5pct.expand.bed | while read line; do echo $line | tr ' ' '\t' | bedtools intersect -a stdin -b Tx-JGI_G.hirsutum_v1.1.geneONLY.gff3 -wo >> QTL.candidateGenes.5pct.expand.named;  done
cut -f4 QTL.candidateGenes.5pct.expand.named | sort | uniq | while read line; do grep $line QTL.candidateGenes.5pct.expand.named > $line.5pct.genes; done





cut -f5 QTL.candidateGenes.5pct.expand.named | sort | uniq | while read line; do grep $line QTL.candidateGenes.5pct.expand.named > $line.5pct.genes; done



grep -f ../QTL.genes.DGE Fiber_quality.5pct.genes > QTL.fiberQuality.5pct.DGE.genes
grep -f ../QTL.genes.DGE Fiber_length.5pct.genes > QTL.fiberLength.5pct.DGE.genes
grep -f ../QTL.genes.DGE Fiber_color.5pct.genes > QTL.fiberColor.5pct.DGE.genes

grep -f ../QTL.JJ.fiberColor.genes QTL.fiberColor.5pct.DGE.genes > QTL.JJ.fiberColor.5pct.DGE.genes
grep -f ../QTL.JJ.fiberQuality.genes QTL.fiberQuality.5pct.DGE.genes > QTL.JJ.fiberQuality.5pct.DGE.genes
grep -f ../QTL.JJ.fiberLength.genes QTL.fiberLength.5pct.DGE.genes > QTL.JJ.fiberLength.5pct.DGE.genes


grep -f ../genes.with.aa.variants Fiber_quality.5pct.genes > QTL.SNPaa.FiberQuality.5pct.genes
grep -f ../genes.with.eff.variants Fiber_quality.5pct.genes > QTL.SNPeff.FiberQuality.5pct.genes
grep -f ../genes.with.aa.variants Fiber_length.5pct.genes > QTL.SNPaa.FiberLength.5pct.genes
grep -f ../genes.with.eff.variants Fiber_length.5pct.genes > QTL.SNPeff.FiberLength.5pct.genes
grep -f ../genes.with.aa.variants Fiber_color.5pct.genes > QTL.SNPaa.FiberColor.5pct.genes
grep -f ../genes.with.eff.variants Fiber_color.5pct.genes > QTL.SNPeff.FiberColor.5pct.genes













