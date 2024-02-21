###### STEP 1: download data

# sra-tools version: 3.0.7

## download genome and annotations
mkdir data/genome cd data/genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Scyliorhinus_canicula/annotation_releases/current/100/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Scyliorhinus_canicula/annotation_releases/current/100/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_genomic.gtf.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Scyliorhinus_canicula/annotation_releases/current/100/GCF_902713615.1_sScyCan1.1/md5checksums.txt
wgat https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_other/Scyliorhinus_canicula/annotation_releases/current/100/GCF_902713615.1_sScyCan1.1/GCF_902713615.1_sScyCan1.1_assembly_report.txt
md5sum --check md5checksums.txt
gunzip *.gz
samtools faidx data/genome/GCF_902713615.1_sScyCan1.1_genomic.fna
cd ../..
## download rnaseq data
mkdir data/rnaseq ; cd data/rnaseq
for i in ERR6338447 ERR6338403 ERR6338409 ERR6338415 ERR6338421 ERR6338427 ERR6338422 ERR6338430 ERR6338431 SRR5179116 SRR5179117 SRR8179289 SRR8179290 SRR8179291 SRR1514131 SRR1514130 SRR1514129 ; do
	# prefetch
	prefetch $i
	# dump fastq
	fasterq-dump $i --split-files
	# gzip
	gzip ${i}_*.fastq
	# remove download dir
	rm -rf $i
done

###### STEP 2: use the nf-core rnavar pipeline for mapping and SNP calling of each sample

# nextflow version: 23.10.0
# rnavar version: 1.1.0dev
# rnavar must be installed in nfcore/rnavar/dev/

mkdir results; mkdir results/var

# run rnavar separately for samples of each read length because the star index is optimised for the correct read length
for(L in 51 100 101 150 151) ; do
	mkdir results/var/L${L}
	nextflow run nfcore/rnavar/dev \
	-profile scw \
	--input data/sample_sheet_L${L}.csv \
	--outdir results/var/L${L} \
	--fasta data/genome/GCF_902713615.1_sScyCan1.1_genomic.fna \
	--gtf data/genome/GCF_902713615.1_sScyCan1.1_genomic.gtf \
	--remove_duplicates true \
	--skip_baserecalibration true \
	--skip_variantannotation true \
	--star_twopass true \
	--generate_gvcf true \
	--read_length $L 
done

# move vcfs into one dir
mkdir results/var/all_gvcf
mv results/var/*/*.haplotypecaller.combined.g.vcf.* results/var/all_gvcf

###### STEP 3: combine gvcfs, run joint genotyping and filter

# GATK version: 4.5.0
# bcftools version: 1.16
# plink version: 2.0
# vcftools version: 0.1.16
# tabix version: 1.15
# bgzip version: 1.15

# combine gvcfs
gatk CombineGVCFs \
--variant results/var/all_gvcf/A_SRR1514129.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/A_SRR1514130.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/A_SRR1514131.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/A_SRR8179289.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/A_SRR8179291.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338403.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338409.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338415.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338421.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338422.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338427.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338430.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338431.haplotypecaller.combined.g.vcf.gz \
--variant results/var/all_gvcf/M_ERR6338447.haplotypecaller.combined.g.vcf.gz \
--reference data/genome/GCF_902713615.1_sScyCan1.1_genomic.fna \
-L NC_052146.1 -L NC_052147.1 -L NC_052148.1 -L NC_052149.1 -L NC_052150.1 -L NC_052151.1 -L NC_052152.1 -L NC_052153.1 -L NC_052154.1 -L NC_052155.1 -L NC_052156.1 -L NC_052157.1 -L NC_052158.1 -L NC_052159.1 -L NC_052160.1 -L NC_052161.1 -L NC_052162.1 -L NC_052163.1 -L NC_052164.1 -L NC_052165.1 -L NC_052166.1 -L NC_052167.1 -L NC_052168.1 -L NC_052169.1 -L NC_052170.1 -L NC_052171.1 -L NC_052172.1 -L NC_052173.1 -L NC_052174.1 -L NC_052175.1 -L NC_052176.1 \
--output results/var/comb.g.vcf.gz &> results/var/CombineGVCFs.log

# joint genotyping to make all-sites vcf for pixy
gatk GenotypeGVCFs \
--variant results/var/comb.g.vcf.gz \
--reference data/genome/GCF_902713615.1_sScyCan1.1_genomic.fna \
-all-sites \
--output results/var/jgen.g.vcf.gz &> results/var/GenotypeGVCFs.log

# filters applied to both invariant and variant sites
vcftools --gzvcf results/var/jgen.g.vcf.gz \
--remove-indels \
--max-missing 0.85 \
--min-meanDP 5 \
--max-meanDP 500 \
--recode --stdout | bgzip -c > results/var/filt1.g.vcf.gz
tabix results/var/filt1.g.vcf.gz

# separate invariant and variant sites and apply additional filters to variants
# inv
vcftools --gzvcf results/var/filt1.g.vcf.gz \
--max-maf 0 \
--recode --stdout | bgzip -c > results/var/filt.inv.vcf.gz
tabix results/var/filt.inv.vcf.gz

# var
# no minor allele frequency filter applied because of the small number of samples
# keep biallelic sites only
vcftools --gzvcf results/var/filt1.g.vcf.gz \
--mac 1 \
--minQ 30 \
--maf 0.1 \
--recode --stdout | bgzip -c > results/var/filt.var.vcf.gz
tabix results/var/filt.var.vcf.gz

# combine vcfs
bcftools concat \
--allow-overlaps \
results/var/filt.inv.vcf.gz results/var/filt.var.vcf.gz \
-O z -o results/var/filt.all.vcf.gz
tabix results/var/filt.all.vcf.gz

# make plink format file for PCadapt
plink2 --vcf results/var/filt.var.vcf.gz \
--max-alleles 2 \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--make-bed --out results/var/filt.var

###### STEP 4: calculate pi, fst & dxy

# pixy version: 1.2.10.beta2
# vcftools version: 0.1.16

# make bed with windows
head -n31 data/genome/GCF_902713615.1_sScyCan1.1_genomic.fna.fai | cut -f1,2 > results/var/genome.tsv
bedtools makewindows -g results/var/genome.tsv -w 1000000 > results/var/1Mb_window.bed 
mkdir results/pixy
# 1 Mb non-overlapping window
pixy --stats pi fst dxy \
--vcf results/var/filt.all.vcf.gz \
--populations data/pops.txt \
--bed_file results/var/1Mb_window.bed \
--n_cores 30 \
--output_folder results/pixy \
--output_prefix 1Mb_window &> results/pixy/pixy_1Mb.log
# per site for variant sites only
mkdir results/vcftools_fst
vcftools --gzvcf results/var/filt.var.vcf.gz \
--weir-fst-pop data/med.txt \
--weir-fst-pop data/atl.txt \
--out results/vcftools_fst/persite_fst

###### STEP 5: run population structure analyses

# admixture version: 1.3.0
# plink version: 2.0

# Thin snps to one every 100 kb to reduce the influence of linkage (there are too few samples for linkage pruning)
plink2 --bfile results/var/filt.var \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--make-bed \
--bp-space 100000 \
--out results/var/thin_100kb

# estimate allele frequency for PCA
mkdir results/pca
plink2 --bfile results/var/thin_100kb \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--freq \
--out results/pca/freq 

# run PCA
plink2 --bfile results/var/thin_100kb \
--double-id --allow-extra-chr --set-missing-var-ids @:# \
--pca \
--read-freq results/pca/freq.afreq  \
--out results/pca/pca

# run ADMIXTURE

mkdir results/admixture
# rename chromosomes to numbers for admixture
cp results/var/thin_100kb* results/admixture/
awk '{$1="0";print $0}' results/admixture/thin_100kb.bim  > results/admixture/thin_100kb.bim.tmp
mv results/admixture/thin_100kb.bim.tmp results/admixture/thin_100kb.bim
cd results/admixture
# run admixture
for i in {2..6} ; do
   admixture --cv thin_100kb.bed $i > log${i}.out
done
# collect CV errors
awk '/CV/ {print $3,$4}' *out | cut -c 4,7-20 > cv.error
cd ../..
# see: pop_structure_plots.R for plotting

###### STEP 6: calculate genome structure statistics
# bedtools version: 2.27.1
# bedops version (gff2bed): 2.4.41
# gc content 
mkdir results/gc
bedtools nuc -fi data/genome/GCF_902713615.1_sScyCan1.1_genomic.fna -bed results/var/1Mb_window.bed > results/gc/GC_1Mb.tsv
# gene density
mkdir results/gene_density
awk '$3=="gene" {print $0}' data/genome/GCF_902713615.1_sScyCan1.1_genomic.gff | gff2bed > results/gene_density/genes.bed
bedtools intersect -a results/var/1Mb_window.bed -b results/gene_density/genes.bed -c > results/gene_density/ngenes_1Mb.tsv 

###### STEP 7: PCadapt
# see: pcadapt.R

# then intersect significantly differentiated SNPs with genes:

# bedtools version: 2.27.1

bedtools intersect -wb -a results/pcadapt/sig_snps.bed -b results/gene_density/genes.bed > results/pcadapt/sig_snp_genes.tsv

###### STEP 8: Manhattan plots
# see: manhattan_plots.R
