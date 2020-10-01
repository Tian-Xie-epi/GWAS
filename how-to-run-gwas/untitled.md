# Quality control

```text
#!/bin/bash
#SBATCH --time=4:30:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=15G

module load PLINK/1.90b3.44
module load R/3.6.1-foss-2018a

#### Name: Tian Xie
#### Date: 20200805
#### Project: SFS quality control and imputation

### The raw data from Behrooz Alizadeh are "SFS_GWS_670.bed" "SFS_GWS_670.fam" "SFS_GWS_670.bim"
### 318237 markers (build 36) (chr 1-23, 25) pseudoautosomal XY (chr25 in plink), X (chr23)

##################################### Part 1 Quality control
## step 1 exclude samples and markers with call rate<95%
## Step 2 exclude markers with MAF<1% or HWE P<10e-6
## step 3 lift markers from build 36 to build 37
## step 4 exclude samples with Het>4SD or <int-4SD(res)+m*ROH
## step 5 exclude gender mismatches and duplicated samples (check twins/relatedness)
## step 6 population stratification (1000G)
## step 7 HRC check before imputation (duplicated markers, strand, AF concordance, Palindromic SNPs with Freq > 0.4)
## step 8 HRC imputation

################### step 1 exclude samples and markers with call rate<95% (first 80% and then 95%) (45 samples and 3054 markers were excluded)
plink --bfile SFS_GWS_670 --missing --out SFS_callrate
awk '{if ($6>=0.2) print $1,$2}' SFS_callrate.imiss > excl80_callrate_sample.txt
awk '{if ($5>=0.2) print $2}' SFS_callrate.lmiss > excl80_callrate_snp.txt
plink --bfile SFS_GWS_670 --remove excl80_callrate_sample.txt --exclude excl80_callrate_snp.txt --make-bed --out SFS_step1a

plink --bfile SFS_step1a --missing --out SFS_step1a_callrate
awk '{if ($6>0.05) print $1,$2}' SFS_step1a_callrate.imiss > excl95_callrate_sample.txt
plink --bfile SFS_step1a --remove excl95_callrate_sample.txt --make-bed --out SFS_step1b

plink --bfile SFS_step1b --missing --out SFS_step1b_callrate
awk '{if ($5>0.05) print $2}' SFS_step1b_callrate.lmiss > excl95_callrate_snp.txt
plink --bfile SFS_step1b --exclude excl95_callrate_snp.txt --make-bed --out SFS_step1c


##################### step 2 -MAF and HWE filtering-

### calculate MAF and HWE
plink --bfile SFS_step1c --freq --out SFS_step1c_maf

plink --bfile SFS_step1c --chr 1-22,25 --make-bed --out SFS_step1c_autosome
plink --bfile SFS_step1c --chr 23 --make-bed --out SFS_step1c_chrX

plink --bfile SFS_step1c_autosome --hardy --out SFS_step1c_autosome_hwe
plink --bfile SFS_step1c_chrX --hardy --filter-females --out SFS_step1c_chrX_hwe

### exclude snps with MAF<1% and hwe<1*10e-6

awk '$5<0.01 {print $2}' SFS_step1c_maf.frq > excl_maf001
awk '$9<0.000001 {print $2}' SFS_step1c_autosome_hwe.hwe > excl_hwe6_autosome
awk '$9<0.000001 {print $2}' SFS_step1c_chrX_hwe.hwe > excl_hwe6_chrX

cat excl_maf001 excl_hwe6_autosome excl_hwe6_chrX  > excl_maf_hwe

plink --bfile SFS_step1c --exclude excl_maf_hwe --make-bed --out SFS_step2


##################### step 3 -lift markers from build 36 to build 37

awk '{print "chr"$1,$4-1,$4,$2}' SFS_step2.bim > SFS_step2_build36.txt
awk '{gsub("chr23","chrX",$1)}1' SFS_step2_build36.txt | awk '{gsub("chr25","chrX",$1)}1' > SFS_step2_build36_forlift.txt

## use Assembly Converter http://www.ensembl.org/Homo_sapiens/Tools/AssemblyConverter?db=core to lift markers from build 36 to build 37.
## 314767 markers were lifted (32 markers cannot be lifted mainly because they are in HG7_PATCH and HG183_PATCH in build37)
## lifted file is "SFS_step2_build37.bed"

awk '{print $4}' SFS_step2_build37.bed > SFS_step2_liftSNPs.txt
awk 'BEGIN {FS=OFS="\t"} {print $1,$3,$4}' SFS_step2_build37.bed > SFS_step2_build37_chrpos.txt

plink --bfile SFS_step2 --extract SFS_step2_liftSNPs.txt --make-bed --out SFS_step3a

paste SFS_step3a.bim SFS_step2_build37_chrpos.txt | awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$8,$5,$6}' > SFS_step3b_build37.bim
cp SFS_step3a.fam SFS_step3b_build37.fam
cp SFS_step3a.bed SFS_step3b_build37.bed


##################### step 4 -Samples Heterozygosity-

### first, select independent SNPs (pruning --indep [SNPwindow] [shift] [LD threshold in 1/(1-r2)])
### separate the HLA SNPs
### then extract these 
### LD-based pruning requires a sorted .bim.  

## sort bim data (sort_bim.R)
## SFS_bim<-read.table("SFS_step3b_build37.bim")
## SFS_bim_sorted<-SFS_bim[order(SFS_bim$V1,SFS_bim$V4),]
## write.table(SFS_bim_sorted,"SFS_sorted_build37.bim",row.names = F,quote = F,col.names = F)
Rscript sort_bim.R
awk '{print $2}' SFS_sorted_build37.bim > SFS_sorted_SNPid.txt
plink --bfile SFS_step3b_build37 --extract SFS_sorted_SNPid.txt --make-bed --out SFS_step4a

plink --bfile SFS_step4a --chr 1-22 --indep 50 5 2.5 --out prune_indep
awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' SFS_step4a.bim > exclHLA.txt
plink --bfile SFS_step4a --extract prune_indep.prune.in --exclude exclHLA.txt --make-bed --out SFS_indep

### perform heterozygosity in these independent SNPs without HLA region

plink --bfile SFS_indep --het --homozyg --out SFS_het
plink --bfile SFS_step4a --missing --out SFS_step4a_miss

### use R software (script "SFS_sampleQC_het.R")to create file with samples to exclude (het>4sd)
### "excl_sample_het.txt" is the list of individuals which should be excluded because of Samples Heterozygosity 
### exclude these (13 samples)

Rscript SFS_sampleQC_het.R
plink --bfile SFS_step4a --remove excl_sample_het.txt --make-bed --out SFS_step4b


#################### - step 5 exclude gender mismatches and duplicated samples (check twins/relatedness)
## identify duplicated sample and siblings 

plink --bfile SFS_indep --genome --out SFS_indep_genome
awk '$10>0.2 {print $1,$2,$3,$4,$10}' SFS_indep_genome.genome> SFS_relatedness.txt

## six pairs of siblings or DZ twins, keep them in the dataset
#FID1 IID1 FID2 IID2 PI_HAT
#1 1394 1 1398 0.5052
#1 2186 1 2190 0.5000
#1 496 1 501 0.5035
#1 596 1 600 0.4996
#1 306 1 310 0.5000
#1 1248 1 1275 0.5000

### sex check 
### exclude 14 samples with mismatched sex

plink --bfile SFS_step4b --check-sex --out checksex
awk 'BEGIN { FS=" "; OFS=" "} /PROBLEM/ {print $1,$2,$3,$4,$5,$6}' checksex.sexcheck >incorrectsex.txt
awk '{print $1,$2}' incorrectsex.txt > excl_incorrectsex.txt

plink --bfile SFS_step4b --remove excl_incorrectsex.txt --make-bed --out SFS_step5


#################### - step 6 population stratification (1000G)

### extract SFS SNPs from 1000G and compare concordance (only for autosomes)
awk '{print $1,$4,$4,$2}' SFS_step5.bim > SFS_step5SNPs_chrpos.txt
awk '{print $2}' SFS_step5.bim > SFS_step5SNPs_ID.txt
plink --bfile SFS_step5 --chr 1-22 --make-bed --out SFS_step5_autosome

###  SNPs are extracted from 1000G (chr1-22) (sbatch extract1000G.sh)
plink --bfile /data/p282717/1000G/1000G_all --extract SFS_step5SNPs_ID.txt --make-bed --out 1000G_SFSallSNPs_ID
plink --bfile /data/p282717/1000G/1000G_all --extract SFS_step5SNPs_chrpos.txt --range --make-bed --out 1000G_SFSallSNPs_chrpos

### use 1000G_SFSallSNPs_chrpos (extract 305051 SNPs)
### assign new SNP ID (Chr:Pos) to SNP and 1000G SNPs because some SNPs in SFS don't have standard rsID
awk '{print $1,$1":"$4,$3,$4,$5,$6}' SFS_step5_autosome.bim > SFS_step5_autosome_newID.bim
cp SFS_step5_autosome.fam SFS_step5_autosome_newID.fam
cp SFS_step5_autosome.bed SFS_step5_autosome_newID.bed

awk '{print $1,$1":"$4,$3,$4,$5,$6}' 1000G_SFSallSNPs_chrpos.bim > 1000G_SFSallSNPs_newID.bim
cp 1000G_SFSallSNPs_chrpos.fam 1000G_SFSallSNPs_newID.fam
cp 1000G_SFSallSNPs_chrpos.bed 1000G_SFSallSNPs_newID.bed

plink --bfile 1000G_SFSallSNPs_newID --bmerge SFS_step5_autosome_newID.bed SFS_step5_autosome_newID.bim SFS_step5_autosome_newID.fam --merge-mode 6 --out SFS_1000G_mismatch 

plink --bfile SFS_step5_autosome_newID --flip SFS_1000G_mismatch.missnp --make-bed --out SFS_step5_autosome_flip 

plink --bfile 1000G_SFSallSNPs_newID --bmerge SFS_step5_autosome_flip.bed SFS_step5_autosome_flip.bim SFS_step5_autosome_flip.fam --merge-mode 6 --out SFS_1000G_mismatch2

plink --bfile SFS_step5_autosome_flip --exclude SFS_1000G_mismatch2.missnp --make-bed --out SFS_1000G_merge
plink --bfile 1000G_SFSallSNPs_newID --exclude SFS_1000G_mismatch2.missnp --make-bed --out 1000G_SFSallSNPs_merge

plink --bfile 1000G_SFSallSNPs_merge --bmerge SFS_1000G_merge.bed SFS_1000G_merge.bim SFS_1000G_merge.fam --out SFS_1000G 

### generate SFS and 1000G independent markers
awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' SFS_1000G.bim > exclHLA_PCA.txt
plink --bfile SFS_1000G --maf 0.1 --geno 0.05 --exclude exclHLA_PCA.txt --make-bed --out SFS_1000G1
plink --bfile SFS_1000G1 --indep-pairwise 1000 5 0.2 --out SFS_1000G_indep
plink --bfile SFS_1000G1 --extract SFS_1000G_indep.prune.in --make-bed --out SFS_1000G_common

### generate SFS and 1000G individual file
awk '{print $1,$2,"SFS"}' SFS_1000G_merge.fam > SFSsample
awk '{print $1,$2,"1000G"}' 1000G_SFSallSNPs_merge.fam >1000Gsample
cat 1000Gsample SFSsample > SFS_1000Gsample

###-make PCA on 1000G and project SFS on them
plink --bfile SFS_1000G_common --within SFS_1000Gsample --pca 10 --pca-cluster-names 1000G --out SFS1000G_PCAresult

### make plots of PCA and identify SFS non-European individuals 
### (output: "Gecko check ethnicity.csv","Gecko_non-European_sample.txt","Gecko and 1000G PCA1.tiff","Gecko and 1000G PCA2-European.tiff",) 
cp /data/p282717/1000G/phase3_orig.psam phase3_orig.psam
Rscript SFS_PCAplot.R

## all samples are European ancestry

#############---calculate Gecko PCAs for further analyses (GWAS)---############################

awk '{if ($1 == 6 && $4 >= 28477797 && $4 <= 35000000) print $2}' SFS_step5.bim > exclHLA_SFS_PCA.txt
plink --bfile SFS_step5 --maf 0.1 --geno 0.05 --exclude exclHLA_SFS_PCA.txt --make-bed --out SFSPCA1
plink --bfile SFSPCA1 --indep-pairwise 1000 5 0.2 --out SFSPCA2
plink --bfile SFSPCA1 --extract SFSPCA2.prune.in --make-bed --out SFS_PCA_common

plink --bfile SFS_PCA_common --pca --out SFS_PCAresult
 

###############################################################################################################
#############---export cleaned data and report---############################

mkdir SFS_cleaned_data #for SFS cleaned data
mkdir SFS_QCreport     #for SFS QC report (figures, text, log)
cp SFS_step5.fam SFS_step5.bed SFS_step5.bim SFS_cleaned_data
cp "SFS and 1000G PCA1.tiff" "SFS and 1000G PCA2-European.tiff" "SFS_sample heterozygosity_call rate.tiff" "SFS sample heterozygosity.tiff" SFS_PCAresult.eigenvec SFS_QCreport
cp -R *.log SFS_QCreport


############## step 7 HRC check before imputation (duplicated markers, strand, AF concordance, Palindromic SNPs with Freq > 0.4)
### recode chromosome 25 to 23 
cd SFS_cleaned_data

awk '{gsub("25","23",$1)}1' SFS_step5.bim > SFS_step7.bim
cp SFS_step5.bed SFS_step7.bed
cp SFS_step5.fam SFS_step7.fam

plink --bfile SFS_step7 --freq --out SFS_clean_frq

unzip HRC-1000G-check-bim-v4.2.13-NoReadKey.zip
gunzip HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz

perl HRC-1000G-check-bim-NoReadKey.pl -b SFS_step7.bim -f SFS_clean_frq.frq -r HRC.r1-1.GRCh37.wgs.mac5.sites.tab -h
sh Run-plink.sh



```



