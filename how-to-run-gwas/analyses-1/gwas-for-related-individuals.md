# SAIGEgds for related individuals

## Background

On this page you will know how to run GWAS in related individuals using the R package `SAIGEgds` .  

`SAIGEgds`  leverages the Genomic Data Structure \(GDS\) format to provide efficient mixed-model analyses which allow for adjusting family relatedness. See more details in below papers and tutorials.

**Paper** 

**SAIGEgds** \(Zheng et al, 2020\) [https://doi.org/10.1093/bioinformatics/btaa731](https://doi.org/10.1093/bioinformatics/btaa731)

**SAIGE** \(Zhou et al, 2018\) [https://www.nature.com/articles/s41588-018-0184-y](https://www.nature.com/articles/s41588-018-0184-y) 

**Tutorial** [http://www.bioconductor.org/packages/devel/bioc/vignettes/SAIGEgds/inst/doc/SAIGEgds.html](http://www.bioconductor.org/packages/devel/bioc/vignettes/SAIGEgds/inst/doc/SAIGEgds.html)

[https://bioconductor.riken.jp/packages/3.10/bioc/manuals/SAIGEgds/man/SAIGEgds.pdf](https://bioconductor.riken.jp/packages/3.10/bioc/manuals/SAIGEgds/man/SAIGEgds.pdf)

## Required Data

<table>
  <thead>
    <tr>
      <th style="text-align:left">File name</th>
      <th style="text-align:left">File type</th>
      <th style="text-align:left">Description</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td style="text-align:left">BP_phenotypes.csv</td>
      <td style="text-align:left">phenotype data</td>
      <td style="text-align:left">can be in other formats (e.g. spss, txt)</td>
    </tr>
    <tr>
      <td style="text-align:left">PCAresult.eigenvec</td>
      <td style="text-align:left">PCA data</td>
      <td style="text-align:left">This file contains the first 20 PCs of samples</td>
    </tr>
    <tr>
      <td style="text-align:left">chr_all_genotype.bed, chr_all_genotype.bim, chr_all_genotype.fam</td>
      <td
      style="text-align:left">genotyped data (Plink file)</td>
        <td style="text-align:left">hard called genotypes (e.g. AG, CT), genotyped by specific array/chip
          (e.g. Illumina GSA array)</td>
    </tr>
    <tr>
      <td style="text-align:left">
        <p>chr1.dose.vcf.gz ...</p>
        <p>chr22.dose.vcf.gz</p>
      </td>
      <td style="text-align:left">imputed data (VCF file)</td>
      <td style="text-align:left">1000G or HRC imputation data, include dosage of markers</td>
    </tr>
    <tr>
      <td style="text-align:left">chr1.info.gz ... chr22.info.gz</td>
      <td style="text-align:left">info data along with imputed data</td>
      <td style="text-align:left">This file contains imputation quality score of markers</td>
    </tr>
  </tbody>
</table>

## Preparation for analyses

### Individual exclusions and stratified analyses

Always check guidelines from analyses plan before running GWAS. We may need to do some individual exclusions and stratified analyses. 

For example, in our meta-GWAS of child BP project, we should exclude**:**  

1. individuals&gt;=18 years of age. 
2. reported cases of hypertension \(children with a diagnose and thereby treatment\). 

In addition, stratified analyses are required:

1. Ethnicity: If we have data from different ethnicity, then perform all analyses stratified by ethnicity \(Caucasians/Europeans, African Americans, etc.\). If the data is mainly from one ethnicity, the few individuals from other ethnicities can be excluded. 
2. Gender based: Perform analyses stratifying for sex, separately in males and females as well as a combined analyses, pooling all. 

### Prepare phenotype

Before running GWAS, we always need to prepare phenotype \(e.g. trait transformation\). According to the analyses plan of our meta-GWAS project, we should do the following trait transformation:

1. Calculate residuals by running this regression model for each trait separately in males and females: Trait \(SBP/DBP/MAP/PP/HR\) ~ age + height + PCs \(if available\) + \(study specific covariates if required\)
2. Perform rank-based Inverse Normal transformation of these residuals

{% hint style="info" %}
Study specific covariate can be used in below two situations. 

1. If you have genetic data genotyped by different arrays/chips but imputed with the same imputation reference panel \(1000G or HRC\), you can combine them into one analysis and add a covariate indicating different chips \(e.g. 1-chip A, 2-chip B\).
2. If you have genetic data \(same chip and same imputation panel\) from different cohorts, you can also combine them into one analysis and add a covariate indicating different cohorts \(e.g. 1-cohort A, 2-cohort B\).  
{% endhint %}

Here I present a example of R code for preparing phenotype. The R code follows steps: 

1. import **phenotype data** and **PCA data** 
2. merge phenotype and PCA, by this we can also select individuals with both phenotype and genotype
3. here we can describe phenotype \(mean, sd\) for supplementary table
4. trait transformation \(rank-based Inverse Normal transformation residuals\) in males and females separately
5. export data for GWAS \(**"pheno\_pooled\_invBP.txt", "pheno\_males\_invBP.txt", "pheno\_females\_invBP.txt"**\)

Three R packages are used to increase the efficiency:

* package `data.table`, function `fread` and `fwrite` can import and export data faster and more convenient. [https://www.rdocumentation.org/packages/data.table/versions/1.13.0/topics/fread](https://www.rdocumentation.org/packages/data.table/versions/1.13.0/topics/fread) 
* package `dplyr`, which is useful in data manipulation \(e.g. function `mutate` for creating new variables\). [https://dplyr.tidyverse.org/reference/](https://dplyr.tidyverse.org/reference/)
* package `tableone`,  which eases the construction of “Table 1”, _i.e._, individuals characteristics table commonly found in supplementary tables of meta-GWAS paper \(used in step 3\). [https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html](https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html)

{% code title="Rscript prepare\_pheno.R" %}
```r
library(data.table) 
library(dplyr)
library(tableone)

##### Step 1 data input
pheno<-fread("BP_phenotypes.csv")  #ID,Age,Sex(1-males,2-females),Height,Weight,SBP,DBP,HR
colnames(pheno)<-c("IID","Age","Sex","Height","Weight","SBP","DBP","HR")  
pheno<-pheno%>% mutate(PP=SBP-DBP,MAP=1/3*SBP+2/3*DBP)

PCA<-fread("PCAresult.eigenvec")
colnames(PCA)<-c("FID","IID",paste0("PC",1:20))
PCA<-PCA[,1:12]

##### Step 2 merge phenotype and PCA
pheno_forGWAS<-merge(pheno,PCA,by="IID")

##### Step 3 describe phenotype (mean, sd) for supplementary table
tableone_pheno<-cbind(print(CreateTableOne(data = pheno_forGWAS[,c("Age","Height","Sex","SBP","DBP","PP","MAP","HR")])),print(CreateTableOne(strata = "Sex", data = pheno_forGWAS[,c("Age","Height","Sex","SBP","DBP","PP","MAP","HR")])))
write.csv(tableone_pheno,"phenotype_tableone.csv")

##### Step 4 Trait transformation (Inverse Normal transformation of residuals)
pheno_males<-filter(pheno_forGWAS,Sex==1)
pheno_females<-filter(pheno_forGWAS,Sex==2)

compute_INR<-function(pheno,dataset){
resid_linear=residuals(lm(pheno~Age+Height+PC1+PC2+PC3+PC4+PC5+PC6+PC7+PC8+PC9+PC10,data=dataset))
inv_pheno=qnorm((rank(resid_linear,na.last="keep")-0.5)/sum(!is.na(resid_linear)))
return(inv_pheno)}

pheno_males<- pheno_males %>% mutate(inv_SBP=compute_INR(SBP,pheno_males),inv_DBP=compute_INR(DBP,pheno_males),inv_MAP=compute_INR(MAP,pheno_males),inv_PP=compute_INR(PP,pheno_males),inv_HR=compute_INR(HR,pheno_males))
pheno_females<- pheno_females %>% mutate(inv_SBP=compute_INR(SBP,pheno_females),inv_DBP=compute_INR(DBP,pheno_females),inv_MAP=compute_INR(MAP,pheno_females),inv_PP=compute_INR(PP,pheno_females),inv_HR=compute_INR(HR,pheno_females))

pheno_males_invBP<-pheno_males %>% select(FID,IID,inv_SBP:inv_HR)
pheno_females_invBP<-pheno_females %>% select(FID,IID,inv_SBP:inv_HR)
pheno_pooled_invBP<-rbind(pheno_males_invBP,pheno_females_invBP)

##### Step 5 export data for GWAS  
fwrite(pheno_pooled_invBP,file="pheno_pooled_invBP.txt",sep="\t")
fwrite(pheno_males_invBP,file="pheno_males_invBP.txt",sep="\t")
fwrite(pheno_females_invBP,file="pheno_females_invBP.txt",sep="\t")
```
{% endcode %}

## Run GWAS using SAIGEgds 

`SAIGEgds` is used to run GWAS using mixed model which allow for adjusting family relatedness. All analyses are conducted in the computer cluster of  University of Groningen \([https://wiki.hpc.rug.nl/peregrine/start](https://wiki.hpc.rug.nl/peregrine/start)\) based on Linux operating system \(A beginners guide [http://www.ee.surrey.ac.uk/Teaching/Unix/](http://www.ee.surrey.ac.uk/Teaching/Unix/)\). 

### Step 1 Install SAIGEgds and relevant package

### Step 2 convert plink and vcf file to gds file

### Step 3 Preparing SNP data for genetic relationship matrix

### Step 4 Fitting the null model

### Step 5 Fitting the null model



