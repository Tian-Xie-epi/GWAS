# SAIGEgds for related individuals

## Background

On this page you will know how to run GWAS in related individuals using the R package `SAIGEgds` .  

`SAIGEgds`  leverages the Genomic Data Structure \(GDS\) format to provide efficient mixed-model analyses which allow for adjusting family relatedness. See more details in paper \(Zheng et al, 2020\) [https://doi.org/10.1093/bioinformatics/btaa731](https://doi.org/10.1093/bioinformatics/btaa731) and tutorial [https://github.com/AbbVie-ComputationalGenomics/SAIGEgds](https://github.com/AbbVie-ComputationalGenomics/SAIGEgds).

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



