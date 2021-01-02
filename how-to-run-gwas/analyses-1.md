# 3 GWAS analyses

There are several softwares for running GWAS.

For independent \(unrelated\) individuals:

* **PLINK** \(Purcell et al. 2007\)--module available, vcf format, limited details in output. [https://www.cog-genomics.org/plink/2.0/](https://www.cog-genomics.org/plink/2.0/)
* **SNPtest** \(Howie at al. 2012\)--stand-alone executable, vcf format, much information in output. [https://mathgen.stats.ox.ac.uk/genetics\_software/snptest/snptest.html\#introductio](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html#introduction)[ ](https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html)

For related individuals:

* **Bolt-LMM** --module available, mixed model analysis, Genetic relationship matrix \(GRM\), vcf or IMPUTE2 \(gen\) format, only continous traits, sample size should be larger than 5000. [https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/BOLT-LMM\_manual.html](https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/BOLT-LMM_manual.html)
* **SAIGE\(gds\)** --Uses R, GDS format, continuous and binary traits. [https://github.com/AbbVie-ComputationalGenomics/SAIGEgds](https://github.com/AbbVie-ComputationalGenomics/SAIGEgds)



