This directory contains all of the scripts for the manuscript 'Identification of 5mC within heterogenous tissue using de-novo somatic mutations'


Scripts are divided into subdirectories based on the analyses for which they were used, as below:

**Annotations** contains all scripts used to obtain the reference genome and annotate the features that it contains

**General** contains all of the general scripts for obtaining, preprocessing, and processing reads. This includes all scripts for finding the proportion blood contamination 'P', methylation calling, and deconvolution.

**GermlineSimulations** contains all of the scripts used to run simulations of germline heterozygote sites

**SomaticSimulations** contains all of the scripts used to simulate somatic mutations

The README.md files within each subdirectory contains the details of the purpose for which each script was used.

All bash files were run on LiDO3 computing cluster with the general specifications detailed in: LiDO3_first_contact_handout.pdf

R scripts were performed in R (R Core Team, 2023) using the RStudio Build 494 (Posit, Boston MA, USA) console on Windows 10, running R version 4.3.2 "Eye Holes" (2023-10-31 ucrt).
