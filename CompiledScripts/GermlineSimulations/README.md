This subdirectoy contains the scripts used to simulate germline-heterozyous sites, align them to the reference genome and subsequent analyses.


The contained files were used for the purposes listed below:

**BismarkAlign100X-ShermanSim-GreatTitAbel.bash** was used to align simulated bisulfite reads to the reference genome with Bismark

**GermlineHalfBinom.R** was used to simulate a binomial distirbution of germline allele frequecies at half the depth of the SNVs aliged simulated reads

**GermlineSomaticHalfBinom.R** was used to simulate a binomial distirbution of somatic and germline allele frequecies at half the depth of the SNVs aliged simulated reads

**GermlineSummarizeFilteredSNVs.R** was used to assess distributions of simulated allele frequency calls at sites SNV filtered for repetitive content

**GermlineSummarizeSNVs.R** was used to perform basic summary analyses of simulated germline heterozygous sites

**afFiltSiteSpectra-BismarkAlign100X-ShermanSim.bash** was used to obtain site/allele frequency information for simulated germline-heterozygous SNVs that were not associated with repetitive content

**afSiteSpectra-BismarkAlign100X-ShermanSim.bash** was used to obtain site/allele frequency information for simulated germline-heterozygous SNVs

**bisulphiteShermanSimHaplotype1-GreatTitAbel.bash** was used to simulate bisulfite reads from the first haplotype of the simulated individual genome

**bisulphiteShermanSimHaplotype2-GreatTitAbel.bash** was used to simulate bisulfite reads from the second haplotype of the simulated individual genome

**bsSnper-BismarkAlign100X-ShermanSim.bash** was used to call SNVs from Sherman simulated reads using BS-SNper

**gzipMERGEDread-ShermanSim-GreatTitAbel.bash** was used to gzip the merged simulated bisulfite reads

**gzipRawRead-ShermanSim-GreatTitAbel.bash** was used to gzip the simulated bisulfite reads

**mergeReads-ShermanSim-GreatTitAbel.bash** was used to merge the simulated bisulfite reads

**mutSimHeterozygote-GreatTitAbel.bash** was used to introduce 2 haplotypes of germline-heterygous mutations onto reference genome

**sortDedupBismarkAlign100X-ShermanSim.bash** was used to sort the deduplicated bisulfite reads aligned to the reference genome. This is a necessary step for variant calling with BS-SNper.

**subsample100X-ShermanSim-GreatTitAbel.bash** was used to subsample simulated reads to 100X for alignment

**trim100X-ShermanSim-GreatTitAbel.bash.bash** was used to trim simulated subsampled reads prior to alignment using Trimmomatic

