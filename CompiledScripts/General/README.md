This subdirectoy contains the general scripts used to obtain, preprocess, and process reads. This includes all scripts for finding the proportion blood contamination 'P', methylation calling, and deconvolution.


The contained files were used for the purposes listed below:

**Abel-bsReadCounter.bash** was used to count the reads in each subsample at each stage of preprocssing

**BS-SNper_Subsamp-Decon.bash** was used to call SNV variants for estimating 'P' from bisulfite reads using BS-SNper

**BismarkBisulphiteGenome.bash** was used to create bisulfite converted genomes from the reference genome using Bismark. This is necessary to align bisulfite genes to the genome with Bismark.

**BismarkMethylation-Subsampler-Decon.bash** was used to make methylation calls on the reference genome with the bisulfite reads

**DepthAnalysis.R** was used to analyze outputs relating to depth of reads following Bismark alignment

**DownloadDeconvolutionData.bash** was used to obtain the bisulfite reads from the reference individual use in this analysis. See--Derks, M.F., Schachtschneider, K.M., Madsen, O., Schijlen, E., Verhoeven, K.J. and van Oers, K., 2016. Gene and transposable element methylation in great tit (Parus major) brain and blood. Bmc Genomics, 17, pp.1-13.

**GenomeFeature-MethylationCaller.R** was used to run Fisher's Exact Tests to calculate differences between methylation levels between the blood and brain for genomic features larger than one site, i.e. CpG islands and gene bodies

**P-BasicsAnalyis.R** was used to calculate P and variation in somatic sites across samples

**ReadSummaryGraph.R** was used to summarize read information at each stage of preprocessing

**SensitivitySpecificity.R** was used to analyze data on the sensitivity and specificity of calls between the blood and the brain of the reference bird before and following deconvolution

**SiteSpecificValidations.R** was used for site-specific validations of deconvolution

**SubsampleTrimmer-Decon.bash** was used to trim subsampled bisulfite reads with trimmomatic

**Subsampler-Decon.bash** was used to subsample bisulfite reads from the blood and brain of the reference individual

**SubsamplerSortBismarkDedup-Decon.bash** was used to sort Bisamrk alignments of subsampled reads to the reference genome--a necessary step for SNV calling with BS-SNper

**WG-MethylationCaller.R** was used to call methylation differences between the blood and brain at individual sites using Fisher's Exact Tests

**alignSubsamplerBismarkDecon.bash** was used to align bisulfite reads to the bisulfite converted reference genome using Bismark

**bsBloodSomatic-x-DerksETAL2016.bash** was used to compare somatic sites found in blood reads with bisulfite sequecning to those found independent non-converted blood reads from the same individual See--Derks, M.F., Schachtschneider, K.M., Madsen, O., Schijlen, E., Verhoeven, K.J. and van Oers, K., 2016. Gene and transposable element methylation in great tit (Parus major) brain and blood. Bmc Genomics, 17, pp.1-13.

**corrB-GeneBody-MethyC.ultFiltMixed.bash** was used to create tables of methylation levels of gene bodies in subsamples with and without deconvolution

**corrB-MethyC.ultFiltMixed.bash** was used to compile methylation levels for all cytosines in subsamples with and without deconvolution

**corrB-Subset-MethyC.ultFiltMixed.bash** was used to subsample cytosines from all subsamples for validation analyses

**corrB-iCpG-MethyC.ultFiltMixed.bash** was used to create tables of methylation levels of CpG islands in subsamples with and without deconvolution

**decon-MethyC.ultFiltMixed.bash** was used to deconvolute all cytosines in all subsamples

**dedupSubsamplerBismarkDecon.bash** was used to deduplicate bisulfite reads aligned to the reference genome with Bismark. This is a necessary step to call methylation levels with Bismark.

**makeMethylationBed-BismarkSubsampler-Decon.bash** was used to process methylation data for deconvoltuion analysis. This preprocessing allowed improved parallelization of deconvolution.

**pspDeduplicateExtractor-Subsampler-Decon.bash** was used to extract informaiton on aligned reads following deduplication

**pspDepthFinder-Subsampler-Decon.bash** was used to obtain informaiton on depth of aligned reads at various stages of preprocessing.

**ultFiltMixed-Boostrap-P.bash** was used to bootstrap analyses for calculation of P using subsamples of 100 SNVs

**ultFiltMixed-EstimateP.bash** was used to find somatic sites in blood and estimate levels of P

**validation-BasicStatsP.ultFiltMixed.bash** was used to prepare fiels for further validation and assessment of somatic sites used to estimate P