This subdirectory contains all of the scripts used to obtain the great tit reference genome GCF_001522545.3 and annotate features contained within it for the manuscript 'Identification of 5mC within heterogenous tissue using de-novo somatic mutations'


The contained files were used for the purposes listed below:

**Abel-BloodBAM2Variants.bash** was used to call variants on the refernce genome with mpileup command using a bam file generated from the non-bisulfite blood reads from 'Derks, M.F., Schachtschneider, K.M., Madsen, O., Schijlen, E., Verhoeven, K.J. and van Oers, K., 2016. Gene and transposable element methylation in great tit (Parus major) brain and blood. Bmc Genomics, 17, pp.1-13.'

**Abel-CpG-IslandsEMBOSS.bash** was used to annotate CpG islands on the reference genome using EMBOSS

**Abel-ReferenceDownload.bash** was used to obtain the reference genome

**Abel-RepeatMaskerDFAM3.7Aves.bash** was used to annotate repeats on the reference genome using the DFAM3.7 database and nhmmer

**Abel-RepeatMaskerRM2.0.4DFAM3.7Aves.bash** was used to annotate repeats on the reference genome using a *de novo* library created with RepeatModeler2

**Abel-RepeatModel2DFAM3.7.bash** was used to generate the *de novo* repeat library with RepeatModeler2

**Abel-SegDupBiser.bash** was used to annotate segmental duplications on the reference genome with Biser
