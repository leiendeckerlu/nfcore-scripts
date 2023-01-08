# Command collection for calling nf-core pipelines
# singularity is pre-loaded on HPC CBE
# generally genomes hg38 for human samples and mm38 for mouse samples (BL6)

# Run nf-core/atacseq pipeline on sample list provided in ../samples/samples.csv using the reference genome in ../genomes/genome.fasta
# ATACseq QC and peak calling pipelinee
# https://nf-co.re/atacseq
# with HPC CBE-specific settings

module load nextflow/21.10.6
nextflow run nf-core/atacseq -profile cbe --input samples.csv --design ./design.csv --genome hg38 -resume --save_reference --save_trimmed --narrow_peak --outdir hg38_results_narrowPeak

# Run nf-core/viralrecon pipeline on sample list provided in ../samples/samples.csv using the reference genome in ../genomes/genome.fasta
# Rapid analysis of targeted/non-targeted virus NGS data (QC, variants, consensus genome)
# https://nf-co.re/viralrecon
# with HPC CBE-specific settings

module load nextflow/21.10.6
nextflow run nf-core/viralrecon -latest --input ../samples/samples.csv --platform illumina --skip_assembly --spades_mode meta --save_trimmed_fail --max_memory 400.GB --skip_markduplicates 0 --skip_snpeff --skip_bandage --max_cpus 20 --max_time 8h --protocol metagenomic --fasta ../genomes/genome.fasta --skip_plasmidid --skip_nextclade --skip_pangolin --save_mpileup --skip_assembly_quast -profile cbe -bg -resume --outdir ./ --skip_kraken2 -r 2.5 --skip_fastp


# Run nf-core/fetchngs pipeline on sample list provided in id.txt
# Rapid download of NGS data from GEO repository
# https://nf-co.re/fetchngs
# with HPC CBE-specific settings

module load nextflow/21.10.6
nextflow run nf-core/fetchngs --input id.txt -profile cbe --force_sratools_download

# Run nf-core/hlatyping pipeline on samples provided as FASTQ.GZ in /RNAseq/FASTQ/. Input format is here RNAseq.
# Identification of HLA type
# https://nf-co.re/hlatyping
# with HPC CBE-specific settings

module load nextflow/21.10.6
nextflow run nf-core/hlatyping -profile cbe --input './RNAseq/FASTQ/*_R{1,2}.fastq.gz' --outdir ./ --max_cpus 20 --max_time 8h

# Run nf-core/sarek pipeline on samples provided as ./samples.csv 
# Genome mapping and variant calling pipeline
# https://nf-co.re/sarek
# with HPC CBE-specific settings

module load nextflow/21.10.6

# genome alignment: BWA-MEM
nextflow run nf-core/sarek -profile cbe -latest --input ./samples.csv --outdir ./ --genome hg38 -bg -resume --max_memory 400.GB --step mapping --only_paired_variant_calling --skip_tools baserecalibrator --max_cpus 20 --max_time 12h

# variant calling: Strelka (SNPs); Manta (SVs); MSISensor (Microsatallite instability)
nextflow run nf-core/sarek -profile cbe -latest --input ./samples.csv --outdir ./ --genome hg38 -bg -resume --max_memory 400.GB --step variant_calling --only_paired_variant_calling --max_cpus 20 --max_time 12h --tools strelka,manta,msisensorpro
