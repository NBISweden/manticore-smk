rule all_qc:
    input: multiqc = __REPORTS__ / "qc/multiqc.html",
           fastqc = all_fastqc,
           jellyfish = all_jellyfish

rule multiqc:
    output: "{reports}/qc/multiqc.html"
    input: unpack(all_multiqc)
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['multiqc']['runtime']
    params: config["qc"]["multiqc"]["options"]
    log: "logs/{reports}/qc/multiqc.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/multiqc"


rule fastqc:
    output: html = "{results}/qc/fastqc/{prefix}{bamfastq}{gz}.html",
            zip = "{results}/qc/fastqc/{prefix}{bamfastq}{gz}_fastqc.zip"
    input:  "{prefix}{bamfastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['fastqc']['runtime']
    params: config["qc"]["fastqc"]["options"]
    threads: config["qc"]["fastqc"]["threads"]
    log: "logs/{results}/qc/fastqc/{prefix}{bamfastq}{gz}.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/fastqc"


rule jellyfish_count:
    output: jf = temp("{results}/qc/jellyfish/{sample}{trimmed}{ext}{gz}.{kmer}mer_counts.jf")
    input: unpack(jellyfish_count)
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['jellyfish']['count']['runtime']
    wildcard_constraints:
        ext = "(.fq|.fastq|.fa|.fasta|.txt|)",
        trimmed = "(|.trimmed)"
    params:
        options = config['qc']['jellyfish']['count']['options']
    threads: config['qc']['jellyfish']['count']['threads']
    log: "logs/{results}/qc/jellyfish/{sample}{trimmed}{ext}{gz}.{kmer}mer_counts.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/count"


rule jellyfish_histo:
    output: hist = "{prefix}.{kmer}_jf.hist"
    input: counts = "{prefix}.{kmer}mer_counts.jf"
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['jellyfish']['histo']['runtime']
    params:
        options = config['qc']['jellyfish']['histo']['options']
    threads: config['qc']['jellyfish']['histo']['threads']
    log: "logs/{prefix}.{kmer}mer_counts.jf.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/histo"


rule picard_collect_alignment_summary_metrics:
    output: "{results}/qc/align/{prefix}{bam}.align_metrics.txt"
    input: bam = "{prefix}{bam}",
           ref = config['db']['ref']
    log: "logs/{results}/qc/align/{prefix}{bam}.align_metrics.log"
    params:
        config['qc']['picard']['options']
    threads: config['qc']['picard']['collect_alignment_summary_metrics']['threads']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/collectalignmentsummarymetrics"


rule picard_collect_insert_size_metrics:
    output: txt = "{results}/qc/align/{prefix}{bam}.insert_metrics.txt",
            pdf = "{results}/qc/align/{prefix}{bam}.insert_metrics.pdf",
    input: "{prefix}{bam}"
    log: "logs/{results}/qc/align/{prefix}{bam}.insert_metrics.log"
    params:
        config['qc']['picard']['options']
    threads: config['qc']['picard']['collect_insert_size_metrics']['threads']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/collectinsertsizemetrics"


rule picard_mark_duplicates:
    output: metrics = "{interim}/map/bwa/dedup/{sample}{bam}.dup_metrics.txt",
            bam = "{interim}/map/bwa//dedup/{sample}{bam}"
    input: "{interim}/map/bwa/{sample}{bam}"
    log: "logs/{interim}/qc/align/{sample}{bam}.dup_metrics.log"
    params:
        config['qc']['picard']['mark_duplicates']['options']
    threads: config['qc']['picard']['mark_duplicates']['threads']
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/markduplicates"


rule qualimap_bamqc_pe:
    output:
        genome_results = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/genome_results.txt",
        html = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/qualimapReport.html",
        coverage_across_reference = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/coverage_across_reference.txt",
        coverage_histogram = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/coverage_histogram.txt",
        duplication_rate_histogram = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/duplication_rate_histogram.txt",
        genome_fraction_coverage = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/genome_fraction_coverage.txt",
        homopolymer_indels = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/homopolymer_indels.txt",
        mapped_reads_clipping_profile = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/mapped_reads_clipping_profile.txt",
        mapped_reads_gc_content_distribution = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        mapped_reads_nucleotide_content = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/mapped_reads_nucleotide_content.txt",
        mapping_quality_across_reference = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/mapping_quality_across_reference.txt",
        mapping_quality_histogram = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/mapping_quality_histogram.txt",
        insert_size_across_reference = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/insert_size_across_reference.txt",
        insert_size_histogram = "{results}/qc/qualimap/{sample}{bam}.pe.qualimap/raw_data_qualimapReport/insert_size_histogram.txt"
    input: bam = bwa_mem_sample
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['qualimap']['runtime'],
        mem_mb = lambda wildcards, attempt: attempt * config['qc']['qualimap']['mem_mb']
    threads: config['qc']['qualimap']['threads']
    log: "logs/{results}/qc/qualimap/{sample}{bam}.pe.qualimap.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/qualimap/bamqc"


rule qualimap_bamqc_se:
    output:
        genome_results = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/genome_results.txt",
        html = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/qualimapReport.html",
        coverage_across_reference = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/coverage_across_reference.txt",
        coverage_histogram = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/coverage_histogram.txt",
        duplication_rate_histogram = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/duplication_rate_histogram.txt",
        genome_fraction_coverage = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/genome_fraction_coverage.txt",
        homopolymer_indels = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/homopolymer_indels.txt",
        mapped_reads_clipping_profile = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/mapped_reads_clipping_profile.txt",
        mapped_reads_gc_content_distribution = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt",
        mapped_reads_nucleotide_content = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/mapped_reads_nucleotide_content.txt",
        mapping_quality_across_reference = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/mapping_quality_across_reference.txt",
        mapping_quality_histogram = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/mapping_quality_histogram.txt",
        insert_size_across_reference = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/insert_size_across_reference.txt",
        insert_size_histogram = "{results}/qc/qualimap/{sample}{bam}.se.qualimap/raw_data_qualimapReport/insert_size_histogram.txt"
    input: bam = bwa_mem_sample
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['qualimap']['runtime'],
        mem_mb = lambda wildcards, attempt: attempt * config['qc']['qualimap']['mem_mb']
    threads: config['qc']['qualimap']['threads']
    log: "logs/{results}/qc/qualimap/{sample}{bam}.se.qualimap.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/qualimap/bamqc"


rule sambamba_depth:
    output: bed = "{results}/qc/sambamba/{sample}.depth.w{window_size}.bed{gz}"
    input: bam = bwa_mem_sample,
           bai = bwa_mem_sample_bai
    resources:
        runtime = lambda wildcards, attempt: attempt * config['qc']['sambamba']['runtime'],
        mem_mb = lambda wildcards, attempt: attempt * config['qc']['sambamba']['mem_mb']
    params: window_size = lambda wildcards: wildcards.window_size
    threads: config['qc']['sambamba']['threads']
    log: "logs/{results}/qc/sambamba/{sample}.depth.w{window_size}.bed{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/sambamba/depth"
