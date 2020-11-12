rule all_qc:
    input: __REPORTS__ / "qc/multiqc.html",
           unpack(all_multiqc)

rule qc_multiqc:
    output: "{reports}/qc/multiqc.html"
    input: unpack(all_multiqc)
    resources:
        runtime = lambda wildcards, attempt: resources("qc_multiqc", "runtime", attempt)
    params: get_params("qc_multiqc", "options")
    log: "logs/{reports}/qc/multiqc.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/multiqc"


rule qc_fastqc:
    output: html = "{results}/qc/fastqc/{prefix}{bamfastq}{gz}.html",
            zip = "{results}/qc/fastqc/{prefix}{bamfastq}{gz}_fastqc.zip"
    input:  "{prefix}{bamfastq}{gz}"
    resources:
        runtime = lambda wildcards, attempt: resources("qc_fastqc", "runtime", attempt)
    params: get_params("qc_fastqc", "options")
    threads: get_params("qc_fastqc", "threads")
    log: "logs/{results}/qc/fastqc/{prefix}{bamfastq}{gz}.log"
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/fastqc"


rule qc_jellyfish_count:
    output: jf = temp("{results}/qc/jellyfish/{sample}{trimmed}{ext}{gz}.{kmer}mer_counts.jf")
    input: unpack(jellyfish_count)
    resources:
        runtime = lambda wildcards, attempt: resources("qc_jellyfish_count", "runtime", attempt)
    wildcard_constraints:
        ext = "(.fq|.fastq|.fa|.fasta|.txt|)",
        trimmed = "(|.trimmed)"
    params:
        options = get_params("qc_jellyfish_count", "options")
    threads: get_params("qc_jellyfish_count", "threads")
    log: "logs/{results}/qc/jellyfish/{sample}{trimmed}{ext}{gz}.{kmer}mer_counts.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/count"


rule qc_jellyfish_histo:
    output: hist = "{prefix}.{kmer}_jf.hist"
    input: counts = "{prefix}.{kmer}mer_counts.jf"
    resources:
        runtime = lambda wildcards, attempt: resources("qc_jellyfish_histo", "runtime", attempt)
    params:
        options = get_params("qc_jellyfish_histo", "options")
    threads: get_params("qc_jellyfish_histo", "threads")
    log: "logs/{prefix}.{kmer}mer_counts.jf.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/jellyfish/histo"


rule qc_picard_collect_alignment_summary_metrics:
    output: "{results}/qc/align/{prefix}{bam}.align_metrics.txt"
    input: bam = "{prefix}{bam}",
           ref = config['db']['ref']
    log: "logs/{results}/qc/align/{prefix}{bam}.align_metrics.log"
    resources:
        runtime = lambda wildcards, attempt: resources("qc_picard_collect_alignment_summary_metrics", "runtime", attempt)
    params:
        get_params('qc_picard_collect_alignment_summary_metrics', 'options')
    threads: get_params('qc_picard_collect_alignment_summary_metrics', 'threads')
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/collectalignmentsummarymetrics"


rule qc_picard_collect_insert_size_metrics:
    output: txt = "{results}/qc/align/{prefix}{bam}.insert_metrics.txt",
            pdf = "{results}/qc/align/{prefix}{bam}.insert_metrics.pdf",
    input: "{prefix}{bam}"
    log: "logs/{results}/qc/align/{prefix}{bam}.insert_metrics.log"
    resources:
        runtime = lambda wildcards, attempt: resources("qc_picard_collect_insert_size_metrics", "runtime", attempt)
    params:
        get_params('qc_picard_collect_insert_size_metrics', 'options')
    threads: get_params('qc_picard_collect_insert_size_metrics', 'threads')
    wrapper: f"{SMK_WRAPPER_PREFIX}/bio/picard/collectinsertsizemetrics"


rule qc_picard_mark_duplicates:
    output: metrics = "{interim}/map/bwa/dedup/{sample}{bam}.dup_metrics.txt",
            bam = "{interim}/map/bwa/dedup/{sample}{bam}"
    input: "{interim}/map/bwa/{sample}{bam}"
    log: "logs/{interim}/qc/align/{sample}{bam}.dup_metrics.log"
    resources:
        runtime = lambda wildcards, attempt: resources("qc_picard_mark_duplicates", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("qc_picard_mark_duplicates", "mem_mb", attempt)
    params:
        get_params('qc_picard_mark_duplicates', 'options')
    threads: get_params('qc_picard_mark_duplicates', 'threads')
    wrapper: f"{WRAPPER_PREFIX}/bio/picard/markduplicates"


rule qc_qualimap_bamqc_pe:
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
    input: bam = map_sample_target
    params:
        options = get_params("qc_qualimap_bamqc_pe", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("qc_qualimap_bamqc_pe", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("qc_qualimap_bamqc_pe", "mem_mb", attempt)
    threads: get_params("qc_qualimap_bamqc_pe", "threads")
    log: "logs/{results}/qc/qualimap/{sample}{bam}.pe.qualimap.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/qualimap/bamqc"


rule qc_qualimap_bamqc_se:
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
    input: bam = map_sample_target
    params:
        options = get_params("qc_qualimap_bamqc_se", "options")
    resources:
        runtime = lambda wildcards, attempt: resources("qc_qualimap_bamqc_se", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("qc_qualimap_bamqc_se", "mem_mb", attempt)
    threads: get_params("qc_qualimap_bamqc_se", "threads")
    log: "logs/{results}/qc/qualimap/{sample}{bam}.se.qualimap.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/qualimap/bamqc"


rule qc_sambamba_depth:
    output: bed = "{results}/qc/sambamba/{sample}.depth.w{window_size}.bed{gz}"
    input: bam = map_sample_target,
           bai = map_sample_target_bai
    resources:
        runtime = lambda wildcards, attempt: resources("qc_sambamba_depth", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("qc_sambamba_depth", "mem_mb", attempt)
    params: window_size = lambda wildcards: wildcards.window_size
    threads: get_params("qc_sambamba_depth", "threads")
    log: "logs/{results}/qc/sambamba/{sample}.depth.w{window_size}.bed{gz}.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/sambamba/depth"


rule qc_bcftools_stats:
    """Run bcftools stats to generate caller statistics on single files"""
    output: stats = "{results}/qc/variantstats/{group}/{callset}/{caller}/{stage}/{region}{vartype}.vcf{gz}.stats"
    input: vcf = "{results}/{group}/{callset}/{caller}/{stage}/{region}{vartype}.vcf{gz}",
           ref = config["db"]["ref"]
    wildcard_constraints:
        stage = "(unfiltered|filter|select)",
        vartype = "(|.indel|.snp)"
    resources:
        runtime = lambda wildcards, attempt: resources("qc_bcftools_stats", "runtime", attempt),
        mem_mb = lambda wildcards, attempt: resources("qc_bcftools_stats", "mem_mb", attempt)
    threads: 1
    log: "logs/{results}/qc/bcftools/{group}/{callset}/{caller}/{stage}/{region}{vartype}.vcf{gz}.stats.log"
    wrapper: f"{WRAPPER_PREFIX}/bio/bcftools/stats"
