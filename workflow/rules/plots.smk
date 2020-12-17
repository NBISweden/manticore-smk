rule manticore_plotvcf:
    """Plot vcf file in various ways"""
    output:
        report("{results}/{group}/analysis/{analysis}/plots/{label}/plot.{type}.{region}.{ext}",
               caption="../report/plotvcf.rst", category="Variant plots", subcategory="{type}")
    input:
        unpack(manticore_plotvcf_input)
    wildcard_constraints:
        type = "(pca)"
    params:
        options = lambda wildcards: get_plot_options(wildcards)
    log:
        "logs/{results}/{group}/analysis/{analysis}/plots/{label}/plot.{type}.{region}.{ext}"
    script:
        "../scripts/manticore-plotvcf.py"
