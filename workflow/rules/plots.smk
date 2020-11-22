rule manticore_plotvcf:
    """Plot vcf file in various ways"""
    output:
        report("{results}/{group}/analysis/{analysis}/plots/{plotdir}/plot.{type}.{region}.png", caption="../reports/plotvcf.rst")
    input:
        unpack(get_plot_input)
    wildcard_constraints:
        type = "(pca)"
    log:
        "logs/{results}/{group}/analysis/{analysis}/plots/{plotdir}/plot.{type}.{region}.log"
    script:
        "../scripts/manticore-plotvcf.py"
