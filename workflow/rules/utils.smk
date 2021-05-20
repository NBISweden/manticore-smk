rule manticore_save_config:
    """Save manticore configuration"""
    output: report("results/config/manticore.config.yaml", caption="../report/config.rst", category="Configuration")
    log: "logs/manticore/manticore_save_config.log"
    script: "../scripts/manticore_save_config.py"


rule manticore_save_readme:
    """Save manticore readme"""
    output: report("README.md", caption="../report/readme.rst", category="Documentation")
    log: "logs/manticore_save_readme.log"
    script: "../scripts/manticore_save_readme.py"


localrules: manticore_save_config, manticore_save_readme
