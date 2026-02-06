############################################
# Include rules
############################################

include:
    "workflow/rules/all_sites_vcf.smk"
    # later:
    # "workflow/rules/pixy.smk"

############################################
# Default target
############################################

rule all:
    input:
        "results/all_sites/all_sites.vcf.gz.tbi"

