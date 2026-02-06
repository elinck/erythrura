rule all_sites_vcf:
    """
    Call variant + invariant sites from BAMs using bcftools mpileup.
    Produces an all-sites VCF suitable for pixy.
    """
    input:
        bamlist="/home/k14m234/erythrura/config/bamlist.txt",
        ref="/home/k14m234/erythrura_assembly/results/GCF_005870125.1/data/genome/GCF_005870125.1.fna",
        fai="/home/k14m234/erythrura_assembly/results/GCF_005870125.1/data/genome/GCF_005870125.1.fna.fai",
    output:
        vcf="results/all_sites/all_sites.vcf.gz",
        tbi="results/all_sites/all_sites.vcf.gz.tbi",
    params:
        min_bq=20,
        min_mq=30,
        max_depth=1000,   # global cap, adjust if needed
    threads: 8
    conda:
        "envs/bcftools.yaml"
    shell:
        r"""
        set -euo pipefail
        mkdir -p results/all_sites

        bcftools mpileup \
          -f {input.ref} \
          -b {input.bamlist} \
          -a AD,DP \
          -q {params.min_mq} \
          -Q {params.min_bq} \
          --max-depth {params.max_depth} \
          --threads {threads} \
          -Ou |

        bcftools call \
          -m \
          --keep-alts \
          --threads {threads} \
          -Ou |

        bcftools filter \
          -e 'INFO/DP==0 || QUAL<20' \
          -Oz \
          -o {output.vcf}

        tabix -f -p vcf {output.vcf}
        """

