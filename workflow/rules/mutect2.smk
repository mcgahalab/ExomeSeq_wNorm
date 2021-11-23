rule MuTect2:
  input:
    bam = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
  params:
    intervals = get_intervals
  output: "results/MuTect2/{sample}/{sample}_{interval}.mut2.vcf"
  params:
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  threads: 2
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    gatk3 -Xmx12g \
    -T MuTect2  \
    -R {input.ref} \
    -L {params.intervals} \
    --input_file:tumor {input.bam} \
    -o {output} \
    """

rule MuTect2Merge:
  input: get_MuTect2_output,
  params:
    script="scripts/concatvcfs",
    out="results/MuTect2Merge/{sample}/",
    samp="{sample}"
  output: "results/MuTect2Merge/{sample}/{sample}_merged_mut2.vcf"
  threads: 2
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    ff=$(find {input} -type f -name '*.mut2.vcf')
    echo 'sh {params.script}' \
    $ff$' > {output}' > {params.out}/'merge_{params.samp}_VCFs.sh'
    sh {params.out}/merge_{params.samp}_VCFs.sh
    """

rule filterMuTect2:
  input:
    vcf = "results/MuTect2Merge/{sample}/{sample}_merged_mut2.vcf",
  params:
    outdirsnv="results/MuTect2Merge/{sample}/{sample}.snvs",
    outdirindel="results/MuTect2Merge/{sample}/{sample}.indels",
  output:
    snv="results/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    indel="results/MuTect2Merge/{sample}/{sample}.indels.recode.vcf"
  threads: 3
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    vcftools --vcf {input.vcf} --remove-indels --recode --recode-INFO-all --out {params.outdirsnv} --remove-filtered-all
    vcftools --vcf {input.vcf} --keep-only-indels --recode --recode-INFO-all --out {params.outdirindel} --remove-filtered-all
    """
