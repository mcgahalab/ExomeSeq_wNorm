
rule haploGATK:
  input:
    bam="results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
  output: "results/Haplotype/{sample}/{sample}.raw.snps.indels.vcf"
  params:
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml"
  shell:
    """
    gatk3 -Xmx12g -T HaplotypeCaller -nct 4 -R {input.ref} -I {input.bam} -o {output}
    """
