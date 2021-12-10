rule varscanCopyNumber:
  input:
    ref = 'ref/genome.fa',
    normal =lambda w: expand("results/alignment/{ctrl}/{ctrl}.realigned.recal.bam", ctrl=get_sample_control(w)),
    tumor = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    bed = region
  params:
    outdir="results/Varscan/cnv/{sample}/{sample}.vscn",
    varscan="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/e19f0037a250972dbd4b858a29e1bd01/share/varscan-2.4.4-1",
    control=has_a_control,
  output: "results/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
  threads: 3
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        samtools mpileup -B -q 1 -d 1000000 -l {input.bed} -f {input.ref} {input.normal} {input.tumor} | \
        java -Xmx12g -jar {params.varscan}/VarScan.jar copynumber - {params.outdir} --mpileup 1
    else
        touch {output}
    fi
    """
rule varscanSomatic:
  input:
    ref = 'ref/genome.fa',
    normal = lambda w: expand("results/alignment/{ctrl}/{ctrl}.realigned.recal.bam", ctrl=get_sample_control(w)),
    tumor = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    bed = region
  output:
    snp="results/Varscan/snv/{sample}/{sample}.snp",
    indel="results/Varscan/snv/{sample}/{sample}.indel"
  params:
    varscan="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/e19f0037a250972dbd4b858a29e1bd01/share/varscan-2.4.4-1",
    control=has_a_control,
  threads: 3
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        samtools mpileup -B -q 1 -d 1000000 -l {input.bed} -f {input.ref} {input.normal} {input.tumor} | \
        java -Xmx12g -jar {params.varscan}/VarScan.jar somatic --mpileup 1 --output-snp {output.snp} --output-indel {output.indel}
    else
        touch {output.snp}
        touch {outpt.indel}
    fi
    """

rule varscanProcessSomatic:
  input:
    snp="results/Varscan/snv/{sample}/{sample}.snp",
    indel="results/Varscan/snv/{sample}/{sample}.indel"
  output:
    snp="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc",
    indel="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc"
  params:
    varscan="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/e19f0037a250972dbd4b858a29e1bd01/share/varscan-2.4.4-1",
    control=has_a_control,
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        java -Xmx8g -jar {params.varscan}/VarScan.jar processSomatic {input.snp} {output.snp}
        java -Xmx8g -jar {params.varscan}/VarScan.jar processSomatic {input.indel} {output.indel}
    else
        touch {output.snp}
        touch {output.indel}
    fi
    """

rule varscanToVCF:
  input:
    snp="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc",
    indel="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc"
  output:
    snp="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc.vcf",
    indel="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc.vcf"
  params:
    control=has_a_control,
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        python scripts/varscan2vcf.py {input.snp} > {output.snp}
        python scripts/varscan2vcf.py {input.indel} > {output.indel}
    else
        touch {output.snp}
        touch {output.indel}
    fi
    """
