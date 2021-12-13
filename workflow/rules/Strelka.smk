rule Strelka:
  input:
    tumor="results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
    conf="config/strelka_config_bwa.ini",
    normal=lambda w: expand("results/alignment/{ctrl}/{ctrl}.realigned.recal.bam", ctrl=get_sample_control(w)),
  params:
    strelka="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/236aa367b1347b9561439ce4facd36c0/share/strelka-2.9.10-1/bin",
    outdir="results/Strelka/{sample}/{sample}",
    control=has_a_control,
  output:
    dir=directory("results/Strelka/{sample}/{sample}.myAnalysis"),
    indel="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/{sample}_Slk_somatic.indels.vcf.gz",
    snv="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/{sample}_Slk_somatic.snvs.vcf.gz",
    indeltbi="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/{sample}_Slk_somatic.indels.vcf.gz.tbi",
    snvtbi="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/{sample}_Slk_somatic.snvs.vcf.gz.tbi",
    tdir="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/strelka.yaml",
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        ## configuration
        {params.strelka}/configureStrelkaSomaticWorkflow.py \
        --normal={input.normal}  \
        --tumor={input.tumor} \
        --ref={input.ref} \
        --config={input.conf} \
        --runDir={output.dir}

        ## running pipeline
        {output.dir}/runWorkflow.py -m local -j 20

        echo "Current Dir: "$(pdir)
        mv {output.dir}/results/variants/somatic.indels.vcf.gz {output.indel}
        mv {output.dir}/results/variants/somatic.snvs.vcf.gz {output.snv}
        mv {output.dir}/results/variants/somatic.indels.vcf.gz.tbi {output.indeltbi}
        mv {output.dir}/results/variants/somatic.snvs.vcf.gz.tbi {output.snvtbi}
    else
        mkdr -p {output.tdir}
        touch {output.indel}
        touch {output.snv}
        touch {output.indeltbi}
        touch {output.snvtbi}
    fi
    """
