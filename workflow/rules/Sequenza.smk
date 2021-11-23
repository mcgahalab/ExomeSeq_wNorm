
rule Sequenza:
  input:
    bam="results/alignment/{sample}/{sample}.realigned.recal.bam",
    normal=norm,
  params:
    refgc='ref/GCgenome.wig',
    ref='ref/genome.fa'
  output: "results/Sequenza/{sample}.gz"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/sequenza.yaml",
  shell:
    """
    sequenza-utils \
    bam2seqz -t {input.bam} \
    -n {input.normal} \
    -gc {params.refgc} \
    -F {params.ref} \
    -o {output}
    """
