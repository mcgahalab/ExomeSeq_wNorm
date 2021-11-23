rule MuTect1:
  input:
    bam = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
    interval = region
  output:
    vcf="results/MuTect1/{sample}/{sample}.mut1.vcf",
    stats="results/MuTect1/{sample}/{sample}.call_stats",
    coverage="results/MuTect1/{sample}/{sample}.wig.txt"
  params:
    mutect="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/ce5c78a09c5970de54b46955738bc233/share/mutect=1.1.6-1"
  threads: 2
#  conda:
#    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/mutect1.yaml",
  shell:
    """
    module load mutect/1.1.5
    java -Xmx12g -jar $mutect_dir/muTect.jar \
    --analysis_type MuTect \
    -R {input.ref} \
    -L {input.interval} \
    --input_file:tumor {input.bam} \
    --vcf {output.vcf} \
    --out {output.stats} \
    --coverage_file {output.coverage} \
    --downsampling_type NONE \
    --fraction_contamination 0.02
    """


#    java -Xmx12g -jar {params.mutect}/muTect-1.1.6.jar -T
