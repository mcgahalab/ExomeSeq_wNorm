rule MuTect1:
  input:
    tumor = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    normal=lambda w: expand("results/alignment/{ctrl}/{ctrl}.realigned.recal.bam", ctrl=get_sample_control(w)),
    ref = 'ref/genome.fa',
    interval = region
  output:
    vcf="results/MuTect1/{sample}/{sample}.mut1.vcf",
    stats="results/MuTect1/{sample}/{sample}.call_stats",
    coverage="results/MuTect1/{sample}/{sample}.wig.txt"
  params:
    mutect="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/ce5c78a09c5970de54b46955738bc233/share/mutect=1.1.6-1",
    control=has_a_control,
  threads: 2
#  conda:
#    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/mutect1.yaml",
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        module load mutect/1.1.5
        java -Xmx12g -jar $mutect_dir/muTect.jar \
        --analysis_type MuTect \
        -R {input.ref} \
        -L {input.interval} \
        --input_file:tumor {input.tumor} \
        --input_file:normal {input.normal} \
        --vcf {output.vcf} \
        --out {output.stats} \
        --coverage_file {output.coverage} \
        --downsampling_type NONE \
        --fraction_contamination 0.02
    else
        mkdir -p results/MuTect1/{sample}/
        touch {output.vcf}
        touch {output.stats}
        touch {output.stats}
    fi
    """


#    java -Xmx12g -jar {params.mutect}/muTect-1.1.6.jar -T
