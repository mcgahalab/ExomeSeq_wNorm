# adapted from: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912
# https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl

rule mapFASTQ:
  input:
    f1=get_map_reads_input_R1,
    f2=get_map_reads_input_R2
  output: temp("results/alignment/{sample}/{sample}.sam")
  params:
    ref = config['ref_index']['genome'],
    conda=config['env']['conda_shell'],
    env=directory(config['env']['r41'])
  threads: 16
  conda:
    # "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    module load bwa/0.7.15
    
    bwa mem -K 100000000 \
    -p -v 3 -t 16 -Y \
    -R "@RG\\tID:{wildcards.sample}\\tLB:Exome\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
    {params.ref} \
    {input.f1} \
    {input.f2} > {output}
    """

rule samtoolsSORT:
  input: "results/alignment/{sample}/{sample}.sam"
  output: "results/alignment/{sample}/{sample}_sorted.bam"
  threads: 4
  conda:
    # "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    module load samtools/1.17
    
    samtools sort -@4 {input} > {output}
    """

rule samtoolsINDEX:
  input: "results/alignment/{sample}/{sample}_sorted.bam"
  output: "results/alignment/{sample}/{sample}_sorted.bam.bai"
  threads: 2
  conda:
    # "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    module load samtools/1.17
    
    samtools index {input} > {output}
    """

rule MarkDuplicates:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.bam",
    bai="results/alignment/{sample}/{sample}_sorted.bam.bai",
  output:
    dedup="results/alignment/{sample}/{sample}_sorted.dedup.bam",
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
  threads: 4
  conda:
    # "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    module load gatk/4.2.5.0
    
     gatk --java-options "-Xmx12g" \
      MarkDuplicates \
      --INPUT {input.bam} \
      --OUTPUT {output.dedup} \
      --METRICS_FILE {output.metrics} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --CREATE_MD5_FILE true
    """

rule SortAndFixTags:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.dedup.bam",
  output:
    bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam",
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
    ref = config['ref_index']['genome'],
  threads: 4
  conda:
    # "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    module load gatk/4.2.5.0
    
    gatk --java-options "-Xmx12g" \
      SortSam \
      --INPUT {input.bam} \
      --OUTPUT /dev/stdout \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --CREATE_MD5_FILE false \
      | \
      gatk4 --java-options "-Xmx12g" \
      SetNmMdAndUqTags \
      --INPUT /dev/stdin \
      --OUTPUT {output.bam} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE {params.ref}
    """

rule gatkRealignerTargetCreator:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.dedup.bam",
    ref= 'ref/genome.fa',
    region=region,
    known1=known_mills,
    known2=known_1000G,
  output: "results/alignment/{sample}/{sample}.IndelRealigner.intervals"
  params:
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    module load gatk/4.2.5.0
    
    gatk3 -Xmx8g -T RealignerTargetCreator \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -nt 4 \
    -I {input.bam} \
    -R {input.ref} \
    --intervals {input.region} \
    --interval_padding 100 \
    -known {input.known1} \
    -known {input.known2} \
    -dt None \
    -o {output}
    """

rule gatkIndelRealigner:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.dedup.bam",
    ref= 'ref/genome.fa',
    interval="results/alignment/{sample}/{sample}.IndelRealigner.intervals",
    known1=known_mills,
    known2=known_1000G,
  output: "results/alignment/{sample}/{sample}.realigned.bam"
  params:
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  threads: 2
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    gatk3 -Xmx12g -T IndelRealigner \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -I {input.bam} \
    -o {output} \
    -R {input.ref} \
    -targetIntervals {input.interval} \
    -known {input.known1} \
    -known {input.known2} \
    -dt None \
    -compress 0
    """

rule gatkBaseRecalibrator:
  input:
    bam="results/alignment/{sample}/{sample}.realigned.bam",
    ref= 'ref/genome.fa',
    dbsnp=dbsnp,
    region=region,
  output: "results/alignment/{sample}/{sample}.recal_data.grp"
  params:
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    gatk3 -Xmx12g -T BaseRecalibrator \
    -nct 4 \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -I {input.bam} \
    -o {output} \
    -R {input.ref} \
    -knownSites {input.dbsnp} \
    -knownSites {input.region} \
    -rf BadCigar \
    -cov ReadGroupCovariate \
    -cov ContextCovariate \
    -cov CycleCovariate \
    -cov QualityScoreCovariate \
    -dt None
    """

rule gatkPrintReads:
  input:
    bam="results/alignment/{sample}/{sample}.realigned.bam",
    ref= 'ref/genome.fa',
    recal="results/alignment/{sample}/{sample}.recal_data.grp"
  output: "results/alignment/{sample}/{sample}.realigned.recal.bam"
  params:
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    gatk3 -Xmx12g -T PrintReads \
    --disable_auto_index_creation_and_locking_when_reading_rods \
    -nct 4 \
    -I {input.bam} \
    -o {output} \
    -R {input.ref} \
    -BQSR {input.recal} \
    -rf BadCigar \
    -dt None
    """
