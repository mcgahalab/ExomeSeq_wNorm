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
  shell:
    """
    module load samtools/1.17
    
    samtools sort -@4 {input} > {output}
    """

rule samtoolsINDEX:
  input: "results/alignment/{sample}/{sample}_sorted.bam"
  output: "results/alignment/{sample}/{sample}_sorted.bam.bai"
  threads: 2
  shell:
    """
    module load samtools/1.17
    
    samtools index {input} > {output}
    """

rule MarkDuplicates:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.bam",
    bai="results/alignment/{sample}/{sample}_sorted.bam.bai"
  output:
    dedup="results/alignment/{sample}/{sample}_sorted.dedup.bam",
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
  threads: 4
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
    bam="results/alignment/{sample}/{sample}_sorted.dedup.bam"
  output:
    bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam",
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
    ref = config['ref_index']['genome'],
  threads: 4
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

rule BaseRecalibrator:
  input:
    bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam"
  output:
    # bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam",
    metrics="results/alignment/{sample}/{sample}.recal_data.csv"
  params:
    ref = config['ref_index']['genome'],
    dbsnp = config['snvdb']['dbsnp'],
    dbmills = config['snvdb']['indels'],
    dbindel = config['snvdb']['indels2'],
  threads: 4
  shell:
    """
    module load gatk/4.2.5.0
    
    gatk --java-options "-Xmx12g" \
      BaseRecalibrator \
      -R {params.ref} \
      -I {input.bam} \
      --use-original-qualities \
      -O {output.metrics} \
      --known-sites {params.dbsnp} \
      --known-sites {params.dbmills} \
      --known-sites {params.dbindel} \
      -rf BadCigar \
      -cov ReadGroupCovariate \
      -cov ContextCovariate \
      -cov CycleCovariate \
      -cov QualityScoreCovariate \
      -dt None
    """

rule ApplyBQSR:
  input:
    bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam",
    recal="results/alignment/{sample}/{sample}.recal_data.csv"
  output:
    bam="results/alignment/{sample}/{sample}..aligned.duplicates_marked.recalibrated.bam"
  params:
    ref = config['ref_index']['genome'],
    dbsnp = config['snvdb']['dbsnp'],
    dbmills = config['snvdb']['indels'],
    dbindel = config['snvdb']['indels2'],
  threads: 4
  shell:
    """
    module load gatk/4.2.5.0
    
    gatk --java-options "-Xmx12g" \
      ApplyBQSR \
      -R {params.ref} \
      -I {input.bam} \
      -O {output.bam} \
      -bqsr {input.recal} \
      --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
      --add-output-sam-program-record \
      --create-output-bam-md5 \
      --use-original-qualities
    """
