rule mapFASTQ:
  input:
    f1 =  get_r1,
    f2 =  get_r2,
    ref = 'ref/BWAgenome.fa'
  output: temp("results/alignment/{sample}/{sample}.sam")
  params:
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    bwa mem -p -t4 -R "@RG\\tID:{wildcards.sample}\\tLB:Exome\\tSM:{wildcards.sample}\\tPL:ILLUMINA" {input.ref} {input.f1} {input.f2} > {output}
    """
rule samtoolsSORT:
  input: "results/alignment/{sample}/{sample}.sam"
  output: "results/alignment/{sample}/{sample}_sorted.bam"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    samtools sort -@4 {input} > {output}
    """

rule samtoolsINDEX:
  input: "results/alignment/{sample}/{sample}_sorted.bam"
  output: "results/alignment/{sample}/{sample}_sorted.bam.bai"
  threads: 2
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    samtools index {input} > {output}
    """

rule picardMarkDuplicates:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.bam",
    bai="results/alignment/{sample}/{sample}_sorted.bam.bai",
  output:
    dedup="results/alignment/{sample}/{sample}_sorted.dedup.bam",
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
    picard="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9b770440ff173434e53ee101c7452a0a/share/picard-2.26.0-0"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    java -Xmx12g -jar {params.picard}/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output.dedup} METRICS_FILE={output.metrics} ASSUME_SORTED=true MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
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
