# adapted from: https://gatk.broadinstitute.org/hc/en-us/articles/360035535912
# https://github.com/gatk-workflows/gatk4-data-processing/blob/master/processing-for-variant-discovery-gatk4.wdl

rule mapFASTQ:
  input:
    f1=get_map_reads_input_R1,
    f2=get_map_reads_input_R2
  output:
    bam=temp("results/alignment/{sample}/{sample}_sorted.bam"),
    bai=temp("results/alignment/{sample}/{sample}_sorted.bam.bai")
  params:
    ref = config['ref_index']['genome'],
    conda=config['env']['conda_shell'],
    env=directory(config['env']['r41'])
  threads: 14
  shell:
    """
    module load bwa/0.7.15
    module load samtools/1.17
    
    bwa mem -K 100000000 \
    -v 3 -t 14 -Y \
    -R "@RG\\tID:{wildcards.sample}\\tLB:Genome\\tSM:{wildcards.sample}\\tPL:ILLUMINA" \
    {params.ref} \
    {input.f1} \
    {input.f2} | \
    samtools sort -@14 - > {output.bam};
    
    samtools index {output.bam} > {output.bai}
    """

rule extractUnmapped:
  input:
    bam="results/alignment/{sample}.aligned.duplicate_marked.recalibrated.bam",
  output:
    unmappedbam=temp("results/alignment/{sample}/{sample}_unmapped.bam"),
    unpaired="results/alignment/{sample}/{sample}_unmapped.unpaired.fastq",
    fq1="results/alignment/{sample}/{sample}_unmapped.R1.fastq",
    fq2="results/alignment/{sample}/{sample}_unmapped.R2.fastq"
  threads: 2
  params:
    seqtkdir=config['env']['seqtkdir']
  shell:
    """
    module load samtools/1.17
    module load picard/2.10.9
    
    ## Extract unmapped reads
    samtools view -b -f4  {input.bam} | \
     samtools collate - -O > {output.unmappedbam}
    
    ## Extract fastq reads from the unmapped bam
    java -jar $picard_dir/picard.jar \
         SamToFastq \
         INPUT={output.unmappedbam} \
         FASTQ={output.fq1} \
         SECOND_END_FASTQ={output.fq2} \
         UNPAIRED_FASTQ={output.unpaired} \
         VALIDATION_STRINGENCY=SILENT
    
    """
    #/cluster/projects/mcgahalab/bin/seqtk/seqtk subseq $(zcat PANX_1213.R1.merged.fastq.gz) unmapped_ids.lst > PANX_1213.unmapped.R1.fastq

rule AlignUnmappedToCHM13:
  input:
    fq1="results/alignment/{sample}/{sample}_unmapped.R1.fastq",
    fq2="results/alignment/{sample}/{sample}_unmapped.R2.fastq"
  output:
    bam="results/alignment/{sample}/{sample}_chm13.bam",
    unmappedbam=temp("results/alignment/{sample}/{sample}_chm13.unmapped.bam"),
    unpaired="results/alignment/{sample}/{sample}_chm13.unmapped.unpaired.fastq",
    fq1="results/alignment/{sample}/{sample}_chm13.unmapped.R1.fastq",
    fq2="results/alignment/{sample}/{sample}_chm13.unmapped.R2.fastq"
  threads: 2
  params:
    chm13genome=config['ref_index']['chm13genome']
  shell:
    """
    module load bowtie2/2.4.5
    module load samtools/1.17
    module load picard/2.10.9
    
    bowtie2 -p 2 -q --local \
        -x {params.chm13genome} \
        -1 {input.fq1} \
        -2 {input.fq2} | \
        samtools view -bS - > {output.bam}
    
    ## Extract unmapped reads
    samtools view -b -f4  {output.bam} | \
     samtools collate - -O > {output.unmappedbam}
    
    ## Extract fastq reads from the unmapped bam
    java -jar $picard_dir/picard.jar \
         SamToFastq \
         INPUT={output.unmappedbam} \
         FASTQ={output.fq1} \
         SECOND_END_FASTQ={output.fq2} \
         UNPAIRED_FASTQ={output.unpaired} \
         VALIDATION_STRINGENCY=SILENT
    """

rule KrakenUniq:
  input:
    fq1="results/alignment/{sample}/{sample}_chm13.unmapped.R1.fastq",
    fq2="results/alignment/{sample}/{sample}_chm13.unmapped.R2.fastq"
  output:
    reportfile="results/metagenomics/reportfiles/reportfile.{sample}.tsv",
    readclassification="results/metagenomics/classifications/readclassifcation.{sample}.tsv",
    classified="results/metagenomics/classifications/classified.{sample}.txt",
  threads: 2
  params:
    krakencontainer=config['env']['krakencontainer']
    krakendb=config['db']['kraken']
  shell:
    """
    module load apptainer/1.0.2
    
    apptainer exec {params.krakencontainer}/krakenuniq.sif \
    krakenuniq --db {params.krakendb} \
    --paired \
    --classified-out {output.classified} \
    --threads {threads} \
    --report-file {output.reportfile} \
    {input.fq1} {input.fq2} \
    > {output.readclassification}
    """

rule MarkDuplicates:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.bam",
    bai="results/alignment/{sample}/{sample}_sorted.bam.bai"
  output:
    dedup=temp("results/alignment/{sample}/{sample}_sorted.dedup.bam"),
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
    tmpdir = "results/alignment/{sample}/markdup"
  threads: 4
  shell:
    """
    module load gatk/4.2.5.0
    
    mkdir -p {params.tmpdir}
    
    gatk --java-options "-Xmx25g" \
      MarkDuplicates \
      --INPUT {input.bam} \
      --OUTPUT {output.dedup} \
      --METRICS_FILE {output.metrics} \
      --VALIDATION_STRINGENCY SILENT \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
      --ASSUME_SORT_ORDER "queryname" \
      --TMP_DIR {params.tmpdir} \
      --CREATE_MD5_FILE true
    """

rule SortAndFixTags:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.dedup.bam"
  output:
    bam=temp("results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam"),
    tmp=temp("results/alignment/{sample}/{sample}.tmp.bam")
  params:
    ref = config['ref_index']['genome'],
    conda = config['env']['conda_shell'],
    env = config['env']['picard'],
    tmpdir = "results/alignment/{sample}"
  threads: 4
  shell:
    """
    source {params.conda} && conda activate {params.env};
    module load gatk/4.2.5.0
    
    gatk --java-options "-Xmx25g" \
      SortSam \
      --INPUT {input.bam} \
      --OUTPUT {output.tmp} \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX false \
      --TMP_DIR {params.tmpdir} \
      --CREATE_MD5_FILE false;
           
    gatk --java-options "-Xmx12g" \
      SetNmMdAndUqTags \
      --INPUT {output.tmp} \
      --OUTPUT {output.bam} \
      --CREATE_INDEX true \
      --CREATE_MD5_FILE true \
      --REFERENCE_SEQUENCE {params.ref}
    """

rule BaseRecalibrator:
  input:
    bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam"
  output:
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
    
    gatk --java-options "-Xmx25g" \
      BaseRecalibrator \
      -R {params.ref} \
      -I {input.bam} \
      --use-original-qualities \
      -O {output.metrics} \
      --known-sites {params.dbsnp} \
      --known-sites {params.dbmills} \
      --known-sites {params.dbindel}
    """

rule ApplyBQSR:
  input:
    bam="results/alignment/{sample}/{sample}.aligned.duplicate_marked.sorted.bam",
    recal="results/alignment/{sample}/{sample}.recal_data.csv"
  output:
    bam="results/alignment/{sample}.aligned.duplicate_marked.recalibrated.bam"
  params:
    ref = config['ref_index']['genome'],
    dbsnp = config['snvdb']['dbsnp'],
    dbmills = config['snvdb']['indels'],
    dbindel = config['snvdb']['indels2'],
  threads: 4
  shell:
    """
    module load gatk/4.2.5.0
    
    gatk --java-options "-Xmx25g" \
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
