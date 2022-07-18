
rule vcftoMAFsnv:
  input:
    ref = 'ref/genome.fa',
    vcf_inter = "results/vcfIntersect/snvs/{sample}/{snv}.vcf",
  params:
    samp="{sample}",
    snvs = "{snv}",
#    snv = get_snvs,
#    intersect = get_snv_intersects,
  output:
    vcf_fil = "results/vcfIntersect/snvs/{sample}/fil_{snv}.vcf",
    maf = "results/MAF_38_final/snv/{sample}/{snv}.maf",
  threads: 4
  conda:
    "../envs/VCFtoMAF.yaml",
  shell:
    """
    if [ {params.snvs} == '0003' ]; then
      bcftools view -f PASS {input.vcf_inter} > {output.vcf_fil};
      perl scripts/vcf2maf.pl \
        --input-vcf {output.vcf_fil} \
        --output-maf {output.maf} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 100 \
        --ref-fasta={input.ref} \
        --filter-vcf ref/VEP_cache/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=ref/98 \
        --vep-data=ref/98 \
        --normal-id unmatched \
        --vcf-tumor-id TUMOR \
        --vcf-normal-id NORMAL

    elif [ {params.snvs} == '0001' ]; then
      bcftools view -f PASS {input.vcf_inter} > {output.vcf_fil};
      perl scripts/vcf2maf.pl \
        --input-vcf {output.vcf_fil} \
        --output-maf {output.maf} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 100 \
        --ref-fasta={input.ref} \
        --filter-vcf ref/VEP_cache/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=ref/98 \
        --vep-data=ref/98 \
        --normal-id unmatched \
        --vcf-tumor-id TUMOR \
        --vcf-normal-id NORMAL

    else
      bcftools view -f PASS {input.vcf_inter} > {output.vcf_fil};
      perl scripts/vcf2maf.pl \
        --input-vcf {output.vcf_fil} \
        --output-maf {output.maf} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 100 \
        --ref-fasta={input.ref} \
        --filter-vcf ref/VEP_cache/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=ref/98 \
        --vep-data=ref/98

    fi
    """
