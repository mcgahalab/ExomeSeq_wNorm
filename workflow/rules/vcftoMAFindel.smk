rule vcftoMAFindel:
  input:
    ref = 'ref/genome.fa',
    vcf_inter = "results/vcfIntersect/indels/{sample}/{indel}.vcf",
  params:
    samp="{sample}",
    indel = "{indel}",
    control=has_a_control,
  output:
    vcf_fil = "results/vcfIntersect/indels/{sample}/fil_{indel}.vcf",
    maf = "results/MAF_38_final/indel/{sample}/{indel}.maf",
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/VCFtoMAF.yaml",
  shell:
    """
    if [ '{params.control}' == 'True' ]; then
        if [ {params.indel} != '0001' ]; then
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
            --vep-path=ref/98 \
            --vep-data=ref/98 \
            --ncbi-build GRCh38
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
                --vep-path=ref/98 \
                --vep-data=ref/98 \
                --ncbi-build GRCh38 \
                --tumor-id={params.samp} \
                --normal-id unmatched \
                --vcf-tumor-id TUMOR \
                --vcf-normal-id NORMAL
        fi
    else
        mkdir -p results/MAF_38_final/indel/{params.samp}
        touch results/vcfIntersect/indels/{params.samp}/fil_{params.indel}.vcf
        touch results/MAF_38_final/indel/{params.samp}/{params.indel}.maf
    fi
    """
