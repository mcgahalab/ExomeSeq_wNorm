rule Sequenza:
    input:
        snp="results/Varscan/snv/{sample}/{sample}.snp",
        copynum="results/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
    params:
        script="scripts/SequenzaSingleSample_v2.1_hg38.R",
    output:
        dir=directory("results/Sequenza/{sample}"),
        complete="results/Sequenza/{sample}_complete"
    threads: 4
    shell:
     """
     module load R/3.3.0
     Rscript {params.script} -s {input.snp} -c {input.copynum} -o {output.dir} -f 10
     touch {output.complete}
     """
