1. For S2788/89 (Arabidopsis): featureCounts_GENE failed
work dir: /home/kai.hu-umw/pi/hira.goel-umw/OneStopRNAseq/users/318/At_GSE169302_07.48.13-04.25.2025/analysis_1

    output: feature_count_gene_level/counts.s2.strict.txt, feature_count_gene_level/counts.s2.strict.txt.summary
    log: feature_count_gene_level/log/counts.s2.strict.log (check log file(s) for error message)
    conda-env: /pi/mccb-umw/shared/OneStopRNAseq/conda/osr_envs/e47ec768116a51673cc018986be20d0f_
    shell:

        # gene level count
        featureCounts -a /pi/mccb-umw/shared/OneStopRNAseq/genome/arabidopsis_thaliana/tair10/Arabidopsis_thaliana.TAIR10.57.gtf -o feature_count_gene_level/counts.s2.strict.txt -T 4         -t gene -g gene_id -s 2         --minOverlap 1 --fracOverlap 0           -Q 20  mapped_reads/GSE169302_ColOmock1.bam mapped_reads/GSE169302_ColOmock2.bam mapped_reads/GSE169302_ColOmock3.bam mapped_reads/GSE169302_ColOIAA1.bam mapped_reads/GSE169302_ColOIAA2.bam mapped_reads/GSE169302_ColOIAA3.bam         > feature_count_gene_level/log/counts.s2.strict.log 2>&1

        # get larger count (exon vs gene) for each gene
        python workflow/script/get_max_of_featureCountsTable.py         feature_count/counts.s2.strict.txt feature_count_gene_level/counts.s2.strict.txt          >> feature_count_gene_level/log/counts.s2.strict.log 2>&1

        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

conda activate /pi/mccb-umw/shared/OneStopRNAseq/conda/osr_envs/e47ec768116a51673cc018986be20d0f_



        """
        featureCounts -a {input.gtf} -o {output.COUNT} \
        -T {threads} -g gene_id -s {wildcards.strand} \
        {params.pe}  {params.mode} \
        {input.bams} > {log} 2>&1
        """

        """
        # gene level count
        featureCounts -a {input.gtf} -o {output.ct} -T {threads} \
        -t gene -g gene_id -s {wildcards.strand} \
        --minOverlap 1 --fracOverlap 0 \
        {params.pe} {params.mode} {input.bams} \
        > {log} 2>&1

        # get larger count (exon vs gene) for each gene
        python workflow/script/get_max_of_featureCountsTable.py \
        {input.exon_counts} {output.ct}  \
        >> {log} 2>&1
        """