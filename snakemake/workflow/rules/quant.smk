"""
Gene Expression Quantification:
- featureCounts
- SalmonTE
- Salmon
"""


rule featureCounts_EXON:
    '''
    exon read count
    recognizes config.yaml: 
        - MODE: strict, liberal
        - config['PAIR_END']: True, False
    '''
    input:
        bams=mapped_bam_inputs(config,SAMPLES),
        gtf=config["GTF"]
    output:
        COUNT="feature_count/counts.s{strand}.{mode}.txt",
        summary="feature_count/counts.s{strand}.{mode}.txt.summary"
    params:
        pe='-p' if config['PAIR_END'] else ' ',
        mode='-Q 20 --minOverlap 1 --fracOverlap 0 -B -C' \
            if config['MODE'] == 'strict' else \
            '-M -Q 0 --primary --minOverlap 1 --fracOverlap 0'
    conda:
        "../envs/subread.yaml"
    priority:
        100
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    threads:
        4
    log:
        "feature_count/log/counts.s{strand}.{mode}.txt.log"
    benchmark:
        "feature_count/log/counts.s{strand}.{mode}.txt.benchmark"
    shell:
        """
        featureCounts -a {input.gtf} -o {output.COUNT} \
        -T {threads} -g gene_id -s {wildcards.strand} \
        {params.pe}  {params.mode} \
        {input.bams} > {log} 2>&1
        """
        

rule featureCounts_EXON_multiqc:
    input:
        "feature_count/counts.s{strand}.{mode}.txt.summary"
    output:
        "feature_count/counts.s{strand}.{mode}.txt.summary.multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "feature_count/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.log"
    benchmark:
        "feature_count/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.benchmark"
    shell:
        """
        multiqc {input} -f -n {output} > {log} 2>&1;
        """


rule featureCounts_GENE:
    '''
    gene read count, including exon and intron reads
    the larger one of gene-count/exon-count for each gene, is used as gene-count
    recognizes config.yaml: 
        - config['MODE']: strict, liberal
        - config['PAIR_END']: True, False
    '''
    input:
        bams=mapped_bam_inputs(config,SAMPLES),
        exon_counts="feature_count/counts.s{strand}.{mode}.txt",
        gtf=config["GTF"]
    output:
        ct='feature_count_gene_level/counts.s{strand}.{mode}.txt',
        summary='feature_count_gene_level/counts.s{strand}.{mode}.txt.summary'
    params:
        pe='-p -B -C ' if config['PAIR_END'] else ' ',
        mode='-Q 20 ' if config['MODE'] == 'strict' else '-M --primary -Q 0'
    conda:
        "../envs/subread.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4000,
    threads:
        4
    log:
        'feature_count_gene_level/log/counts.s{strand}.{mode}.log'
    benchmark:
        'feature_count_gene_level/log/counts.s{strand}.{mode}.benchmark'
    shell:
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


rule featureCounts_GENE_multiqc:
    input:
        'feature_count_gene_level/counts.s{strand}.{mode}.txt.summary'
    output:
        "feature_count_gene_level/counts.s{strand}.{mode}.txt.summary.multiqc_report.html"
    conda:
        "../envs/fastqc.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    threads:
        1
    log:
        "feature_count_gene_level/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.log"
    benchmark:
        "feature_count_gene_level/log/counts.s{strand}.{mode}.txt.summary.multiqc_report.html.benchmark"
    shell:
        """
        multiqc {input} -f -n {output} > {log} 2>&1;
        """


rule Strand_Detection:
    input:
        strand_detection_input(config)
    output:
        "feature_count/counts.strict.txt" if config['MODE'] == 'strict' else "feature_count/counts.liberal.txt",
        "meta/strandness.detected.txt"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    log:
        "meta/log/strandness.detected.txt.log"
    script:
        "../script/strandness_detection.py"


rule Strand_Detection_GeneLevel:
    input:
        ["feature_count_gene_level/counts.s{}.{}.txt.summary". \
             format(strand,config['MODE']) for strand in config['STRAND']]
    output:
        "feature_count_gene_level/counts.{}.txt".format(config['MODE']),
        "feature_count_gene_level/strandness.detected.txt"
    threads:
        1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000
    log:
        "feature_count_gene_level/log/strandness.detected.txt.log"
    script:
        "../script/strandness_detection.ge.py"


# rule SalmonTE_prep
if config['PAIR_END']:
    rule SalmonTE_Prep:
        input:
            r1=expand("trimmed/{sample}.R1.fastq.gz",sample=SAMPLES),
            r2=expand("trimmed/{sample}.R2.fastq.gz",sample=SAMPLES)
        output:
            r1=temp(expand("fastq_salmon/{sample}_1.fastq.gz",sample=SAMPLES)),
            r2=temp(expand("fastq_salmon/{sample}_2.fastq.gz",sample=SAMPLES))
        log:
            "fastq_salmon/SalmonTE_prep.log"
        resources:
            mem_mb=lambda wildcards, attempt: attempt * 1000
        threads:
            1
        run:
            for i, f in enumerate(input):
                os.system("ln -s " + "../" + input[i] + " " + output[i])

rule SalmonTE:
    # only start from Fastq
    input:
        reads=expand("fastq_salmon/{sample}_{n}.fastq.gz",n=[1, 2],sample=SAMPLES) \
            if config['PAIR_END'] else \
            expand("trimmed/{sample}.fastq.gz",sample=SAMPLES),
        raw_reads=expand("trimmed/{sample}.R1.fastq.gz",sample=SAMPLES) \
            if config['PAIR_END'] else \
            expand("trimmed/{sample}.fastq.gz",sample=SAMPLES),
        raw_reads2=expand("trimmed/{sample}.R2.fastq.gz",sample=SAMPLES) \
            if config['PAIR_END'] else \
            expand("trimmed/{sample}.fastq.gz",sample=SAMPLES),
        meta="meta/meta.csv",
        contrast="meta/contrast.de.csv"
    output:
        "SalmonTE_output/EXPR.csv"
    conda:
        "../envs/deseq2_salmonte.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1000,
    params:
        ref=config['TE_REFERENCE'],
        fdr=config['MAX_FDR'],
        lfc=config['MIN_LFC'],
        independentFilter=config["independentFilter"],
        cooksCutoff=config["cooksCutoff"],
        blackSamples=config['blackSamples'] if 'blackSamples' in config else "",
        anno_tab=config['ANNO_TAB']
    threads:
        16
    log:
        "log/SalmonTE_output/SalmonTE.log"
    benchmark:
        "log/SalmonTE_output/SalmonTE.benchmark"
    shell:
        """
        rm -rf SalmonTE_output/
        python workflow/envs/SalmonTE/SalmonTE.py --version >> {log}

        # if custom repeat library is provided, use it
        if [ -f custom_repeat_lib.fa ]; then
            python workflow/envs/SalmonTE/SalmonTE.py index \
                --ref_name=custom \
                --input_fasta=custom_repeat_lib.fa > {log} 2>&1
            python workflow/envs/SalmonTE/SalmonTE.py quant \
                --reference=custom --exprtype=count \
                --num_threads={threads} \
                {input.reads} > {log} 2>&1
        else
            # if custom repeat library is not provided, use default
            python workflow/envs/SalmonTE/SalmonTE.py quant \
                --reference={params.ref} --exprtype=count \
                --num_threads={threads} \
                {input.reads} > {log} 2>&1
        fi

        # Perform statistical test using SalmonTE built-in test, this does not perform batch correction
        num_cols=$(awk -F, 'NR==1 {{print NF}}' {input.contrast})
        mv SalmonTE_output/condition.csv SalmonTE_output/original_condition.csv
        mv SalmonTE_output/EXPR.csv SalmonTE_output/original_EXPR.csv

        for ((i=1; i<=${{num_cols}}; i++)); do
            python workflow/envs/SalmonTE/scripts/prepare_condition.py \
                {input.meta} \
                {input.contrast} \
                ${{i}} \
                SalmonTE_output/original_condition.csv \
                SalmonTE_output/original_EXPR.csv \
                SalmonTE_output/condition.csv \
                SalmonTE_output/EXPR.csv

            python workflow/envs/SalmonTE/SalmonTE.py test --inpath=SalmonTE_output --outpath=SalmonTE_output/DET_contrast_${{i}} --conditions=control,treatment
            mv SalmonTE_output/condition.csv SalmonTE_output/DET_contrast_${{i}}
            mv SalmonTE_output/EXPR.csv SalmonTE_output/DET_contrast_${{i}}
        done

        mv SalmonTE_output/original_condition.csv SalmonTE_output/condition.csv
        mv SalmonTE_output/original_EXPR.csv SalmonTE_output/EXPR.csv

        # Perform statistical test using DESeq2, this performs batch correction
        mkdir SalmonTE_output/DET_batch_corrected
        Rscript workflow/envs/SalmonTE/scripts/DESeq2.R \
            {input.meta} \
            {input.contrast} \
            SalmonTE_output/EXPR.csv \
            SalmonTE_output/DET_batch_corrected \
            {params.fdr} \
            {params.lfc} \
            {params.independentFilter} \
            {params.cooksCutoff} \
            {params.blackSamples} \
            {params.anno_tab} > SalmonTE_output/DET_batch_corrected/DESeq2_log.txt 2>&1
        """

rule Merge_TE_and_GE:
    input:
        gene=Merge_TE_and_Gene_input(config),
        te="SalmonTE_output/EXPR.csv"
    output:
        "feature_count_gene_level/TE_included.txt" \
            if config['INTRON'] \
            else "feature_count/TE_included.txt"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
    threads:
        1
    log:
        "feature_count_gene_level/log/TE_included.txt.log" \
            if config['INTRON'] \
            else "feature_count/log/TE_included.txt.log"
    benchmark:
        "feature_count_gene_level/log/TE_included.txt.benchmark" \
            if config['INTRON'] \
            else "feature_count/log/TE_included.txt.benchmark"
    shell:
        """
        python workflow/script/merge_featureCount_and_SalmonTE.py \
        {input.gene} {input.te} {output} > {log} 2>&1;
        """

# Below are added by Kai
# I would like to add a rule that quantify the reads using Salmon if the genome is human or mouse
rule salmon:
    input:
        reads=["trimmed/{sample}.R1.fastq.gz", "trimmed/{sample}.R2.fastq.gz"] \
            if config["PAIR_END"] else \
            "trimmed/{sample}.fastq.gz"
    output:
        quant="feature_count/salmon/{sample}/quant.sf"
    params:
        salmon_index = lambda wildcards: config["SALMON_INDEX"],
        libtype = "A",
        outdir = lambda wildcards, output: "Salmon_output/{sample}".format(
            sample=wildcards.sample,
        ),
        input_args = lambda wildcards, input: (
            f"-1 {input.reads[0]} -2 {input.reads[1]}" if config["PAIR_END"] 
            else f"-r {input.reads[0]}"
        )
    conda:
        "../envs/salmon.yaml"
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 16000,
    threads:
        4
    log:
        "Salmon_output/{sample}/log/quant.log"
    benchmark:
        "Salmon_output/{sample}/log/quant.benchmark"
    shell:
        """
        salmon quant \
            -i {params.salmon_index} \
            -l {params.libtype} \
            {params.input_args} \
            --validateMappings \
            -p {threads} \
            -o {params.outdir} \
            > {log} 2>&1
        """