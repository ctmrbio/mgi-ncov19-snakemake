# vim: syntax=python expandtab ts=4
from snakemake.exceptions import WorkflowError
SAMPLES = glob_wildcards("input/{sample}_1.fq.gz").sample
print(SAMPLES)

rule all:
    input:
        expand("output_dir/assembly/{sample}/contigs.fasta", sample=SAMPLES),
        expand("output_dir/bwa/Ref_nCoV_{sample}.sam", sample=SAMPLES),
        #"output_dir/samtools/Ref_nCoV_GZMU0014-Meta-Microbiome-smd.matrics"

rule extract_cov19:
    """extract reads matching CoV_19 database using Kraken2"""
    input:
        read1="input/{sample}_1.fq.gz",
        read2="input/{sample}_2.fq.gz",
    output:
        read1="output_dir/kraken/{sample}.cov19_1.fq.gz",
        read2="output_dir/kraken/{sample}.cov19_2.fq.gz",
        rest1="output_dir/kraken/{sample}.rest_1.fq.gz",
        rest2="output_dir/kraken/{sample}.rest_2.fq.gz",
        kraken="output_dir/kraken/{sample}.kraken",
        kreport="output_dir/kraken/{sample}.kreport",
    params:
        db="database_cov19",
        classified=lambda w: f"output_dir/kraken/{w.sample}.cov19#.fq",
        unclassified=lambda w: f"output_dir/kraken/{w.sample}.rest#.fq",
        fq_to_compress=lambda w: f"output_dir/kraken/{w.sample}*.fq",
    threads:
        8 
    conda:
        "envs/qc_assembly.yaml",
    log:
        stderr="output_dir/kraken/{sample}.kraken2.log",
    shadow:
        "shallow"
    shell:
        """
        kraken2 \
            --db {params.db} \
            --threads {threads} \
            --output {output.kraken}\
            --report {output.kreport}\
            --classified-out {params.classified} \
            --unclassified-out {params.unclassified}\
            --paired \
            {input.read1} {input.read2} 
            2> {log.stderr}
        pigz \
            --processes {threads} \
            --verbose \
            --force \
            {params.fq_to_compress} \
            2> {log.stderr}
        """
   
rule fastp:
    """ quality trimming the cov-19 reads """
    input:
        read1="output_dir/kraken/{sample}.cov19_1.fq.gz",
        read2="output_dir/kraken/{sample}.cov19_2.fq.gz",
    output:
        read1="output_dir/fastp/{sample}.cov19_1.fq.gz",
        read2="output_dir/fastp/{sample}.cov19_2.fq.gz",
        json="output_dir/fastp/{sample}.fastp.json",
        html="output_dir/fastp/{sample}.fastp.html",
    log:
        stdout="output_dir/fastp/{sample}.stdout.log",
        stderr="output_dir/fastp/{sample}.stderr.log",
    shadow:
        "shallow"
    conda:
        "envs/qc_assembly.yaml"
    threads:
        8 
    shell:
        """
        fastp \
        -q 20 -u 20 -n 1 -l 50 \
        -i {input.read1} \
        -I {input.read2}\
        -o {output.read1}\
        -O {output.read2}\
        -j {output.json}\
        -h {output.html}\
        --adapter_sequence AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        --adapter_sequence_r2 AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        --detect_adapter_for_pe \
        --disable_trim_poly_g \
        --thread {threads}\
        --low_complexity_filter \
        --complexity_threshold 7 \
        > {log.stdout} \
        2> {log.stderr}
        """

rule soapnuke:
    """ furthur remove adaptor contaminant from fastq output reads"""
    """ github:https://github.com/BGI-flexlab/SOAPnuke"""
    input:
        read1="output_dir/fastp/{sample}.cov19_1.fq.gz",
        read2="output_dir/fastp/{sample}.cov19_2.fq.gz",
    output:
        read1="output_dir/soapnuke/{sample}/{sample}.cov19_1.fq.gz",
        read2="output_dir/soapnuke/{sample}/{sample}.cov19_2.fq.gz",
    log:
        stdout="output_dir/soapnuke/{sample}/{sample}.stdout.log",
        stderr="output_dir/soapnuke/{sample}/{sample}.stderr.log",
    params:
        outdir=lambda w: f"output_dir/soapnuke/{w.sample}",
        clean1=lambda w: f"{w.sample}.cov19_1.fq.gz",
        clean2=lambda w: f"{w.sample}.cov19_2.fq.gz",
    shadow:
        "shallow"
    conda:
        "envs/qc_assembly.yaml"
    threads:
        8
    shell:
        """
        SOAPnuke filter \
        -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        -l 20 -q 0.2 -n 0.02 -5 0 -Q 2 -G 2 \
        --fq1 {input.read1} \
        --fq2 {input.read2}\
        --cleanFq1 {params.clean1} \
        --cleanFq2 {params.clean2} \
        --outDir {params.outdir} \
        -T {threads} \
        > {log.stdout} \
        2> {log.stderr} 
        """

rule assembly:
        """ de novo assembly """
        input:
            read1="output_dir/soapnuke/{sample}/{sample}.cov19_1.fq.gz",
            read2="output_dir/soapnuke/{sample}/{sample}.cov19_2.fq.gz",
        output:
            contigs="output_dir/assembly/{sample}/contigs.fasta"
        params:
            outdir=lambda w: f"output_dir/assembly/{w.sample}"
        conda:
            "envs/qc_assembly.yaml"
        shadow:
            "shallow"
        log:
            stdout="output_dir/soapnuke/{sample}.stdout.log",
            stderr="output_dir/soapnuke/{sample}.stderr.log",
        threads:
            8
        shell:
            """
            spades.py \
            -1 {input.read1} \
            -2 {input.read2} \
            -o {params.outdir} \
            -t {threads} \
            > {log.stdout} \
            2> {log.stderr}
            """

rule coverage_bwa_sai:
        input:
            read1="output_dir/soapnuke/{sample}/{sample}.cov19_1.fq.gz",
            read2="output_dir/soapnuke/{sample}/{sample}.cov19_2.fq.gz",
        output:
            bwa1="output_dir/bwa/Ref_nCov_{sample}_1.sai",
            bwa2="output_dir/bwa/Ref_nCov_{sample}_2.sai",
        params:
            db="database_cov19/HCoV-19"
        conda:
            "envs/coverage_mapping.yaml"
        shadow:
            "shallow"
        log:
            stderr1="output_dir/bwa/Ref_nCov_{sample}_1.sai.stderr.log",
            stderr2="output_dir/bwa/Ref_nCov_{sample}_2.sai.stderr.log",
        shell:
            """
            bwa aln -t 4 \
            {params.db} {input.read1} > {output.bwa1}\
            2>{log.stderr1}
            
            bwa aln -t 4 \
            {params.db} {input.read2} > {output.bwa2}\
            2>{log.stderr2}
            """

rule coverage_bwa_sam:
        input:
            read1="output_dir/soapnuke/{sample}/{sample}.cov19_1.fq.gz",
            read2="output_dir/soapnuke/{sample}/{sample}.cov19_2.fq.gz",
            bwa1="output_dir/bwa/Ref_nCov_{sample}_1.sai",
            bwa2="output_dir/bwa/Ref_nCov_{sample}_2.sai",
        output:
            sam="output_dir/bwa/Ref_nCoV_{sample}.ori.sam",
        params:
            db="database_cov19/HCoV-19"
        conda:
            "envs/coverage_mapping.yaml"
        shadow:
            "shallow"
        log:
            stderr="output_dir/bwa/Ref_nCoV_{sample}.ori.sam.stderr.log"
        shell:
            """
            bwa sampe {params.db} \
            {input.bwa1} {input.bwa2} {input.read1} {input.read2}\
            >{output.sam} \
            2>{log.stderr}
            """

rule coverage_bwa_filter:
        """filter the sam file, coverage 0.95, identity 0.90"""
        input:
            "output_dir/bwa/Ref_nCoV_{sample}.ori.sam"
        output:
            "output_dir/bwa/Ref_nCoV_{sample}.sam"
        shadow:
            "shallow"
        log:
            stdout="output_dir/bwa/Ref_nCoV_{sample}.sam.FLT.stdout.log",
            stderr="output_dir/bwa/Ref_nCoV_{sample}.sam.FLT.stderr.log",
        shell:
            """
            perl scripts/BWA_sam_Filter_identity_cvg.pl \
            -i {input} \
            -o {output}\
            -m 0.95 \
            -s 0.90 \
            > {log.stdout}\
            2> {log.stderr}
            """

rule coverage_samtools:
        """
        https://github.com/BGI-IORI/nCoV_Meta/blob/master/nCoV_Finder.pl        
        the {output.Lxxx} indicate the line number of the sbove script where the output file were specified
        """
        input:
            "output_dir/bwa/Ref_nCoV_{sample}.sam",
        output:
            L178="output_dir/samtools/Ref_nCoV_{sample}-uF.bam",
            L179="output_dir/samtools/Ref_nCoV_{sample}.bam",
            L180="output_dir/samtools/Ref_nCoV_{sample}.bam.flagstat",
            L181="output_dir/samtools/Ref_nCoV_{sample}-s.bam",
            L183a="output_dir/samtools/Ref_nCoV_{sample}-smd.bam",
            L183b="output_dir/samtools/Ref_nCoV_{sample}-smd.matrics",
        params: 
            db="envs/database_cov19/HCoV-19.fai"
        shadow:
            "shallow"
        log:
            stdout="output_dir/samtools/Ref_nCoV_{sample}-smd.matrics.stdout.log",
            stderr="output_dir/samtools/Ref_nCoV_{sample}-smd.matrics.stderr.log"
        shell:
            """
            samtools view -bt $DB.fai {input} > {output.L178} 
         
            samtools sort -n {output.L178} | samtools fixmate: - {output.L179}
         
            samtools flagstat {output.L179} > {output.L180}
         
            samtools sort {output.L179} -o {output.L181} --reference {params.db}
         
            samtools index {output.L181} 
         
            java -jar bin/picard.jar \
            MarkDuplicates AS=TRUE \
            VALIDATION_STRINGENCY=LENIENT \
            MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 \
            REMOVE_DUPLICATES=TRUE INPUT={output.L181} \
            OUTPUT={output.L183a} \
            METRICS_FILE={output.L183b} \
            > {log.stdout} \
            2> {log.stderr}
            """ 
