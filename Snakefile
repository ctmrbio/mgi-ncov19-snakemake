# vim: syntax=python expandtab ts=4
from snakemake.exceptions import WorkflowError

rule all:
    input:
        "output_dir/kraken/GZMU0014-Meta-Microbiome.kraken",
        "output_dir/kraken/GZMU0014-Meta-Microbiome.kreport",
        "output_dir/kraken/GZMU0014-Meta-Microbiome.cov19.1.fq.gz",
        "output_dir/kraken/GZMU0014-Meta-Microbiome.cov19.2.fq.gz",
        "output_dir/kraken/GZMU0014-Meta-Microbiome.rest.1.fq.gz",
        "output_dir/kraken/GZMU0014-Meta-Microbiome.rest.2.fq.gz",

rule extract_cov19:
    """extract reads matching CoV_19 database using Kraken2"""
    input:
        read1="input/{sample}.1.fq.gz",
        read2="input/{sample}.2.fq.gz",
    output:
        read1="output_dir/kraken/{sample}.cov19.1.fq.gz",
        read2="output_dir/kraken/{sample}.cov19.2.fq.gz",
        rest1="output_dir/kraken/{sample}.rest.1.fq.gz",
        rest2="output_dir/kraken/{sample}.rest.2.fq.gz",
        kraken="output_dir/kraken/{sample}.kraken",
        kreport="output_dir/kraken/{sample}.kreport",
    params:
        db="/ctmr/projects/MGI_ncov-19/db_cov19_kraken2",
        classified=lambda w: f"output_dir/kraken/{w.sample}.cov19.#.fq",
        unclassified=lambda w: f"output_dir/kraken/{w.sample}.rest.#.fq",
        fq_to_compress=lambda w: f"output_dir/kraken/{w.sample}.*.fq",
    threads:
       16
    conda:
        "../envs/qc_assembly.yaml",
    log:
        stderr="output_dir/kraken/kraken2.log",
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
            2>> {log.stderr}
        """
   
rule fastp:
    """ quality trimming the cov-19 reads """
    input:
        read1="output_dir/kraken/{sample}.cov19.1.fq.gz",
        read2="output_dir/kraken/{sample}.cov19.2.fq.gz",
    output:
        read1="output_dir/fastp/{sample}.cov19.1.fq.gz",
        read2="output_dir/fastp/{sample}.cov19.2.fq.gz",
        json="output_dir/fastp/{sample}.fastp.json",
        html="output_dir/fastp/{sample}.fastp.html",
    log:
        stdout="output_dir/fastp/{sample}.stdout.log",
        stderr="output_dir/fastp/{sample}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../envs/qc_assembly.yaml"
    threads:
        4
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
        --disable_trim_poly_g -w 4\
        > {log.stdout} \
        2> {log.stderr}
        """

rule soapnuke:
    """ 
    test with raw code 
    working dir: /ctmr/projects/MGI_ncov-19/mgi_raw_codes
    error message when running the below code: check the options
    SOAPnuke filter -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG  -l 20 -q 0.2  -n 0.02 -5 0 -Q 2 -G 2 -d --fq1 output_dir/fastp/cl.reformatHeader.1.fq.gz --fq2 output_dir/fastp/cl.reformatHeader.2.fq.gz --cleanFq1 cl.1.fq.gz --cleanFq2 cl.2.fq.gz --outDir output_dir/soapnuke/ -T 4
    """
    """ furthur remove adaptor contaminant from fastq output reads"""
    input:
        read1="output_dir/fastp/{sample}.cov19.1.fq.gz",
        read2="output_dir/fastp/{sample}.cov19.2.fq.gz",
    output:
        read1="output_dir/soapnuke/{sample}.cov19.1.fq.gz",
        read2="output_dir/soapnuke/{sample}.cov19.2.fq.gz",
    log:
        stdout="output_dir/soapnuke/{sample}.stdout.log",
        stderr="output_dir/soapnuke/{sample}.stderr.log",
    shadow:
        "shallow"
    conda:
        "../envs/qc_assembly.yaml"
    threads:
        4
    shell:
        """
        SOAPnuke filter \
        -f AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA \
        -r AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG \
        -l 20 -q 0.2  -n 0.02 -5 0 -Q 2 -G 2 -d \
        --fq1 {input.read1} \
        --fq2 {input.read2}\
        --cleanFq1 {output.read1} \
        --cleanFq2 {output.read2} \
        --outDir directory(output_dir/soapnuke/) \
        -T {threads} \
        > {log.stdout} \
        2> {log.stderr} 
        """

rule prinseq:
    """
    test with raw code  
    working dir: /ctmr/projects/MGI_ncov-19/mgi_raw_codes
    
    code1 (not work with paird-end reads): prinseq-lite.pl -lc_method dust -lc_threshold 7 -fastq <(gzip -dc output_dir/fastp/cl._1.fq.gz) -fastq2 <(gzip -d    c output_dir/fastp/cl._2.fq.gz) -out_good output_dir/prinseq/cl.good -out_bad output_dir/prinseq/cl.bad
    running error: the number of beases and quality scores are not the same for sequence

    code2 (WORKED with single-end reads): gzip -dc output_dir/fastp/cl._1.fq.gz | prinseq-lite.pl -fastq stdin -out_good output_dir/prinseq/output_good.1 -lc_method dust -lc_threshold 7 | gzip > output_dir/prinseq/ output_dir/prinseq/output_good.1.fq.gz
    """
    
    """ further filtering of the low complexity reads"""
    input:
        read1="output_dir/soapnuke/{sample}.cov19.1.fq.gz",
        read2="output_dir/soapnuke/{sample}.cov19.2.fq.gz",
    output:
        read1="output_dir/prinseq/{sample}.good.1.fq.gz",
        read2="output_dir/prinseq/{sample}.good.2.fq.gz",
        fastq1=temp("output_dir/prinseq/{sapmle}.good.1.fastq"),
        fastq2=temp("output_dir/prinseq/{sapmle}.good.2.fastq"),
        good=lambda w: f"output_dir/prinseq/{w.sample}.good",
        bad=lambda w: f"output_dir/prinseq/{w.sample}.bad",
    log:
        stdout1="output_dir/prinseq/{sample}.1.stdout.log",
        stderr1="output_dir/prinseq/{sample}.1.stderr.log",
        stdout2="output_dir/prinseq/{sample}.2.stdout.log",
        stderr2="output_dir/prinseq/{sample}.2.stderr.log",
    shadow:
        "shallow"
    conda:
        "../envs/qc_assembly.yaml"
    threads:
        4
    shell:
        """
        gzip -dc {input.read1}| \
        prinseq-lite.pl \
        -fastq stdin \
        -lc_method dust -lc_threshold 7 \
        -out_good {output.fastq1}|\
        gzip > {output.read1} \
        > {log.stdout1} \
        2> {log.stderr1} 

        gzip -dc {input.read2}| \
        prinseq-lite.pl \
        -fastq stdin \
        -lc_method dust -lc_threshold 7 \
        -out_good {output.fastq2}|\
        gzip > {output.read2} \ 
        > {log.stdout2} \
        2> {log.stderr2} 
        """

rule assembly:
        """ raw code testing
        test directory "/ctmr/projects/MGI_ncov-19/mgi_raw_codes"
        error: no equal amount of input reads
        """
        input:
            read1="output_dir/prinseq/{sample}.good.1.fq.gz",
            read2="output_dir/prinseq/{sample}.good.2.fq.gz",
        output:
            directory("output_dir/assembly"),
        threads:
            20
        shell:
        """
        spades.py \
        -1 {input.read1} \
        -2 {input.read2} \
        -o {output}
        """
