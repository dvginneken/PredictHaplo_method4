import subprocess

configfile: "config.json"

SAMPLES, = glob_wildcards(config['haplohiv_folder']+"ec/{id}.fastq")
REPLICATES, = glob_wildcards(config['output_dir']+"/consensus/consensusses_{id}.txt")

rule all:
    input:          
        expand("{prefix}/haplotypes/{samples}_haplotypes.fa", prefix=config['output_dir'], samples=SAMPLES)


rule concatenate:
    input:
        "{prefix}/consensus/consensusses_{name}.txt"
    output:
        "{prefix}/consensus/{name}_concat.fa"
    shell:
        """
        files=`cat {input}`
        cat $files > {output}
        """ 

rule align:
    input:
        "{prefix}/consensus/{name}_concat.fa"
    output:
        "{prefix}/consensus/{name}_aligned.fa"
    shell:
        "mafft {input} > {output}"

rule consensus:
    input:
        "{prefix}/consensus/{name}_aligned.fa"
    output:
        "{prefix}/consensus/{name}_consensus.fa"
    shell:
        """
        if [ `grep -c ">" {input}` -gt 1 ]
        then
            cons -sequence {input} -outseq {output} -name {wildcards.name}
        else
            cat {input} > {output}
            sed -i "1s/.*/>{wildcards.name}/" {output}
        fi
        bwa index {output}
        """

rule map_to_consensus:
    input:
        reads = config['haplohiv_folder']+"ec/{samples}.fastq",
        consensus = expand("{prefix}/consensus/{name}_consensus.fa", prefix=config['output_dir'], name=REPLICATES)
    output:
        "{prefix}/mapped/{samples}.sam"
    shell:
        """
        name=`grep {wildcards.samples} {config[table]} | awk '{{print $3}}'`
        bwa mem -p {config[output_dir]}/consensus/${{name::-1}}_consensus.fa {input.reads} > {output}
        """

rule configure_predicthaplo:
    input:
        dummy = "PredictHaplo.conf",
        sample = "{prefix}/mapped/{samples}.sam"
    output:
        "{prefix}/predicthaplo/{samples}.conf"
    shell:
        """
        cp {input.dummy} {output}
        sed -i "s|prefix_line|{config[output_dir]}/predicthaplo/output_{wildcards.samples}/{wildcards.samples}|" {output}
        name=`grep {wildcards.samples} {config[table]} | awk '{{print $3}}'`
        sed -i "s|reference_line|{config[output_dir]}/consensus/${{name::-1}}_consensus.fa|" {output}
        sed -i "s|reads_line|{input.sample}|" {output}
        sed -i "s|true_haplotypes_line|{config[output_dir]}/consensus/${{name::-1}}_aligned.fa|" {output}
        """

rule predicthaplo:
    input:
        "{prefix}/predicthaplo/{samples}.conf"
    output:
        "{prefix}/predicthaplo/output_{samples}/{samples}.conf"
    shell:
        """
        PredictHaplo-Paired {input}
        mv {input} {output}
        """

rule gather_haplotype:
    input:
        "{prefix}/predicthaplo/output_{samples}/{samples}.conf"
    output:
        "{prefix}/haplotypes/{samples}_haplotypes.fa"
    shell:
        """
        mkdir -p {config[output_dir]}/haplotypes
        python scripts/gather_haplotypes.py {config[output_dir]}/predicthaplo/output_{wildcards.samples}/{wildcards.samples}global_*.fas \
        {output} {wildcards.samples}
        """


