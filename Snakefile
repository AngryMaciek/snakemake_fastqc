##############################################################################
#
#   Snakemake pipeline:
#   FastQC
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 21-11-2019
#   LICENSE: GPL v3.0
#
##############################################################################

# imports
import sys
import os

# local rules
localrules: create_output_dir, all

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        TXT_final_results = \
            expand(os.path.join("{output_dir}", "results.txt"),
                output_dir=config["output_dir"])

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_random_samples = os.path.join("{output_dir}", "random_samples"),
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_random_samples}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Sample some random data
##############################################################################

rule generate_files:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created"),
        SCRIPT = \
            os.path.join(config["src_dir"], "mb_random_sample.py")
    output:
        TXT_random_sample = \
            os.path.join("{output_dir}", "random_samples", "{file}")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "generate_files_{file}.log"),
        queue = "30min",
        time = "0:05:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "generate_files_{file}.log"),
    resources:
        threads = 1,
        mem = 5000
    benchmark:
        os.path.join("{output_dir}",
            "cluster_log", "generate_files_{file}_benchmark.log")
    conda:
        "packages.yaml"
    singularity:
        ""
    shell:
        """
        python {input.SCRIPT} \
        --outfile {output.TXT_random_sample} \
        &> {log.LOG_local_log};
        """

##############################################################################
### Merge the results
##############################################################################

rule merge_results:
    input:
        TXT_result_files = \
            lambda wildcards: [os.path.join(wildcards.output_dir,
                "random_samples", f) for f in config["samples_filenames"]]
    output:
        TXT_final_results = os.path.join("{output_dir}", "results.txt")
    params:
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log/merge_results.log"),
        queue = "30min",
        time = "00:05:00"
    resources:
        threads = 1,
        mem = 5000
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", "merge_results.log")
    run:
        # read all the sampled numbers:
        numbers = []
        for i in input.TXT_result_files:
            with open(i) as f:
                numbers.append(f.read().splitlines()[0])
        # save into one file:
        with open(output.TXT_final_results, "w") as outfile:
                outfile.write("\n".join(numbers))









'''
Small pipeline to get quality scores by FastQC.

Maciek Bak
'''

import os
import sys
import pandas as pd

localrules: final, create_out_dir, extract_adapters

def get_all_fastq():
    design_table = pd.read_csv(config["design_file"],sep="\t", index_col=0)
    x = list(design_table["fq1"])+list(design_table["fq2"])
    x = [i.split("/")[-1].split(".")[0] for i in x if str(i)!="nan"]
    return x

def get_full_path(sample_name):
    design_table = pd.read_csv(config["design_file"],sep="\t", index_col=0)
    for i,row in design_table.iterrows():
        if str(row["fq1"])!="nan":
            if row["fq1"].find(sample_name)>-1:
                return row["fq1"]
        if str(row["fq2"])!="nan":
            if row["fq2"].find(sample_name)>-1:
                return row["fq2"]

#################################################################################
### Target rule with final outfiles
#################################################################################

rule final:
    input:
        reports = expand("{output_dir}/{sample}_fastqc.html", output_dir=config["output_dir"], sample=get_all_fastq())

#################################################################################
### Create directories for the result and copy config immediately
#################################################################################

rule create_out_dir:
    input:
        configfile = config["this_file"]
    output:
        configfile = "{output_dir}/config.yaml"
    params:
        main_dir = "{output_dir}",
        cluster_log = "{output_dir}/cluster_log",
    log:
        local_log = "{output_dir}/local_log"
    shell:
        """
        mkdir -p {params.main_dir}; \
        mkdir -p {params.cluster_log}; \
        mkdir -p {log.local_log}; \
        cp {input.configfile} {output.configfile}
        """

#################################################################################
### Extract adapter per fastq file
#################################################################################

rule extract_adapters:
    input:
        configfile = "{output_dir}/config.yaml"
    output:
        adapter_dir = directory("{output_dir}/adapters")
    params:
        cluster_log = "{output_dir}/cluster_log/extract_adapters.log",
        queue = "30min",
        time = "00:05:00"
    log:
        local_log = "{output_dir}/local_log/extract_adapters.log",
    resources:
        threads = 1,
        mem = 5000
    run:
        os.mkdir(output.adapter_dir)
        design_table = pd.read_csv(config["design_file"],sep="\t", index_col=0)
        for i,row in design_table.iterrows():
            if str(row["fq1"])!="nan" and str(row["adapter1"])!="nan":
                fname = row["fq1"].split("/")[-1].split(".")[0]+".txt"
                with open(os.path.join(output.adapter_dir,fname),"w") as f:
                    f.write(row["adapter1"]+"\t"+row["adapter1"]+"\n")
            if str(row["fq2"])!="nan" and str(row["adapter2"])!="":
                fname = row["fq2"].split("/")[-1].split(".")[0]+".txt"
                with open(os.path.join(output.adapter_dir,fname),"w") as f:
                    f.write(row["adapter2"]+"\t"+row["adapter2"]+"\n")

#################################################################################
### Reads alignment
#################################################################################

rule run_FastQC:
    input:
        adapter_dir = "{output_dir}/adapters",
        fastq_path = lambda wildcards: get_full_path(wildcards.sample)
    output:
        fastqc_report = "{output_dir}/{sample}_fastqc.html"
    params:
        outdir = "{output_dir}",
        sample_adapter = os.path.join("{output_dir}","adapters","{sample}.txt"),
        cluster_log = "{output_dir}/cluster_log/run_FastQC_{sample}.log",
        queue = "30min",
        time = "00:30:00"
    log:
        local_log = "{output_dir}/local_log/run_FastQC_{sample}.log",
    resources:
        threads = 4,
        mem = 10000
    benchmark:
        "{output_dir}/cluster_log/run_FastQC_{sample}.benchmark.log"
    conda:
        "yaml_files/fastqc.yaml"
    shell:
        """
        fastqc {input.fastq_path} \
        --outdir {params.outdir} \
        --format fastq \
        --nogroup \
        --extract \
        --adapters {params.sample_adapter} \
        --threads {resources.threads} \
        --kmers 7 \
        &> {log.local_log}
        """