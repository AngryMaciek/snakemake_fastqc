##############################################################################
#
#   Snakemake pipeline:
#   Run FastQC on RNA-Seq samples (quality analysis)
#
#   AUTHOR: Maciej_Bak
#   AFFILIATION: Swiss_Institute_of_Bioinformatics
#   CONTACT: maciej.bak@unibas.ch
#   CREATED: 21-11-2019
#   LICENSE: Apache_2.0
#
##############################################################################

# imports
import sys
import os
import pandas as pd

# local rules
localrules: all, create_output_dir, extract_adapters,

# get all fastq files
def get_all_fastq():
    design_table = pd.read_csv(config["design_file"], sep="\t", index_col=0)
    x = list(design_table["fq1"])+list(design_table["fq2"])
    x = [i.split("/")[-1].split(".")[0] for i in x if str(i)!="nan"]
    return x

# get full path for a given RNA-Seq sample
def get_full_path(sample_name):
    design_table = pd.read_csv(config["design_file"], sep="\t", index_col=0)
    for i,row in design_table.iterrows():
        if str(row["fq1"])!="nan":
            if row["fq1"].find(sample_name)>-1:
                return row["fq1"]
        if str(row["fq2"])!="nan":
            if row["fq2"].find(sample_name)>-1:
                return row["fq2"]

##############################################################################
### Target rule with final output of the pipeline
##############################################################################

rule all:
    input:
        HTML_fastqc_report = \
            expand(os.path.join("{output_dir}","{sample}_fastqc.html"), \
                output_dir=config["output_dir"], \
                sample=get_all_fastq())

##############################################################################
### Create directories for the result
##############################################################################

rule create_output_dir:
    output:
        TMP_output = temp(os.path.join("{output_dir}", "dir_created"))
    params:
        DIR_results_dir = "{output_dir}",
        DIR_cluster_log = os.path.join("{output_dir}", "cluster_log"),
    log:
        DIR_local_log = os.path.join("{output_dir}", "local_log"),
    shell:
        """
        mkdir -p {params.DIR_results_dir}; \
        mkdir -p {params.DIR_cluster_log}; \
        mkdir -p {log.DIR_local_log}; \
        touch {output.TMP_output}
        """

##############################################################################
### Extract adapter per fastq file
##############################################################################

rule extract_adapters:
    input:
        TMP_output = os.path.join("{output_dir}", "dir_created")
    output:
        DIR_adapter_dir = \
            directory(os.path.join("{output_dir}", "adapters"))
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "extract_adapters.log"),
    run:
        # extract adapter sequence per each fastq file
        # and save it to a separate text file
        os.mkdir(output.DIR_adapter_dir)
        design_table = \
            pd.read_csv(config["design_file"],sep="\t", index_col=0)
        for i,row in design_table.iterrows():
            if str(row["fq1"])!="nan" and str(row["adapter1"])!="nan":
                fname = row["fq1"].split("/")[-1].split(".")[0]+".txt"
                path_adapter = os.path.join(output.DIR_adapter_dir,fname)
                with open(path_adapter,"w") as f:
                    f.write(row["adapter1"]+"\t"+row["adapter1"]+"\n")
            if str(row["fq2"])!="nan" and str(row["adapter2"])!="":
                fname = row["fq2"].split("/")[-1].split(".")[0]+".txt"
                path_adapter = os.path.join(output.DIR_adapter_dir,fname)
                with open(path_adapter,"w") as f:
                    f.write(row["adapter2"]+"\t"+row["adapter2"]+"\n")

##############################################################################
### Reads alignment
##############################################################################

rule run_FastQC:
    input:
        DIR_adapter_dir = os.path.join("{output_dir}", "adapters"),
        STRING_fastq_path = \
            lambda wildcards: get_full_path(wildcards.sample)
    output:
        HTML_fastqc_report = \
            os.path.join("{output_dir}", "{sample}_fastqc.html")
    params:
        DIR_outdir = "{output_dir}",
        TXT_sample_adapter = \
            os.path.join("{output_dir}", "adapters", "{sample}.txt"),
        LOG_cluster_log = \
            os.path.join("{output_dir}", "cluster_log", \
                "run_FastQC_{sample}.log"),
        queue = "30min",
        time = "00:30:00"
    log:
        LOG_local_log = \
            os.path.join("{output_dir}", "local_log", \
                "run_FastQC_{sample}.log"),
    resources:
        threads = 4,
        mem = 10000
    benchmark:
        os.path.join("{output_dir}",
            "cluster_log", "run_FastQC_{sample}.benchmark.log")
    conda:
        "packages.yaml"
    shell:
        """
        fastqc {input.STRING_fastq_path} \
        --outdir {params.DIR_outdir} \
        --format fastq \
        --nogroup \
        --extract \
        --adapters {params.TXT_sample_adapter} \
        --threads {resources.threads} \
        --kmers 7 \
        &> {log.LOG_local_log}
        """
