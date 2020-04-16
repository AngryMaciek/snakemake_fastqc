# Snakemake pipeline for FastQC
*Maciej_Bak  
Swiss_Institute_of_Bioinformatics*

[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is a nice tool to inspect the quality of RNA-Seq samples.  
This repository is a very small snakemake workflow that I use for automated and reproducible quality analyses of sequencing samples in my reseach.

## Snakemake pipeline execution
Snakemake is a workflow management system that helps to create and execute data processing pipelines. It requires Python 3 and can be most easily installed via the bioconda package from the anaconda cloud service.

### Step 1: Download and installation of Miniconda3
Linux:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
  bash Miniconda3-latest-Linux-x86_64.sh
  source .bashrc
  ```

macOS:
  ```bash
  wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
  bash Miniconda3-latest-MacOSX-x86_64.sh
  source .bashrc
  ```

### Step 2: Pandas and Snakemake installation

To execute the workflow one would require pandas python library and snakemake workflow menager.  
Unless a  specific snakemake version is specified explicitly it is most likely the best choice to install the latest versions:
  ```bash
  conda install -c conda-forge pandas
  conda install -c bioconda snakemake
  ```

In case you are missing some dependancy packages please install them first (with `conda install ...` as well).

### Step 3: Pipeline execution
Specify all the required information (input/output/parameters) in the config.yaml.  
The main input to the pipeline is the design table which has to have the following format:

sample  fq1 fq2 adapter1  adapter2  
[sample_name] [path_to_fq1] [path_to_fq2] [adapter1_sequence] [adapter2_sequence]  
[sample_name] [path_to_fq1] [path_to_fq2] [adapter1_sequence] [adapter2_sequence]  
...

Where:  
* Each row is a sequencing sample.
* All the fastq files need to have a different name, regardless of their location/directory.
* fq1 stands for forward-strand reads, fq2 for reverse-strand reads and adapter1/adapter2 stand for matching adapter sequences.
* In case we have single-end sequencing data please leave fq2 and adapter2 columns with empty strings.
* Design table might contain more columns than these specified above.

Once the metadata are ready write a DAG (directed acyclic graph) into dag.pdf:
  ```bash
  bash snakemake_dag_run.sh
  ```

There are two scripts to start the pipeline, depending on whether you want to run locally or on a SLURM computational cluster. In order to execute FastQC snakemake automatically creates an internal conda virtual environment and installs the tool from anaconda cloud service. For the cluster execution it might be required to adapt the 'cluster_config.json' and submission scripts before starting the run.
  ```bash
  bash snakemake_local_run_conda_env.sh
  bash snakemake_cluster_run_conda_env.sh
  ```

## License

Apache 2.0
