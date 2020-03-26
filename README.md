# exSEEK

exSEEK is an integrated computational framework to discover and evaluate exRNA biomarkers for liquid biopsy.

The exSEEK framework consists of:
+ Pre_processing:
   
   + Building index with various types of genomes and annotations.
   + Quality control and removing adaptors. [`exseek.py quality_control`] [`exseek.py cutadapt`]
   + Sequential mapping for small/long RNA-seq. [`exseek.py mapping`]
   + Counting expression matrix. [`exseek.py count_matrix`]
+ Main function:
   
   + Peak calling for fragments of long RNAs. [`exseek.py call_domains`]
   + Normalization and batch removal. [`exseek.py normalization`]
   + Feature selection and machine learning. [`exseek.py feature_selection`]
   + Evaluation for selected biomarker panels. [`exseek.py evaluate_features`]

Table of Contents:

* [Installation](#istallation)
* [Usage](#Usage)
  * [Index preparing](#Index_preparing)
  * [Small RNA-seq mapping](#Small_RNA-seq_mapping)
  * [Long RNA-seq mapping](#Long_RNA-seq_mapping)
  * [Peak calling for fragments](#Peak_calling_for_fragments)
  * [Normalization and batch removal](Normalization_and_batch_removal)
  * [Feature selection](Feature_selection)
  * [Evaluation for biomarkers](Evaluation_for_biomarkers)
* [Copyright and License Information](#copyright-and-license-information)
* [Citation](#citation)

---


## Installation

For easy installation, you can use the [exSEEK image](https://hub.docker.com/r/ltbyshi/exseek) of [docker](https://www.docker.com) with all dependencies installed:

```bash
docker pull ltbyshi/exseek
```
All required software and packages are already installed in docker, so there is no more requirements. 

To test the installation and get information about the command-line interface of exSEEK, you can execute:

```bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py -h
```
A helper message is shown:

```bash
usage: exseek.py [-h] --dataset DATASET [--config-dir CONFIG_DIR] [--cluster]
                 [--cluster-config CLUSTER_CONFIG]
                 [--cluster-command CLUSTER_COMMAND] [--singularity]  
                 {quality_control,quality_control_clean,cutadapt,rename_fastq,fastq_to_fasta,prepare_genome,bigwig,mapping,
                 filter,count_matrix,call_domains,merge_domains,combine_domains,normalization,feature_selection,
                 differential_expression,evaluate_features,igv,update_sequential_mapping,update_singularity_wrappers}

exseek main program

positional arguments:
  {quality_control,quality_control_clean,cutadapt,rename_fastq,fastq_to_fasta,prepare_genome,bigwig,mapping,
  filter,count_matrix,call_domains,merge_domains,combine_domains,normalization,feature_selection,
  differential_expression,evaluate_features,igv,update_sequential_mapping,update_singularity_wrappers}

optional arguments:
  -h, --help                                    show this help message and exit
  --dataset DATASET, -d DATASET                 dataset name
  --config-dir CONFIG_DIR, -c CONFIG_DIR        directory for configuration files
  --cluster                                     submit to cluster
  --cluster-config CLUSTER_CONFIG               cluster configuration file ({config_dir}/cluster.yaml by default)
  --cluster-command CLUSTER_COMMAND             command for submitting job to cluster (default read from
                                                {config_dir}/cluster_command.txt
  --singularity                                 use singularity
```

The exSEEK directory was cloned to /apps/exseek in the docker.
You can create a bash script named `exseek` and set the script executable:

```bash
#! /bin/bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py "$@"
```
After adding the file to one of the directory in the `$PATH` variable, you can simply run: `exseek`.

### Input files

An example of input files can be found in `example_data` directory with the following structure:
```
example_data/
├── config
|   ├── example.yaml
│   └── default_config.yaml
├── data
│   └── example
|       ├── fastq
│       ├── batch_info.txt
│       ├── compare_groups.yaml
│       ├── sample_classes.txt
│       └── sample_ids.txt
└── output
    └── example
        └── ...
```

> **Note:**
> * `config/example.yaml`: configuration file with frequently changed parameters, such as file paths and mapping parameters.
> * `config/default_config.yaml`: configuration file with additional detailed parameters for each step.
> * `data/example/batch_info.txt`: table of batch information.
> * `data/example/compare_groups.yaml`: configuration file for definition of positive and negative samples.
> * `data/example/sample_classes.txt`: table of sample labels.
> * `output/example/`: input matrix of read counts

## Usage
### Index preparing

exSEEK docker contains a varitey of commonly used genomes and annotations. Besides of RNA types extracted from GENCODE V27, exSEEK can also analyze rRNA from NCBI refSeq 109, miRNA from miRBase, piRNA from piRNABank, circRNA from circBase, lncRNA and TUCP from mitranscriptome, repeats from UCSC Genome Browser (rmsk) and promoter and enhancer from ChromHMM tracks. You can use these .fa and .gtf files to generate the index you wanted:

### Small RNA-seq mapping

Run:

```bash
exseek normalization -d ${dataset}
```

This will generate normalized expression matrix for every combination of methods with the following file name pattern:

`output/${dataset}/matrix_processing/filter.${imputation_method}.Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}.txt`

You can specify normalization methods by setting the value of `normalization_method` and the batch removal methods
by setting the value of `batch_removal_method` in in `config/${dataset}.yaml`.

Supported normalization methods: TMM, RLE, CPM, CPM_top, UQ, null

Supported batch removal methods: limma, ComBat, RUV, null

When the method name is set to "null", the step is skipped.

`${batch_index}` is the column number (start from 1) in `config/${dataset}/batch_info.txt` to be used to remove batch effects.

### Feature selection

Run:

```bash
exseek feature_selection -d ${dataset}
```

This will evaluate all combinations of feature selection methods and classifiers by cross-validation.

Three summary files will be generated:

* `output/${dataset}/summary/cross_validation/metrics.test.txt`
* `output/${dataset}/summary/cross_validation/metrics.train.txt`
* `output/${dataset}/summary/cross_validation/feature_stability.txt`

Cross-validation results and trained models for individual combinations are in this directory:

`output/${dataset}/feature_selection/filter.${imputation_method}.Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}/${compare_group}/${classifier}.${n_select}.${selector}.${fold_change_filter_direction}`

Selected list of features are in `features.txt`.

> **Note:**
> More information about output files can be found on [File format](docs/file_format.md) page. Detailed parameters of feature selection and classifiers can be found in [config/machine_learning.yaml](config/machine_learning.yaml).


### Advanced Usage

* [Click here](docs/README.md) to see details

## Copyright and License Information

Copyright (C) 2019 Tsinghua University, Beijing, China 

This program is licensed with commercial restriction use license. Please see the [LICENSE](https://github.com/lulab/exSEEK_docs/blob/master/LICENSE) file for details.

## Citation

Binbin Shi, Jingyi Cao, Xupeng Chen and Zhi John Lu (2019) exSEEK: an integrative computational framework for identifying extracellular RNA biomarkers in liquid biopsy

