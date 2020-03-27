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
---

Table of Contents:

* [Installation](#istallation)
* [Usage](#Usage)
  * [1. Index preparing](#1._Index_preparing)
  * [2. Small RNA-seq mapping](#2._Small_RNA-seq_mapping)
  * [3. Peak calling](#3._Peak_calling)
  * [4. Long RNA-seq mapping](#4._Long_RNA-seq_mapping)
  * [5. Counting expression matrix](#5._Counting expression matrix)
  * [6. Normalization and batch removal](#6._Normalization_and_batch_removal)
  * [7. Feature selection](#7._Feature_selection)
  * [8. Evaluation for biomarkers](#8._Evaluation_for_biomarkers)
* [Copyright and License Information](#copyright-and-license-information)
* [Citation](#citation)

---


## Installation

For easy installation, you can use the [exSEEK image](https://hub.docker.com/r/ltbyshi/exseek) of [docker](https://www.docker.com) with all dependencies installed:

```bash
docker pull ltbyshi/exseek
```
All required software and packages are already installed in docker, so there are no more requirements. To test the installation and get information about the command-line interface of exSEEK, you can execute:

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
After adding the file to one of the directories in the `$PATH` variable, you can simply run: `exseek`.

---

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
> * `output/example/`: input matrix of read counts.

---

## Usage

You can use the provided `example_data` to run exSEEK.

### 1.Index preparing

exSEEK docker contains a variety of commonly used genomes and annotations. Besides of RNA types extracted from GENCODE V27, exSEEK can also analyze rRNA from NCBI refSeq 109, miRNA from miRBase, piRNA from piRNABank, circRNA from circBase, lncRNA and TUCP from mitranscriptome, repeats from UCSC Genome Browser (rmsk) and promoter and enhancer from ChromHMM tracks. You can use these `.fa` and `.gtf` files to generate the index you needed:

---

### 2.Small RNA-seq mapping
#### 2.1 Update sequential mapping order

The default mapping order is set as `rna_types` variable in `config/default_config.yaml`:

```yaml
rna_types: [rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, 
  snRNA, srpRNA, tRNA, tucpRNA, Y_RNA]
```

You can change the mapping order by add a `rna_types` variable in `config/${dataset}.yaml`. For example, add spike-in sequences as the first RNA type:

```yaml
rna_types: [spikein, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, 
  snRNA, srpRNA, tRNA, tucpRNA, Y_RNA]
```
```bash
exseek.py update_sequential_mapping -d example
```

#### 2.2 Add new reference sequence

If a new RNA type is added, you should also add a sequence file in FASTA format: `${genome_dir}/fasta/${rna_type}.fa`. Then build a FASTA index \(`${genome_dir}/fasta/${rna_type}.fa.fai`\):

```bash
samtools faidx ${genome_dir}/fasta/${rna_type}.fa
```

Then build a bowtie2 index \(`${genome_dir}/index/bowtie2/${rna_type}`\):

```bash
bowtie2-build ${genome_dir}/fasta/${rna_type}.fa ${genome_dir}/index/bowtie2/${rna_type}
```

#### 2.3 Quality control \(before adaptor removal\)

```bash
exseek.py quality_control -d example
```
> **Note:**
> * The detailed results for each sample are in `example_data/output/example/fastqc`. 
> * You can quickly check the summary results with the `fastqc.txt` file in `example_data/output/example/summary`.

#### 2.4 Remove adapter

```bash
exseek.py cutadapt -d example
```
> **Note:**
> * Make sure that you have added your adaptor information in `example_data/config/example.yaml` file. 
> * You can check the adaptor revmoval summary with `example_data/output/example/summary/cutadapt.txt` file.

#### 2.5 Quality control \(after adapter removal\)

```bash
exseek.py quality_control_clean -d example
```

#### 2.6 Mapping

```bash
exseek.py mapping -d example
```
> **Note:**
> * Make sure that the parameter `small_rna` is `True` in `example_data/config/example.yaml`.
> * The output folder `example_data/output/example/gbam` contains genome bam files.
> * The output folder `example_data/output/example/tbam` contains transcriptome bam files for all types of RNA.
> * You can check the read length distribution for each type of RNA in folder `example_data/output/example/stats/mapped_read_length/`.
> * You can also check the summary of read counts mapped to all RNA types for all smaples with the file `example_data/output/example/summary/read_counts.txt`.

#### 2.7 Generate BigWig files

```bash
exseek.py bigwig -d example
```

---

### 3. Peak (domains) Calling
exSEEK provides peak calling methods for identifying conserved fragments (domains) of long exRNAs. These domains can be used to conduct differential expression analysis and combined into the following expression matrix and serve as potential biomarker candidates.

```bash
exseek.py call_domains -d example
```

**Notes:**
* Domain calling parameters in `example_data/config/example.yaml`:
> * `call_domain_pvalue: "05"`: adjusted p-value threshold for defining peaks.
> * `bin_size: 20`: size of bins for calculating read coverage
> * `cov_threshold: 0.05`: The fraction of samples that have the called peak. Peaks with cov_threshold above 0.05 are dedined as domains.
* Output files:
> * `example_data/output/example/domains_localmax_recurrence/recurrence.bed` contains .
> * `example_data/output/example/domains_localmax/domains.bed` contains .

---

### 4. Long RNA-seq mapping

The methods for long RNA-seq mapping are very similar to **2. Small RNA-seq mapping**. You can use the above commandlines for long RNA-seq by setting `small_rna` to `False` in file `example_data/config/example.yaml`. There is no peak calling step for long RNA-seq, beacuse there are no significantly conserved domains detected in long RNA-seq datasets. 

---

### 5. Counting expression matrix

```bash
exseek.py count_matrix -d example
```
```bash
exseek.py combine_domains -d example
```

**Notes:**
* The default methods for counting expression matrix includes:
> * `mirna_and_domains`:
> * `domains_long`:
> * `transcript`:
> * `transcript_mirna`:
> * `domains_long`:

---

### 6. Normalization and batch removal

exSEEK supports 5 kinds of normalization methods and 4 kinds of batch removal methods:

```yaml
normalization_method: ["TMM", "RLE", "CPM", "CPM_top", "null"]
batch_removal_method: ["ComBat", "limma", "RUV", "null"]
```

You can get normalized expression matrix generated by any combination of normalization and batch removal methods by executing:

```bash
exseek normalization -d example
```

> **Notes:**
> * When the method name is set to "null", the step is skipped.
> * `${batch_index}` is the column number (start from 1) in `config/example/batch_info.txt` to be used to remove batch effects.
> * The name pattern of output files in folder `example_data/output/example/matrix_processing` is:   `Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}.txt`.

You can choose the best combination based on `UCA` score and `mKNN` score, which is summarized in folder `example_data/output/example/select_preprocess_method/uca_score` and `example_data/output/example/select_preprocess_method/knn_score`.

The `UCA` score files look like this:

And the `mKNN` score files look like this:

After deciding the most proper combination of normalization and batch removal methods, you can specify the exact normalization method by setting the value of `normalization_method` and the batch removal method by setting the value of `batch_removal_method` in `config/sample.yaml`.

---

### Feature selection and biomarker evaluation
This step identifies and evaluates exRNA biomarker panels generated from various feature selection methods and machine learning classifiers. 

exSEEK supported feature selection methods include:

```bash
[DiffExp_TTest, MaxFeatures_RandomForest, MaxFeatures_LogRegL1, MaxFeatures_LogRegL2, MaxFeatures_ElasticNet, RandomSubset_RandomForest, RandomSubset_LogRegL1, RandomSubset_LogRegL2, SIS, ReliefF, SURF, MultiSURF]
```

exSEEK supported classifiers include:

```bash
[LogRegL2, RandomForest, RBFSVM, DecisionTree, MLP]
```

You can evaluate all combinations of feature selection methods and classifiers by cross-validation:

```bash
exseek feature_selection -d example
```
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

Binbin Shi, Jingyi Cao, Xupeng Chen and Zhi John Lu (2019) exSEEK: an integrative computational framework for identifying extracellular
