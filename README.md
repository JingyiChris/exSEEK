# exSEEK

exSEEK is an integrated computational framework to discover and evaluate exRNA biomarkers for liquid biopsy.

The exSEEK framework consists of:
+ Pre_processing:
   
   + Building index with various types of genomes and annotations. []
   + Quality control and removing adaptors. [`exseek.py quality_control`] [`exseek.py cutadapt`]
   + Sequential mapping for small/long RNA-seq. [`exseek.py mapping`]
   + Counting expression matrix. [`exseek.py count_matrix`]

+ Main function:
   
   + Peak calling for recurring fragments of long RNAs. [`exseek.py call_domains`]
   + Normalization and batch removal. [`exseek.py normalization`]
   + Feature selection and classification. [`exseek.py feature_selection`]
   + Biomarker evaluation. [`exseek.py feature_selection`]

Table of Contents:

* [Installation](#installation)
* [Usage](#usage)
  * [Index preparing](#index-preparing)
  * [Small RNA-seq mapping](#small-rna-seq-mapping)
  * [Peak calling](#peak-(domain)-calling)
  * [Long RNA-seq mapping](#long-rna-seq-mapping)
  * [Counting expression matrix](#counting-expression-matrix)
  * [Normalization and batch removal](#normalization-and-batch-removal)
  * [Feature selection and biomarker evaluation](#feature-selection-and-biomarker-evaluation)
* [Copyright and License Information](#copyright-and-license-information)
* [Citation](#citation)


## Installation

For easy installation, you can use the [exSEEK image](https://hub.docker.com/r/ltbyshi/exseek) of [docker](https://www.docker.com) with all dependencies installed:
```bash
docker pull ltbyshi/exseek
```

All required software and packages are already installed in docker, so there are no more requirements. To test the installation and get information about the command-line interface of exSEEK, you can execute:
```bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py -h
```
The -v flag mounts the current working directory `$PWD` into the `/workspace` in docker image, so you can easily check the output files in `/workspace` directory after exiting docker.

You can create a bash script named `exseek` and set the script executable: 
```bash
#! /bin/bash
docker run --rm -it -v $PWD:/workspace -w /workspace ltbyshi/exseek exseek.py "$@"
```
After adding the file to one of the directories in the `$PATH` variable, you can simply run: `exseek`.


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

The basic usage of exSEEK is:
```bash
exseek ${step_name} -d ${dataset}
```

> **Note:**
> * `${step_name}` is one of the step listed in 'positional arguments'.
> * `${dataset}` is the name of your dataset that should match the prefix of your configuration file described in the following section.



## Usage

You can use the provided `example_data` to run exSEEK:
```bash
cp /apps/example_data /workspace
```

The `example_data` folder has the following structure:
```
example_data/
├── config
|   ├── example.yaml
|   ├── default_config.yaml
│   └── machine_learning.yaml

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
> * `config/example.yaml`: configuration file with frequently adjusted parameters, such as file paths and mapping parameters.
> * `config/default_config.yaml`: configuration file with additional detailed parameters for each step.
> * `config/machine_learning.yaml`: configuration file with parameters used for feature selection and classification step.
> * `data/example/batch_info.txt`: table of batch information.
> * `data/example/compare_groups.yaml`: table for definition of positive and negative samples.
> * `data/example/sample_classes.txt`: table of sample labels.
> * `output/example/`: output folder.


### Index preparing

exSEEK docker contains a variety of commonly used genomes and annotations. Besides of RNA types extracted from GENCODE V27, exSEEK can also analyze rRNA from NCBI refSeq 109, miRNA from miRBase, piRNA from piRNABank, circRNA from circBase, lncRNA and TUCP from mitranscriptome, repeats from UCSC Genome Browser (rmsk) and promoter and enhancer from ChromHMM tracks. You can use these `.fa` and `.gtf` files to generate the index you needed:


### Small RNA-seq mapping

#### Update sequential mapping order

The default mapping order is set as `rna_types` variable in `config/default_config.yaml`:
```yaml
rna_types: [rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, 
  snRNA, srpRNA, tRNA, tucpRNA, Y_RNA]
```

You can change the mapping order by add a `rna_types` variable in `config/example.yaml`. For example, add spike-in sequences as the first RNA type:
```yaml
rna_types: [spikein, rRNA, lncRNA, miRNA, mRNA, piRNA, snoRNA, 
  snRNA, srpRNA, tRNA, tucpRNA, Y_RNA]
```
```bash
exseek.py update_sequential_mapping -d example
```

#### Add new reference sequence

If a new RNA type is added, you should also add a sequence file in FASTA format: `${genome_dir}/fasta/${rna_type}.fa`. Then build a FASTA index \(`${genome_dir}/fasta/${rna_type}.fa.fai`\):
```bash
samtools faidx ${genome_dir}/fasta/${rna_type}.fa
```

Then build a bowtie2 index \(`${genome_dir}/index/bowtie2/${rna_type}`\):
```bash
bowtie2-build ${genome_dir}/fasta/${rna_type}.fa ${genome_dir}/index/bowtie2/${rna_type}
```

#### Quality control \(before adaptor removal\)

```bash
exseek.py quality_control -d example
```
> **Note:**
> * The detailed results for each sample are in folder `example_data/output/example/fastqc`. 
> * You can quickly check the summary results with the `fastqc.txt` file in `example_data/output/example/summary`.

#### Remove adapter

```bash
exseek.py cutadapt -d example
```
> **Note:**
> * Make sure that you have added your adaptor information in `example_data/config/example.yaml` file. 
> * You can check the adaptor revmoval summary with `example_data/output/example/summary/cutadapt.txt` file.

#### Quality control \(after adapter removal\)

```bash
exseek.py quality_control_clean -d example
```

#### Mapping

```bash
exseek.py mapping -d example
```

> **Note:**
> * Make sure that the parameter `small_rna` is `True` in `example_data/config/example.yaml`.
> * The output folder `example_data/output/example/gbam` contains genome bam files.
> * The output folder `example_data/output/example/tbam` contains transcriptome bam files for all types of RNA.
> * You can check the read length distribution for each type of RNA in folder `example_data/output/example/stats/mapped_read_length/`.
> * You can also check the summary of read counts mapped to all RNA types for all samples with the file `example_data/output/example/summary/read_counts.txt`.

#### Generate BigWig files

```bash
exseek.py bigwig -d example
```


### Peak (domains) Calling

exSEEK provides local maximum-based peak calling methods for identifying recurring fragments (domains) of long exRNAs, such as mRNA and lncRNA. These called domains can be used to conduct differential expression analysis and combined into the following expression matrix and serve as potential biomarker candidates.
```bash
exseek.py call_domains -d example
```

**Notes:**
* Domain calling parameters in `example_data/config/example.yaml`:
> * `call_domain_pvalue: "05"`: adjusted p-value threshold for defining peaks.
> * `bin_size: 20`: size of bins for calculating read coverage.
> * `cov_threshold: 0.05`: The fraction of samples that have the called peak. Peaks with cov_threshold above 0.05 are dedined as domains.

* Output files:
> * `example_data/output/example/domains_localmax_recurrence/recurrence.bed` contains .
> * `example_data/output/example/domains_localmax/domains.bed` contains .

The `recurrence.bed` file looks like:
| Transcript ID | start | end | X| frequency | strand |
| :--- | :--- | :--- | :--- | :--- |
| ENST00000365118.2 | 0 | 30 | X | 8 | + |
| ENST00000365223.1 | 0 | 61 | X | 12 | + |
| ENST00000365436.1 | 69 | 92 | X | 2 | + |
| ENST00000366365.2 | 236 | 261 | X | 1 | + |

The `domains.bed` file looks like:
| Transcript ID | start | end | | frequency | strand |
| :--- | :--- | :--- | :--- | :--- |
| ENST00000365118.2 | 0 | 30 | X | 8 | + |
| ENST00000365223.1 | 0 | 61 | X | 12 | + |
| ENST00000365436.1 | 69 | 92 | X | 2 | + |
| ENST00000366365.2 | 236 | 261 | X | 1 | + |

### Long RNA-seq mapping

The methods for long RNA-seq mapping are very similar to **2. Small RNA-seq mapping**. You can use the above command lines for long RNA-seq by setting `small_rna` to `False` in file `example_data/config/example.yaml`. There is no peak calling step for long RNA-seq, because there are no significantly recurring fragments (domains) detected in long RNA-seq datasets. 


### Counting expression matrix

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


### 6. Normalization and batch removal

exSEEK supports 5 kinds of normalization methods and 4 kinds of batch removal methods:

```yaml
normalization_method: ["TMM", "RLE", "CPM", "CPM_top", "null"]
batch_removal_method: ["ComBat", "limma", "RUV", "null"]
```

You can get normalized expression matrix generated by any combinations of normalization and batch removal methods by executing:

```bash
exseek normalization -d example
```

> **Notes:**
> * When the method name is set to "null", the step is skipped.
> * `${batch_index}` is the column number (start from 1) in `config/example/batch_info.txt` to be used to remove batch effects.
> * The name pattern of output files in folder `example_data/output/example/matrix_processing` is:   `Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}.txt`.

You can choose the best combination based on `UCA` score and `mKNN` score, which is summarized in folder `example_data/output/example/select_preprocess_method/uca_score` and `example_data/output/example/select_preprocess_method/knn_score`.
For a perfectly corrected expression matrix, both UCA and mKNN scores approach 1.

The `UCA` score files look like this:
| preprocess_method | uca_score |
| :--- | :--- |
| filter.null.Norm_CPM_top.Batch_limma_1 | 0.578 |
| filter.null.Norm_CPM.Batch_limma_1 | 0.563 |
| filter.null.Norm_CPM_top.Batch_ComBat_1 | 0.563 |
| filter.null.Norm_CPM_top.Batch_RUV_1 | 0.564 | 

And the `mKNN` score files look like this:
| preprocess_method | knn_score |
| :--- | :--- |
| filter.null.Norm_CPM_top.Batch_limma_1 | 0.940 |
| filter.null.Norm_CPM.Batch_limma_1 | 0.936 |
| filter.null.Norm_CPM_top.Batch_ComBat_1 | 0.936 |
| filter.null.Norm_CPM_top.Batch_RUV_1 | 0.927 | 

After deciding the most proper combination of normalization and batch removal methods, you can specify the exact normalization method by setting the value of `normalization_method` and the batch removal method by setting the value of `batch_removal_method` in `config/sample.yaml`.


### 7. Feature selection and biomarker evaluation

This step identifies and evaluates exRNA biomarker panels selected by various combinations of feature selection methods and machine learning classifiers. 

exSEEK supported feature selection methods:
```bash
[DiffExp_TTest, MaxFeatures_RandomForest, MaxFeatures_LogRegL1, 
MaxFeatures_LogRegL2, MaxFeatures_ElasticNet, RandomSubset_RandomForest, 
RandomSubset_LogRegL1, RandomSubset_LogRegL2, SIS, ReliefF, SURF, MultiSURF]
```

exSEEK supported classifiers:
```bash
[LogRegL2, RandomForest, RBFSVM, DecisionTree, MLP]
```

You can evaluate all combinations of feature selection methods and classifiers by cross-validation:
```bash
exseek feature_selection -d example
```

> **Note:**
> * The detailed parameters of feature selection and machine learning can be found in file:
`config/machine_learning.yaml`. 
> * The cross-validation results and trained models for individual combinations are in this directory:
`output/example/feature_selection/.Norm_${normalization_method}.Batch_${batch_removal_method}_${batch_index}.${count_method}/${compare_group}/${classifier}.${n_select}.${selector}.${fold_change_filter_direction}`.
> * Selected features (biomarker panels) for each combination can be found in `features.txt`.

Three summary files will be generated in this step:

* `output/example/summary/cross_validation/metrics.test.txt`
* `output/example/summary/cross_validation/metrics.train.txt`
* `output/sxample/summary/cross_validation/feature_stability.txt`

The `metrics.*.txt` file looks like:
| classifier | n_features | selector | fold_change_direction | compare_group | filter_method | imputation | normalization | batch_removal | count_method | preprocess_method | split | accuracy | average_precision | f1_score | precision | recall | roc_auc |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| LogRegL2 | 10 | MaxFeatures_RandomForest | any | Normal-HCC | filter | null | Norm_RLE | Batch_limma_1 | mirna_and_domains_rna | filter.null.Norm_RLE.Batch_limma_1 | 48 | 0.928 | 0.916 | 0.800 | 1.000 | 0.666 | 0.969 |
| LogRegL2 | 5 | DiffExp_TTest| any | Normal-HCC | filter | null | Norm_RLE | Batch_limma_1 | mirna_and_domains_rna | filter.null.Norm_RLE.Batch_limma_10 | 0 | 0.928 | 0.743 | 0.800 | 1.000 | 0.666 | 0.696 |

The `feature_stability.txt` file looks like:
| classifier | n_features | selector | fold_change_direction | compare_group | filter_method | imputation | normalization | batch_removal | count_method | preprocess_method | feature_stability
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| LogRegL2 | 5 | DiffExp_TTest | any | Normal-HCC | filter | null | Norm_RLE | Batch_limma_1 | mirna_and_domains_rna | filter.null.Norm_RLE.Batch_limma_1 | 0.450 |
| RBFSVM | 5 | DiffExp_TTest | any | Normal-stage_A | filter | null | Norm_RLE | Batch_limma_1 | mirna_and_domains_rna | filter.null.Norm_RLE.Batch_limma_1 | 0.473 |

You can choose the most proper combination and its identified features (biomarker panel) base on ROC_AUC and feature stability score summarized in the above three files. 


## Copyright and License Information

Copyright (C) 2019 Tsinghua University, Beijing, China 

This program is licensed with commercial restriction use license. Please see the [LICENSE](https://github.com/lulab/exSEEK_docs/blob/master/LICENSE) file for details.

## Citation

Binbin Shi, Jingyi Cao, Xupeng Chen and Zhi John Lu (2019) exSEEK: an integrative computational framework for identifying extracellular
