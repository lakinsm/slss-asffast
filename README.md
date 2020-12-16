# slss-asffast
ASF-FAST Nextflow pipeline for USDA APHIS SLSS

## Authors:
 - Roger Barrette (original methods)
 - Oleg Osipenko (Python development)
 - Steven Lakin (optimization)
 
## Description:

Takes as input the root data directory for an ONT run and performs alignment of reads against a database of choice.  Options 
are available for Oxford Nanopore platforms to produce time-series coverage information, regarding the rate of coverage
gain as sequencing is performed.  Reports are generated for overall coverage data in PDF and CSV format.  A consensus
sequence is produced using available data when the sequencing run is terminated or when the user specifies.

## Installation and Requirements

**Runs on:**
- Debian or CentOS Linux
- MacOS

**Depends on:**
 - Python3
 - Python matplotlib
 - Python numpy
 - Nextflow v20+ (Java 8+)
 - BWA v0.7.17
 - Samtools 1.10
 - Bcftools 1.10.2
 - (Optional) Singularity v3+
 
**Installation:**

Install nextflow:

```shell script
curl -s https://get.nextflow.io | bash
```

Clone this repository:

```shell script
git clone https://github.com/lakinsm/slss-asffast
```

## Quickstart Usage

**Example usage:**

```shell script
slss-asffast/bin/asffast.py -i data/ -o desired_output_dir -r slss-asffast/containers/data/databases/asfv_refseq_db.fasta
```

## Documentation

**Command-line options:**

```
Usage: bin/asffast.py [-h] -i INPUT -o OUTPUT -r REFERENCE

Input/output options:

    -i --input      STR      path to ONT data run directory
    -o --output     STR      path to output directory
    -r --reference  STR      path to multi-FASTA file of reference genomes

Algorithm options:

    --wait          INT      number of seconds to watch for new data files before finalizing [1, inf)

Help options:

    -h --help                display this message
``` 

**Nextflow.config file:**

In the root directory of this repository is a file called nextflow.config.  You can modify this file to change certain
Nextflow parameters.  The most important of which is the "maxForks" option for each process in the pipeline.  This is
in the process section in that file.  MaxForks determines how many *processes* can be run in parallel.  To then determine 
maximum thread usage, you must multiply the maxForks * threads, since threads determine per-process thread utilization, 
but multiple process instances can run simultaneously. Nextflow handles the scheduling of these tasks.

Singularity options can also be configured in the nextflow.config file.  If Singularity is not desired, it can be disabled here.

## Description of Output
