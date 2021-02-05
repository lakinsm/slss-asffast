# slss-asffast
ASF-FAST Nextflow pipeline for USDA APHIS SLSS

## Authors:
 - Roger Barrette (original methods)
 - Oleg Osipenko (Python development)
 - Steven Lakin (optimization, containerization)
 
## Description:

Takes as input the root data directory for an ONT run and performs alignment of reads against a database of choice.  Options 
are available for Oxford Nanopore platforms to produce time-series coverage information, regarding the rate of coverage
gain as sequencing is performed.  Reports are generated for overall coverage data in PDF and CSV format.  A consensus
sequence is produced using available data when the sequencing run is terminated or when the user specifies.

## Installation and Requirements

**Runs on:**
- Debian or CentOS Linux

**Depends on:**

All of the following *must* be in the user's $PATH:
 - Python3.5+ with the numpy package
 - Nextflow v20+ (Java 8+)
 - Singularity v3+
 
**Installation:**

Install nextflow:

```shell script
curl -s https://get.nextflow.io | bash
```

[Install Singularity](https://sylabs.io/guides/3.7/user-guide/quick_start.html#quick-installation-steps) and make the 
following changes to the file singularity.conf (usually /usr/local/etc/singularity/singularity.conf):
 
 - [enable overlay](https://singularity-admindoc.readthedocs.io/en/latest/the_singularity_config_file.html#user-bind-control-boolean-default-yes)
 - add /run (and /scratch if on a distributed cluster) as a default bind path in singularity.conf, for example as below: 
 
```
# BIND PATH: [STRING]
# DEFAULT: Undefined
# Define a list of files/directories that should be made available from within
# the container. The file or directory must exist within the container on
# which to attach to. you can specify a different source and destination
# path (respectively) with a colon; otherwise source and dest are the same.
# NOTE: these are ignored if singularity is invoked with --contain except
# for /etc/hosts and /etc/localtime. When invoked with --contain and --net,
# /etc/hosts would contain a default generated content for localhost resolution.
#bind path = /etc/singularity/default-nsswitch.conf:/etc/nsswitch.conf
#bind path = /opt
#bind path = /scratch
bind path = /etc/localtime
bind path = /etc/hosts
bind path = /run
```

Clone this repository:

```shell script
git clone https://github.com/lakinsm/slss-asffast
```

## Quickstart Usage

**Standard (local executor) usage:**

```shell script
slss-asffast/bin/asffast.py -i data/ -o desired_output_dir -r slss-asffast/containers/data/databases/asfv_refseq_db.fasta
```

**SLURM (MPI executor) usage:**

```shell script
slss-asffast/bin/asffast.py -i data/ -o desired_output_dir -r slss-asffast/containers/data/databases/asfv_refseq_db.fasta --slurm
```

## Documentation

**Command-line options:**

```
Usage: bin/asffast.py [-h] -i INPUT -o OUTPUT -r REFERENCE

Input/output options:

    -i --input          STR      path to ONT data run directory
    -o --output         STR      path to output directory
    -r --reference      STR      path to multi-FASTA file of reference genomes

Algorithm and execution options:

    -s --singularity    STR      path to Singularity image (uses cloud by default if this isn't passed)
    --slurm             FLAG     Flag to execute pipeline using the SLURM executor for HPC clusters
    --wait              INT      number of seconds to watch for new data files before finalizing [1, inf)

Help options:

    -h --help                display this message
``` 

**Nextflow.config file:**

In the root directory of this repository is a file called nextflow.config.  You can modify this file to change certain
Nextflow parameters.  The most important of which is the "maxForks" option for each process in the pipeline.  This is
in the process section in that file.  MaxForks determines how many *processes* can be run in parallel.  To then determine 
maximum thread usage, you must multiply the maxForks * threads, since threads determine per-process thread utilization, 
but multiple process instances can run simultaneously. Nextflow handles the scheduling of these tasks.

**slss_hpc_slurm.conf file:**

This is the config file for SLURM execution across MPI-based distributed clusters.  You can modify SLURM options and 
change how Nextflow interfaces with SLURM in this file.  This file is independent of nextflow.config, so any changes 
made to nextflow.config must also be made to this file for equivalent effect.
