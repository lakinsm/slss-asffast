# slss-viral-direct
ASF-FAST Nextflow pipeline for USDA APHIS SLSS

## Authors:
 - Roger Barrette (original methods)
 - Oleg Osipenko (Python development)
 - Steven Lakin (optimization)
 
## Description:

Takes as input single-end FASTQ or FASTA data and performs nucleotide BLAST against a database of choice.  Options 
are available for Oxford Nanopore platforms to produce time-series coverage information, regarding the rate of coverage
gain as sequencing is performed.  Reports are generated for overall coverage data in PDF and CSV format.

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
 - (Optional) Singularity v3+
 
**Installation:**

Install nextflow:

```shell script
curl -s https://get.nextflow.io | bash
```

Clone this repository:

```shell script
git clone https://github.com/lakinsm/slss-viral-direct
```

## Quickstart Usage

**Create a BLAST database using the option -parse_seqids:**

```shell script
makeblastdb -parse_seqids -dbtype nucl -in /path/to/reference_genomes.fasta -out /path/to/db/my_db_name
```

**Run example pipeline:**

```shell script
cd slss-viral-direct
nextflow viral-direct.nf --reads "containers/data/raw/*.fastq" --db containers/data/databases/tempdb --output temp_out --timeseries
```

**Notes:**

Input can be FASTQ or FASTA or a mixture thereof.  Use of globstar (*) must be surrounded by double quotes. 
The database must be formatted to work with your current version of BLAST.  You can optionally provide a list of 
sequence headers from the reference database to exclude from analysis (these must exactly match the headers in the 
reference database), one header per line in a .txt file.


## Documentation

**Command-line options:**

```
Usage: viral-direct.nf [options]

Input/output options:

    --reads         STR      path to single-end FASTQ input files (if using globstar, wrap in double quotes)
    --output        STR      path to output directory
    --db            STR      path to pre-made BLAST database to query against
    --screen        STR      path to .txt file of reference headers to exclude, one per line

Algorithm options:

    --threads       INT      number of process threads, default 2 (max thread use = maxForks * threads)
    --timeseries    FLAG     produce coverage curves over time

Important Nextflow options (note single dashes for arguments):

    -resume         FLAG     use cached work to avoid recompute (this option is recommended for all runs)
    -work-dir       STR      path to directory for intermediate data (default current directory ./work)
    -help           FLAG     display help for Nextflow-specific options

Help options:

    --help                   display this message
``` 

**Nextflow.config file:**

In the root directory of this repository is a file called nextflow.config.  You can modify this file to change certain
Nextflow parameters.  The most important of which is the "maxForks" option for each process in the pipeline.  This is
in the process section in that file.  MaxForks determines how many *processes* can be run in parallel.  To then determine 
maximum thread usage, you must multiply the maxForks * threads, since threads determine per-process thread utilization, 
but multiple process instances can run simultaneously. Nextflow handles the scheduling of these tasks.

Singularity options can also be configured in the nextflow.config file.  If Singularity is not desired, it can be disabled here.


## Description of Output

**Directory BlastResults:**
 - **[Sample_name].sam** - These files are the Sequence Alignment Map (SAM) files output by BLAST.  They are symlinked 
 by default and actually reside in the Nextflow working directory in physical form, so if you delete your Nextflow 
 working directory, these files will also be deleted.  They are symlinked due to their size (to avoid copies).
 
**Directory OverallResults:**
 - **aligned_reads.fasta** - FASTA format file of reads aligning to the reference database, excluding reads screened out, 
 if a screen file was provided with the screen option.
 - **overall_coverage_report.csv** - Comma-separated value format file with coverage information for each reference genome
 in the database, using all reads that aligned and that were not screened out.  The following fields are present:
    1. **target** - The header of the reference genome in the reference database, as determined by makeblastdb.
    2. **avg_coverage** - The total genome average coverage (# of bases aligned / genome length).
    3. **percent_coverage** - The percent (0-100%) of bases in the genome at >= 1x coverage.
    4. **genome_length** - The length of the reference genome.
 - **top_genomic_targets.pdf** - A PDF file containing a rolling average coverage graph for each reference genome, one per page.
 Only the 5 most covered genomes are displayed.  The x-axis is the genomic position (location), and the y-axis is raw coverage 
 calculated on a rolling average.  The rolling average window size varies depending on genome size: it is calculated 
 to produce 500 windows over the genome with 50% window overlap.

**Directory TimeSeriesResults:**
 - **time_series_results.csv** -  Comma-separated value format file with coverage information for each reference genome
 in the database at each time point.  Time points are determined by the filename produced by the Oxford Nanopore sequencer. 
 The option for time series analysis will not work if the sequence filenames do not include this metadata. The following 
 fields are present:
    1. **timepoint** - The order in which the file was produced, starting with 1.
    2. **target** - The header of the reference genome in the reference database, as determined by makeblastdb.
    3. **percent_coverage** - The percent (0-100%) of bases in the genome at >= 1x coverage, up until this time point.  
    This is a cumulative measure, so time point 1 will have less than or equal coverage to time point 2, and so on.
    4. **genome_length** - The length of the reference genome.
 - **time_series_results.pdf** - A PDF file containing a multi-line graph of the top 5 most covered genomes over time. 
 The x-axis is the time point, starting with 0 coverage at 0 time for all genomes.  The y-axis is the percent coverage 
 (0-100%) for that genome cumulatively over time.  The color, labelled in the legend, represents the differnet genomes.
 By default, a maximum of 100 time points are displayed.
