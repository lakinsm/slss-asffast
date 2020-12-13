#!/usr/bin/env nextflow

if( params.help ) {
	return help()
}

blastdb_name = file(params.db).name
blastdb_dir = file(params.db).parent
threads = params.threads
timeseries = params.timeseries

/*
screen_file = file(params.screen)
*/

Channel
    .fromFilePairs( params.reads, size: 1 ) { file -> file.getSimpleName() }
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
    .into { long_reads_fasta; long_reads_fastq }


process FastqToFasta {
	tag {file_id}

	input:
		set file_id, file(reads) from long_reads_fastq
	output:
		set file_id, file("${file_id}.fasta") into (fasta_from_fastq)

	when:
		isFastq(reads)

	"""
	fastq_to_fasta.py $reads > ${file_id}.fasta
	"""
}


process FastaToFasta {
	tag {file_id}

	input:
		set file_id, file(reads) from long_reads_fasta
	output:
		set file_id, file("${file_id}.fasta") into (fasta_from_fasta)

	when:
		isFasta(reads)

	"""
	if [ ! -f ${file_id}.fasta ]; then
		cp $reads ${file_id}.fasta
	fi
	"""
}


process NucBlast {
    tag {file_id}
    
    publishDir "${params.output}/BlastResults", mode: "symlink"
    
    input:
        set file_id, file(reads) from fasta_from_fastq.mix(fasta_from_fasta)
		path db_dir from blastdb_dir
    output:
        file("${file_id}_blast.sam") into (overall_data, timeseries_data)
    
    """
	blastn -db $db_dir/$blastdb_name -query $reads -out ${file_id}_blast.sam -num_threads $threads \
	-num_alignments 100000 -outfmt "17 SQ SR" -parse_deflines
    """
}


process CombineSamFiles {
	input:
		val(sam_list) from overall_data.toList()
	output:
		file("aligned_reads.sam") into (all_coverage)

	"""
	merge_sam_files.py ${sam_list} > aligned_reads.sam
	"""
}


process OverallAnalysis {
	publishDir "${params.output}/OverallResults", mode: "copy"

	input:
		file(all_sam) from all_coverage
	output:
		file("overall_coverage_report.csv")
		file("top_genomic_targets.pdf")
		file("aligned_reads.fasta")

	"""
	sam_parser.py -i $all_sam -oc overall_coverage_report.csv -op top_genomic_targets.pdf -r > aligned_reads.fasta
	"""
}


process TimeSeriesAnalysis {
	publishDir "${params.output}/TimeSeriesResults", mode: "copy"

	input:
		val(file_list) from timeseries_data.toList()
	output:
		file("time_series_results.pdf")
		file("time_series_results.csv")

	when:
		params.timeseries

	"""
	sam_parser.py -i '$file_list' -ot time_series_results.csv
	"""
}


def isFastq(def filename) {
	def pattern = ~/(\.fastq$)|(\.fq$)/
	return filename =~ pattern
}


def isFasta(def filename) {
	def pattern = ~/(\.fa$)|(\.fasta$)/
	return filename =~ pattern
}


def help() {
    println ""
    println "Program: viral-direct"
    println "Developer: Steven Lakin, USDA APHIS"
    println "Description: BLAST screening pipeline for single-end, long read FASTQ files from ONT platforms."
    println "             Files can optionally also be FASTA (automatically detected based on extension)."
    println ""
    println "Usage: viral-direct.nf [options]"
    println "Example: nextflow viral-direct.nf --reads \"containers/data/raw*\" --db containers/databases/tempdb --output temp_out --timeseries"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to single-end FASTQ input files (if using globstar, wrap in double quotes)"
    println "    --output        STR      path to output directory"
    println "    --db            STR      path to pre-made BLAST database to query against"
    println "    --screen        STR      path to .txt file of reference headers to exclude, one per line"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of process threads, default 2 (max thread use = maxForks * threads)"
    println "    --timeseries    FLAG     produce coverage curves over time"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}