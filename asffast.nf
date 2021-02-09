#!/usr/bin/env nextflow

if( params.help ) {
	return help()
}

reference_db = file(params.db)
threads = params.threads
forks = params.forks
final_flag = params.final_flag


Channel
    .fromPath( params.reads )
    .map{ in -> tuple("$in", in.getSimpleName(), file(in)) }
    .into{ fastq_barcodes; fastq_nobarcodes }


process BwaIndexReference {
	input:
		path db_dir from reference_db
	output:
		file '*' into (asfv_reference)

	"""
	bwa index $db_dir
	"""
}


process BwaAlignNoBarcodes {
    tag {file_id}

    input:
        set full_name, file_id, path(reads) from fastq_nobarcodes
        file(index) from asfv_reference
		file(reference_db)
	output:
		file("*${file_id}.sam") into (final_data1, coverage_data1, final_coverage_data1)

	when:
		!params.barcodes

	script:
		def barcode_match = "barcode00"
    """
	bwa mem -t $threads $reference_db $reads > ${barcode_match}_${file_id}.sam
    """
}


process BwaAlignWithBarcodes {
    tag {file_id}

    input:
        set full_name, file_id, path(reads) from fastq_barcodes
        file(index) from asfv_reference
		file(reference_db)
	output:
		file("*${file_id}.sam") into (final_data2, coverage_data2, final_coverage_data2)

	when:
		params.barcodes

    script:
        def barcode_match = getBarcode(full_name).findAll().first()[1]
    """
	bwa mem -t $threads $reference_db $reads > ${barcode_match}_${file_id}.sam
    """
}


process CoverageAnalysisIntermediate {
	publishDir "${params.output}/CoverageAnalysis", mode: "copy"

	input:
		val(file_list) from coverage_data1.mix(coverage_data2).toList()
	output:
		file("coverage_results.csv")

	when:
		!final_flag

	"""
	echo "$file_list" > list_of_inputs.txt
	sam_parser_parallel.py -i list_of_inputs.txt -oc coverage_results.csv --threads $forks
	"""
}


process CoverageAnalysisFinal {
	publishDir "${params.output}/CoverageAnalysis", mode: "copy"

	input:
		val(file_list) from final_coverage_data1.mix(final_coverage_data2).toList()
	output:
		file("coverage_results.csv") into (cov_res)
		file("*coverage_results.pdf")
		file("*timeseries_results.csv")
		file("*timeseries_results.pdf")
		file("*aligned_reads.fasta")

	when:
		final_flag

	"""
	echo "$file_list" > list_of_inputs.txt

	sam_parser_parallel.py \
	-i list_of_inputs.txt \
	-oc coverage_results.csv \
	-op coverage_results.pdf \
	-ot timeseries_results.csv \
	-r \
	--threads $forks \
	--final

	merge_fasta_files.py .
	"""
}


process MergeAlignedSamFiles {
	input:
		val(file_list) from final_data1.mix(final_data2).toList()
	output:
		file("*aligned_reads.sam") into (consensus_sam)

	when:
		final_flag

	"""
	echo "$file_list" > list_of_inputs.txt
	merge_sam_files.py list_of_inputs.txt
	"""
}


process ProduceConsensus {
	publishDir "${params.output}/ConsensusSequences", mode: "copy"
	input:
		each file(sam) from consensus_sam
		file(ref) from reference_db
	output:
		file("*consensus.fasta")

	when:
		final_flag

	"""
	samtools view -bS $sam | samtools sort - -o ${sam}_sorted.bam

	bcftools mpileup -d 500 -f $ref ${sam}_sorted.bam | bcftools call -c | vcfutils.pl vcf2fq > ${sam}.fastq
	fastq_to_fasta.py ${sam}.fastq > ${sam}_consensus.fasta
	"""
}


def getBarcode(def filename) {
	def pattern = ~/(barcodes?[0-9]{1,2})/
	return filename =~ pattern
}


def help() {
    println ""
    println "Program: asffast.nf"
    println "Developer: Steven Lakin, USDA APHIS"
    println "Description: Alignment screening pipeline for single-end, long read FASTQ files from ONT platforms."
    println ""
    println "Usage: asffast.nf [options]"
    println "Example: asffast.nf"
    println ""
    println "Input/output options:"
    println ""
    println "    --reads         STR      path to single-end FASTQ input files (if using globstar, wrap in double quotes)"
    println "    --output        STR      path to output directory"
    println "    --db            STR      path to pre-made BLAST database to query against"
    println "    --barcodes      FLAG     input directory structure includes barcode folders"
    println ""
    println "Algorithm options:"
    println ""
    println "    --threads       INT      number of process threads, default 1 (max thread use = maxForks * threads)"
    println "    --final_flag    FLAG     indicate that this is the final run"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
