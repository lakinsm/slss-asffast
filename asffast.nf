#!/usr/bin/env nextflow

if( params.help ) {
	return help()
}

reference_db = file(params.db)
threads = params.threads
forks = params.forks
final_flag = params.final_flag
out_prefix = params.out_prefix


if( params.throughput != 'NONE_T' ) {
	throughput = Channel.fromPath("$params.throughput")
}
else {
	throughput = params.throughput
}


if( params.sequencing != 'NONE_S' ) {
	sequencing = Channel.fromPath("$params.sequencing")
}
else {
	sequencing = params.sequencing
}

Channel
    .fromPath( params.reads )
    .map{ in -> tuple("$in", in.getSimpleName(), file(in)) }
    .into{ fastq_barcodes; fastq_nobarcodes }


process PlotMinknowMetadata {
	publishDir "${params.output}/FlowcellRunMetadata", mode: "copy"

	input:
		file thrfile from throughput.collect()
		file seqfile from sequencing.collect()
		each file(sam) from alignment_curves
	output:
		file("nanopore_filecounts.csv")
		file("nanopore_throughput.csv")
		file("nanopore_stats_overall.txt")
		file("filecount_timeseries_graph.pdf")
		file("throughput_timeseries_graph.pdf")
		file("aligned_timeseries_graph.pdf")
		file("alignment_timeseries.csv")

	when:
		final_flag

	script:
		if( (thrfile.name != 'NONE_T') && (seqfile.name != 'NONE_S') )
			"""
			plot_nanopore_metadata.py $seqfile $thrfile $sam .
			"""
		else
			"""
			touch nanopore_filecounts.csv
			touch nanopore_throughput.csv
			touch nanopore_stats_overall.txt
			touch filecount_timeseries_graph.pdf
			touch throughput_timeseries_graph.pdf
			touch aligned_timeseries_graph.pdf
			touch alignment_timeseries.csv
			"""
}


process BwaIndexReference {
	input:
		path db_dir from reference_db
	output:
		file '*' into (asfv_reference)

	"""
	bwa index $db_dir
	samtools faidx $db_dir
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
		file("*aligned_reads.sam") into (consensus_sam, alignment_curves)

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
		file("${out_prefix}_${this_barcode}.vcf.gz")
		file("${out_prefix}_${this_barcode}_consensus.fasta")

	when:
		final_flag

	script:
		def this_barcode = getBarcode(sam.name).findAll().first()[1]
		"""
		samtools view -bS $sam | samtools sort - -o ${out_prefix}_${this_barcode}.bam
		freebayes -p 1 --standard-filters --min-coverage 10 -f $ref ${out_prefix}_${this_barcode}.bam | vcffilter -f "QUAL > 20" | bcftools view -Oz -o ${out_prefix}_${this_barcode}.vcf.gz
		bcftools index ${out_prefix}_${this_barcode}.vcf.gz
		cat $ref | bcftools consensus ${out_prefix}_${this_barcode}.vcf.gz > ${out_prefix}_${this_barcode}_consensus.fasta
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
    println "    --out_prefix    STR      output file/consensus prefix, e.g. flowcell or sample ID"
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
