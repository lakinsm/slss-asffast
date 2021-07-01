#!/usr/bin/env nextflow

if( params.help ) {
	return help()
}

reference_db = file(params.db)
threads = params.threads
forks = params.forks
out_prefix = params.out_prefix


if( params.throughput != 'NONE_T' ) {
	throughput = Channel.fromPath("$params.throughput")
}


if( params.sequencing != 'NONE_S' ) {
	sequencing_qc = Channel.fromPath("$params.sequencing")
	sequencing_cov = Channel.fromPath("$params.sequencing")
}

if( params.final_info != 'NONE_F') {
	final_info_cov = Channel.fromPath("$params.final_info")
	final_info_merge = Channel.fromPath("$params.final_info")
	final_info_consensus = Channel.fromPath("$params.final_info")
}

Channel
    .fromPath( params.reads )
    .map{ in -> tuple("$in", in.getSimpleName(), file(in)) }
    .into{ fastq_barcodes; fastq_nobarcodes }


process PlotMinknowMetadata {
	publishDir "${params.output}/FlowcellRunMetadata", mode: "copy"

	input:
		file thrfile from throughput.collect()
		file seqfile from sequencing_qc.collect()
	output:
		file("nanopore_filecounts.csv")
		file("nanopore_throughput.csv")
		file("nanopore_stats_overall.txt")
		file("filecount_timeseries_graph.pdf")
		file("throughput_timeseries_graph.pdf")

	when:
		params.final_info != 'NONE_F'

	script:
		if( (params.throughput != 'NONE_T') && (params.sequencing != 'NONE_S') )
			"""
			plot_nanopore_metadata.py $seqfile $thrfile .
			"""
		else
			"""
			touch nanopore_filecounts.csv
			touch nanopore_throughput.csv
			touch nanopore_stats_overall.txt
			touch filecount_timeseries_graph.pdf
			touch throughput_timeseries_graph.pdf
			"""
}


process MinimapIndexReference {
	input:
		path db_fpath from reference_db
	output:
		file 'slss_asffast2_reference.mmi' into (asfv_reference)

	"""
	minimap2 -x map-ont -d slss_asffast2_reference.mmi $db_fpath
	"""
}


process MinimapAlignNoBarcodes {
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
	minimap2 -N 1000 -a --eqx -x map-ont -t $threads -R "@RG\\tID:${barcode_match}\\tSM:${barcode_match}" $index $reads > ${barcode_match}_${file_id}.sam
    """
}


process MinimapAlignWithBarcodes {
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
	minimap2 -N 1000 -a --eqx -x map-ont -t $threads -R "@RG\\tID:${barcode_match}\\tSM:${barcode_match}" $index $reads > ${barcode_match}_${file_id}.sam
    """
}


process CoverageAnalysisIntermediate {
	publishDir "${params.output}/CoverageAnalysis", mode: "copy"

	input:
		val(file_list) from coverage_data1.mix(coverage_data2).toList()
	output:
		file("coverage_results.csv")

	when:
		params.final_info == 'NONE_F'

	"""
	echo "$file_list" > list_of_inputs.txt
	sam_parser_parallel.py -i list_of_inputs.txt -oc coverage_results.csv --threads $forks
	"""
}


process CoverageAnalysisFinal {
	publishDir "${params.output}/CoverageAnalysis", mode: "copy"

	input:
		val(file_list) from final_coverage_data1.mix(final_coverage_data2).toList()
		file final_file from final_info_cov.collect()
	output:
		file("coverage_results.csv") into (cov_res)
		file("*coverage_plots.pdf")

	when:
		params.final_info != 'NONE_F'

	"""
	echo "$file_list" > list_of_inputs.txt

	sam_parser_parallel.py \
	-i list_of_inputs.txt \
	-oc coverage_results.csv \
	-op coverage_plots.pdf \
	-r \
	--threads $forks \
	--final $final_file

	merge_fasta_files.py .
	"""
}


process MergeAlignedSamFiles {
	publishDir "${params.output}/AlignmentFiles", mode: "copy"

	input:
		val(file_list) from final_data1.mix(final_data2).toList()
		file final_file from final_info_merge.collect()
	output:
		file("*aligned_reads.sam") into (consensus_sam, alignment_curves)

	when:
		params.final_info != 'NONE_F'

	"""
	echo "$file_list" > list_of_inputs.txt
	merge_sam_files.py list_of_inputs.txt $final_file
	"""
}


process PlotAlignmentCurves {
	publishDir "${params.output}/AlignmentTimeseries", mode: "copy"

	input:
		file seqfile from sequencing_cov.collect()
		each file(sam) from alignment_curves
	output:
		file("*alignment_timeseries_graph.pdf")
		file("*alignment_timeseries_data.csv")

	when:
		(params.final_info != 'NONE_F') && (params.sequencing != 'NONE_S')

	script:
		"""
		plot_alignment_curves.py $sam $seqfile .
		"""
}


process ProduceConsensus {
	publishDir "${params.output}/ConsensusSequences", mode: "copy"
	input:
		each file(sam) from consensus_sam
		file(ref) from reference_db
		file final_file from final_info_consensus.collect()
	output:
		file("*.vcf.gz")
		file("*_consensus.fasta")

	when:
		params.final_info != 'NONE_F'

	script:
		def this_barcode = getBarcode(sam.name).findAll().first()[1]
		"""
		samtools faidx $ref
		samtools view -bS $sam | samtools sort - -o ${out_prefix}_${this_barcode}.bam
		freebayes -p 1 --standard-filters --min-coverage 10 -f $ref ${out_prefix}_${this_barcode}.bam | vcffilter -f "QUAL > 20" | bcftools view -Oz -o ${out_prefix}_${this_barcode}.vcf.gz
		bcftools index ${out_prefix}_${this_barcode}.vcf.gz
		extract_fasta_record.py $ref $final_file ${this_barcode} | bcftools consensus ${out_prefix}_${this_barcode}.vcf.gz > ${out_prefix}_${this_barcode}_consensus.fasta
		"""
}


def getBarcode(def filename) {
	def pattern = ~/(barcode[0-9]{1,2}[a-z]?)/
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
    println "    --final_info    STR      if specified, do final analysis, value is .tsv with best genome for each barcode"
    println ""
    println "Help options:"
    println ""
    println "    --help                   display this message"
    println ""
    return 1
}
