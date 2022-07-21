#!/usr/bin/env nextflow

/* --------------------------
  Set up variables
 ----------------------------*/

// user input bed file
my_bed_ch = Channel
            .fromPath(params.input_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.input_bed}" }

// aggV2 bed chunks
aggv2_bed_ch = Channel
            .fromPath(params.aggv2_chunks_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.aggv2_chunks_bed}" }

// VEP severity scale
severity_scale_ch = params.severity_scale ? Channel.fromPath(params.severity_scale, checkIfExists: true) : Channel.empty()


/*---------------------
  Start pipeline
 ----------------------*/

/*---------------------
  Find the genomic and annotated aggV2 vcf chunk 
 ----------------------*/

process find_chunk {
    
    publishDir "${params.outdir}", mode: 'copy'

    input:
    file(my_bed) from my_bed_ch
    file(aggv2_bed) from aggv2_bed_ch

    output:
    file(geno_files) into geno_vcf_list_ch
    file(anno_files) into vep_vcf_list_ch

    shell:

    '''
    
    while read -r line; do
    gene="$(echo "${line}"| awk '{print $4}')";
    printf "${line}" > ${gene}.bed;
    gvcf="$(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 10)";
    gvcf_index="$(echo $(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 10).csi)";
    avcf="$(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 11)";
    avcf_index="$(echo $(bedtools intersect -wo -a ${gene}.bed -b !{aggv2_bed} |cut -f 11).csi)";
    echo "$gene,$gvcf,$gvcf_index" >> geno_files
    echo "$gene,$avcf,$avcf_index" >> anno_files
    done < !{my_bed}

    '''

}

/*
 * Modify channels
 */

vep_vcf_list_ch
    	.splitCsv()
    	.map {row -> [row[0], file(row[1]), file(row[2])] }
    	.set {vep_vcf_ch}

geno_vcf_list_ch
    	.splitCsv()
    	.map {row -> [row[0], file(row[1]), file(row[2])] }
    	.set {geno_vcf_ch}


/*---------------------
  Extract variants in the functional annotation VCF 
 ----------------------*/

process extract_variant_vep {

    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(gene), file(avcf), file(avcf_index) from vep_vcf_ch
    each file(severity_scale) from severity_scale_ch

    output:
    tuple val(gene), file("${gene}_annotation.vcf.gz"), file("${gene}_annotation.vcf.gz.csi") into annotation_vcf_ch

    script:
    """
    if [ ${params.worst_consequence} == yes ]
    then
    	bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL -s worst:${params.severity}+ -S ${severity_scale} ${avcf} -O z -o ${gene}_annotation.vcf.gz
    else
    	bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL ${avcf} -O z -o ${gene}_annotation.vcf.gz
    fi

    bcftools index ${gene}_annotation.vcf.gz
    """

}

/*
 * Join channel geno_vcf_ch and annotation_vcf_ch
 */

 geno_vcf_ch
		.join(annotation_vcf_ch)
		.set {intersect_input_ch}


/*---------------------
Intersect functional annotation VCF with genotype VCF
----------------------*/

process intersect_annotation_genotype_vcf {

    publishDir "${params.outdir}/intersect", mode: 'copy'

    input:
    tuple val(gene), file(gvcf), file(gvcf_index), file(avcf_subset), file(avcf_subset_index) from intersect_input_ch

    output:
    tuple val(gene), file("${gene}_intersect/0000.vcf.gz"), file("${gene}_intersect/0000.vcf.gz.tbi") into intersect_out_vcf_ch

    script:

    """
    bcftools isec -i ${params.expression} -e- -p ${gene}_intersect -n=2 -O z ${gvcf} ${avcf_subset}
    """

}

/*---------------------
Process results
----------------------*/

process find_samples {

    publishDir "${params.outdir}/combined_queries", mode: 'copy'

    input:
    tuple val(gene), file(int_vcf), file(int_vcf_index) from intersect_out_vcf_ch

    output:
    file("${gene}_results.tsv") into query_result_ch

    script:

    """
    bcftools query -i ${params.expression} -f ${params.format} ${int_vcf} > ${gene}_results.tsv
    """
}

/*----------------------
Create summary tables 
-----------------------*/

process summarise_output {

    publishDir "${params.outdir}/combined_queries", mode: 'copy'

    input:
    file(query_result) from query_result_ch

    output:
    file("*_summary.tsv") 

    script:

    """
    summarise.R ${query_result}
    """

}