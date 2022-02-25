#!/usr/bin/env nextflow

/* 
 * Set up variables
 */

my_bed_ch = Channel
            .fromPath(params.input_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.input_bed}" }

aggv2_bed_ch = Channel
            .fromPath(params.aggv2_chunks_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.aggv2_chunks_bed}" }

Channel
      .fromPath(params.severity_scale)
      .ifEmpty { exit 1, "Cannot find severity scale : ${params.severity_scale}" }
      .set {severity_scale_ch}




/*
 * Start pipeline
 */


process find_chunk {
    
    publishDir "/Users/marcocrotti/Desktop/nextflow_pipelines/testing_bedtools/bed_files", mode: 'copy'

    input:
    file(my_bed) from my_bed_ch
    file(aggv2_bed) from aggv2_bed_ch

    output:
    file(geno_files) into geno_vcf_ch
    file(anno_files) into vep_vcf_ch

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


process extract_variant_vep {

    input:
    tuple val(gene), val(avcf), val(avcf_index) from vep_vcf_ch.splitCSV().map {row -> tuple(row.gene, val(row.avcf), val(row.avcf_index)) }
    file(severity_scale) from severity_scale_ch

    output:
    file("*.vcf.gz"), file("*.vcf.gz.tbi") into annotation_vcf_ch

    script:

    """
    bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL -s worst:missense+ -S ${severity_scale} ${vcf} -O z -o ${gene}_annotation.vcf.gz
    tabix -p vcf ${gene}_annotation.vcf.gz
    """

}


process intersect_annotation_genotype_vcf {

    input:
    tuple val(gene), val(gvcf), val(gvcf_index) from geno_vcf_ch.splitCSV().map {row -> tuple(row.gene, val(row.gvcf), val(row.gvcf_index)) }
    file(avcf_subset), file(avcf_subset_index) from annotation_vcf_ch

    output:
    tuple val(gene), file("${gene}_intersect/0000.vcf.gz"), file("${gene}_intersect/0000.vcf.gz.tbi") into intersect_out_vcf_ch

    script:

    """
    bcftools isec -i 'GT="AA" & INFO/AF<=0.05' -e- -p ${gene}_intersect -n=2 -O z ${vcf} ${anno_vcf}
    """

}


process find_samples {

    publishDir "${params.outdir}/combined_queries", mode: 'copy'

    input:
    tuple val(gene), file(int_vcf), file(int_vcf_index) from intersect_out_vcf_ch

    output:
    file("*.tsv")

    script:

    """
    bcftools query -i 'GT="AA"' -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/OLD_MULTIALLELIC\t%INFO/OLD_CLUMPED\t%FILTER\t%GT\n]' ${int_vcf} > ${gene}_results.tsv
    """
}
