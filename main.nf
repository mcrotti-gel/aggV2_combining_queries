#!/usr/bin/env nextflow

/* 
 * Set up variables
 */

vep_vcf = Channel
            .fromFilePairs(params.vep_vcf)
            .ifEmpty { exit 1, "Cannot find input file : ${params.vep_vcf}" }
vep_index = Channel
              .fromFilePairs(params.vep_vcf_index)
              .ifEmpty { exit 1, "Cannot find input file : ${params.vep_vcf_index}" }

vep_vcf
  .combine(vep_index, by:0)
  .set {vep_vcf_ch}

Channel
      .value(params.expression)
      .ifEmpty { exit 1, "Cannot find expression : ${params.expression}" }
      .into { expression1_ch; expression2_ch}

Channel
      .value(params.gene_symbol)
      .ifEmpty { exit 1, "Cannot find input gene : ${params.gene_symbol}" }
      .into { gene_symbol_ch; gene_name_ch}

Channel
      .fromPath(params.severity_scale)
      .ifEmpty { exit 1, "Cannot find severity scale : ${params.severity_scale}" }
      .set {severity_scale_ch}


geno_vcf = Channel
            .fromFilePairs(params.geno_vcf)
            .ifEmpty { exit 1, "Cannot find input file : ${params.geno_vcf}" }
geno_index = Channel
              .fromFilePairs(params.geno_vcf_index)
              .ifEmpty { exit 1, "Cannot find input file : ${params.geno_vcf_index}" }

geno_vcf
  .combine(geno_index, by:0)
  .set {geno_vcf_ch}




/*
 * Start pipeline
 */

process extract_variant_vep {

    input:
    val(gene) from gene_symbol_ch
    tuple val(vcf_name), file(vcf), file(vcf_index) from vep_vcf_ch
    file(severity_scale) from severity_scale_ch

    output:
    tuple val(gene), file("*.vcf.gz"), file("*.vcf.gz.tbi") into annotation_vcf_ch

    script:

    """
    bcftools +split-vep -i 'SYMBOL="'"${gene}"'"' -c SYMBOL -s worst:missense+ -S ${severity_scale} ${vcf} -O z -o ${gene}_annotation.vcf.gz
    tabix -p vcf ${gene}_annotation.vcf.gz
    """

}


process intersect_annotation_genotype_vcf {

    input:
    tuple val(vcf_name), file(vcf), file(vcf_index) from geno_vcf_ch
    tuple val(gene_name), file(anno_vcf), file(anno_vcf_index) from annotation_vcf_ch
	val(expression) from expression1_ch

    output:
    tuple file("intersect/0000.vcf.gz"), file("intersect/0000.vcf.gz.tbi") into intersect_out_vcf_ch

    script:

    """
    bcftools isec -i ${expression} -e- -p intersect -n=2 -O z ${vcf} ${anno_vcf}
    """

}


process find_samples {

    publishDir "${params.outdir}/combined_queries", mode: 'copy'

    input:
    tuple file(int_vcf), file(int_vcf_index) from intersect_out_vcf_ch
    val(gene) from gene_name_ch
	val(expression) from expression2_ch

    output:
    file("*.tsv")

    script:

    """
    bcftools query -i ${expression} -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/OLD_MULTIALLELIC\t%INFO/OLD_CLUMPED\t%FILTER\t%GT\n]' ${int_vcf} > ${gene}_results.tsv
    """
}
