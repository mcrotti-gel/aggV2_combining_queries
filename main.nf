#!/usr/bin/env nextflow

/* 
 * Set up variables
 */

Channel
      .fromPath( ["'params.vep_vcf, params.vep_vcf_index'"] )
      .ifEmpty { exit 1, "Cannot find input file : ${params.vep_vcf}" }
      .set {vep_vcf_ch}

Channel
      .value(params.gene_symbol)
      .ifEmpty { exit 1, "Cannot find input gene : ${params.gene_symbol}" }
      .set {gene_symbol_ch}

Channel
      .fromPath(params.severity_scale)
      .ifEmpty { exit 1, "Cannot find severity scale : ${params.severity_scale}" }
      .set {severity_scale_ch}

Channel
      .fromPath( ["'params.geno_vcf, params.geno_vcf_index'"] )
      .ifEmpty { exit 1, "Cannot find input geno file : ${params.geno_vcf}" }
      .set {geno_vcf_ch}


/*
 * Start pipeline
 */

process extract_variant_vep {

    input:
    val(gene) from gene_symbol_ch
    tuple file(vcf), file(vcf_index) from vep_vcf_ch
    file(severity_scale) from severity_scale_ch

    output:
    tuple val(gene), file("*.vcf.gz"), file("*.vcf.gz.tbi") into annotation_vcf_ch

    script:

    """
    bcftools +split-vep -i "'SYMBOL="IKZF1"'" -c SYMBOL -s worst:missense+ -S ${severity_scale} ${vcf} -O z -o ${gene}_annotation.vcf.gz
    tabix -p vcf ${gene}_annotation.vcf.gz
    """

}


process intersect_annotation_genotype_vcf {

    input:
    tuple file(vcf), file(vcf_index) from geno_vcf_ch
    tuple file(anno_vcf), file(anno_vcf_index) from annotation_vcf_ch

    output:
    tuple file("0000.vcf.gz"), file("0000.vcf.gz.tbi") into intersect_out_vcf_ch

    script:

    """
    bcftools isec -i 'GT="AA" & INFO/AF<=0.05' -e- -p intersect -n=2 -O z ${vcf} ${anno_vcf}
    """

}


process find_samples {

    input:
    tuple file(int_vcf), file(int_vcf_index) from intersect_out_vcf_ch

    output:
    file("*.tsv")

    script:

    """
    bcftools query -i "'GT="AA"'" -f "'[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/OLD_MULTIALLELIC\t%INFO/OLD_CLUMPED\t%FILTER\t%GT\n]'" ${int_vcf} > ${gene}_results.tsv
    """
}