#!/usr/bin/env nextflow

/* 
 * Set up variables
 */

my_bed_ch = Channel
            .fromPath(params.input_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.input_bed}" }
            .splitCsv(sep='\t')
            .map { row -> tuple(row[0], file(row[1]), file(row[2]), file(row[3]))}

aggv2_bed_ch = Channel
            .fromPath(params.aggv2_chunks_bed)
            .ifEmpty { exit 1, "Cannot find input file : ${params.aggv2_chunks_bed}" }


Channel
      .value(params.expression)
      .ifEmpty { exit 1, "Cannot find expression : ${params.expression}" }
      .into { expression1; expression2}



Channel
      .fromPath(params.severity_scale)
      .ifEmpty { exit 1, "Cannot find severity scale : ${params.severity_scale}" }
      .set {severity_scale_ch}




/*
 * Start pipeline
 */


process find_chunk {

    input:
    file(bed) from my_bed
    file(bed) from aggv2_bed_ch

    output:
    tuple file(geno_vcf), file(geno_vcf_index)  into geno_vcf_ch
    tuple file(vep_vcf), file(vep_vcf_index) into vep_vcf_ch
    val(gene) into gene_symbol_ch

    script:

    """
    bedtools intersect -wo -a my_regions.bed -b aggV2_chunk_names.bed | cut -f 4 > gene
    bedtools intersect -wo -a my_regions.bed -b aggV2_chunk_names.bed | cut -f 10 > geno_vcf
    bedtools intersect -wo -a my_regions.bed -b aggV2_chunk_names.bed | cut -f 11 > vep_vcf
    echo $(cat geno_vcf).csi > geno_vcf_index
    echo $(cat vep_vcf).csi > vep_vcf_index
    """

}


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

    output:
    tuple file("intersect/0000.vcf.gz"), file("intersect/0000.vcf.gz.tbi") into intersect_out_vcf_ch

    script:

    """
    bcftools isec -i "'"${expression1}"'" -e- -p intersect -n=2 -O z ${vcf} ${anno_vcf}
    """

}


process find_samples {

    publishDir "${params.outdir}/combined_queries", mode: 'copy'

    input:
    tuple file(int_vcf), file(int_vcf_index) from intersect_out_vcf_ch
    val(gene) from gene_name_ch

    output:
    file("*.tsv")

    script:

    """
    bcftools query -i "'"${expression2}"'" -f '[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\t%INFO/OLD_MULTIALLELIC\t%INFO/OLD_CLUMPED\t%FILTER\t%GT\n]' ${int_vcf} > ${gene}_results.tsv
    """
}
