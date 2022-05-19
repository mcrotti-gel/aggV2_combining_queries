# Combine genotype and functional annotation queries

This workflow allows you to extract variants and samples that comply to both a set of genotype and functional annotation filters, by 
intersecting the genotype VCFs with the functional annotation VCFs.

## Table of contents
- [Pipeline overview](#pipeline-overview)
- [Required inputs](#required-inputs)
  * [input_bed](#inputbed)
  * [aggv2_chunks_bed](#aggv2chunksbed)
  * [worst_consequence](#worstconsequence)
  * [severity_scale](#severityscale)
  * [severity](#severity)
  * [expression](#expression)
- [Outputs](#outputs)

## Pipeline overview
The pipeline has the following main processes:
* find_chunk: finds the genomic and functional annotation aggV2 chunks of interest.
* extract_variant_vep: filters the annotation aggV2 vcf.
* intersect_annotation_genotype_vcf: intersects the genomic vcf with the filtered annotation vcf.
* find_samples: finds samples of interest.
* summarise_output: produces summary tables.

## Required inputs

### input_bed

This is a region file of your genes of interest. This must be a three or column tab-delimited file of chromosome, start, and stop (with an option fourth column of an identifier - i.e. a gene name). The file should have the .bed extension. 

Example of `input_bed` file:

```bash
chr2	213005363	213151603	IKZF2
chr7	50304716	50405101	IKZF1
```

