# Combine genotype and functional annotation queries

This workflow allows you to extract variants and samples that comply to both a set of genotype and functional annotation filters, by 
intersecting the genotype VCFs with the functional annotation VCFs.

## Table of contents
- [Pipeline overview](#pipeline-overview)
- [Required inputs](#required-inputs)
  * [input_bed](#input_bed)
  * [aggv2_chunks_bed](#aggv2_chunks_bed)
  * [worst_consequence](#worst_consequence)
  * [severity_scale](#severity_scale)
  * [severity](#severity)
  * [expression](#expression)
  * [format](#format)
- [Outputs](#outputs)
- [Examples](#examples)

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

### aggv2_chunks_bed

This is the list of chunk names and full file paths to both the genotype and functional annotation VCFs. This can be found under
`GEL data resources > aggregations > gel_mainProgramme > aggv2 > genomic > additional data > chunk names > aggv2_chunk_names.bed`

### worst_consequence

In the pipeline we use the bcftools plugin `+split-vep` to run annotation queries on vcfs. This parameter allows us to limit the output of a query to transcripts with the most severe consequence. It takes as value `yes` or `no`.

### severity_scale

This file lists the severity of variants. It can be found under
`GEL data resources > aggregations > gel_mainProgramme > somAgg > v0.2 > additional data > vep severity scale > VEP_severity_scale_2020.txt`.

### severity

With this parameter we choose the severity of variants we are interested in for our query. For example, if you want look only at missense variants or worse, the input value would be `missense`. It can be omitted if the parameter `worst_consequence` is set to `no`.

### expression

This parameter defines the bcftools filter of your query. See bcftools `EXPRESSIONS` for accepted filters https://samtools.github.io/bcftools/bcftools.html#expressions.

### format

This parameter defines the format of the query, see https://samtools.github.io/bcftools/bcftools.html#query for details. For the process to run, you should add the following fields `'[%SAMPLE\t%CHROM\t%POS\t%REF\t%ALT\n]'`, but you can also specify additional fields after the initial list.


## Outputs

This workflows produces three ouputs for each gene in your input bed file. 
* `*_result.tsv` file: this is a tab-delimited output from `bcftools query` command.
* `*_platekey_summary.tsv` file: this is a two-column tab-delimited file, where one column is the list of platekeys recovered by the query, and the second column is the number of variants per each participant that satisfied the query.
* `*_variant_summary.tsv` file: this is a two-column tab-delimited file, where one column is the list of variants that satisfied the query, and the second column is the number of participants that have that query.


## Examples

### Example 1

An example question would be: "I want to extract the samples in aggV2 who are homozygous alt for missense (or worse) rare variants within the gene IKZF1". 

The final command would look like this:

![Example 1](images/example1.png)

### Example 2

An example question would be: "I want to extract the samples in aggV2 who are homozygous alt for any type of variant within the gene IKZF1". 

The final command would look like this:

![Example 2](images/example2.png)
