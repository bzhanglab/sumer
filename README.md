# SUMER

An R package for summarizing multiple enrichment analysis results


# Installation
```r
# install.packages("devtools")
devtools::install_github("bzhanglab/sumer")
```

# Usage 
The main function of the package is `sumer()`. This function takes two
parameters. The first parameter is the name of a file in json format
that, among other things, specifies the location of input data files for
each omics platform. The second parameter specifies the location of
output files. Currently, enrichment results from up to 7 platforms can
be summarized with `sumer()`.


```r
sumer("/path/to/config.json", "/path/to/output_dir")
```

A sample configuration file is shown below:

```js
{
  "project":"My study",
  "top_num": 50,
  "similarity": "Jaccard",
  "data": [
    {
      "platform_name":"rna-seq",
      "platform_abbr":"rna",
      "gmt_file":"/data/rna_seq.gmt",
      "score_file":"/data/rna_seq_score.txt"
    },
    {
      "platform_name":"proteomics",
      "platform_abbr": "pro",
      "gmt_file":"/data/proteomics.gmt",
      "score_file":"/data/proteomics_score.txt"
    },
    {
      "platform_name":"phospho-proteomics",
      "platform_abbr": "phospho",
      "gmt_file":"/data/phospho_proteomics.gmt",
      "score_file":"/data/phospho_proteomics_score.txt"
    }
  ]
}
```

The keys for the configuration file are explained in the following
table.

| name          | description                                                                                                  |
|---------------|--------------------------------------------------------------------------------------------------------------|
| project       | a brief description for your project                                                                         |
| top_num       | maximum number of top sets selected by set cover algorithm                                                   |
| similarity    | (optional) Similarity measure for gene sets, default is "Jaccard". The other supported measures are "Simpson", "Dice". |
| data          | an array of objects that describe the data for each platform                                                 |
| platform_name | name of platform                                                                                             |
| platform_abbr | abbreviation of the platform name                                                                            |
| gmt_file      | path to the gmt file                                                                                         |
| score_file    | path to the score file                                                                                       |


The `gmt_file` specifies a tab delimited file that provides gene set
information. Each row describes a gene set where the first column lists
the name of gene set and the optional second column a brief description
of the set. The rest of columns list the genes included in the set. This
is the same format as described
[here](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29).

The `score_file` lists the measurement of significance for each gene
set. There is one row for each gene set. The first column is the gene
set name and second column is the score for the set. Typical examples of
the score value are signed <img src="https://latex.codecogs.com/gif.latex?-\log\{p-value\}" title="-\log\{p-value\}" />  or signed
<img src="https://latex.codecogs.com/gif.latex?-\log\{FDR\}" title="-\log\{FDR\}" />. Columns are separated by tab.

The gene set names in `gmt_file` and `score_file` should match.

To run the sample data provided by the package, first make sure the
working directory (`/path/to/your_work_dir`) and output directory
(`/path/to/your_output_dir`) exist in your system.

```r
> setwd("/path/to/your_work_dir")
> file.copy(system.file("data", "sample.tgz", package="sumer"), ".")
> untar("sample.tgz")
> sumer("config.json", "/path/to/your_output_dir")
```

The final output files are in `your_output_dir`. You can open
`index.html` to explore the summarized results.

