# LocusCompareR

## 1. Installation
LocusCompareR is an R package for visualization of GWAS-eQTL colocalization events. 


Use the following commands to install LocusCompareR. If you don't have `devtools`, uncomment the first line to install it. 

```
# install.packages("devtools")
devtools::install_github("boxiangliu/locuscomparer")
```

## 2. Example

To illustrate the use of locuscompare, provide two data frames with chromosome, position, and -log10 p-value:

```
library(locuscomparer)
gwas_df = data.frame(chromosome = '1', position = 1:1000, check.names = FALSE)
gwas_df[['-log10 p-value']] = runif(1000, 0, 12)
eqtl_df = data.frame(chromosome = '1', position = 1:1000, check.names = FALSE)
eqtl_df[['-log10 p-value']] = runif(1000, 0, 12)
lead_ld = data.frame(chromosome = '1', position = 1:1000, r2 = runif(1000, 0, 1))
locuscompare(in_fn1 = gwas_df, in_fn2 = eqtl_df, lead_ld = lead_ld,
             title1 = 'CAD GWAS', title2 = 'Coronary Artery eQTL')
```

The output from the `main` function is a figure like the following:

![](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/fig/locuscompare.png)

The labeled SNP is the lead SNP (in this case for both studies), and other SNPs are colored according to their LD $r^2$ with the lead SNP.

## 3. Using your own dataset:

The input to `locuscompare()` is two data frames with three required columns:

1. chromosome
2. position
3. -log10 p-value

Here is an example data frame:

```
chromosome	position	-log10 p-value
1	12345	3.21
1	12500	1.87
1	13000	5.42
1	14000	0.90
```

You can download the example files below:  [GWAS](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/inst/extdata/gwas.tsv) and [eQTL](https://raw.githubusercontent.com/boxiangliu/locuscomparer/master/inst/extdata/eqtl.tsv) datasets. 

Then run the following commands: 
```
library(locuscomparer)
gwas_df = read.table('path/to/gwas.tsv', header = TRUE, check.names = FALSE)
eqtl_df = read.table('path/to/eqtl.tsv', header = TRUE, check.names = FALSE)
lead_ld = read.table('path/to/lead_ld.tsv', header = TRUE)
locuscompare(in_fn1 = gwas_df, in_fn2 = eqtl_df, lead_ld = lead_ld,
             title1 = 'GWAS', title2 = 'eQTL')
```

## 4. Documentations

To view documentation for each function, type ?[function name] in the R console. 

LocusCompareR current export the following functions:

**Data munging**

- `assign_color`: Assign color to each SNP according to LD. 
- `get_lead_snp`: Add a column of SNP labels to input data.frame.
- `get_position`: Append two columns, chromosome (chr) and position (pos), to the input data.frame.

**Plotting**

- `locuscompare`: Make a locuscompare plot.
- `make_combined_plot`: Generated a combined plot with two locuszoom plots and a locuscompare plot.
- `make_locuszoom`: Make a simple locuszoom plot.
- `make_scatterplot`: Make a scatter plot (called the LocusCompare plot)

**Data loading**

- `read_metal`" Read association summary statistics from file. 
- `retrieve_LD`: Retrive SNP pairwise LD from database.

## 5. Citation

If you use locuscompare, please cite the following paper: https://www.nature.com/articles/s41588-019-0404-0


Boxiang Liu, Michael J. Gloudemans, Abhiram S. Rao, Erik Ingelsson & Stephen B. Montgomery (2019) Abundant associations with gene expression complicate GWAS follow-up, *Nature Genetics*
