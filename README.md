# Shiny Reports

Shiny App to explore results from different analyses, such as Differential Expression analysis



## Shiny Apps available

- [Differential Expression Analysis](###differential-expression)

### Differential Expression

[Shiny App to explore Differential Expression analysis results](expression_shiny.R). The input must be an RData or RDS file storing an R list object. Each element of the list must be a data.frame storing the differential expression results for a dataset. The code is designed to work with `limma::topTable()` output format, with the gene names added as a new column. If the results come from, for example, `DESeq2` package, the column names (or the code) have to be manually edited. The essential colums are: logFC, se.coef, P.Value, adj.P.Val and geneNames. 

The se.coef column is computed by default by `DESeq2`, but not by `limma`. Following [Gordon Smith's advice](https://support.bioconductor.org/p/70175/ "Bioconductor post"), it can be computed like this:


```r
se.coef <- sqrt(fit$s2.post) * fit$stdev.unscaled
```

#### Demo

![](example.gif)
