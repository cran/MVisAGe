% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/perm.significance.R
\name{perm.significance}
\alias{perm.significance}
\title{A Function for Computing a Vector of Pearson Correlation Coefficients}
\usage{
perm.significance(exp.mat, cn.mat, gene.annot, method = "pearson",
  digits = 5, num.perms = 100, random.seed = NULL,
  alternative = "greater")
}
\arguments{
\item{exp.mat}{A matrix of gene-level expression data (rows = genes, columns = samples).  Missing values are not permitted.}

\item{cn.mat}{A matrix of gene-level DNA copy number data (rows = genes, columns = samples).  Both genes and samples should
appear in the same order as exp.mat.  Missing values are not permitted.}

\item{gene.annot}{A three-column matrix containing gene position information.  Column 1 = chromosome number written in 
the form 'chr1' (note that chrX and chrY should be written chr23 and chr24), Column 2 = position (in base pairs), Column 3 = cytoband.
Genes should appear in the same order as exp.mat and cn.mat.}

\item{method}{A character string (either "pearson" or "spearman") specifying the method used to calculate the correlation coefficient
(default = "pearson").}

\item{digits}{Used with signif() to specify the number of significant digits (default = 5).}

\item{num.perms}{Number of permutations used to assess significance (default = 1e2).}

\item{random.seed}{Random seed (default = NULL).}

\item{alternative}{A character string ("greater" or "less") that specifies the direction of the alternative hypothesis, 
either rho > 0 or rho < 0 (default = "greater").}
}
\value{
Returns a five-column matrix.  The first three columns are the same as gene.annot.  The fourth column contains
gene-specific Pearson or Spearman correlation coefficients based on the entries in each row of exp.mat and cn.mat, 
respectively (column name = "R").  The fifth column contains squared Pearson correlation coefficients (column name = "R^2").  
The sixth column contains the permutation-based right-tailed p-value of the correlation coefficient (column name = "perm_pValue").
The seventh column contains Benjamini-Hochberg q-values corresponding to the p-values.  Genes with constant gene expression
or DNA copy number are removed because they have zero variance.
}
\description{
This function computes Pearson correlation coefficients on a row-by-row basis for two numerical input matrices of the same size.
}
\examples{
exp.mat = tcga.exp.convert(exp.mat)

 cn.mat = tcga.cn.convert(cn.mat)

 prepped.data = data.prep(exp.mat, cn.mat, gene.annot, sample.annot, log.exp = FALSE)

 perm.significance(prepped.data[["exp"]], prepped.data[["cn"]], prepped.data[["gene.annot"]])

}
