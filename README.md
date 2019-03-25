# make9mer

### What is make9mer

**make9mer** is a 9-mer amino acid generator from protein sequence of length at least 18 amino acids. The 9-mers can then be used to query immunogenicity scores from Immune Epitope Data Base (IEDB.org)

# Installation

The package can be installed with

```r
install_github("nahtan/make9mer")
```

After installation, the package can be loaded into R.

    library(make9mer)

# Using make9mer

The main function in the **make9mer** package is `make9mer()`.

See example below.

# Parameters

**protseq** Sequence of the protein. Data.frame with: codon number in codon column, and amino acid for each codon in aminoacid column.


# Example
```r
protseq <- data.frame(  codon = seq(1:20)
                      , aminoacid = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
                      )
make9mer(protseq = protseq)
```
