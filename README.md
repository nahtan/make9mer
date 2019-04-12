# make9mer

### What is make9mer

**make9mer** is a 9-mer amino acid generator from protein sequence of length at least 18 amino acids. The 9-mers can then be used to query immunogenicity scores from Immune Epitope Data Base (IEDB.org)

**make9mer** generates a list which length is equal to the number of codons and includes a data.frame of all possible 9-mers by replacing the amino acid of the specific codon by any possible amino acid at each possible position in the 9-mer. 

# Installation

The package can be installed and loaded with

```r
install_github("nahtan/make9mer")
library(make9mer)
```

# Using make9mer

The main function in the **make9mer** package is `make9mer(protseq=protseq)`.

See example below.

# Parameters

**protseq** Sequence of the protein. Data.frame with: codon number in codon column, and amino acid for each codon in aminoacid column.

# Output

**list9mers** List of length equal of protein length. Includes data frames with 9-mers.

# Example
```r
protseq <- data.frame(  codon = seq(1:20)
                      , aminoacid = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y")
                      )
make9mer(protseq = protseq)
```
