# In silico TFBS analysis

The *in silico* study of regulatory elements depends on the correct screening of the regions of the genome where these elements are located, and on the identification of characteristic and recurrent motifs in these regions, which are detected from a set of related sequences. The *in silico* detection of TFBSs can be divided into two approaches: i) detection of binding sites of TFs with known biding motifs, and ii) *de novo* motif discovery.

## 1. Detection of TFBSs with known PWMs
This approach consists of mapping a set of sequences against a databse of TFBSs, and then identifying the over-represented motifs in these sequences, in order to find those that possibly regulate the corresponding gene.
For this analysis, we use the software [FIMO](https://meme-suite.org/meme/doc/fimo.html) to search for instances of PWMs of plant TFs obtained from the [JASPAR CORE non-redundant collection](http://jaspar.genereg.net/). The promoter region was defined as a stretch of 800 bp upstream of the TSS, and TSSs were predicted using [TSSFinder](https://tssfinder.github.io/).
To run the script

## 2. *De novo* motif discovery
This approach consists of finding the most over-represented motif in a given set of sequences.
