# HBV_TREE

```R
library(Biostrings)
library(msa)
library(phangorn)
library(ggtree)

for (genotype in c("A_Genomes", "B_Genomes", "C_Genomes")) {

    # Load in the libraries and sequences, and make an alignment
    seqs <- readAAStringSet(file.path(getwd(), "data", paste0(genotype, ".fas")))
    aln <- msa::msa(seqs, method=c("ClustalOmega"))

    # FIXME
    # WARNING: Overriding automatically determined seq-type DNA to Protein as requested

    # Convert the alignment to the phyat object
    aln <- as.phyDat(aln, type = "AA")

    # Make UPGMA and neighbor-joining trees from a distance matrix
    dist_mat <- dist.ml(aln)
    upgma_tree <- upgma(dist_mat)
    ggtree(upgma_tree)

    # Save as PNG
    ggsave(paste0("assets/", genotype, ".png"))
}
```

![Genotype A](https://raw.githubusercontent.com/kojix2/HBV_TREE/master/assets/A_Genomes.png)
![Genotype B](https://raw.githubusercontent.com/kojix2/HBV_TREE/master/assets/B_Genomes.png)
![Genotype C](https://raw.githubusercontent.com/kojix2/HBV_TREE/master/assets/C_Genomes.png)

## References

* [R Bioinformatics Cookbook](https://github.com/PacktPublishing/R-Bioinformatics-Cookbook) Chapter4
* [HBVdb::Index - Hepatitis B Virus Database](https://hbvdb.lyon.inserm.fr/)
