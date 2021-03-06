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
