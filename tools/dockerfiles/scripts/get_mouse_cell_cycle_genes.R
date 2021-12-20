# https://www.r-bloggers.com/2016/10/converting-mouse-to-human-gene-names-with-biomart-package/
# https://github.com/satijalab/seurat/issues/2493

library(Seurat)
library(biomaRt)

convert_human_to_mouse <- function(x){
    human <- useMart("ensembl", dataset="hsapiens_gene_ensembl")
    mouse <- useMart("ensembl", dataset="mmusculus_gene_ensembl")
    genes <- getLDS(
        attributes=c("hgnc_symbol"),
        filters="hgnc_symbol",
        values=x,
        mart=human,
        attributesL=c("mgi_symbol"),
        martL=mouse,
        uniqueRows=T
    )
    return(unique(genes[, 2]))
}

g2m_genes <- convert_human_to_mouse(cc.genes$g2m.genes)
s_genes <- convert_human_to_mouse(cc.genes$s.genes)

phase <- c(
    rep("G2/M", length(g2m_genes)),
    rep("S", length(s_genes))
)
gene_id <- c(g2m_genes, s_genes)

write.table(
    data.frame(phase, gene_id),
    file="mouse_cell_cycle_genes.csv",
    sep=",",
    row.names=FALSE,
    col.names=TRUE,
    quote=FALSE
)