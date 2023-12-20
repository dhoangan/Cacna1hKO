library("tximeta")
library("DESeq2")
library("org.Mm.eg.db")
library("clusterProfiler")
library("enrichplot")
library("ensembldb")
library("apeglm")
library("utils")

#generate own annotation data base file from gtf file
#gtffile <- "gtf/Mus_musculus.GRCm39.104.gtf.gz"
#DB <- ensDbFromGtf(gtf=gtffile, outfile="Mus_musculus.GRCm39.104.DB")

EDB <- EnsDb("Mus_musculus.GRCm39.104.DB")


# Read sample infromation and quantification files
samples <- read.csv("sample_table.csv", header=TRUE, sep=";", check.names=FALSE)
files <- file.path("quants_filtered_GRCm39", "filtered_data", samples$quant_name, "quant.sf")
coldata <- samples
coldata$files <- files
coldata$names <- samples$quant_name
file.exists(coldata$files)

# change animal and genotype to factors
coldata$animal <- factor(coldata$animal)
coldata$genotype <- factor(coldata$genotype)
coldata$genotype <- relevel(coldata$genotype, ref="WT")
coldata$sex <- factor(coldata$sex)


# Create summarizedExperiment object from quantification files
se <- tximeta(coldata)

# Create summarizedExperiment object from quantification files on gene level 
gse <- summarizeToGene(se)

# Create DESeqDataSet object from summarizedExperiment
#ddsSE <- DESeqDataSet(se, design = ~ genotype)
ddsgeneSE <- DESeqDataSet(gse, design = ~ genotype)


# Differential expression analysis
#dds <- DESeq(ddsSE)
gene_dds <- DESeq(ddsgeneSE)


# define the levels to be compared
contrast <- c("genotype", "KO", "WT")

# Extract results from a DESeq analysis
#res <- results(dds)
gene_res <- results(gene_dds, contrast=contrast, alpha=0.05)
summary(gene_res)
mcols(gene_res)$description


#lfc shrinkage
resLFCshrink <- lfcShrink(gene_dds, coef="genotype_KO_vs_WT", type="apeglm")
summary(resLFCshrink)

#Annotation
#resAnno <- res
gene_resAnno <- gene_res
shrink_Anno <- resLFCshrink


#### annnotate gene-level results
gene_resAnno$GENENAME <- mapIds(EDB,
                           keys=rownames(gene_res),
                           column="GENENAME",
                           keytype="GENEID",
                           multiVals="first")

gene_resAnno$SYMBOL <- mapIds(EDB,
                                keys=rownames(gene_res),
                                column="SYMBOL",
                                keytype="GENEID",
                                multiVals="first")

gene_resAnno$ENTREZID <- mapIds(org.Mm.eg.db,
                              keys=rownames(gene_res),
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")


shrink_Anno$GENENAME <- mapIds(EDB,
                                keys=rownames(resLFCshrink),
                                column="GENENAME",
                                keytype="GENEID",
                                multiVals="first")

shrink_Anno$SYMBOL <- mapIds(EDB,
                              keys=rownames(resLFCshrink),
                              column="SYMBOL",
                              keytype="GENEID",
                              multiVals="first")


shrink_Anno$ENTREZID <- mapIds(org.Mm.eg.db,
                              keys=rownames(resLFCshrink),
                              column="ENTREZID",
                              keytype="ENSEMBL",
                              multiVals="first")



# Order annotated results
#resAnnoOrdered <- resAnno[order(resAnno$padj),]
gene_resAnnoOrdered <- gene_resAnno[order(gene_resAnno$padj),]
gene_resAnnoOrdered_l2fc <- gene_resAnno[order(gene_resAnno$log2FoldChange,
                                               decreasing=TRUE),]

# Export annotate results
write.csv(as.data.frame(shrink_Anno),file="shrink_Anno.csv")

# Export count matrix
count_gene_matrix_normalized = counts(gene_dds, normalized=TRUE)
count_gene_matrix_raw = counts(gene_dds, normalized=FALSE)

write.csv(count_gene_matrix_raw, file = "results_gene_counts_raw_221011.csv")
write.csv(count_gene_matrix_normalized, file = "results_gene_counts_normalized_221011.csv")


# subset genes 
padj_threshold = 0.05
lfc_threshold = 1
#bM_threshold = 20
res_subset_sig <- gene_resAnno[which(gene_resAnno$padj < padj_threshold &
                                          abs(gene_resAnno$log2FoldChange) >= lfc_threshold), ]

shrink_subset_sig <- shrink_Anno[which(shrink_Anno$padj < padj_threshold &
                                       abs(shrink_Anno$log2FoldChange) >= lfc_threshold), ]


res_subset_sig$ENSEMBL <- row.names(res_subset_sig)
shrink_subset_sig$ENSEMBL <- row.names(shrink_subset_sig)


## GO Analysis

# ORA:
# set background genes (all assayed genes)
genes_all <- rownames(gene_res)
genes_all_entrezid <- gene_resAnno$ENTREZID
# set differentially expressed genes(threshold for padj and log2fc)
genes_de <- rownames(gene_res)[which(abs(gene_resAnno$log2FoldChange) >= lfc_threshold &
                                       gene_resAnno$padj < padj_threshold) ]
genes_de_entrezid <- gene_resAnno$ENTREZID[which(abs(gene_resAnno$log2FoldChange) >= lfc_threshold &
                                                   gene_resAnno$padj < padj_threshold) ]

# default values for enrichGO arguments
# pAdjustMethod = "BH",
# qvalueCutoff = 0.2,
# minGSSize = 10,
# maxGSSize = 500,

results_ora_go <- enrichGO(gene = genes_de, # interested genes
                universe = genes_all,       # background genes
                OrgDb = org.Mm.eg.db,
                keyType = 'ENSEMBL',
                minGSSize = 1,
                ont = 'ALL',
                readable = TRUE)


head(results_ora_go)

#KEGG pathway over-representation analysis
ekegg <- enrichKEGG(gene = genes_de_entrezid,
                    universe = genes_all_entrezid,
                 organism     = 'mmu',
                 minGSSize = 1,
                 pvalueCutoff = 0.1)
head(ekegg)
