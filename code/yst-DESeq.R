library(DESeq2)

url <- "https://raw.githubusercontent.com/ggloor/datasets/main/transcriptome.tsv"
yst <- read.table(url, header=T, row.names=1)

# remove the one gene with 0 reads

yst <- yst[rownames(yst) != "YOR072W-B",]

# Gierlinski:2015aa
yst[,c('SNF2.6', 'SNF2.13','SNF2.25','SNF2.35')] <- NULL 
yst[,c('WT.21','WT.22','WT.25','WT.28','WT.34','WT.36')] <- NULL  

conds <- c(rep('S', 44), rep('W', 42))
coldata <- data.frame(conds)
dds <- DESeqDataSetFromMatrix(countData = yst,
          colData = coldata, design= ~ conds)
dds <- DESeq(dds)
res <- results(dds, name="conds_W_vs_S")
save(res, file='analysis/res.Rda')

