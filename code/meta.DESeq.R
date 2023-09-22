
# get the dataset
load(url('https://raw.githubusercontent.com/ggloor/datasets/main/ko.both.all.Rda'))

conds.K0 <- c(rep('H',8), rep('B',14), rep('B',14), rep('H', 8)) 

# DESeq2
library(DESeq2)

# make a vector of conditions
coldata <- data.frame(conds.K0)

dds <- DESeqDataSetFromMatrix(countData = ko.both.all,
          colData = coldata, design= ~ conds.K0)
dds <- DESeq(dds)
meta.DES.res <- results(dds, name="conds.K0_H_vs_B")
save(meta.DES.res, file='analysis/meta.DES.res.Rda')

# edgeR
# qlf is recommended for bulk RNA seq 
library(edgeR)
group=factor(conds.K0)
y <- DGEList(counts=ko.both.all,group=group)
y <- normLibSizes(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)
fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
meta.edg.qlf <- qlf$table
save(meta.edg.qlf, file='analysis/meta.edge.qlf.Rda')


