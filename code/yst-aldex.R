library(DESeq2)
devtools::load_all('~/Documents/0_git/ALDEx_bioc')
#library(ALDEx2)

url <- "https://raw.githubusercontent.com/ggloor/datasets/main/transcriptome.tsv"
yst <- read.table(url, header=T, row.names=1)
# remove the one gene with 0 reads

yst <- yst[rownames(yst) != "YOR072W-B",]

# Gierlinski:2015aa
yst[,c('SNF2.6', 'SNF2.13','SNF2.25','SNF2.35')] <- NULL 
yst[,c('WT.21','WT.22','WT.25','WT.28','WT.34','WT.36')] <- NULL  

conds <- c(rep('S', 44), rep('W', 42))
coldata <- data.frame(conds)

set.seed(2023)
#ALDEx2
x <- aldex.clr(yst, conds, gamma=1e-3,verbose=F)
x.e <- aldex.effect(x, include.sample.summary=T, verbose=F)
x.t <- aldex.ttest(x, verbose=F)
x.all <- cbind(x.e, x.t)
save(x.all, file='analysis/x.all.Rda')

set.seed(2023)
x.s <- aldex.clr(yst, conds, gamma=0.3, verbose=F)
x.s.e <- aldex.effect(x.s, include.sample.summary=T, verbose=F)
x.s.t <- aldex.ttest(x.s, verbose=F)
x.s.all <- cbind(x.s.e, x.s.t)
save(x.s.all, file='analysis/x.s.all.Rda')

