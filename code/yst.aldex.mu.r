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

mu.dat <- aldex.makeScaleMatrix(0.5,c(1,1.04),conds,mc.samples=128)

set.seed(2023)
x.s.mu.all <- aldex(yst, conds, gamma=mu.dat, verbose=F)
save(x.s.mu.all, file='analysis/x.s.mu.all.Rda')
