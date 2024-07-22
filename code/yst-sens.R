# this gives a sensitivity analysis for transcriptome
# maybe supplement?
# x.sens <- aldex.senAnalysis(x, gamma=c(1e-3, 0.1, 0.2, 0.3, 0.4,0.5, 0.7,1))

library(ALDEx2)

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
x <- aldex.clr(yst, conds, gamma=1e-3, verbose=F)
x.sens <- aldex.senAnalysis(x, gamma=c(1e-3, 0.1, 0.2, 0.3, 0.4,0.5, 0.7,1))
save(x.sens, file="analysis/yst.sens.Rda")