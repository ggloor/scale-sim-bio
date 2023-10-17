# mu.vec = c(log2(rep(1,7)), log2(rep(1.2,7)))
# scale_samples <- t(sapply(mu.vec, FUN = function(mu) rlnorm(128, mu, 1e-3)))
library(edgeR)

devtools::load_all('~/Documents/0_git/ALDEx_bioc')
data(selex)

load(url('https://raw.githubusercontent.com/ggloor/datasets/main/ko.both.Rda'))

url <- "https://raw.githubusercontent.com/ggloor/datasets/main/transcriptome.tsv"
yst <- read.table(url, header=T, row.names=1)
# remove the one gene with 0 reads
yst <- yst[rownames(yst) != "YOR072W-B", ]
# Gierlinski:2015aa
yst[,c('SNF2.6', 'SNF2.13','SNF2.25','SNF2.35')] <- NULL 
yst[,c('WT.21','WT.22','WT.25','WT.28','WT.34','WT.36')] <- NULL  

url <- "https://raw.githubusercontent.com/ggloor/datasets/main/meta16S.tsv"
rRNA <- read.table(url, header=T, row.names=1, sep='\t')
url <- "https://raw.githubusercontent.com/ggloor/datasets/main/singleCell.tsv"
ss <- read.table(url, header=T, row.names=1, sep='\t')
ss <- ss[,c(1:100,1502:1601)]

conds <- list()
conds$rconds <- c(rep(1,198), rep(2,161))
conds$sconds <- c(rep(1,7), rep(2,7))
conds$yconds <- c(rep(1,44), rep(2,42))
conds$mconds <- c(rep(1,8), rep(2,28), rep(1,8))
conds$ssconds <- c(rep(1,100), rep(2,100))

datasets <- list()
datasets$rRNA <- rRNA
datasets$selex <- selex
datasets$yst <- yst
datasets$meta <- ko.both
datasets$ss <- ss

norm.method <- c('RLE', 'TMM', 'TMMwsp', 'upperquartile')

# samples by column, offsets by row
# using sweep '+'
norm.data.out <- matrix(data=NA, ncol=6, nrow=5)
norm.nl.data.out <- matrix(data=NA, ncol=6, nrow=5)
for(i in 1:5){
	x <- aldex.clr(datasets[[5]], conds[[5]], gamma=0.3)
	x.e <- aldex.effect(x)
	norm.data.out[i,6] <- mean(x.e$diff.btw)
	norm.nl.data.out[i,6] <- mean(x.e$diff.btw)

	group=factor(conds[[i]])
	y <- DGEList(counts=datasets[[i]],group=group)
		
	norm.factor <- rep(1, length(conds[[i]]))
	
scale.matrix <- aldex.makeScaleMatrix(gamma=0.5, mu=norm.factor, conditions =conds[[i]])

	x <- aldex.clr(datasets[[i]], conds[[i]], gamma=scale.matrix )
	x.e <- aldex.effect(x, verbose=F)
	norm.data.out[i,5] <- mean(x.e$diff.btw)

scale.matrix <- aldex.makeScaleMatrix(gamma=0.5, mu=norm.factor, conditions =conds[[i]], log=F)

	x <- aldex.clr(datasets[[i]], conds[[i]], gamma=scale.matrix )
	x.e <- aldex.effect(x, verbose=F)
	norm.nl.data.out[i,5] <- mean(x.e$diff.btw)


for(j in 1:4){

# edgeR
# qlf is recommended for bulk RNA seq 
	group=factor(conds[[i]])
	y <- DGEList(counts=datasets[[i]],group=group)
	
	y <- normLibSizes(y, method=norm.method[j])
	
	norm.factor <- 2^(-1*log2(y[[2]]$norm.factors))
	
scale.matrix <- aldex.makeScaleMatrix(gamma=0.5, mu=norm.factor, conditions =conds[[i]])

	x <- aldex.clr(datasets[[i]], conds[[i]], gamma=scale.matrix )
	x.e <- aldex.effect(x, verbose=F)
	norm.data.out[i,j] <- mean(x.e$diff.btw)
scale.matrix <- aldex.makeScaleMatrix(gamma=0.5, mu=norm.factor, conditions =conds[[i]], log=F)

	x <- aldex.clr(datasets[[i]], conds[[i]], gamma=scale.matrix )
	x.e <- aldex.effect(x, verbose=F)
	norm.nl.data.out[i,j] <- mean(x.e$diff.btw)

}
}

rownames(norm.data.out) <- names(datasets)
colnames(norm.data.out) <- c(norm.method, 'iso', 'GM')

rownames(norm.nl.data.out) <- names(datasets)
colnames(norm.nl.data.out) <- c(norm.method, 'iso', 'GM')

save(norm.data.out, file='analysis/normalizations.Rda')
save(norm.nl.data.out, file='analysis/normalizations.nl.Rda')

#### build a scale model of the yeast transcriptome dataset that recapitulates the
# result but with a different denominator
