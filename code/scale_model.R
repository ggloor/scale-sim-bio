# mu.vec = c(log2(rep(1,7)), log2(rep(1.2,7)))
# scale_samples <- t(sapply(mu.vec, FUN = function(mu) rlnorm(128, mu, 1e-3)))

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

scale <- c(-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1)

# samples by column, offsets by row
# using sweep '+'
data.out <- matrix(data=NA, ncol=5, nrow=20)

for(i in 1:5){
  x <- aldex.clr(datasets[[i]], conds[[i]])
  x.e <- aldex.effect(x, verbose=F)
  data.out[20,i] <- mean(x.e$diff.btw)
for(j in 1:19){
  mu <- c(1, 1+scale[j])
  g.scale <- aldex.makeScaleMatrix(0.5, mu, conds[[i]], log=FALSE)
  x <- aldex.clr(datasets[[i]], conds[[i]], gamma=g.scale, verbose=F)
  x.e <- aldex.effect(x)
  data.out[j,i] <- mean(x.e$diff.btw)
}
}

data.out <- cbind(data.out, c(scale, 0))

save(data.out, file='analysis/data.out.Rda')

#### build a scale model of the yeast transcriptome dataset that recapitulates the
# result but with a different denominator
