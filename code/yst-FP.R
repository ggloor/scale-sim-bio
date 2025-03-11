#devtools::load_all('~/Documents/0_git/ALDEx_bioc')
library(ALDEx2)
### permutations below se li Li:2022aa
url <- "https://raw.githubusercontent.com/ggloor/datasets/main/transcriptome.tsv"
yst <- read.table(url, header=T, row.names=1)
# remove the one gene with 0 reads

yst <- yst[rownames(yst) != "YOR072W-B",]

# Gierlinski:2015aa
yst[,c('SNF2.6', 'SNF2.13','SNF2.25','SNF2.35')] <- NULL 
yst[,c('WT.21','WT.22','WT.25','WT.28','WT.34','WT.36')] <- NULL  

conds <- c(rep('S', 44), rep('W', 42))
coldata <- data.frame(conds)

# what I want is an analysis of the FP with a threshold
# subset vs remainder of set
nreps=25
ss.vec <- c(3,5,7,10,14,20)
counter <- 0

load(file='analysis/x.all.Rda')
load(file='analysis/x.s.all.Rda')

mat <- matrix(data=NA, ncol=20, nrow=nreps * length(ss.vec))

# note ref, sub  are the numbers in the reference and subset
# the reference is over-powered
# FP is the number in the subset not in the reference
FP.df <- data.frame(mat)
names(FP.df) <- c('count', 'ssize','ref.unscale', 'sub.unscale', 'ref.scale',
  'sub.scale', 'sig.s.t.ref1', 'sig.t.ref1','sig.s.t.ref1.4', 'sig.t.ref1.4','sig.s.t.ref2', 'sig.t.ref2','FP.scale', 'FP.unscale','FP.scale1', 'FP.unscale1', 'FP.scale1.4', 
  'FP.unscale1.4', 'FP.scale2', 'FP.unscale2')
# keep
# ssvec[j], nrep, sig.ref, sig.sub, sig.sc, sig.sc.sub, FP.sig.sub, FP.sig.sub1, FP.sig.sub2,
#   FP.sig.sc.sub, FP.sig.sc.sub1, FP.sig.sc.sub2
for(j in 1:length(ss.vec)){
  for(i in 1:nreps){
  	counter <- counter+1
  	print(counter)
    sample.size=ss.vec[j]
    
    FP.df[counter, 'count'] <- i
    FP.df[counter, 'ssize'] <- ss.vec[j]
    
	yst.sub <- cbind(yst[,sample(1:44,sample.size)],
	   yst[,sample(45:86,sample.size)])
	sub.conds <-  c(rep('S', sample.size), rep('W', sample.size))
	yst.ref <- yst[,!(colnames(yst) %in% colnames(yst.sub))]
	ref.conds <- c(rep('S', 44-sample.size), rep('W', 42-sample.size))
	
	# unscaled
	x.sub <- aldex(yst.sub, sub.conds)
	x.ref <- aldex(yst.ref, ref.conds)
	#scaled
	xs.sub <- aldex(yst.sub, sub.conds, gamma=0.5)
	xs.ref <- aldex(yst.ref, ref.conds, gamma=0.5)
	
	sig.ref <- rownames(x.ref)[x.ref$we.eBH < 0.05]
	sig.s.ref <- rownames(xs.ref)[xs.ref$we.eBH < 0.05]
	
	sig.sub <- rownames(x.sub)[x.sub$we.eBH < 0.05]
	sig.s.sub <- rownames(xs.sub)[xs.sub$we.eBH < 0.05]
	
	sig.t.ref1 <- rownames(x.ref)[x.ref$we.eBH < 0.05 & abs(x.ref$diff.btw) >1]
	sig.t.sub1 <- rownames(x.sub)[x.sub$we.eBH < 0.05 & abs(x.sub$diff.btw) >1]
	sig.t.ref1.4 <- rownames(x.ref)[x.ref$we.eBH < 0.05 & abs(x.ref$diff.btw) >1.4]
	sig.t.sub1.4 <- rownames(x.sub)[x.sub$we.eBH < 0.05 & abs(x.sub$diff.btw) >1.4]
	sig.t.ref2 <- rownames(x.ref)[x.ref$we.eBH < 0.05 & abs(x.ref$diff.btw) >2]
	sig.t.sub2 <- rownames(x.sub)[x.sub$we.eBH < 0.05 & abs(x.sub$diff.btw) >2]
	
	sig.s.t.ref1 <- rownames(xs.ref)[xs.ref$we.eBH < 0.05 & abs(xs.ref$diff.btw) >1]
	sig.s.t.sub1 <- rownames(xs.sub)[xs.sub$we.eBH < 0.05 & abs(xs.sub$diff.btw) >1]
	sig.s.t.ref1.4 <- rownames(xs.ref)[xs.ref$we.eBH < 0.05 & abs(xs.ref$diff.btw) >1.4]
	sig.s.t.sub1.4 <- rownames(xs.sub)[xs.sub$we.eBH < 0.05 & abs(xs.sub$diff.btw) >1.4]
	sig.s.t.ref2 <- rownames(xs.ref)[xs.ref$we.eBH < 0.05 & abs(xs.ref$diff.btw) >2]
	sig.s.t.sub2 <- rownames(xs.sub)[xs.sub$we.eBH < 0.05 & abs(xs.sub$diff.btw) >2]
	
	FP.df[counter,'ref.unscale'] <-  length(sig.ref)
	FP.df[counter,'sub.unscale'] <-  length(sig.sub)
	FP.df[counter,'ref.scale'] <-  length(sig.s.ref)
	FP.df[counter,'sub.scale'] <-  length(sig.s.sub)
	FP.df[counter,'sig.s.t.ref1'] <-  length(sig.s.t.ref1)
	FP.df[counter,'sig.t.ref1'] <-  length(sig.t.ref1)
	FP.df[counter,'sig.s.t.ref1.4'] <-  length(sig.s.t.ref1.4)
	FP.df[counter,'sig.t.ref1.4'] <-  length(sig.t.ref1.4)
	FP.df[counter,'sig.s.t.ref2'] <-  length(sig.s.t.ref2)
	FP.df[counter,'sig.t.ref2'] <-  length(sig.t.ref2)
	FP.df[counter,'FP.scale'] <- sum( !(sig.s.sub %in% sig.s.ref) )
	FP.df[counter,'FP.unscale'] <- sum( !(sig.sub %in% sig.ref) )
	FP.df[counter,'FP.scale1'] <- sum( !(sig.s.t.sub1 %in% sig.s.t.ref1) )
	FP.df[counter,'FP.unscale1'] <- sum( !(sig.t.sub1 %in% sig.t.ref1) )
	FP.df[counter,'FP.scale1.4'] <- sum( !(sig.s.t.sub1.4 %in% sig.s.t.ref1.4) )
	FP.df[counter,'FP.unscale1.4'] <- sum( !(sig.t.sub1.4 %in% sig.t.ref1.4) )
	FP.df[counter,'FP.scale2'] <- sum( !(sig.s.t.sub2 %in% sig.s.t.ref2) )
	FP.df[counter,'FP.unscale2'] <- sum( !(sig.t.sub2 %in% sig.t.ref2) )
	print(FP.df[counter,])
  }
}

# this is the summary file
save(FP.df, file="analysis/FP.df.Rda")

# these are final examples of the subset and ref data, 20/group in sub, 
# 24/20 in ref so this is really about concordance
save(x.sub, file="analysis/x.sub.Rda")
save(x.ref, file="analysis/x.ref.Rda")
save(xs.sub, file="analysis/xs.sub.Rda")
save(xs.ref, file="analysis/xs.ref.Rda")


#### Try some permutations
#### neither ALDEx2 nor DESeq2 returned FP on permutation
nperms <- 25
tally.BH <- vector(mode='numeric', length=nperms)
tally.p <- vector(mode='numeric', length=nperms)
tally.DES <- vector(mode='numeric', length=nperms)

for(i in 1:nperms){
yst.perm <- yst[,sample(colnames(yst))]
x.p <- aldex.clr(yst.perm, conds, gamma=1e-3)
x.pt <- aldex.ttest(x.p)
aldex.plot(x.pall)
tally.BH[i] <- sum(x.pt$we.eBH < 0.05)
tally.p[i] <- sum(x.pt$we.ep < 0.05)

dds <- DESeqDataSetFromMatrix(countData = yst.perm,
          colData = coldata, design= ~ conds)
dds <- DESeq(dds)
meta.DES.res <- results(dds, name="conds_W_vs_S")
tally.DES[i] <- sum(meta.DES.res@listData$p.adj < 0.05)
}


dds <- DESeqDataSetFromMatrix(countData = sample(yst),
          colData = coldata, design= ~ conds)
dds <- DESeq(dds)
meta.DES.res <- results(dds, name="conds.W_vs_S")
sum(meta.DES.res@listData$p.adj < 0.05)

