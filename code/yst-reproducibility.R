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

# what I want is an analysis of the reproducibility with sample size of scaled and unscaled data
nreps=25
ss.vec <- c(0,3,5,7,10,14,20)

load(file='analysis/x.all.Rda')
load(file='analysis/x.s.all.Rda')

out.mat.scale <- matrix(data=NA, ncol=nreps, nrow=length(ss.vec))
out.mat.raw <- matrix(data=NA, ncol=nreps, nrow=length(ss.vec))
for(j in 2:length(ss.vec)){
  for(i in 1:nreps){
    sample.size=ss.vec[j]
	yst.sub <- cbind(yst[,sample(1:44,sample.size)],
	   yst[,sample(45:86,sample.size)])
	sub.conds <-  c(rep('S', sample.size), rep('W', sample.size))
	
	x.sub <- aldex.clr(yst.sub, sub.conds)
	x.sub.t <- aldex.ttest(x.sub)
	
	xs.sub <- aldex.clr(yst.sub, sub.conds, gamma=0.5)
	xs.sub.t <- aldex.ttest(xs.sub)
	
	sig.all <- rownames(x.all)[x.all$we.eBH < 0.05]
	sig.s.all <- rownames(x.s.all)[x.s.all$we.eBH < 0.05]
	sig.sub.t <- rownames(x.sub.t)[x.sub.t$we.eBH < 0.05]
	sig.s.sub.t <- rownames(xs.sub.t)[xs.sub.t$we.eBH < 0.05]
	out.mat.scale[1,i] <- i
	out.mat.scale[j,i] <- round(length(intersect(sig.s.sub.t, sig.s.all))/length(sig.s.all),2)
	out.mat.raw[1,i] <- i
	out.mat.raw[j,i] <- round(length(intersect(sig.sub.t, sig.all))/length(sig.all),2)
  }
}

save(out.mat.raw, file="analysis/out.mat.raw.Rda")
save(out.mat.scale, file="analysis/out.mat.scale.Rda")