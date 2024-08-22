library(ALDEx2)


set.seed(2023)
load(url('https://raw.githubusercontent.com/ggloor/datasets/main/ko.both.Rda'))

# make a vector of conditions
conds.K0 <- c(rep('H',8), rep('B',14), rep('B',14), rep('H', 8)) 

# practically is identical to 0, and fills scale info slot
xt <- aldex.clr(ko.both, conds.K0, gamma=1e-3)
xt.e <- aldex.effect(xt, include.sample.summary=T)
xt.t <- aldex.ttest(xt)
xt.all <- cbind(xt.e, xt.t)

##### find scale estimate for hk functions
# commented out since just need values
# hk.off <- xt.all$diff.win < 2.5 & xt.all$diff.btw > 1 & xt.all$diff.btw < 3
# xhk <- aldex.clr(ko.both[hk.off,], conds.K0, gamma=1e-3)
# round(2^(mean(xhk@scaleSamps[c(1:8,37:44),]) -  mean(xhk@scaleSamps[c(9:36),])),2) # 1.14
# 
# x.lvha <- aldex.clr(ko.both, conds.K0, denom="lvha")
# denom <- x.lvha@denom
# 
# xd <- aldex.clr(ko.both[denom,], conds.K0, gamma=1e-3)
# round(2^(mean(xd@scaleSamps[c(1:8,37:44),]) -  mean(xd@scaleSamps[c(9:36),])),2) # 1.05

#####


save(xt.all, file='analysis/xt.Rda')

## single gamma model
xg <- aldex.clr(ko.both, conds.K0, gamma=0.5)
xg.e <- aldex.effect(xg)
xg.t <- aldex.ttest(xg)
xg.all <- cbind(xg.e, xg.t)

# can get the mean scale values
#  mean( xt@scaleSamps[c(1:8,37:44),]) #H
#  mean( xt@scaleSamps[c(9:36),]) #BV
for(i in 1:length(xt@dirichletData)){xt@dirichletData[[1]] <- NULL}
save(xt, file='analysis/xt.clr.Rda')
for(i in 1:length(xg@dirichletData)){xg@dirichletData[[1]] <- NULL}
save(xg, file='analysis/xg.clr.Rda')

save(xg.all, file='analysis/xg.Rda')
# this is the difference in housekeeping
mu <- c(1,1.14) # 14% difference in 
mu.mod <- aldex.makeScaleMatrix(gamma=0.5, mu=mu, conds.K0)
 
xt.m.all <- aldex(ko.both, conds.K0, gamma=mu.mod)

save(xt.m.all, file='analysis/xt.m.Rda')

# this is the difference from the lvha
mu <- c(1,1.05) # 14% difference in 
mu.mod <- aldex.makeScaleMatrix(gamma=0.5, mu=mu, conds.K0)
 
xt.lv.all <- aldex(ko.both, conds.K0, gamma=mu.mod)

save(xt.lv.all, file='analysis/xt.lv.Rda')




