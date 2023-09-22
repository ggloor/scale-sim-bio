devtools::load_all('~/Documents/0_git/ALDEx_bioc')


set.seed(2023)
load(url('https://raw.githubusercontent.com/ggloor/datasets/main/ko.both.Rda'))

# make a vector of conditions
conds.K0 <- c(rep('H',8), rep('B',14), rep('B',14), rep('H', 8)) 

xt <- aldex.clr(ko.both, conds.K0)
xt.e <- aldex.effect(xt, include.sample.summary=T)
xt.t <- aldex.ttest(xt)
xt.all <- cbind(xt.e, xt.t)

save(xt.all, file='analysis/xt.Rda')

## single gamma model
xg <- aldex.clr(ko.both, conds.K0, gamma=0.5)
xg.e <- aldex.effect(xg)
xg.t <- aldex.ttest(xg)
xg.all <- cbind(xg.e, xg.t)

save(xg.all, file='analysis/xg.Rda')
mu <- c(1,1.15) # 15% difference in 
mu.mod <- aldex.makeScaleMatrix(gamma=0.5, mu=mu, conds.K0)
 
xt.m <- aldex.clr(ko.both, conds.K0, gamma=mu.mod)
xt.m.e <- aldex.effect(xt.m, include.sample.summary=T)
xt.m.t <- aldex.ttest(xt.m)
xt.m.all <- cbind(xt.m.e, xt.m.t)

save(xt.m.all, file='analysis/xt.m.Rda')



