library(ALDEx2)
data(selex)

# this is pretty counter-intuitive but here we go
# we start by getting our prior values from the point clr
# we need to do this to log(e)  because the rlnorm function
# uses this base.  
# tldr: the makeScaleMatrix function is messed up and needs fixing


conds <- c(rep('NS',7), rep('S',7))

sel.all <- aldex(selex, conditions=conds, CI=T, gamma=0.5)

sel.mu <- aldex.makeScaleMatrix(gamma=0.5, mu=c(10,10), conds)
sel.1.all <- aldex(selex, conditions=conds, CI=T, gamma=sel.mu)

# setting ratio as 1:50 explicitly
#sel.mu <- aldex.makeScaleMatrix(gamma=0.5, mu=c(1,50), conds, log=FALSE)
scale <- c(rep(1,7), rep(50,7)) 
sel.mu <- aldex.makeScaleMatrix(gamma=0.5, mu=scale, conds, log=FALSE)
sel.5.all <- aldex(selex, conditions=conds, CI=T, gamma=sel.mu)

# setting ratio in log space
# 2^-5.64 : 2^0  = 0.02 : 1
# this also works, but reversed because of factor()  
# sel.mu <- aldex.makeScaleMatrix(gamma=0.5, mu=c(0,-5.64), conds, log=TRUE)
scale <- c(rep(-5.6,7), rep(0,7))
sel.mu <- aldex.makeScaleMatrix(gamma=0.5, mu=scale, conds, log=TRUE)

sel.2.all <- aldex(selex,conditions=conds, CI=T, gamma=sel.mu)

save(sel.2.all, file="analysis/sel.2.all.Rda")
save(sel.5.all, file="analysis/sel.5.all.Rda")
save(sel.1.all, file="analysis/sel.1.all.Rda")
save(sel.all, file="analysis/sel.all.Rda")
