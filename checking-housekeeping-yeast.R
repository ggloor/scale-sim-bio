# note: need to manually pull changed genes and then get the equivalent in KEGG search

up.path <-  read.table('up.path.txt', row.names=1, sep='\t')
dn.path <-  read.table('dn.path.txt', row.names=1, sep='\t')

both.path <- intersect(rownames(up.path), rownames(dn.path))
up.only.path <- setdiff(rownames(up.path), rownames(dn.path))
dn.only.path <- setdiff(rownames(dn.path), rownames(up.path))
all.path <- union(rownames(up.path), rownames(dn.path))

all.path.df <- as.data.frame(matrix(data=NA, ncol=3, nrow=length(all.path))

colnames(all.path.df) <- c('pathway', 'N.up', 'N.dn')
 
rownames(all.path.df) <- all.path

rownames(all.path.df) %in% up.only.path

all.path.df[rownames(all.path.df) %in% up.only.path,1:2] <- up.path[up.only.path,1:2]
all.path.df[rownames(all.path.df) %in% dn.only.path,c(1,3)] <- dn.path[dn.only.path,c(1,2)]
all.path.df[rownames(all.path.df) %in% both.path,c(1,3)] <- dn.path[both.path,c(1,2)]
all.path.df[rownames(all.path.df) %in% both.path,2] <- up.path[both.path,2]

up.path[up.only.path,1]

ribo <- c('YBL027W','YBL072C','YBR146W','YBR181C','YBR189W','YDL075W','YDL082W','YDL136W','YDR064W','YDR237W','YDR382W','YDR418W','YDR447C','YEL054C','YER102W','YER131W','YGL068W','YGL123W','YGR085C','YHL004W','YHL015W','YIL069C','YIL133C','YJL177W','YKL180W','YKR094C','YLL045C','YLR048W','YLR333C','YLR441C','YML026C','YML063W','YMR142C','YMR242C','YNL069C','YNL178W','YOL040C','YOR096W','YOR312C','YPL131W','YPR102C')

ribo.unch <- c('YBR084C-A','YBR251W','YDL081C','YDL202W','YER050C','YER117W','YGL103W','YGR027C','YGR118W','YGR220C','YHR147C','YLR287C-A','YML024W','YML025C','YMR143W','YMR194W','YMR286W','YNL306W','YOR167C','YPL173W','YPL249C-A','YPR132W')