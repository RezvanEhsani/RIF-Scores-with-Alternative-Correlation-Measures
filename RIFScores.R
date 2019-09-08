######################################################################################################################################################################
#
#
#   - This code is to analysis Regulatory Factors(RF) e.g. Transcription Factor(TF), Epigenetic Factor(EF), and long-non-coding RNA (lncRNA) 
#	based on Regulatory Impact Factor(RIF) scores.
# 
#   - Please first download related data from supplementary and change below addresses according to your system location of these files.
#   - Main function is 'RIF_Analysis' that integrates three steps: Normalization and identify DE genes, compute RIF scores and identify k.top RIF scores, and
#     using 'PCIT' package to find significant interaction between RFs and plot the network.
#
#             Input: 1- ExprData, a dataframe of microarray counts or expressions in two conditions. 
#                    2- Method, different correlation metrics to compute IRF score. Options: Pearson, Spearman, Sobolev, and Fisher. 
#                    3- k_top, k top RFs with the best absolute RIF scores. 
#		     4- adjp, adjusted p-value to identify most significant DE genes. 
#		     5- study, Options:prostate and NULL. prostate (The data we used in our study)
#
#             Output: A list include four objects:
#                     1- DEGs, A dataframe of most significant DE genes.
#                     2- RIF.Scores, A dataframe of all RFs and thier RIF scores.
#                     3- k.Top.RIF, A dataframe of k the best RIF scores.
#                     4- Plot, A network graph that shows most significant interaction between RFs.
#    
#
########################################################################################################################################################################

options(stringsAsFactors = FALSE)
library(edgeR)
library(data.table)
library(PCIT)
library(igraph)

####################################################################
#
# Edit below values and addresses accordingly
#
####################################################################

N.Sample.columns = 43
Len.Group1 = 21
Len.Group2 = 21

data_CAT <- as.data.frame(fread('/home/rezvane/IRF_LINUX/Data/genes.csv'))
All.RFs <- read.table('/home/rezvane/IRF_LINUX/Data/All_Regulator_Factors.txt', sep='\t', head=TRUE, stringsAsFactors=FALSE)
Prostate.IRF.Scores <- read.table('/home/rezvane/test_2019/IRF_Scores_Prostate.txt', sep='\t', head=TRUE, stringsAsFactors=FALSE)
ExprData <- '/home/rezvane/IRF_LINUX/Data/ALL_DATA.txt'


TFs.Ensembl.ID <- All.RFs[(All.RFs$Regulatory.Type == 'TF'),1]
EFs.Ensembl.ID <- All.RFs[(All.RFs$Regulatory.Type == 'EF'),1]
lncRNAs.Ensembl.ID <- All.RFs[(All.RFs$Regulatory.Type == 'lncRNA'),1]

TFs.Gene.ID <- All.RFs[(All.RFs$Regulatory.Type == 'TF'),2]
EFs.Gene.ID <- All.RFs[(All.RFs$Regulatory.Type == 'EF'),2]
lncRNAs.Gene.ID <- All.RFs[(All.RFs$Regulatory.Type == 'lncRNA'),2]

####################################################################
#
# Step 1:
#	Normalization and identify DE genes
#
####################################################################

edgeR_Normalization_DEGs <- function(ExprData){

	raw.data <- read.delim(ExprData)
	d <- raw.data[,2:N.Sample.columns]  
	rownames(d) <- raw.data[,1]
	group <- c(rep('A',Len.Group1),rep('B',Len.Group2)) 
	d <- DGEList(counts = d, group=group)
	nc <- cpm(d, normalized.lib.sizes = FALSE)

	d = calcNormFactors(d)
	d = estimateCommonDisp(d, verbose = TRUE)
	et <- exactTest(d)
	er <- data.frame(et$table)
	FDR <- p.adjust(et$table$PValue, method='BH')
	rd <- cbind(er,FDR)
	rd <- rd[order(rd$FDR),]

	return(list(NormalizedExpr=nc, DEGs=rd))
	}

####################################################################
#
# Step 2:
#	Computate RIF scores and identify k.top RIF scores
#
####################################################################
Sobolev_func <- function(x, y) {

	x1 <- x**2 / (sum(x**2))
	y1 <- y**2 / (sum(y**2))
	z1 <- x1-y1
	FT <- fft(z1)
	w <- 2*pi*(1:length(FT))/(length(FT))
	s <- sum((1+w)*abs(FT)**2)**(1/2)
	if(!is.nan(s)) return(s)
	else return(100000)
	}

Fisher_func <- function(x, y) {

	x1 <- x**2 / (sum(x**2))
	y1 <- y**2 / (sum(y**2))
	t <- x1 * y1
	s <- acos(sum(sqrt(t)))
	if(!is.nan(s)) return(s)
	else return(100000)
	}

RIF1 <- function(e1, e2, r1, r2, l){

	preRIF <- (1/2)*(e1^2-e2^2)*(r1-r2)^2
	preRIF <- na.omit(preRIF)
	preRIF <- apply(preRIF,2,sum)
	preRIFdividlength <- preRIF / l
	preRIFdividlength[is.na(preRIFdividlength)] <- 0
	return(preRIFdividlength)
	}

RIF2 <- function(e1, e2, r1, r2, l){

	preRIF <- (e1*r1)^2-(e2*r2)^2
	preRIF <- na.omit(preRIF)
	preRIF <- apply(preRIF,2,sum)
	preRIFdividlength <- preRIF /l 
	preRIFdividlength[is.na(preRIFdividlength)] <- 0
	return(preRIFdividlength)
	}

RIF <- function(Exprs, Exprs1, Exprs2, Regulatory.Type, Method, adjp){

	if(Regulatory.Type == 'TF') RFs <- TFs.Ensembl.ID
	if(Regulatory.Type == 'EF') RFs <- EFs.Ensembl.ID
	if(Regulatory.Type == 'lncRNA') RFs <- lncRNAs.Ensembl.ID

	RFs <- data.frame(RF = unique(as.data.frame(RFs)))
	Reg.Exprs1 <- unique(merge(Exprs1, RFs, by.x='row.names', by.y='RFs'))
	rownames(Reg.Exprs1) <- Reg.Exprs1[,1]
	RIF.Reg <- data.frame(RF = Reg.Exprs1[,'Row.names'], RIF1=FALSE, RIF2=FALSE)

	Reg.Exprs1 <- Reg.Exprs1[,2:ncol(Reg.Exprs1)]
	Reg.Exprs2 <- unique(merge(Exprs2, RFs, by.x="row.names", by.y='RFs'))
	rownames(Reg.Exprs2) <- Reg.Exprs2[,1]
	Reg.Exprs2 <- Reg.Exprs2[,2:ncol(Reg.Exprs2)]

	cutoff <- min(500, nrow(Exprs1))
	ncol.DERes1 <- ncol(Reg.Exprs1)
	DERes.Exprs1 <- Exprs1[c(1:cutoff), 1:ncol.DERes1]
	mean.DERes.Exprs1 <- apply(DERes.Exprs1, 1, mean)
	tmp1 <- names(mean.DERes.Exprs1)
	mean.DERes.Exprs1.dataframe <- data.frame(tmp1, mean.DERes.Exprs1)

	ncol.DERes2 <- ncol(Reg.Exprs2)
	DERes.Exprs2 <- Exprs2[c(1:cutoff), 1:ncol.DERes2]
	mean.DERes.Exprs2 <- apply(DERes.Exprs2, 1, mean)
	tmp2 <- names(mean.DERes.Exprs2)
	mean.DERes.Exprs2.dataframe <- data.frame(tmp2, mean.DERes.Exprs2)

	l <- nrow(mean.DERes.Exprs1.dataframe)

	nReg.Exprs <- nrow(Reg.Exprs1)
	e1 <- mean.DERes.Exprs1.dataframe
	e1 <- e1[,'mean.DERes.Exprs1']
	e2 <- mean.DERes.Exprs2.dataframe
	e2 <- e2[,'mean.DERes.Exprs2']
	rownames(RIF.Reg)<-RIF.Reg[,1]

	for( i in 1:nReg.Exprs){
		n <- nrow(DERes.Exprs1)
		m <- nrow(DERes.Exprs1[rownames(RIF.Reg)[i],])

		if(Method == 'Pearson'){
			r1 <- cor(t(DERes.Exprs1),t(Reg.Exprs1[rownames(RIF.Reg)[i],]),method="pearson",use="pairwise.complete.obs")
			r2 <- cor(t(DERes.Exprs2),t(Reg.Exprs2[rownames(RIF.Reg)[i],]),method="pearson",use="pairwise.complete.obs")
			RIF.Reg[i,2] <- RIF1(e1, e2, r1, r2, l)
			RIF.Reg[i,3] <- RIF2(e1, e2, r1, r2, l)

		}
		if(Method == 'Spearman'){
			r1 <- cor(t(DERes.Exprs1),t(Reg.Exprs1[rownames(RIF.Reg)[i],]),method="spearman",use="pairwise.complete.obs")
			r2 <- cor(t(DERes.Exprs2),t(Reg.Exprs2[rownames(RIF.Reg)[i],]),method="spearman",use="pairwise.complete.obs")
			RIF.Reg[i,2] <- RIF1(e1, e2, r1, r2, l)
			RIF.Reg[i,3] <- RIF2(e1, e2, r1, r2, l)
		}
		if(Method == 'Sobolev'){
			Sobolev1 <- sapply(1:n, function(x) sapply(1:m, function(y) Sobolev_func(t(DERes.Exprs1[x,]),t(Reg.Exprs1[rownames(RIF.Reg)[i],][y,]))))
			Sobolev2 <- sapply(1:n, function(x) sapply(1:m, function(y) Sobolev_func(t(DERes.Exprs2[x,]),t(Reg.Exprs2[rownames(RIF.Reg)[i],][y,]))))
			r1 <-data.frame(1/(1+Sobolev1))
			r2 <-data.frame(1/(1+Sobolev2))

			RIF.Reg[i,2] <- RIF1(e1, e2, r1, r2, l)
			RIF.Reg[i,3] <- RIF2(e1, e2, r1, r2, l)
		}
		if(Method == 'Fisher'){
			Fisher1 <- sapply(1:n, function(x) sapply(1:m, function(y) Fisher_func(t(DERes.Exprs1[x,]),t(Reg.Exprs1[rownames(RIF.Reg)[i],][y,]))))
			Fisher2 <- sapply(1:n, function(x) sapply(1:m, function(y) Fisher_func(t(DERes.Exprs2[x,]),t(Reg.Exprs2[rownames(RIF.Reg)[i],][y,]))))
			r1 <- data.frame(1/(1+Fisher1))
			r2 <- data.frame(1/(1+Fisher2))
			RIF.Reg[i,2] <- RIF1(e1, e2, r1, r2, l)
			RIF.Reg[i,3] <- RIF2(e1, e2, r1, r2, l)
		}
	}
	df <- RIF.Reg[,-c(1)]
	cols <- colnames(df)
	df[cols] <- round(scale(df[cols]), 2)
	RIF.Reg <- cbind(RFs=RIF.Reg[,1], df)
	rownames(RIF.Reg) <- c()
	return(RIF.Reg)
	}

k_Top_RF <- function(RIF.Scores, k, x, y){

	RIF.x <- RIF.Scores[rev(order(abs(RIF.Scores[,x]))), ]
	RIF.y <- RIF.Scores[rev(order(abs(RIF.Scores[,y]))), ]
	RIF.xy_k <- rbind(RIF.x[c(1:k), c(1,2,3)], RIF.y[c(1:k), c(1,2,3)])
	return(unique(RIF.xy_k))
	}

get_k_top_all_RFs <- function(RIF.Scores, k, x, y){

	RIF.TF <- RIF.Scores[(RIF.Scores$Regulatory.Type == 'TF'),]
	Top.TFs <- k_Top_RF(RIF.TF, k, x, y)
	Top.TFs <- data.frame(cbind(Top.TFs, Gene.type=rep(1, nrow(Top.TFs))))

	RIF.EF <- RIF.Scores[(RIF.Scores$Regulatory.Type == 'EF'),]
	RIF.EF <- RIF.EF[!(RIF.EF[,1] %in% Top.TFs[,1]),]
	Top.EFs <- k_Top_RF(RIF.EF, k, x, y)
	Top.EFs <- data.frame(cbind(Top.EFs, Gene.type=rep(2, nrow(Top.EFs))))

	RIF.lncRNA <- RIF.Scores[(RIF.Scores$Regulatory.Type == 'lncRNA'),]
	lncrna_type <- c("lncRNA, divergent", "lncRNA, antisense", "lncRNA, intergenic", "lncRNA, sense intronic")
	lncrna_data <- data_CAT[(data_CAT[,4] %in% lncrna_type),1]
	lncrna_data <- substr(lncrna_data,1,15)
	RIF.lncRNA <- RIF.lncRNA[(RIF.lncRNA[,1] %in% lncrna_data),]
	Top.lncRNAs <- k_Top_RF(RIF.lncRNA, k, x, y)
	Top.lncRNAs <- data.frame(cbind(Top.lncRNAs, Gene.type=rep(3, nrow(Top.lncRNAs))))
	Top.lncRNAs$Gene.ID <- substr(Top.lncRNAs$Gene.ID,10,15)

	k.Top.All.RFs <- rbind(Top.TFs, Top.EFs, Top.lncRNAs)
	return(k.Top.All.RFs)
	}

####################################################################
#
# Step 3:
#	using 'PCIT' package to find significant interaction between
#	RFs and plot the network
#
####################################################################

get_links <- function(ExprsN, k_top_all_RFs){

	Top.RFs.Exprs <- ExprsN[(rownames(ExprsN) %in% k_top_all_RFs$RFs),]

	cor.Top.RFs.Exprs <- cor(t(Top.RFs.Exprs), t(Top.RFs.Exprs), method='pearson', use='pairwise.complete.obs')
	diag(cor.Top.RFs.Exprs) <- 0
	W.Edge <- getEdgeList(cor.Top.RFs.Exprs)
	Top.links <- W.Edge[(abs(W.Edge[,3]) > .90), ]
	link1 <- as.character(Top.links[,1])
	link2 <- as.character(Top.links[,2])
	uniq.links <- unique(c(link1, link2))

	nodes <- k_top_all_RFs[(k_top_all_RFs[,1] %in% uniq.links),]
	links <- aggregate(Top.links[,3], Top.links[,-3], sum)
	links <- links[order(links$From, links$To),]
	colnames(links)[3] <- 'Weight'
	rownames(links) <- NULL
	
	return(list(nodes = nodes, links = links))
	}

get_net <- function(nodes, links){

	TF_Nodes <- nodes[(nodes[,3]=='TF'),1]
	EF_Nodes <- nodes[(nodes[,3]=='EF'),1]
	lncRNA_Nodes <- nodes[(nodes[,3]=='lncRNA'),1]
	net <- graph_from_data_frame(d=links, vertices=nodes, directed=F) 
	colrs <- c('blue', 'red', 'black')
	V(net)$color <- colrs[as.numeric(V(net)$Gene.type)]

	E(net)$color <- apply(as.data.frame(get.edgelist(net)), 1, 
               function(x) {
			if(x[1] %in% TF_Nodes & x[2] %in% TF_Nodes) return('lightblue')
			if(x[1] %in% EF_Nodes & x[2] %in% EF_Nodes) return('red')
			if(x[1] %in% lncRNA_Nodes & x[2] %in% lncRNA_Nodes) return('black')
			return('yellow')
			})
	return(net)
	}

####################################################################
#
# 
#	Main function to integrate three steps
#	
#
####################################################################

RIF_Analysis <- function(ExprData, Method, k_top, adjp, study){
		
	ExprsND <- edgeR_Normalization_DEGs(ExprData)

	ExprsN <- as.data.frame(ExprsND$NormalizedExpr)
	rownames(ExprsN) <- substring(rownames(ExprsN), 1, 15)

	DEGs <- as.data.frame(ExprsND$DEGs)
	rownames(DEGs) <- substring(rownames(DEGs), 1, 15)
	DEGs <- DEGs[(DEGs$FDR < adjp),]

	if(!is.null(study)) {
				if(Method == 'Pearson') RIF.All <- Prostate.IRF.Scores[,c(1:3, 4,8)]
				if(Method == 'Spearman') RIF.All <- Prostate.IRF.Scores[,c(1:3, 5,9)]
				if(Method == 'Sobolev') RIF.All <- Prostate.IRF.Scores[,c(1:3, 6,10)]
				if(Method == 'Fisher') RIF.All <- Prostate.IRF.Scores[,c(1:3, 7,11)]
				RIF.All <- unique(RIF.All)
				}
	else{
		Exprs <- ExprsN[(rownames(ExprsN) %in% rownames(DEGs)),]
		Exprs1 <- Exprs[,1:Len.Group1]
		Exprs2 <- Exprs[,(Len.Group1+1):(N.Sample.columns-1)]

		RIF.TF <- RIF(Exprs, Exprs1, Exprs2, Regulatory.Type='TF', Method, adjp)

		RIF.EF <- RIF(Exprs, Exprs1, Exprs2, Regulatory.Type='EF', Method, adjp)
		RIF.EF <- RIF.EF[!(RIF.EF[,1] %in% RIF.TF[,1]),]
		RIF.lncRNA <- RIF(Exprs, Exprs1, Exprs2, Regulatory.Type='lncRNA', Method, adjp)

		RIF.All <- rbind(RIF.TF, RIF.EF, RIF.lncRNA)
		RIF.All <- merge(RIF.All, All.RFs,  by.x=c('RFs'), by.y=c('Ensembl.ID'))
		RIF.All <- RIF.All[,c(1,4,5,2,3)]
		}
	k_top_all_RFs <- get_k_top_all_RFs(RIF.All, k_top, 4, 5)
	k_top_all_RFs <- k_top_all_RFs[order(k_top_all_RFs[,1] ),] 

	k_top_all_RFs_ <- merge(k_top_all_RFs, RIF.All, by.x='RFs', by.y='RFs')
	k_top_all_RFs_ <- k_top_all_RFs_[order(k_top_all_RFs_$Gene.type),]
	k_top_all_RFs_ <- k_top_all_RFs_[,c(1:3,7,8)]

	nodes_links <- get_links(ExprsN, k_top_all_RFs)
	nodes <- nodes_links$nodes
	links <- nodes_links$links
	net <- get_net(nodes, links)
	
	colrs <- c('blue', 'red', 'black')
	Plot <- plot(net, vertex.label = as.character(V(net)$Gene.ID), vertex.label.dist=.75, vertex.label.font=1, vertex.label.color=V(net)$color, vertex.label.cex=.6,
	edge.arrow.size=.1, vertex.size=4, main=sprintf('The RIF Scores Computed by %s Metric', Method))
	legend(x=-1.5, y=-1, c('TF','EF', 'lncRNA'), pch=21,col="#777777", pt.bg=colrs, pt.cex=1, cex=.6, bty='n', ncol=1)
	
	return(list(DEGs=DEGs, RIF.Scores=RIF.All, k.Top.RIF.Scores=k_top_all_RFs_, Plot=Plot))
	}

# Example
Res <- RIF_Analysis(ExprData, Method='Fisher', k_top=10, adjp=0.01, study='prostate')
#Res <- RIF_Analysis(ExprData, Method='Fisher', k_top=10, adjp=0.01, study=NULL)
head(Res$DEGs)
head(Res$RIF.Scores)
head(Res$k.Top.RIF.Scores)
Res$Plot
