require(fastcluster);require(RColorBrewer);require(gplots)
#davidSimple outputs a crude matrix of clustered, matching, unfiltered enrichment results  
# in the process, the function saves a nice heatmap that can be used for identifying important biological activity in your experiment 
# delist is a list of gene lists of DE genes from experimental conditions; it must have two components for each sample: upregulated genes placed first, and downregulated genes placed second
# genesets is the functional enrichment data, in the form of a list of lists of gene sets
# universe is the set of all possible genes. 
# overlaps is an optional time saver if you need to run the code multiple times; you can run 'dave1<-davidSimple(delist,genesets,universe)' the first time, then
# 'davidSimple(delist,genesets,universe,overlaps=dave1$overlaps)' on subsequent runs.
# This allows you to reuse the results of the statistical tests for enrichment.
# clustThresh is the threshold for clustering the "DAVID" output results, using gene overlaps:
#		higher clustThresh; more but smaller clusters
#		lower clustThresh; fewer but larger clusters
# matchThresh is the threshold for matching clustering results across multiple conditions, works the same way
davidSimple<-function(delist,genesets,universe,overlaps=NULL,clustThresh=.3,matchThresh=.3,doClusterMatching=T,assmat=NULL,multiTestCorrect=F,doClustering=T){
	genesets<-checkInput(delist,genesets,universe)
	if(is.null(genesets)) return(NULL)
	if(is.null(overlaps)){
		cat("getting core stats...\n")
		overlaps<-coreStats(delist,genesets,universe,pvalcutoff=.1,multiTestCorrect=multiTestCorrect)
	} else if(!is.null(multiTestCorrect))
		cat("multiTestCorrect status already determined by preloaded overlaps.\n")
	overlapGenes<-overlaps$genes
	testRes<-overlaps$values
	clusterGenes<-vector("list",length=length(delist)) # we return both clustered results, and
	clusterVals<-vector("list",length=length(delist)) 
	names(clusterGenes)<-names(delist)
	names(clusterVals)<-names(delist)
	nonclusterGenes<-vector("list",length=length(delist)) # nonclustered results
	nonclusterVals<-vector("list",length=length(delist))
	for(i in 1:length(delist)){
		cat("processing condition ",i,"\n")
		overlap<-overlapGenes[[i]]
		testR<-testRes[[i]]
		if(length(overlap)>0){
			if(doClustering){
				cl<-kappaCluster(overlap,as.numeric(length(universe)),clustThresh)
				clusterInds<-cl$clustered;nonclustInds<-cl$nonclustered
			}else {
				nonclustInds<-1:length(overlap)
				clusterInds<-numeric()
			}
			####### deal with nonclustered results first
			nonclusterGenes[[i]]<-compileNonClusteredDataObject(nonclustInds,overlap)
			nonclusterVals[[i]]<-testR[nonclustInds]
			#######  now deal with the clustered results
			if(length(clusterInds)>0){
				temp<-compileClusteredDataObject(clusterInds,overlap,testR)
				clusterGenes[[i]]<-temp$genes
				clusterVals[[i]]<-temp$values
			} else{ 
				clusterGenes[[i]]<-list() # if there are no clusters, send an empty list
				clusterVals[[i]]<-list()
			}
		} else{
			nonclusterGenes[[i]]<-list()
			nonclusterVals[[i]]<-list()
			clusterGenes[[i]]<-list() 
			clusterVals[[i]]<-list()
		}
	}
	clustering<-list(clustered=list(genes=clusterGenes,values=clusterVals),nonclustered=list(genes=nonclusterGenes,values=nonclusterVals)) 
	count<-unlist(lapply(clusterGenes,length))
	if(length(which(count>0))>1){
		if(doClusterMatching==F) return(list(clustering=clustering,overlaps=overlaps))
		matching<-processClusteredData(clustering$clustered,genesets,matchThresh,assmat)
		return(list(matching=matching,clustering=clustering,overlaps=overlaps,clustThresh=clustThresh,matchThresh=matchThresh))
	}
	return(list(matching=list(),clustering=clustering,overlaps=overlaps,clustThresh=clustThresh,matchThresh=matchThresh))
}

checkInput<-function(delist,genesets,universe){
	setGenes<-unique(unlist(genesets))
	myGenes<-unique(unlist(delist))
	temp<-length(intersect(universe,myGenes))
	if(temp!=length(myGenes)){
		cat("warning: only",temp/length(myGenes)*100,"% of DE genes agree with universe.\n")
		cat("Please adjust gene list and/or universe.\n")
		return(NULL)
	}
	temp<-length(intersect(universe,setGenes))
	if(temp!=length(setGenes)){
		cat("Only",temp/length(setGenes)*100,"% of set genes agree with universe. Adjusting gene sets...\n")
		newGeneSets<-lapply(genesets,function(x) intersect(x,universe))
		names(newGeneSets)<-names(genesets)
	}else newGeneSets<-genesets
	setGenes<-unique(unlist(newGeneSets))
	frac<-length(intersect(myGenes,setGenes))/length(myGenes)
	cat(frac*100,"percent of DE genes exist in gene sets. Continue?\n")
	if(readline()=="n") return(NULL)
	return(newGeneSets)
}

# this function combines the results from several rounds of enrichment result clustering, and performs the cross-experiment matching to produce
# a functional enrichment matrix. Yields same data structure as davidSimple(), such that the output is ready for finish().
# study names is a character vector consisting of the names of enrichment runs from the experiments to combine. each enrichment run needs to 
# be the output of davidSimple() when doClusterMatching is set to FALSE.
# genesets is the same as the genesets input for davidSimple()
# matchThresh is the same as the matchThresh for davidSimple()
combineStudies<-function(studyNames,genesets,matchThresh=.25,assmat=NULL){
	first<-get(studyNames[1])
	clustering<-first$clustering
	overlaps<-first$overlaps
	cg<-clustering$clustered$genes
	cv<-clustering$clustered$values
	ncg<-clustering$nonclustered$genes
	ncv<-clustering$nonclustered$values
	og<-overlaps$genes
	ov<-overlaps$values
	if(length(studyNames)>1)
		for(i in 2:length(studyNames)){
			this<-get(studyNames[i])
			clustering<-this$clustering
			overlaps<-this$overlaps
			cg<-c(cg,clustering$clustered$genes)
			cv<-c(cv,clustering$clustered$values)
			ncg<-c(ncg,clustering$nonclustered$genes)
			ncv<-c(ncv,clustering$nonclustered$values)
			og<-c(og,overlaps$genes)
			ov<-c(ov,overlaps$values)
		}
	allDave_clustering<-list(clustered=list(genes=cg,values=cv),nonclustered=list(genes=ncg,values=ncv))
	allDave_overlaps<-list(genes=og,values=ov)
	matching<-processClusteredData(allDave_clustering$clustered,genesets,matchThresh,assmat)
	return(list(matching=matching,clustering=allDave_clustering,overlaps=allDave_overlaps,clustThresh=-1,matchThresh=matchThresh))
}

# conduct DAVID style fisher test on all DE genes lists and gene sets
# returns hit genes and p-values
coreStats<-function(delist,genesets,universe,pvalcutoff=.1,multiTestCorrect){
	pvalcutoff <- -log10(pvalcutoff)
	if(is.null(multiTestCorrect)){
		cat("Multiple test correction = ",multiTestCorrect,".\n",sep="")
		multiTestCorrect<-FALSE
	}
	overlaps<-vector("list",length=length(delist)) # the data structure is set up as two data objects (genes and scores) combined into a list of 2.
	testRes<-vector("list",length=length(delist))
	names(overlaps)<-names(delist)
	for(i in 1:length(overlaps)){
		cat("\tcondition",i,"\n")
		overlaps[[i]]<-vector("list",length=length(genesets))
		test<-numeric()
		for(j in 1:length(genesets)){
			#cat("\rtest",j," ")
			test[j]<-fish2(delist[[i]],genesets[[j]],universe,ease=T,quiet=T) # here's the score storage, 
			ol<-overlaps[[i]][[j]]<-intersect(delist[[i]],genesets[[j]]) # and we store the gene hits
			if(length(ol)>0) overlaps[[i]][[j]]<-ol
			else overlaps[[i]][[j]]<-character()
		}
		if(multiTestCorrect) test<-p.adjust(test,method="BH")
		test[which(test < 1e-320)] <- 1e-320 # less than this and the log goes to inf
		test<--log10(test)
		names(overlaps[[i]])<-names(genesets)
		overlaps[[i]]<-overlaps[[i]][which(test>pvalcutoff)]
		test<-test[which(test>pvalcutoff)] # get rid of everything that doesn't meet the threshold
		testRes[[i]]<-test
	}
	return(list(genes=overlaps,values=testRes,multiTestCorrect=multiTestCorrect))
}

# calculate exact fisher test for a list of important genes, a relevant geneset, and a gene universe
# can also enter a confusion matrix
# ease is the DAVID EASE adjustment to penalize hits with small number of genes
fish2<-function(genelist=NULL,geneset=NULL,total=NULL,mat=NULL,quiet=F,ease=T){
	if(is.null(mat)){
		intersection<-length(intersect(genelist,geneset))
		if(!quiet) cat("hits:",intersection,"\n")
		#listnotset<-length(setdiff(genelist,geneset))
		listnotset<-length(genelist)-intersection
		#setnotlist<-length(setdiff(geneset,genelist))
		setnotlist<-length(geneset)-intersection
		#tn<-length(setdiff(total,union(genelist,geneset)))
		tn<-length(total)-length(unique(c(genelist,geneset)))
		mat<-matrix(c(intersection,listnotset,setnotlist,tn),nrow=2)
	}
	if(ease){ # this is the EASE score adjustment that DAVID uses, penalizes gene sets with very few matches
		if(mat[1,1]>0) mat[1,1]<-mat[1,1]-1 #intersection<-intersection-1
		if(mat[1,2]>0) mat[1,2]<-mat[1,2]+1 #setnotlist<-setnotlist+1
	}
	test<-fisher.test(mat,alt="greater")
	if(!quiet) cat("fisher test: p-value = ",test$p.value,", odds ratio = ",test$estimate,"\n")
	return(test$p.value)
}

# kappa clustering uses the kappa statistic, and merges by making a union list of the two merging nodes, thus giving each new merged node a new kappa score with all other nodes
kappaCluster<-function(hitlist,universeLen,thresh,assmat=NULL){
	if(length(hitlist)==1) return(list(list(1),1))
	record<-list() # this preserves a record of the merges that occurred 
	for(i in 1:length(hitlist)) record[[i]]<-i # to start, just start with the index number of each nonclustered hit as each list item  
	originalLen<-length(hitlist) # this preserves the structure before we do any merging
	cat("\tgetting matrix...\n")
	if(is.null(assmat))
		assmat<-getKappaMat(hitlist,universeLen) # get a matrix of the similarities between hits
	cat("\titerating...\n")
	while(TRUE){
		mx<-max(assmat)
		if(mx<=thresh) break # we're done if the max similarity doesn't exceed the threshold
		len<-length(hitlist)
		merge<-which(assmat==mx)[1] # find which pair has the max similarity (take the first if there are ties)
		ceil<-ceiling(merge/len) # identify the pair associated with the matrix index
		if(ceil==merge/len){
			pair<-c(ceil,len)
		}else pair<-c(ceil,merge %% len)
		pair<-sort(pair)
		record[[pair[1]]]<-c(record[[pair[1]]],record[[pair[2]]]) # now update the record and merge the pair;
		record[[pair[2]]]<-NULL # we're preserving the original hit indices in record
		hitlist[[pair[1]]]<-unique(c(hitlist[[pair[1]]],hitlist[[pair[2]]])) # but in hitlist the original structure is destroyed
		hitlist[[pair[2]]]<-NULL
		assmat<-assmat[-pair[2],] # remove the absorbed node's spot in the distance matrix
		if(is.null(dim(assmat))) break
		assmat<-assmat[,-pair[2]]
		updates<-numeric(dim(assmat)[1]) # now update the matrix so it contains accurate scores based on the new node's content and its relationship to all other nodes.
		for(i in 1:length(updates)){
			if(i==pair[1]) updates[i]<-0
			else updates[i]<-kappa.stat(hitlist[[pair[1]]],hitlist[[i]],universeLen)
		}
		assmat[pair[1],]<-updates
		assmat[,pair[1]]<-updates
	}
	len<-lapply(record,length) # we're only interested in hits that got merged
	record<-record[which(len>1)]
	nonclustered<-setdiff(1:originalLen,unlist(record))# identify stuff that didn't get clustered
	return(list(clustered=record,nonclustered=nonclustered,assmat=assmat))
}

# make a matrix of kappa statistics based on all pairs in a list of gene annotations (functional categories)
getKappaMat<-function(hits,universeLen){
	cat(length(hits),"hits\n")
	perms<-combn(length(hits),2) # generate all pair combos
	data<-numeric(dim(perms)[2])
	cat("filling",length(data),"values...\n")
	mat<-matrix(0,nrow=length(hits),ncol=length(hits))
	for(i in 1:length(data)) # first store the data in a vector
		data[i]<-kappa.stat(hits[[perms[1,i]]],hits[[perms[2,i]]],universeLen)
	len<-length(hits)
	spot<-1
	for(i in 2:len){ # fill in a matrix with the results, column by column(just fill the half-matrix)
		mat[i:len,(i-1)]<-data[spot:(spot+len-i)]
		spot<-spot+len-i+1
	}
	#mat<-as.matrix(as.dist(mat)) 
	return(mat)
}

# calculate Kappa statistic for similarity between 2 raters, or in this case, two gene annotations using their gene overlap
kappa.stat<-function(list1,list2,all){
	posagree<- length(intersect(list1,list2))
	negagree<- all-length(unique(c(list1,list2)))
	O<-(posagree+negagree)/all
	len1<-length(list1);len2<-length(list2)
	E<-(len1*len2+(all-len1)*(all-len2))/all^2
	return((O-E)/(1-E))
}

# here we transform hits that did not get clustering into a similar format as the clustered hits
compileNonClusteredDataObject<-function(nonclustInds,overlap){
	if(length(nonclustInds)==0) return(NULL)
	nonclust<-vector("list",length=length(nonclustInds))
	for(j in 1:length(nonclustInds)) # attach hit genes and scores to nonclustered results
		nonclust[[j]]<-overlap[[nonclustInds[j]]]
	names(nonclust)<-names(overlap)[nonclustInds]
	return(nonclust)
}

# the clustering process keeps a record of what nodes got clustered, now we need to transform that record into a data object that contains all the necessary info
compileClusteredDataObject<-function(clusterInds,overlap,testR){
	clustGenes<-vector("list",length=length(clusterInds))
	clustVals<-vector("list",length=length(clusterInds))
	#clustScores<-numeric(length(clusterInds)) # score each cluster based on the mean of all the cluster member scores
	for(j in 1:length(clustGenes)){
		clusterInd<-clusterInds[[j]]
		clustGene<-vector("list",length(clusterInd))
		for(k in 1:length(clusterInd))
			clustGene[[k]]<-overlap[[clusterInd[k]]] # assign hit genes 
		clustVals[[j]]<-testR[clusterInd]
		#clustScores[j]<-mean(clustVals[[j]]) # take the mean and store as the cluster score
		names(clustGene)<-names(overlap)[clusterInd] # name the cluster members with their respective names
		clustGenes[[j]]<-clustGene
	}
	#names(clustVals)<-clustScores # assign cluster scores
	return(list(genes=clustGenes,values=clustVals))
}

# here's where we do the work of clustering items across conditions and merging them into single rows
processClusteredData<-function(clustered,genesets,matchThresh,assmat){
	universe<-length(genesets)
	clusterGenes<-clustered$genes
	lens<-unlist(lapply(clusterGenes,length))
	hits<-vector("list",length=sum(lens))
	counter<-0
	cat("initializing...\n") # first thing to do is to transform the cluster data into a flat structure that can be accessed with a single index
	# this will contain data down to the function level, not all the way to the gene level
	for(i in 1:length(clusterGenes)){
		if(length(clusterGenes[[i]])>0)
			for(j in 1:length(clusterGenes[[i]])){
				counter<-counter+1
				hits[[counter]]<-names(clusterGenes[[i]][[j]]) # functions and values are stored in separate structures, but accessible with the same index
			}			
	}
	kc<-kappaCluster(hits,universe,matchThresh)
	clus<-kc$clustered;nonclus<-kc$nonclustered
	nonclus<-as.list(nonclus)
	all<-c(clus,nonclus)
	#this thing takes the cluster assignments and converts them to the same form that cut() does,
	#which is what we were using before. It's an assignment to a cluster for each elements
	#The outer 'loop' feeds in all the different elements, while the inner loop
	#checks each cluster to see if it contains that element, so that at the end
	#each element has a cluster assignment
	clustassign<-sapply(1:length(hits), function(hit) 
											which(unlist(lapply(all, function(elements)
																		hit  %in% elements
																)
														 )
												  )
						)
	return(list(assmat=assmat,cut=clustassign,hits=hits))
}

# dave is output from davidSimple() 
# sigThresh is the threshold for including results in the final heatmap in -log10(pval) format, defaults to 2 (equivalent to .01). At least
#       one value in a row must meet the threshold for inclusion
# maxShow is an alternative filtering strategy to sigThresh. It is simply the number of top-ranking rows to include, with ranks calculated 
# 	    using the sum of -log10 pvalues for each row.
# mergeThresh is the threshold of when to acknowledge a difference between up and down-regulated functions/clusters. If a function/cluster 
#       is only above the threshold for one direction, the other direction is ignored. If both are above the threshold, then both are acknowledged. 
# bidirectional refers to whether up- and down-regulated genes were included separately in the input data. If bidirectional is TRUE, columns will 
#       be treated with the assumption that the columns are in up/down pairs and expression directions will be represented as blue/red. If 
#       bidirectional is false, significance is colored by intensity only and no assumptions are made regarding column pairing.
# output is the name of the file where the breakdown of the heatmap rows is stored, defaults to clusterInfo.txt
# makeGeneRefFile determines whether to make gene reference files; i.e. files that show what genes were identified from each functional enrichment
#       hit for each condition
# geneInfoMap is a table of gene annotations that can be optionally included to provide more information in the gene reference files
# mapMatchCol is the column in geneInfoMap that matches the IDs used in the analysis
# 	it can be used in conjunction with the clusterInfo.txt file to drill down into what genes are behind each hit or set of hits.
# selcol is the column where the new info of interest is
# seps is the locations of division lines in the heatmap 
finish<-function(dave,sigThresh=2,maxShow=NULL,mergeThresh=1.30103,bidirectional=T,output="clusterInfo.txt",
		makeGeneRefFile=T,geneInfoMap=NULL,mapMatchCol=NULL,selcol=NULL,seps=NULL,mx=6,hmw=1000,hmh=2000){
	clustering<-dave$clustering
	matching<-dave$matching
	overlaps<-dave$overlaps
	if(length(matching)>0){
		clusMatrix<-buildHMmatrix2(matching,overlaps,sigThresh)
		clus<-clusMatrix$FEmatrix
		allTerms<-clusMatrix$allTerms
	} else{	
		clus<-matrix(nrow=0,ncol=length(overlaps$genes))
		colnames(clus)<-names(overlaps$genes)
		allTerms<-list()	
	}
	nonclus<-processNonClusteredData2(clustering$nonclustered,overlaps,allTerms)
	allTerms<-nonclus$allTerms;nonclus<-nonclus$ma
	mat<-cleanUp(clus,nonclus,allTerms)
	temp<-apply(mat,1,function(x) sum(abs(x)))
	temp<--sort(-temp)
	temp<-temp/max(temp)
	if(!is.null(maxShow))
		select<-names(temp[1:maxShow])
	else{
		temp<-apply(mat,1,max)
		select<-names(temp)[which(temp>sigThresh)]
	}
	select<-which(rownames(mat) %in% select) # this prevents the potential problem with matching only the first instance when duplicates exist
	selmat<-mat[select,]
#	return (selmat)
	hm<-process4BidirectionalHM(selmat,mergeThresh,bidirectional,seps,mx=mx,width=hmw,height=hmh) # processes data for useful display in bidirectional 
																	   # heatmap and returns the order of rows in the heatmap
	
	ord<-hm$ord;selmat<-hm$mergedMatrix
	ord<-ord[length(ord):1] # now we use that order to format a cluster reference file that's in the same order as the heatmap rows
	selmat<-selmat[ord,]
	select<-select[ord]
	writeFinalClusters(allTerms,select,output)
	if(makeGeneRefFile) makeGeneReferenceFile(overlaps,map=geneInfoMap,matchcol=mapMatchCol,selcol=selcol)
	dave$fullMat<-mat;dave$finalMat<-selmat;dave$allTerms<-allTerms;dave$sigThresh=sigThresh
	return(dave)
}

cleanUp<-function(clus,nonclus,allTerms){
	singles<-which(unlist(lapply(allTerms,length))==1)
	singleTerms<-unlist(allTerms[singles])
	nonsingles<-setdiff(1:length(allTerms),singles)
	clusterTerms<-sapply(allTerms[nonsingles],"[",1:2)
	strike<-which(singleTerms %in% clusterTerms)
	mat<-rbind(clus,nonclus[-strike,])
	return(mat)
}

# buildMatrix takes the output from processClusteredData() and uses sigThresh to calculate values for heatmap cells, then builds the 
# 	heatmap matrix for the clustered items. 
# matching is the matching results for processClusteredData() consists of the association matrix (assmat), hits (flat list of all FE hits), 
#	and cut (cluster assignments for each FE cluster).
# clustered is the results of the initial FE clustering within each condition
# sigThresh is the limit for inclusion in the final FE heatmap 
buildHMmatrix2<-function(matching,overlaps,sigThresh){
	cut<-matching$cut
	hits<-matching$hits
	resultmat<-matrix(0,nrow=max(cut),ncol=length(overlaps$genes)) # here is where we store the results
	colnames(resultmat)<-names(overlaps$genes)
	allTerms<-list()
	for(i in 1:max(cut)){
		nameBin<-character()
		these<-which(cut==i)
		terms<-unique(unlist(hits[these]))
		for(j in 1:length(overlaps$genes)){
			tms<-names(overlaps$genes[[j]])
			scrs<-overlaps$values[[j]]
			sel<-which(tms %in% terms)
			theseTerms<-tms[sel]
			theseScores<-scrs[sel]
			good<-which(theseScores>=sigThresh)
			if(length(good)>0)
				resultmat[i,j]<-mean(theseScores[good])
			names(theseScores)<-theseTerms
			nameBin<-c(nameBin,theseScores)
		}
		nameBin<-nameBin[order(nameBin,decreasing=T)]
		nameBin<-unique(names(nameBin))
		allTerms[[i]]<-nameBin
	}
	rownames(resultmat)<-makeUpNames(allTerms,"topTerms",3)# get names for the rows based on terms incorporated in that row
	cat(dim(resultmat)[1]," collapsed lines.\n")
	return(list(FEmatrix=resultmat,allTerms=allTerms))
}

# here's where we get the non-clustered items across conditions and insert them appropriately into the matrix
processNonClusteredData2<-function(nonclustered,overlaps,allTerms){
	nms<-nonclustered$genes
	allnames<-character()
	for(i in 1:length(nms))
		allnames<-unique(c(allnames,names(nms[[i]])))
	terms<-overlaps$genes
	values<-overlaps$values
	mat<-matrix(0,nrow=length(allnames),ncol=length(terms))
	rownames(mat)<-allnames;colnames(mat)<-names(terms)
	for(i in 1:length(terms)){
		tms<-names(terms[[i]])
		scrs<-values[[i]]
		sel<-match(allnames,tms)
		theseTerms<-tms[sel]
		theseScores<-scrs[sel]
		mat[,i]<-theseScores
	}
	mat[which(is.na(mat))]<-0
	cat(dim(mat)[1],"non-clustered lines.\n")
	allTerms<-c(allTerms,as.list(rownames(mat)))
	return(list(mat=mat,allTerms=allTerms))
}

# generate rownames for collapsed data matrix based on terms incorporated in each row
# this can be done 2 ways: by simply making a name from the top-scoring 1,2, or 3 terms, or by using the 3 most 
# popular words that occur in the term names
# for the former, choose mode = "topTerms", for the latter, choose mode = "sorted", or both choose mode = "both"
makeUpNames<-function(allTerms,mode="topTerms",numTerms=3){
	dumb<-c("and","of","by","the","to","in","process","activity","pathway","response","signaling") # we don't want to include these in row names
																									#	if we're using the top terms method
	if(numTerms==1) topTerms<-sapply(allTerms,"[",1) # make names out of top scoring terms
	else if(numTerms==2 | numTerms==3){
		topTerms<-sapply(allTerms,"[",1:numTerms)
		topTerms<-apply(topTerms,2,paste,collapse="*")
		topTerms<-gsub("*NA","",topTerms,fixed=T)
	}
	else{ 
		cat("numTerms must be 1 2 or 3\n")
		return()
	}
	bank<-vector("list",length(allTerms))
	for(i in 1:length(bank))
		bank[[i]]<-unlist(strsplit(allTerms[[i]]," "))
	bank<-lapply(bank,tolower) # bank holds the words used in each cluster
	bank<-lapply(bank, function(x) unlist(strsplit(x,"-")))
	bank<-lapply(bank, function(x){  # take out the 'dumb' words
							dumbs<-which(x %in% dumb)
							if(length(dumbs)==0) return(x)
							return(x[-dumbs])
						}
	)
	sorted<- lapply(bank, function(x) { # grab the top 3 words and combine them to name the row
							temp<- names(sort(-table(x)))
							ind<-min(3,length(temp))
							return(paste(temp[1:ind],collapse=" "))
		}
	)
	if(mode=="both")
		return(paste(topTerms,"/",unlist(sorted),sep=""))
	if(mode=="topTerms") return(topTerms)
	return(unlist(sorted))
}

# this function transforms the matrix into a form that is nicely displayed in a heatmap
# downregulated conditions (even numbered columns) are made negative, 
# and dominating trends in one fold change direction are allowed to obliterate the other 
# direction, unless both meet the significance threshold of .05 (or 1.30103 -log10) in which case both are preserved
# this function returns the order used by the heatmap, so an appropriate reference file can be created that orders the 
# rows in the same way.
process4BidirectionalHM<-function(mat,mergeThresh,bidirectional,seps,width=1000,height=2000,output="heatmap.jpg",mx=0){
	rownames(mat)<-gsub("^ ","",rownames(mat))
	if(bidirectional){
		# make downregs negative
		mat[,seq(2,dim(mat)[2],2)]<-mat[,seq(2,dim(mat)[2],2)]*-1
		# merge up downs into matching columns
		for(i in seq(1,dim(mat)[2],2))
			for(j in 1:dim(mat)[1]){
				temp<-mat[j,i:(i+1)]
				if(min(abs(temp))<mergeThresh)
					mat[j,i:(i+1)]<-temp[which(abs(temp)==max(abs(temp)))]
			}
		if(is.null(colnames(mat)))
			colnames(mat)<-paste("condition ",seq(.5,dim(mat)[2]/2,.5),sep="")
		colnames(mat)<-gsub("_d$","",colnames(mat))
		colnames(mat)[seq(1,dim(mat)[2],2)]<-""
		if(mx==0) mx<-max(mat)
		br=seq(-mx,mx,2*mx/11);
		pal<-brewer.pal(11,"RdYlBu")[11:1]
	}
	else{
		if(is.null(colnames(mat)))
			colnames(mat)<-paste("condition ",1:dim(mat)[2],sep="")
		pallete<-colorRampPalette(brewer.pal(9,"Blues"))
		pal<-pallete(15)
		br<-1:16
	}
	cat(dim(mat),"\n")
	# find breaks in the columns
	tempcols<-colnames(mat)[which(colnames(mat)!="")]
	if(is.null(seps)) seps<-findDiv(tempcols)*2
	else if(seps[1]=="none") seps<-numeric()
	else seps<-seps*2
	jpeg(output,height=height,width=width,pointsize=24)
	a<-heatmap.2(mat,trace="none",col=pal,breaks=br,mar=c(15,20),Colv=F,Rowv=T,dendrogram="none",colsep=seps,sepcolor="black")
	dev.off()
	return(list(ord=a[[1]],mergedMatrix=mat))
}

# guess at good places to draw division lines in the heatmap
findDiv<-function(colnms){
	tempcols<-gsub("_[0-9]+.+{1,3}$","",colnms)
	seps<-numeric()
	for(i in 1:(length(tempcols)-1))
		if(tempcols[i+1]!=tempcols[i]) seps<-c(seps,i)
	return(seps)
}

# writes a row name reference file that discloses what functional terms are included in each row of the heatmap
# The rows are written in the same order as they appear in the heatmap
writeFinalClusters<-function(allTerms,select,output){
	terms<-allTerms[select]
	sink(output)
	for(i in 1:length(terms)){
		cat(i,"\n",sep="")
		for(j in 1:length(terms[[i]])) cat(terms[[i]][j],"\n",sep="")
	}
	sink()
}

# used for printing reference file containing genes matching functional categories with score
# overlaps is the output of davidSimple() with returnOverlapGenes=TRUE
# map is the mapping of gene IDs to annotation data, any number of annotation columns will work
# folder is the folder name for output files
# matchcol is the column with the gene IDs that match to the initial data
# selcol is used in a special case for mapping back to targeting miRNAs, see the embedded comment
makeGeneReferenceFile<-function(overlaps,map=NULL,matchcol=1,selcol=2,dir="geneRef"){
	if(!file.exists(dir)) dir.create(dir)
	genes<-overlaps$genes;scores<-overlaps$values
	for(i in 1:length(genes)){
		cat("writing file",i,"\r")
		sink(paste("geneRef/",names(genes)[i],".txt",sep=""))
		if(length(genes[[i]])==0){
			sink()
			next
		}
		ord<-order(scores[[i]],decreasing=T)
		scores[[i]]<-scores[[i]][ord]
		genes[[i]]<-genes[[i]][ord]
		for(j in 1:length(genes[[i]])){
			cat(names(genes[[i]])[j],"\t",scores[[i]][j],"\n",sep="")
			if(!is.null(map)){
				if(is.null(dim(map))){ # written for a case where we supply an individual map file for each condition
					# specifically, this is intended for miRNA target data, where we want to map targeted genes
					# back to the miRNA(s) that targeted them. For this case, the function expects a vector of
					# strings, which are file names where the additional info is stored for each condition.
					f<-read.delim(map[i],header=T,sep="\t",stringsAsFactors=F)
					tempgenes<-genes[[i]][[j]]
					selected<-character(length(tempgenes))
					for(k in 1:length(tempgenes)){
						selRows<-f[which(f[,matchcol]==tempgenes[k]),selcol]
						selected[k]<-paste(selRows,collapse=",")				
					}
					g<-cbind(tempgenes,selected)
				}else  # otherwise the function expects a single table of gene info data.
					g<-map[match(genes[[i]][[j]],map[,matchcol]),]
			}else	
				g<-matrix(genes[[i]][[j]],ncol=1)
			for(k in 1:length(genes[[i]][[j]])){
				for(l in 1:dim(g)[2]) cat("\t",g[k,l],sep="")
				cat("\n")
			}
		}
		sink()
	}
	cat("\n")
}

# generates input file for davidSimple(). Assumes a table with primary gene IDs as row names,
# and a fold change column and significance column for each condition.
# datatable is the data table to use
# FCkey is a string to use to correctly identify the fold change columns
# pvalkey is a string to use to correctly identify the significance columns
# nm is a character vector representing the names to be used for each condition
# nm should have labels in the format <condition>_<time point> e.g. "VN1203_12h"
# conversionMap is the map for converting to gene symbols, or whatever form the gene sets are in, if the 
# primary IDs are not already in this form
# fcThresh is the threshold representing the minimum (absolute) log2 fold change value for counting a gene as changed
# pvalThresh is the threshold representing the maximum significance value for counting a gene as changed
getDElist<-function(datatable,FCkey,pvalKey,nm=NULL,conversionMap=NULL,fcThresh=1,pvalThresh=.05){
	cat("Make sure data row names are primary IDs.\n")
	if(!is.null(conversionMap)){ 
		cat("Make sure current ID is the first column of conversion table; new ID is the second column.\n")
		newIDs<-conversionMap[match(rownames(datatable),conversionMap[,1]),2]
	} else newIDs<-rownames(datatable)
	fc<-grep(FCkey,colnames(datatable))
	pval<-grep(pvalKey,colnames(datatable))
	if(length(fc)!=length(pval)){
		cat("Error: Different number of FC and p-val columns\n")
		return()
	}
	delist<-list()
	for(i in 1:length(fc)){
		theseup<-which(datatable[,fc[i]]>=fcThresh & datatable[,pval[i]]<=pvalThresh)
		thesedown<-which(datatable[,fc[i]]<=-fcThresh & datatable[,pval[i]]<=pvalThresh)
		theseup<-unique(newIDs[theseup])
		theseup<-theseup[which(theseup!="")]
		theseup<-theseup[which(!is.na(theseup))]
		delist[[length(delist)+1]]<-theseup
		thesedown<-unique(newIDs[thesedown])
		thesedown<-thesedown[which(thesedown!="")]
		thesedown<-thesedown[which(!is.na(thesedown))]
		delist[[length(delist)+1]]<-thesedown	
	}
	nm<-rep(nm,each=2)
	nm[seq(1,length(nm),2)]<-paste(nm[seq(1,length(nm),2)],"_u",sep="")
	nm[seq(2,length(nm),2)]<-paste(nm[seq(2,length(nm),2)],"_d",sep="")
	names(delist)<-nm
	return(delist)
}

# drillDown() outputs a table and heatmap of the expression of the genes identified from rows and columns of the functional heatmap 
# data is the gene expression table (or a file name for one). It must have a single gene symbol column labled with the word 
# "symbol", and fold change columns labled with "_Log2FC"
# rws is the row number(s) of interest from the functional heatmap
# cols is the column(s) of interest from the functional heatmap
# dirs is the direction of change of interest for each column investigated, e.g. cols=c(3,6,9,11), dirs=c("u","d","d","b")
# u=up, d=down, b=both
# if dirs is left null, both up and down files will be used for all columns
# clustInfo is the name of the file where the function names are stored. by default this is "clusterInfo.txt"
# geneRef is the folder where the gene reference files are stored. by default this is "geneRef"
# output is the optional name for the heatmap file, must end in ".jpg", defaults to "xprsnHeatmap.jpg"
drillDownGenes<-function(data,symbolKey="ymbol",fcKey="_Log2FC",rws=1,cols=1,dirs=NULL,clustInfo="clusterInfo.txt",
		geneRef="geneRef",output="xprsnHeatmap.jpg"){
	if(is.character(data))	data<-read.delim(data,header=T,sep="\t",stringsAsFactors=F)
	# first make sure the column headers make sense, and throw away everything we don't need.
	symbols<-grep(symbolKey,colnames(data))
	fc<-grep(fcKey,colnames(data))
	if(length(symbols)!=1){
		cat("can't find symbol column.\n")
		return()
	}
	cat(length(fc),"fold change columns found.\n")
	data<-cbind(data[,symbols],data[,fc],stringsAsFactors=F)
	colnames(data)<-gsub("_Log2FC","",colnames(data))
	# read in heatmap row information
	lines<-readLines(clustInfo)
	numbers<-grep("^[0-9]+$",lines) # identify cluster demarcations
	starts<-numbers+1 
	stops<-numbers-1
	stops<-stops[-1]
	stops<-c(stops,length(lines))
	terms<-character()
	# now collect all the function names in the rows we need
	for(i in 1:length(rws))
		terms<-c(terms,lines[starts[rws[i]]:stops[rws[i]]])
	genepile<-character() # this is where the associated genes will go
	# need to build file names; can't do a list.files() because they come back in a crazy order
	files<-colnames(data)[-1] # the first column is the symbols so we don't want that column name
	files<-rep(files,each=2) # duplicate each one since we have up and down files for each, and then finish building the geneRef file names.
	files[seq(1,length(files),2)]<-paste(files[seq(1,length(files),2)],"_u.txt",sep="")
	files[seq(2,length(files),2)]<-paste(files[seq(2,length(files),2)],"_d.txt",sep="")
	files<-paste("geneRef",files,sep="/")
	# need to make sure all the right files are selected according to cols and dirs
	cols<-cols*2 # column indices selected need to be doubled, again because there are 2 geneRef files for each apparent heatmap column
	if(is.null(dirs))
		cols<-c(cols-1,cols)
	else cols[which(dirs=="u")]<-cols[which(dirs=="u")]-1
	cols<-c(cols,cols[which(dirs=="b")]-1)
	# loop through the appropriate files and grab the genes
	for(i in 1:length(cols)){
		file<-try(read.delim(files[cols[i]],header=F,sep="\t",stringsAsFactors=F),silent=T)
		if(class(file)!="data.frame") next # skip files with nothing in them
		fileterms<-which(file[,1]!="") # identify function names and get demarcations
		starts<-fileterms+1
		stops<-fileterms-1
		stops<-stops[-1]
		stops<-c(stops,length(file[,1]))
		fileterms<-file[fileterms,1]
		termInds<-match(terms,fileterms) # find which ones match the terms we need
		termInds<-termInds[which(!is.na(termInds))] # there may be terms we want that aren't here.
		# if there are matches, loop through the terms and collect each one's associated genes
		if(length(termInds)>0)
			for(j in 1:length(termInds))		
				genepile<-c(genepile,file[starts[termInds[j]]:stops[termInds[j]],2])
	}
	genepile<-unique(genepile)
	# match up the genes with the idenfiers in the expression table
	these<-which(data[,1] %in% genepile)
	data<-data[these,]
	labl<-data[,1] # need to save these so we can label the heatmap appropriately. Can't make them rownames cuz they're not unique
	pal<-brewer.pal(11,"RdYlBu")[11:1]
	br<-seq(-4,4,8/11)
	seps<-findDiv(colnames(data)[-1])# get column separations
	jpeg(output,width=2000,height=3000,pointsize=24)
	# remove the labels and double t() to make it a numeric matrix; symkey setting prevents warning when there are few data points
	heatmap.2(t(t(data[,-1])),trace="none",col=pal,mar=c(10,10),breaks=br,Colv=F,Rowv=T,dendrogram="none",labRow=labl,colsep=seps,sepcolor="black",symkey=F)
	dev.off()
	return(data)
}

drillDownText<-function(dave,txt){
	pthwy<-dave$overlaps$genes
	values<-dave$overlaps$values
	len<-length(pthwy)
	allnames<-unique(unlist(lapply(pthwy,names)))
	allnames<-allnames[grep(txt,allnames)]
	mat<-matrix(0,nrow=length(allnames),ncol=len)
	rownames(mat)<-allnames;colnames(mat)<-names(pthwy)
	for(i in 1:len){
		pth<-names(pthwy[[i]])
		val<-values[[i]]
		sel<-grep(txt,pth)
		mat[pth[sel],i]<-val[sel]
	}
	process4BidirectionalHM(mat,1.30103,bidirectional=T,seps=c(4,8),height=750,width=960,output=paste0(txt,"-pathways_heatmap.jpg"),mx=6)
	return(mat)
}


drillDownRows<-function(dave,rws,clustInfo){
	pthwy<-dave$overlaps$genes
	values<-dave$overlaps$values
	lines<-readLines(clustInfo)
	numbers<-grep("^[0-9]+$",lines) # identify cluster demarcations
	starts<-numbers+1 
	stops<-numbers-1
	stops<-stops[-1]
	stops<-c(stops,length(lines))
	terms<-character()
	# now collect all the function names in the rows we need
	terms<-unlist(sapply(rws,function(x) lines[starts[x]:stops[x]]))
	len<-length(pthwy)
	mat<-matrix(0,nrow=length(terms),ncol=len)
	rownames(mat)<-terms;colnames(mat)<-names(pthwy)
	for(i in 1:len){
		pth<-names(pthwy[[i]])
		val<-values[[i]]
		sel<-which(pth %in% terms)
		mat[pth[sel],i]<-val[sel]
	}
	process4BidirectionalHM(mat,1.30103,bidirectional=T,seps=c(4,8),height=750,width=960,output=paste0("row",rws,"-pathways_heatmap.jpg"),mx=6)
	return(mat)
}




# used for comparing the performance of two sets of clustered data (only compare the clustered part of each matrix)
# overlaps is the second item of the list object outputted from davidSimple()
evalPerformance<-function(mat1,mat2,overlaps){
	cor1<-as.dist(cor(mat1))
	cor2<-as.dist(cor(mat2))
	unclust<-processNonClusteredData(overlaps,0,"dummy.txt")
	unlink("dummy.txt")
	corU<-as.dist(cor(unclust))
	distScore1<-mean(abs(cor1-corU))
	distScore2<-mean(abs(cor2-corU))
	corScore1<-cor(cor1,corU)
	corScore2<-cor(cor2,corU)
	cat("Distance measure:\n\tmatrix 1: ",distScore1,"\n\tmatrix 2: ",distScore2,"\n")
	cat("Correlation measure:\n\tmatrix 1: ",corScore1,"\n\tmatrix 2: ",corScore2,"\n")
}


