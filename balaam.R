# MultivarALAAM.R
# source("/Users/johankoskinen/Documents/manchester admin/bela/MultivarALAAM.R")
simulate.alaam <- function(ALAAMobj,statsvec=NULL,theta,contagion ='simple', thinning = 1, NumIterations = 1, burnin = 3000, DoSave=FALSE, returnNet= FALSE,doGOF=FALSE)
{
	y <- ALAAMobj$y
	n <- length(y)
	directed <- ALAAMobj$directed
	interaction <- ALAAMobj$interaction
	covariates <- ALAAMobj$covariates
	p <- dim(covariates)[2]+1# add 1 for intercept
	p1 <- p
	
	if (is.null(statsvec))
	{
		statsvec <- getALAAMstats(ALAAMobj,contagion =contagion)
	}
	if ('simple' %in% contagion){
		p <- p+1
		
		EdgeList <- ALAAMobj$EdgeList
		RowIn <- ALAAMobj$RowIn
		degree <- ALAAMobj$degree
		covariates <- ALAAMobj$covariates
		if (directed==TRUE){
			EdgeListIn <- ALAAMobj$EdgeListIn 
		 	degreein <- ALAAMobj$degreein
		 	RowInIn <- ALAAMobj$RowInIn
		}
		
		}
		
	if (('simple' %in% contagion)==FALSE){
		stop('internal error number 36: cannot simulate wihtout contagion statistic')
		}
		dorec <- FALSE
	if ('recip' %in% contagion){
		p <- p+1
		dorec <- TRUE
		
		
		RecEdgeList <- ALAAMobj$RecEdgeList
		Recdegree <- ALAAMobj$Recdegree
		RecRowIn <- ALAAMobj$RecRowIn
		}
		doindir <- FALSE
	if ('indirect' %in% contagion)
		{
			p <- p+1
			doindir <- TRUE
		
			
		P2EdgeList <- ALAAMobj$P2EdgeList
		P2degree <- ALAAMobj$P2degree
		P2RowIn <- ALAAMobj$P2RowIn
		if (directed==TRUE){
			P2EdgeListIn <- ALAAMobj$P2EdgeListIn
		P2degreeIn <- ALAAMobj$P2degreeIn
		P2RowInIn <- ALAAMobj$P2RowInIn
		}
		if (is.null(P2EdgeList) | is.null(P2degree) | is.null(P2RowIn))
		{
		error('edgelist, 2-path degree, and index to 2-path edgelist required if indirect contagion requested')	
		}
	}
	doCloseindir <- FALSE
	if ('closedind' %in% contagion)
		{
			p <- p+1
			doCloseindir <- TRUE
		CP2EdgeList <- ALAAMobj$CP2EdgeList
		CP2degree <- ALAAMobj$CP2degree
		CP2RowIn <- ALAAMobj$CP2RowIn
		 			
		if (directed==TRUE){
			CP2EdgeListIn <- ALAAMobj$CP2EdgeListIn
	CP2degreeIn <- ALAAMobj$CP2degreeIn
	CP2RowInIn <- ALAAMobj$CP2RowInIn
		}
	}
	doTrans <- FALSE
	if ('transitive' %in% contagion)
		 	{
		 		p <- p+1
		 		doTrans <- TRUE
		 		translist <- ALAAMobj$translist
		 		transFirst <- ALAAMobj$transFirst
		 		transMiddle <- ALAAMobj$transMiddle
		 		transEndnode <- ALAAMobj$transEndnode
	 
		 		}
		 	
	
	p2 <- 0
	if (!is.null(interaction))
	{
		p2 <- length(interaction)
		p <- p+p2
	}

	out <- simulateALAAM(y=y,
	                     EdgeList=EdgeList,
	                     RowIn=RowIn,
	                     degree = degree,
	                     EdgeListIn=EdgeListIn,
	                     RowInIn=RowInIn,
	                     degreein = degreein,covariates = covariates,NumIterations=NumIterations,theta=theta,statsvec=statsvec ,DoSave=DoSave,returnNet=returnNet, directed=directed,interaction=interaction, thinning = thinning,burnin = burnin, 
	                     contagion =contagion,
		RecEdgeList=RecEdgeList,
		Recdegree=Recdegree,
		RecRowIn=RecRowIn, 
		P2EdgeList=P2EdgeList,
		P2degree=P2degree, 
		P2RowIn=P2RowIn, 
		P2EdgeListIn=P2EdgeListIn, 
		P2degreeIn=P2degreeIn, 
		P2RowInIn=P2RowInIn,
		CP2EdgeList = CP2EdgeList,
		CP2degree = CP2degree,
		CP2RowIn =CP2RowIn,
		CP2EdgeListIn = CP2EdgeListIn,
	CP2degreeIn = CP2degreeIn,
	CP2RowInIn = CP2RowInIn,
	translist = translist,
		 		transFirst = transFirst,
		 		transMiddle = transMiddle,
		 		transEndnode = transEndnode,doGOF=doGOF
 )
		
		
	out
}

simulateALAAM <- function(y,EdgeList,RowIn,degree,covariates,NumIterations,theta,statsvec,directed=FALSE,EdgeListIn=NULL, degreein=NULL,RowInIn=NULL,DoSave=FALSE,thinning = NULL,burnin=NULL,canchange=NULL,returnNet=FALSE,interaction=NULL,checkAlgorithm=FALSE, bias = NULL, contagion ='simple',RecEdgeList=NULL,Recdegree=NULL,RecRowIn=NULL,P2EdgeList=NULL,P2degree=NULL,P2RowIn=NULL, P2EdgeListIn=NULL,P2degreeIn=NULL,P2RowInIn	=NULL, CP2EdgeList=NULL, CP2degree = NULL, CP2RowIn = NULL , CP2EdgeListIn = NULL , CP2degreeIn = NULL, CP2RowInIn = NULL, translist = NULL , transFirst = NULL , transMiddle = NULL , transEndnode = NULL,doGOF=FALSE)
{
	#browser()
	####
		
	if (doGOF==TRUE){
		this.theta <- theta[1,]
		num.thetas <- dim(theta)[1]
	}
	
	#####
	n <- length(y)
	p <- dim(covariates)[2]+1# add 1 for intercept
	p1 <- p
	if ('simple' %in% contagion){
		p <- p+1
		}
		
	if (('simple' %in% contagion)==FALSE){
		error('model without simple contagion not yet specified')
		}
		dorec <- FALSE
	if ('recip' %in% contagion){
		p <- p+1
		dorec <- TRUE
		if (is.null(RecEdgeList) | is.null(Recdegree) | is.null(RecRowIn))
		{
		error('edgelist, reciprocated degree, and index to reciprocated edgelist required if reciprocated contagion requested')	
		}
		
		}
		doindir <- FALSE
	if ('indirect' %in% contagion)
		{
			p <- p+1
			doindir <- TRUE
		
			if (is.null(P2EdgeList) | is.null(P2degree) | is.null(P2RowIn))
		{
		error('edgelist, 2-path degree, and index to 2-path edgelist required if indirect contagion requested')	
		}
		 }
	doCloseindir <- FALSE
	if ('closedind' %in% contagion)
		{
			p <- p+1
			doCloseindir <- TRUE
			}
	doTrans <- FALSE
	if ('transitive' %in% contagion)
		 	{
		 		p <- p+1
		 		doTrans <- TRUE
		 				 		}

	p2 <- 0
	if (!is.null(interaction))
	{
	  
		p2 <- length(interaction)
		p <- p+p2
	}
	ystar <- y

	statsvecstar <- matrix(0,p,1)
	
	if (is.null(canchange))
	{
		# regular sampling
		canchange <- c(1:n)
		
	}
	
	if ( DoSave )# we output statistics
	{
		if (is.null(burnin))
		{
			burnin <- n^2
		}
	BigTotal <- thinning*NumIterations
	saveFreq <- ceiling(BigTotal/NumIterations)
	SimStats <- matrix(0,p,NumIterations)
	if (returnNet=='TRUE')
	{
SampleNet<-matrix(0,n,NumIterations)

}
	NumIterations <- BigTotal+burnin
	saveHere <- 1
	

	}
	
	for (iterations in c(1:NumIterations))
	{
		#actor <- ceiling(runif(1)*n)
		actor <- sample(canchange,size=1)
		statsvecstar[1] <- 1 - 2*ystar[actor]
		neig <- 0 # for contagion count
		statsvecstar[2] <- 0
		Recneig <- 0 # for reciprocated contagion count
		Indineig <- 0 # for indirect contagion count
		CIndineig <- 0 # for closed indirect contagion count
		TIndineig <- 0 # for transitive contagion count
		if (degree[actor] > 0)
		{
		for (alter in c(1:degree[actor]))
		{
			neig <- neig + ystar[ EdgeList[RowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		#neig <- neig + sum(ystar[ EdgeList[RowIn[actor]+c(1:degree[actor])-1,2] ])
		}
		
		if (directed)
		{
			if (degreein[actor]>0)
			{
			for (alter in c(1:degreein[actor]))
		{
			neig <- neig + ystar[ EdgeListIn[RowInIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		#neig <- neig + sum(ystar[ EdgeListIn[RowInIn[actor]+c(1:degreein[actor])-1,2] ])
		}
		}
		
		### reciprocated contagion if requested
		if (dorec)
		{
				#ALAAMobj$RecEdgeList <- RecEdgeList
		 		#ALAAMobj$Recdegree <- Recdegree
		 		#ALAAMobj$RecRowIn <- RecRowIn

		if (Recdegree[actor]>0)
		{
			for (alter in c(1:Recdegree[actor]))
		{
			Recneig <- Recneig + ystar[ RecEdgeList[RecRowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
			
		}

		### indirect contagion if requested
		if (doindir)
		{
		 		#ALAAMobj$P2EdgeList <- P2EdgeList
		 		#ALAAMobj$P2degree <- P2degree
		 		#ALAAMobj$P2RowIn <- P2RowIn
if (is.na(P2degree[actor]))
{
	print('degree distribution is undefined')
	#browser()
	
}
		if (P2degree[actor]>0)
		{
			for (alter in c(1:P2degree[actor]))
		{
			Indineig <- Indineig + ystar[ P2EdgeList[P2RowIn[actor]+alter-1,2] ]*P2EdgeList[P2RowIn[actor]+alter-1,3]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		
		if (directed)
		{
			
			if (P2degreeIn[actor]>0)
		{
			for (alter in c(1:P2degreeIn[actor]))
		{
			Indineig <- Indineig + ystar[ P2EdgeListIn[P2RowInIn[actor]+alter-1,2] ]*P2EdgeListIn[P2RowInIn[actor]+alter-1,3]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
			
			}
			
		}
		
		### indirect + direct if requested
		if (doCloseindir)
		{
			if (CP2degree[actor]>0)
		{
			for (alter in c(1:CP2degree[actor]))
		{
			CIndineig <- CIndineig + ystar[ CP2EdgeList[CP2RowIn[actor]+alter-1,2] ]*CP2EdgeList[CP2RowIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		
		if (directed)
		{
			if (CP2degreeIn[actor]>0)
		{
			for (alter in c(1:CP2degreeIn[actor]))
		{
			CIndineig <- CIndineig + ystar[ CP2EdgeListIn[CP2RowInIn[actor]+alter-1,2] ]*CP2EdgeListIn[CP2RowInIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		}

		}
		
		### transitive contagion
		if (doTrans)
		{
			
			copyTP <- transFirst[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],2]]==1 && ystar[translist[copyTP[rowind],3]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			copyTP <- transMiddle[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],3]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			copyTP <- transEndnode[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],2]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			

			
		}
		
		
		statsvecstar[2] <- statsvecstar[1]*neig
			k <- 3
		if (dorec)
		{
			statsvecstar[k] <- statsvecstar[1]*Recneig
			k <- k+1
		}
		if (doindir)
		{
			statsvecstar[k] <- statsvecstar[1]*Indineig
			k <- k+1
		}
		
		if (doCloseindir)
		{
			statsvecstar[k] <- statsvecstar[1]*CIndineig
			k <- k+1
			
		}
		if (doTrans)
		{
			statsvecstar[k] <- statsvecstar[1]*TIndineig
			k <- k+1
		}
		
		for (t in c(1:(p1 - 1) ))
		{
			statsvecstar[k] <- statsvecstar[1]*covariates[actor,t]
			k <- k+1
		}
		
		if (!is.null(interaction))
	{
		#statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction]
		statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction[1]]
	}
	
	if (doGOF==FALSE)
	{
		hrat <- theta %*%  statsvecstar
		}
		
	if (doGOF==TRUE){
		hrat <- this.theta %*%  statsvecstar
		}	
		
		
	if (!is.null(bias))
	{
		hrat  <- hrat + statsvecstar[1]*bias[actor]
	}
		
		if ( log(runif(1)) <  hrat   )
		{
		
			statsvec <- statsvec+statsvecstar
			ystar[actor] <- 1 - ystar[actor]
		}
		
		
		
		
		if ( (DoSave==TRUE) && (iterations>burnin) && ((iterations %% saveFreq)==0) )
{
	#print(paste('you have done ',iterations,' iterations'))
	SimStats[,saveHere] <- statsvec
	
	if (returnNet==TRUE)
	{
SampleNet[,saveHere] <-  ystar
}

saveHere <- saveHere +1
if (doGOF==TRUE) 
{
if (saveHere<= num.thetas){
		this.theta <- theta[saveHere,]
	}
}
	}
	}
	
	if ( DoSave )# we output statistics
	{
		if (checkAlgorithm==TRUE)
		{
			statsvecKoll <- getALAAMstats(y=ystar,EdgeList=EdgeList,RowIn=RowIn,degree = degree,EdgeListIn=EdgeListIn,RowInIn=RowInIn,degreein = degreein,covariates = covariates,directed=directed,interaction=interaction)
if ( any(statsvec!= statsvecKoll) )
{
	print(paste('diff :',statsvec- statsvecKoll))

}
}

		statsvec <- SimStats
		}
	if ( DoSave==FALSE )# we output statistics
	{
	
		}
		
	if (returnNet==TRUE)
	{
		if (DoSave)
		{
			
			ystar <- SampleNet
		}
		statsvec <-list(y=ystar,statsvec=statsvec)
	}
		
		
	statsvec 

}



simulateALAAMtest <- function(y,EdgeList,RowIn,degree,covariates,NumIterations,theta,statsvec,directed=FALSE,EdgeListIn=NULL, degreein=NULL,RowInIn=NULL,DoSave=FALSE,thinning = NULL,burnin=NULL,canchange=NULL,returnNet=FALSE,interaction=NULL,checkAlgorithm=FALSE, bias = NULL, contagion ='simple')
{
	n <- length(y)
	p <- dim(covariates)[2]+2# add 2 for intercept and contagion
	p1 <- p
	if (!is.null(interaction))
	{
		p <- p+length(interaction)
	}
	ystar <- y

	statsvecstar <- matrix(0,p,1)
	
	if (is.null(canchange))
	{
		# regular sampling
		canchange <- c(1:n)
		
	}
	
	if ( DoSave )# we output statistics
	{
		if (is.null(burnin))
		{
			burnin <- n^2
		}
	BigTotal <- thinning*NumIterations
	saveFreq <- ceiling(BigTotal/NumIterations)
	SimStats <- matrix(0,p,NumIterations)
	if (returnNet=='TRUE')
	{
SampleNet<-matrix(0,n,NumIterations)

}
	NumIterations <- BigTotal+burnin
	saveHere <- 1
	

	}
	
	for (iterations in c(1:NumIterations))
	{
		#actor <- ceiling(runif(1)*n)
		actor <- sample(canchange,size=1)
		statsvecstar[1] <- 1 - 2*ystar[actor]
		neig <- 0 # for contagion count
		statsvecstar[2] <- 0
		
		if (degree[actor] > 0)
		{
		for (alter in c(1:degree[actor]))
		{
			neig <- neig + ystar[ EdgeList[RowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		#neig <- neig + sum(ystar[ EdgeList[RowIn[actor]+c(1:degree[actor])-1,2] ])
		}
		
		if (directed)
		{
			if (degreein[actor]>0)
			{
			for (alter in c(1:degreein[actor]))
		{
			neig <- neig + ystar[ EdgeListIn[RowInIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		#neig <- neig + sum(ystar[ EdgeListIn[RowInIn[actor]+c(1:degreein[actor])-1,2] ])
		}
		}
		statsvecstar[2] <- statsvecstar[1]*neig
		for (k in c(3:p1))
		{
			statsvecstar[k] <- statsvecstar[1]*covariates[actor,k-2]
		}
		if (!is.null(interaction))
	{
		statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction]
	}
		hrat <- theta %*%  statsvecstar
	if (!is.null(bias))
	{
		hrat  <- hrat + statsvecstar[1]*bias[actor]
	}
		
		if ( log(runif(1)) <  hrat   )
		{
			statsvec <- statsvec+statsvecstar
			ystar[actor] <- 1 - ystar[actor]
		}
		
		
		
		
		if ( (DoSave==TRUE) && (iterations>burnin) && ((iterations %% saveFreq)==0) )
{
	#print(paste('you have done ',iterations,' iterations'))
	SimStats[,saveHere] <- statsvec
	
	if (returnNet==TRUE)
	{
SampleNet[,saveHere] <-  ystar
}
saveHere <- saveHere +1
	}
	}
	
	if ( DoSave )# we output statistics
	{
		if (checkAlgorithm==TRUE)
		{
			statsvecKoll <- getALAAMstats(y=ystar,EdgeList=EdgeList,RowIn=RowIn,degree = degree,EdgeListIn=EdgeListIn,RowInIn=RowInIn,degreein = degreein,covariates = covariates,directed=directed,interaction=interaction)
if ( any(statsvec!= statsvecKoll) )
{
	print(paste('diff :',statsvec- statsvecKoll))

}
}

		statsvec <- SimStats
		}
	if ( DoSave==FALSE )# we output statistics
	{
	
		}
		
	if (returnNet==TRUE)
	{
		if (DoSave)
		{
			
			ystar <- SampleNet
		}
		statsvec <-list(y=ystar,statsvec=statsvec)
	}
		
		
	statsvec 

}

prepALAAMdata <- function(y,ADJ,
                          covariates,
                          directed=FALSE,
                          useDegree = FALSE,
                          contagion ='simple',
                          interaction=NULL)
{
	n <- dim(ADJ)[1]
	RowIn <- matrix(NA,n,1)
	degree <- matrix(0,n,1)
	#### set nulls
	
	
	
	####
	
	
	
	if (length(y)!=n)
	{
		warning('response variable not same lenght as rows of adjacency matrix')
	}
	 if (directed)
		 {
		 	RowInIn <- matrix(NA,n,1)
	degreein <- matrix(0,n,1)
	if (sum( t(ADJ)!=ADJ  )==0 )
	{
		warning('you defined network as directed but it is perfectly symmetric')
		
	} 
	}
	if (directed==FALSE)
		 {
		 	if (sum( t(ADJ)!=ADJ  )!=0)
		 	{
		 		warning('you defined network as un-directed but it is assymetric')
		 	}
		 }
	
	if (n != dim(ADJ)[2])
	{
		warning('your adjacency matrix is not square')
	}
	if (sum(diag(ADJ) )!=0 )
	{
		warning('diagonal elements of adjacency matrix non-zero')
	}

	
	
		 EdgeList <- which(ADJ==1,arr.ind=TRUE)
		 EdgeList <- EdgeList[order(EdgeList[,1]),]
		 for (actor in c(1:n))
		 {
		 	TheseRows <- which(EdgeList[,1]==actor)
		 	if (length(TheseRows)>0)
		 	{
		 		RowIn[actor] <- min(TheseRows)
		 		degree[actor] <- length(TheseRows)
		 	}
		 	
		 }
		 
		 if (directed)
		 {
		 	EdgeListIn <- which(t(ADJ)==1,arr.ind=TRUE)
		 EdgeListIn <- EdgeListIn[order(EdgeListIn[,1]),]
		 for (actor in c(1:n))
		 {
		 	TheseRows <- which(EdgeListIn[,1]==actor)
		 	if (length(TheseRows)>0)
		 	{
		 		RowInIn[actor] <- min(TheseRows)
		 		degreein[actor] <- length(TheseRows)
		 	}
		 	
		 }

		 	
		 }
if (useDegree==TRUE)
{
	covariates <- cbind(degree,covariates)
}

if ('recip' %in% contagion){
	if (directed==FALSE)
	{
	#	error('reciprochal contagion not defined for undirected networks')
	}
	# print('contagion through reciprocal ties included')
	# produce edgelist for reciprochal ties
	RecRowIn <- matrix(NA,n,1)
	Recdegree <- matrix(0,n,1)
	RecADJ <- ADJ*t(ADJ)
	RecADJ[RecADJ>0] <- 1
	diag(RecADJ) <- 0
	RecEdgeList <- which(RecADJ==1,arr.ind=TRUE)
	RecEdgeList <- RecEdgeList[order(RecEdgeList[,1]),]
		 for (actor in c(1:n))
		 {
		 	TheseRows <- which(RecEdgeList[,1]==actor)
		 	if (length(TheseRows)>0)
		 	{
		 		RecRowIn[actor] <- min(TheseRows)
		 		Recdegree[actor] <- length(TheseRows)
		 	}
		 	
		 }
		 

	
}

if ('indirect' %in% contagion)
{
	# list neighbours at distance 2
	if (directed==FALSE)
	{
	#	error('indirect contagion is not implemented for undirected networks')
	}
	# create matrix of paths i -> k -> j
	path2 <- matrix(0,n,n)
	for (i in c(1:n))
	{
		for (j in c(1:n))
		{
			if (i != j)
			{
			path2[i,j] <- ADJ[i,] %*% ADJ[,j] 
			}
		}
	}
	# produce edgelist for neighbours at distance 2
	P2RowIn <- matrix(NA,n,1)
	P2degree <- matrix(0,n,1)
	P2EdgeList <- which(path2>=1,arr.ind=TRUE)
	
	P2EdgeList <- cbind(P2EdgeList[order(P2EdgeList[,1]),],matrix(0,dim(P2EdgeList)[1],1) )
		 for (actor in c(1:n))
		 {
		 	TheseRows <- which(P2EdgeList[,1]==actor)
		 	if (length( P2EdgeList[TheseRows,c(1:2)])==2)
		 	{
		 	#	print('pucko')
		 	#	browser()
		 		P2EdgeList[TheseRows,3] <-  path2[P2EdgeList[TheseRows,1],P2EdgeList[TheseRows,2]]
		 	}
		 	if (length( P2EdgeList[TheseRows,c(1:2)])>2)
		 	{
		 	P2EdgeList[TheseRows,3] <- path2[P2EdgeList[TheseRows,c(1:2)]]
		 	}
		 	if (length(TheseRows)>0)
		 	{
		 		P2RowIn[actor] <- min(TheseRows)
		 		P2degree[actor] <- length(TheseRows)
		 	}
		 	
		 }
	 if (directed){
	 	path2 <- matrix(0,n,n)
	for (i in c(1:n))
	{
		for (j in c(1:n))
		{
			if (i != j)
			{
			path2[i,j] <- ADJ[j,] %*% ADJ[,i] 
			}
		}
	}
	 	
	#ALAAMobj$P2EdgeListIn <- P2EdgeListIn
	#ALAAMobj$P2degreeIn <- P2degreeIn
	#ALAAMobj$P2RowInIn <- P2RowInIn	 
	P2RowInIn <- matrix(NA,n,1)
	P2degreeIn <- matrix(0,n,1)
	P2EdgeListIn <- which(path2>=1,arr.ind=TRUE)
	P2EdgeListIn <- P2EdgeListIn[order(P2EdgeListIn[,1]),]
	P2EdgeListIn <- cbind(P2EdgeListIn[order(P2EdgeListIn[,1]),],matrix(0,dim(P2EdgeListIn)[1],1) )
	for (actor in c(1:n))
		 {
		 	TheseRows <- which(P2EdgeListIn[,1]==actor)
		 	if (length( P2EdgeListIn[TheseRows,c(1:2)])==2)
		 	{
		 	#	print('pucko')
		 	#	browser()
		 		P2EdgeListIn[TheseRows,3] <-  path2[P2EdgeListIn[TheseRows,1],P2EdgeListIn[TheseRows,2]]
		 	}
		 	if (length( P2EdgeListIn[TheseRows,c(1:2)])>2)
		 	{

		 	P2EdgeListIn[TheseRows,3] <- path2[P2EdgeListIn[TheseRows,c(1:2)]]
		 	}
		 	if (length(TheseRows)>0)
		 	{
		 		P2RowInIn[actor] <- min(TheseRows)
		 		P2degreeIn[actor] <- length(TheseRows)
		 	}
		 	
		 }

	}
	
}
		
	ALAAMobj <- list(y = y,EdgeList = EdgeList, RowIn =RowIn, degree = degree, covariates= covariates)
	ALAAMobj$EdgeListIn <-NULL
		 	ALAAMobj$degreein <- NULL
		 	ALAAMobj$RowInIn <- NULL
		 	ALAAMobj$RecEdgeList <- NULL
		 		ALAAMobj$Recdegree <- NULL
		 		ALAAMobj$RecRowIn <- NULL
		 		ALAAMobj$P2EdgeList <- NULL
		 		ALAAMobj$P2degree <- NULL
		 		ALAAMobj$P2RowIn <- NULL
		 		ALAAMobj$P2EdgeListIn <- NULL
	ALAAMobj$P2degreeIn <- NULL
	ALAAMobj$P2RowInIn <- NULL
	ALAAMobj$CP2EdgeList <- NULL
		 		ALAAMobj$CP2degree <- NULL
		 		ALAAMobj$CP2RowIn <- NULL
		 		ALAAMobj$CP2EdgeListIn <- NULL
	ALAAMobj$CP2degreeIn <- NULL
	ALAAMobj$CP2RowInIn <- NULL	 
	ALAAMobj$translist <- NULL
		 		ALAAMobj$transFirst <- NULL
		 		ALAAMobj$transMiddle <- NULL
		 		ALAAMobj$transEndnode <- NULL
		 	
		 if (directed)
		 {
		 	ALAAMobj$EdgeListIn <- EdgeListIn
		 	ALAAMobj$degreein <- degreein
		 	ALAAMobj$RowInIn <- RowInIn
		 	ALAAMobj$directed <- TRUE
		 	if ('recip' %in% contagion)
		 	{
		 		ALAAMobj$RecEdgeList <- RecEdgeList
		 		ALAAMobj$Recdegree <- Recdegree
		 		ALAAMobj$RecRowIn <- RecRowIn
		 	}
		 	if ('indirect' %in% contagion)
		 	{
		 		ALAAMobj$P2EdgeList <- P2EdgeList
		 		ALAAMobj$P2degree <- P2degree
		 		ALAAMobj$P2RowIn <- P2RowIn
		 		ALAAMobj$P2EdgeListIn <- P2EdgeListIn
	ALAAMobj$P2degreeIn <- P2degreeIn
	ALAAMobj$P2RowInIn <- P2RowInIn	 
		 	}

		 	}
		 	
		 	if ('closedind' %in% contagion)
		 	{
		 		CP2 <- getlistTransbase(ADJ,directed=directed)
		 		ALAAMobj$CP2EdgeList <- CP2$CP2EdgeList
		 		ALAAMobj$CP2degree <- CP2$CP2degree
		 		ALAAMobj$CP2RowIn <- CP2$CP2RowIn
		 		ALAAMobj$CP2EdgeListIn <- CP2$CP2EdgeListIn
	ALAAMobj$CP2degreeIn <- CP2$CP2degreeIn
	ALAAMobj$CP2RowInIn <- CP2$CP2RowInIn	 
		 		}
		 	if ('transitive' %in% contagion)
		 	{
		 		TP <- getlistTrans(ADJ,directed=directed)
		 		
		 		ALAAMobj$translist <- TP$translist
		 		ALAAMobj$transFirst <- TP$First
		 		ALAAMobj$transMiddle <- TP$Middle
		 		ALAAMobj$transEndnode <- TP$Endnode
	 
		 		}
		 	
		 	if (!directed)
		 	{
		 		ALAAMobj$directed <- FALSE
		 	
		 		if ('indirect' %in% contagion)
		 		{
		 		  ALAAMobj$P2EdgeList <- P2EdgeList
		 		  ALAAMobj$P2degree <- P2degree
		 		  ALAAMobj$P2RowIn <- P2RowIn
		 		  #ALAAMobj$P2EdgeListIn <- P2EdgeListIn
		 		  #ALAAMobj$P2degreeIn <- P2degreeIn
		 		  #ALAAMobj$P2RowInIn <- P2RowInIn	 
		 		}
		 		}
		 		
		 		if (!is.null(interaction))
		 		{
		 		  ALAAMobj$interaction <- interaction
		 		  
		 		}
		 		ALAAMobj$contagion <- contagion
		ALAAMobj
}

getALAAMstats <- function(ALAAMobj,contagion ='simple')
{
	#print('heyho')
	y <- ALAAMobj$y
	n <- length(y)
	directed <- ALAAMobj$directed
	interaction <- ALAAMobj$interaction
	covariates <- ALAAMobj$covariates
	p <- dim(covariates)[2]+1# add 1 for intercept
	p1 <- p
	if ('simple' %in% contagion){
		p <- p+1
		
		EdgeList <- ALAAMobj$EdgeList
		RowIn <- ALAAMobj$RowIn
		degree <- ALAAMobj$degree
		covariates <- ALAAMobj$covariates
		if (directed==TRUE){
			EdgeListIn <- ALAAMobj$EdgeListIn 
		 	degreein <- ALAAMobj$degreein
		 	RowInIn <- ALAAMobj$RowInIn
		}
		
		}
		
	if (('simple' %in% contagion)==FALSE){
		error('model without simple contagion not yet specified')
		}
		dorec <- FALSE
	if ('recip' %in% contagion){
		p <- p+1
		dorec <- TRUE
		
		
		RecEdgeList <- ALAAMobj$RecEdgeList
		Recdegree <- ALAAMobj$Recdegree
		RecRowIn <- ALAAMobj$RecRowIn
		}
		doindir <- FALSE
	if ('indirect' %in% contagion)
		{
			p <- p+1
			doindir <- TRUE
		
			
		P2EdgeList <- ALAAMobj$P2EdgeList
		P2degree <- ALAAMobj$P2degree
		P2RowIn <- ALAAMobj$P2RowIn
		if (directed==TRUE){
			P2EdgeListIn <- ALAAMobj$P2EdgeListIn
		P2degreeIn <- ALAAMobj$P2degreeIn
		P2RowInIn <- ALAAMobj$P2RowInIn
		}
		if (is.null(P2EdgeList) | is.null(P2degree) | is.null(P2RowIn))
		{
		print('edgelist, 2-path degree, and index to 2-path edgelist required if indirect contagion requested')	
		}
	}
	doCloseindir <- FALSE
	if ('closedind' %in% contagion)
		{
			p <- p+1
			doCloseindir <- TRUE
		CP2EdgeList <- ALAAMobj$CP2EdgeList
		CP2degree <- ALAAMobj$CP2degree
		CP2RowIn <- ALAAMobj$CP2RowIn
		 			
		if (directed==TRUE){
			CP2EdgeListIn <- ALAAMobj$CP2EdgeListIn
	CP2degreeIn <- ALAAMobj$CP2degreeIn
	CP2RowInIn <- ALAAMobj$CP2RowInIn
		}
	}
	doTrans <- FALSE
	if ('transitive' %in% contagion)
		 	{
		 		p <- p+1
		 		doTrans <- TRUE
		 		translist <- ALAAMobj$translist
		 		transFirst <- ALAAMobj$transFirst
		 		transMiddle <- ALAAMobj$transMiddle
		 		transEndnode <- ALAAMobj$transEndnode
	 
		 		}
		 	
	
	p2 <- 0
	if (!is.null(interaction))
	{
		p2 <- length(interaction)
		p <- p+p2
	}
	ystar <- matrix(0,n,1)
	statsvec <- matrix(0,p,1)
	statsvecstar <- statsvec
	for (actor in c(1:n))
	{
		
		if (y[actor]==1)
		{
		statsvecstar[1] <- 1 - 2*ystar[actor]
		neig <- 0 # for contagion count
		Recneig <- 0 # for reciprocated contagion count
		Indineig <- 0 # for indirect contagion count
		CIndineig <- 0 # for closed indirect contagion count
		TIndineig <- 0 # for transitive contagion count
		statsvecstar[2] <- neig
		if (degree[actor] > 0)
		{
		for (alter in c(1:degree[actor]))
		{
			neig <- neig + ystar[ EdgeList[RowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		}
		#print('hoho')
		
		#browser()
		if (directed)
		{
		if (degreein[actor]>0)
		{
			for (alter in c(1:degreein[actor]))
		{
			neig <- neig + ystar[ EdgeListIn[RowInIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		
		}
		
		}
		### reciprocated contagion if requested
		if (dorec)
		{
				#ALAAMobj$RecEdgeList <- RecEdgeList
		 		#ALAAMobj$Recdegree <- Recdegree
		 		#ALAAMobj$RecRowIn <- RecRowIn

		if (Recdegree[actor]>0)
		{
			for (alter in c(1:Recdegree[actor]))
		{
			Recneig <- Recneig + ystar[ RecEdgeList[RecRowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
			
		}
		### indirect contagion if requested
		if (doindir)
		{
		 		#ALAAMobj$P2EdgeList <- P2EdgeList
		 		#ALAAMobj$P2degree <- P2degree
		 		#ALAAMobj$P2RowIn <- P2RowIn

		if (P2degree[actor]>0)
		{
			for (alter in c(1:P2degree[actor]))
		{
			Indineig <- Indineig + ystar[ P2EdgeList[P2RowIn[actor]+alter-1,2] ]*P2EdgeList[P2RowIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		
		if (directed)
		{
			if (P2degreeIn[actor]>0)
		{
			for (alter in c(1:P2degreeIn[actor]))
		{
			Indineig <- Indineig + ystar[ P2EdgeListIn[P2RowInIn[actor]+alter-1,2] ]*P2EdgeListIn[P2RowInIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		}
			
		}
		
		### indirect + direct if requested
		if (doCloseindir)
		{
			if (CP2degree[actor]>0)
		{
			for (alter in c(1:CP2degree[actor]))
		{
			CIndineig <- CIndineig + ystar[ CP2EdgeList[CP2RowIn[actor]+alter-1,2] ]*CP2EdgeList[CP2RowIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		
		if (directed)
		{
			if (CP2degreeIn[actor]>0)
		{
			for (alter in c(1:CP2degreeIn[actor]))
		{
			CIndineig <- CIndineig + ystar[ CP2EdgeListIn[CP2RowInIn[actor]+alter-1,2] ]*CP2EdgeListIn[CP2RowInIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		}

		}
		
		### transitive contagion
		if (doTrans)
		{
			# if y[i] is toggled, the change statistic is counts of:
			# y[translist[indexF,2]]*y[translistF[index,3]]
			# y[translist[indexM,1]]*y[translist[indexM,3]]
			# y[translist[indexE,1]]*y[translist[indexE,2]]
			# for indices indexF (index first in triple), indexM (index middle in triple), indexE (index end of triple)
			# 
			copyTP <- transFirst[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],2]]==1 && ystar[translist[copyTP[rowind],3]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			copyTP <- transMiddle[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
				 # if (copyTP[rowind]>dim(translist)[1] )
				  #{
				   # print('check trans')
				    #browser()
				  #}
					if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],3]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			copyTP <- transEndnode[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],2]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			

			
		}

		
		if (is.na(neig))
		{
			warning('missing in contagion')
			browser()
		}
		statsvecstar[2] <- statsvecstar[1]*neig
		k <- 3
		if (dorec)
		{
			statsvecstar[k] <- statsvecstar[1]*Recneig
			k <- k+1
		}
		if (doindir)
		{
			statsvecstar[k] <- statsvecstar[1]*Indineig
			k <- k+1
		}
		if (doCloseindir)
		{
			statsvecstar[k] <- statsvecstar[1]*CIndineig
			k <- k+1
			
		}
		if (doTrans)
		{
			statsvecstar[k] <- statsvecstar[1]*TIndineig
			k <- k+1
			}
		for (t in c(1:(p1 - 1) ))
		{
			statsvecstar[k] <- statsvecstar[1]*covariates[actor,t]
			k <- k+1
		}
		if (!is.null(interaction))
	{
		if (k!=p)
		{
			error('dimension missmatch in statistics and effects')
		}
		statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction[1]]
	}
		
			statsvec <- statsvec+statsvecstar
			ystar[actor] <- 1 - ystar[actor]
		}
	}
	if (any(is.na(statsvec)) || any(statsvec==0))
	{
		warning('there are zero statistics')
		#	browser()
	}
statsvec
}


testOfALLAMStats <- function()
{
	load('/Users/johankoskinen/Documents/manchester admin/SOST71032 Social Network Analysis/network regression/orgTrustMat.RData')
	load('/Users/johankoskinen/Documents/manchester admin/SOST71032 Social Network Analysis/network regression/orgVars.RData')
net2_cnorm<-SmallAdj/(colSums(SmallAdj)+1)
covs <- cbind(matrix(1,dim(nodevars)[1],1),nodevars[,1],nodevars[,2])
require(network)
require(sna)
grp1 <- as.numeric(rownames(nodevars))%in%c(21:30)
grp2 <- as.numeric(rownames(nodevars))%in%c(60:70)
grp3 <- as.numeric(rownames(nodevars))%in%c(71:84)
grp4 <- as.numeric(rownames(nodevars))%in%c(175:201)
net4 <- as.network(SmallAdj[grp4,grp4])
covs <- nodevars[grp4,]
y <- as.numeric(covs[,3] > mean(covs[,3]))
sex <- covs[,1]
net4 %v% 'sex' <- sex
plot(net4,vertex.col='sex')


ADJ <- SmallAdj[grp3,grp3]
net3 <- as.network(ADJ)
covs <- nodevars[grp3,]
y <- as.numeric(covs[,3] > mean(covs[,3]))
sex <- covs[,1]
net3 %v% 'sex' <- sex
plot(net3,vertex.col='sex')

ALAAMobj <- prepALAAMdata(y=y,ADJ=ADJ,covariates=covs[,1],directed=TRUE,useDegree = TRUE)

ALAAMobj$EdgeList[ALAAMobj$RowIn[2]:(ALAAMobj$RowIn[2]+ALAAMobj$degree[2]-1),]
ALAAMobj$EdgeList[ALAAMobj$RowIn[3]:(ALAAMobj$RowIn[3]+ALAAMobj$degree[3]-1),]
#ALAAMobj <- list(EdgeList = EdgeList, RowIn =RowIn, degree = degree, covariates= covariates)
statsvec <- getALAAMstats(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,directed=TRUE)

theta <- matrix(0,1,4)
theta[2] <- .1
theta[1] <- -.5
simSt <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein, covariates = ALAAMobj$covariates,directed=TRUE
,NumIterations=1000,theta,statsvec=statsvec ,DoSave=TRUE,thinning = 100)


simSt <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates
,NumIterations=1000,theta,statsvec=statsvec ,DoSave=FALSE,directed=TRUE)
#### initialise:
ans.log.Int <- glm(y~ ALAAMobj$covariates, family = binomial(link = "logit"))
summary(ans.log.Int)
p <- 4
theta[3:p] <- ans.log.Int$coefficients[2:(p-1)]
theta[1] <- ans.log.Int$coefficients[1]
theta[2] <- 0
simSt <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein, covariates = ALAAMobj$covariates,directed=TRUE
,NumIterations=1000,theta,statsvec=statsvec ,DoSave=TRUE,thinning = 100)
PropSigma <- solve(cov(t(simSt)))
require('mvtnorm')
sigma.epsilon <-(1/sqrt(p))*PropSigma

#### begin estimation
Iterations <- 1000
ThetaCorrent <- matrix(0,Iterations,p)
ThetaCorrent[1,] <- theta
for (iteration in c(2:Iterations))
	{
theta1 <- rmvnorm(1,ThetaCorrent[iteration-1,], sigma = sigma.epsilon)
statStar <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates
,NumIterations=1000,theta=theta1,statsvec=statsvec ,DoSave=FALSE,directed=TRUE)
Hratio <- (ThetaCorrent[iteration-1,] - theta1) %*% ( statStar-statsvec )
ThetaCorrent[iteration,] <- ThetaCorrent[iteration-1,]
if (Hratio >= log(runif(1))) {
                
                ThetaCorrent[iteration,] <- theta1
               
            }
}

####### Try the s50 network:

ADJ <- as.matrix(read.table("/Users/johankoskinen/Documents/manchester admin/bela/s50net.txt"))
ADJ[is.na(ADJ)] <- 0
smoke <- as.matrix(read.table("/Users/johankoskinen/Documents/manchester admin/bela/smoke.txt"))
alco <- as.matrix(read.table("/Users/johankoskinen/Documents/manchester admin/bela/alco.txt"))
drugs <- as.matrix(read.table("/Users/johankoskinen/Documents/manchester admin/bela/drugs.txt"))


Thetas <- BayesALAAM(y=smoke,ADJ=ADJ,directed=TRUE,covariates=as.matrix(colSums(ADJ)),NumIterations=1000,useDegree = TRUE,silent=FALSE)

Thetas2 <- BayesALAAM(y=alco,ADJ=ADJ,directed=TRUE,covariates=as.matrix(colSums(ADJ)),NumIterations=1000,useDegree = TRUE,silent=FALSE)

Thetas3 <- BayesALAAM(y=drugs,ADJ=ADJ,directed=TRUE,covariates=as.matrix(colSums(ADJ)),NumIterations=1000,useDegree = TRUE,silent=FALSE)


effname <- c('density','conagion','sender','receiver')
par(mfrow=c(2,2))
for (k in c(1:4))
{
	hist(Thetas3[200:1000,k],main=effname[k])
}
}

tuneALAAMprop <- function(ALAAMobj,directed=FALSE,theta=NULL,interaction=NULL,logitOnly=FALSE,contagion ='simple')
{
	n <- length(ALAAMobj$y)
	p <- dim(ALAAMobj$covariates)[2]+1# add 1 for intercept
	p1 <- p
	pstart <- 2
	if ('simple' %in% contagion){
		p <- p+1
		pstart <- pstart +1
		}
		
	
		dorec <- FALSE
	if ('recip' %in% contagion){
		p <- p+1
		pstart <- pstart +1
		dorec <- TRUE
				
		}
		doindir <- FALSE
	if ('indirect' %in% contagion)
		{
			p <- p+1
			doindir <- TRUE
		pstart <- pstart +1
		 }
		 
doCloseindir <- FALSE
	if ('closedind' %in% contagion)
		{
			p <- p+1
			doCloseindir <- TRUE
		pstart <- pstart +1
		
	}
	doTrans <- FALSE
	if ('transitive' %in% contagion)
		 	{
		 		p <- p+1
		 		doTrans <- TRUE
		 		pstart <- pstart +1
	 
		 		}		 

	p2 <- 0
	if (!is.null(interaction))
	{
		p2 <- length(interaction)
		p <- p+p2
	}

	ans.log.Int <- glm(ALAAMobj$y~ ALAAMobj$covariates, family = binomial(link = "logit"))
	theta[pstart:(p-p2)] <- ans.log.Int$coefficients[2:length(ans.log.Int$coefficients)]
	if (any(ans.log.Int$fitted.values==1) || ans.log.Int$converged==FALSE)
	{
		ans.log.Int <- glm(ALAAMobj$y~ 1, family = binomial(link = "logit"))
		theta[pstart:(p-p2)] <-0
	}
	
	theta[1] <- ans.log.Int$coefficients[1]
	if (any(is.na(theta)))
	{
		print('you have an ill-defined value of theta - this will be set to 0')
	
		theta[is.na(theta)] <- 0
	}
	if (!is.null(interaction))
	{
		theta[p] <- 0

	}
	
	if (logitOnly==FALSE)
	{
	simSt <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein, covariates = ALAAMobj$covariates,directed=directed
,NumIterations=1000,theta,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,thinning = 100,interaction=interaction, contagion =contagion,
RecEdgeList=ALAAMobj$RecEdgeList,
Recdegree=ALAAMobj$Recdegree,
RecRowIn=ALAAMobj$RecRowIn,
 P2EdgeList=ALAAMobj$P2EdgeList,
 P2degree=ALAAMobj$P2degree,
 P2RowIn=ALAAMobj$P2RowIn, 
 P2EdgeListIn=ALAAMobj$P2EdgeListIn,
 P2degreeIn=ALAAMobj$P2degreeIn,
 P2RowInIn=ALAAMobj$P2RowInIn,
		CP2EdgeList = ALAAMobj$CP2EdgeList,
		CP2degree = ALAAMobj$CP2degree,
		CP2RowIn = ALAAMobj$CP2RowIn,
		CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
	CP2degreeIn = ALAAMobj$CP2degreeIn,
	CP2RowInIn = ALAAMobj$CP2RowInIn,
	translist = ALAAMobj$translist,
		 		transFirst = ALAAMobj$transFirst,
		 		transMiddle = ALAAMobj$transMiddle,
		 		transEndnode = ALAAMobj$transEndnode)

if (det(cov(t(simSt))) <0.0000001)
{
	warning('model near singular - proposal will be set to default diagonal')
PropSigma <- diag(p)*0.01
}
if (det(cov(t(simSt))) >0.0000001)
{
PropSigma <- solve(cov(t(simSt)))
}

if (any( diag(PropSigma)<0 ) )
{
	print('proposal covariance ill-conditioned - default set to diagonal 0.01')
	PropSigma <- diag(p)*0.01
}
}

	if (logitOnly==TRUE)
	{
		PropSigma <- NULL
		}
tuneCons <- list(PropSigma = PropSigma, theta=theta)
tuneCons
}

BayesALAAM <- function(y,ADJ,directed=FALSE,covariates,Iterations=1000,useDegree = FALSE,silent=FALSE, PropSigma=NULL,scaling = 1, initcontagion = NULL, burnin = NULL,interaction=NULL,missingCovs = NULL,missingPhi =NULL,priorSigma=NULL,priorMu=NULL,scalePrior=NULL,contagion ='simple',canchange=NULL,MPLE=FALSE, saveFreq=NULL,missFreq=100)
{
	checkMiss <- FALSE
	require(MASS)
	require('mvtnorm')
	require('coda')
	
	ignore.cont <- FALSE
	if ('none' %in% contagion)
	{
		# run independent 
MPLE <- TRUE
ignore.cont <- TRUE
contagion <- 'simple'
	}

	
if (is.null(saveFreq))
{	saveFreq <- 100
}
	if (is.null(burnin))
	{
	burnin <- 3000
	}
	if (any(is.na(y)))
	{
		# set missing to 0 and add loop
		canchangeMiss <- which(is.na(y))
		y[is.na(y)] <- 0
		print(paste('you have ',length(canchangeMiss),' missing entries in response'))
		print(canchange)
		bias <- matrix(0,length(y),1)
		
		missFreq <- missFreq
		save.miss.imp <- matrix(0,length(y),missFreq)
		
		put.imp.here <- 1
		save.miss.imp[,put.imp.here ] <- y
		put.imp.here <- put.imp.here + 1
	}
	else {canchangeMiss <- NULL}
	
	ALAAMobj <- prepALAAMdata(y=y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =contagion, interaction= interaction)

	# statsvec <- getALAAMstats(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,directed=directed,interaction=interaction, contagion =contagion,
	# RecEdgeList=ALAAMobj$RecEdgeList,
	# Recdegree=ALAAMobj$Recdegree,
	# RecRowIn=ALAAMobj$RecRowIn, 
	# P2EdgeList=ALAAMobj$P2EdgeList,
	# P2degree=ALAAMobj$P2degree,
	# P2RowIn=ALAAMobj$P2RowIn, 
	# P2EdgeListIn=ALAAMobj$P2EdgeListIn,
	# P2degreeIn=ALAAMobj$P2degreeIn,
	# P2RowInIn=ALAAMobj$P2RowInIn)
	
	statsvec <- getALAAMstats(ALAAMobj , contagion =contagion)
	
	ALAAMobj$statsvec <- statsvec
	#ALAAMobj$directed <- directed
	#ALAAMobj$canchangeMiss <- NULL
	#ALAAMobj$canchange <- NULL
	if (!is.null(canchangeMiss))
	{
		ALAAMobj$canchangeMiss <- canchangeMiss
		
	}
	
	if (!is.null(canchange))
	{
	ALAAMobj$non.fix <- canchange
	}
	#if (is.null(canchange))
	#{
	#	tem.obj <- NULL
	#ALAAMobj$canchange <- tem.obj
	#}

	
	
	# clean up
	 if (!is.null(ALAAMobj$canchangeMiss) && checkMiss==TRUE){
		 		 	print('check canchange here at the beginning')
		 		 	browser()
		 		 }
	## === below reported in estimate.alaam instead === ####
#	print(paste('you have chosen directed equals ',ALAAMobj$directed))
#	print(paste('you have chosen degree equals ',useDegree))
	
	p <- length(ALAAMobj$statsvec)
	print(paste('you have p: ', p))
	cat('\nobserved stats: ', round(ALAAMobj$statsvec ,3) ) 
	theta <- matrix(0,1,p)
	cat('\nnumber of covariates: ', dim(ALAAMobj$covariates)[2] ) 

#browser()
	tuneCons <- tuneALAAMprop(ALAAMobj,directed=ALAAMobj$directed,theta=theta,interaction=interaction, contagion =contagion)
	if (is.null(PropSigma))
	{
	sigma.epsilon <-(scaling/sqrt(p))*tuneCons$PropSigma
	}
	if (is.null(PropSigma)==FALSE)
	{
		sigma.epsilon <-(scaling/sqrt(p))*PropSigma
	}
	if (!is.null(initcontagion))
	{
		tuneCons$theta[2] <- initcontagion
	}
	if (!is.null(missingPhi))
	{
		# this used to be detemined by  (!is.null(missingCovs))
		# but note that only missingPhi is needed
		print(paste('MNAR mechanism specified based on response and ', dim(missingCovs)[2],' covariates'))
		# define logit(p_miss_1) = missingPhi[1]*y[i]+missingPhi[2]*missingCovs[i,1]
		# define logit(p_miss_0) = missingPhi[1]*y[i]+missingPhi[2]*missingCovs[i,1]
		#bias <- missingPhi[2]+log(1+exp(missingPhi[1]+missingPhi[3]* colSums( ADJ )))-log(1+exp(missingPhi[1]+missingPhi[2]+missingPhi[3]* colSums( ADJ )))
		
		for (i in c(1:length(y)))
		{
			# missingPhi[1]*y[i]
			# missingPhi[2] is intercept
			# missingPhi[3:dim(missingCovs)[2]] are other covariates
			if (i %in% canchangeMiss)
			{
			bias[i] <-  missingPhi[1]
			}
		}
		
	}
	#if (is.null(missingCovs))
	#{
	#	bias <- NULL
	#	}
	if (!is.null(scalePrior))
	{
		cat('\npseudo-conjugate prior will be used with scaling factor: ',scalePrior)
		if (is.null(priorMu))
		{
			# use 0 prior mean
			priorMu <- matrix(0,1,p)
		}
		zStats <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=100,theta=priorMu,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,directed=ALAAMobj$directed,interaction=interaction, thinning = 1000,burnin = 3000, contagion =contagion,
		RecEdgeList=ALAAMobj$RecEdgeList,
		Recdegree=ALAAMobj$Recdegree,
		RecRowIn=ALAAMobj$RecRowIn, 
		P2EdgeList=ALAAMobj$P2EdgeList,
		P2degree=ALAAMobj$P2degree, 
		P2RowIn=ALAAMobj$P2RowIn, 
		P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
		P2degreeIn=ALAAMobj$P2degreeIn, 
		P2RowInIn=ALAAMobj$P2RowInIn ,
		CP2EdgeList = ALAAMobj$CP2EdgeList,
		CP2degree = ALAAMobj$CP2degree,
		CP2RowIn = ALAAMobj$CP2RowIn,
		CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
	CP2degreeIn = ALAAMobj$CP2degreeIn,
	CP2RowInIn = ALAAMobj$CP2RowInIn,
	translist = ALAAMobj$translist,
		 		transFirst = ALAAMobj$transFirst,
		 		transMiddle = ALAAMobj$transMiddle,
		 		transEndnode = ALAAMobj$transEndnode)
 priorSigma <- solve(cov(t(zStats)))*scalePrior
	}
	
	if (!is.null(priorSigma))
{
	require('mvtnorm')
}
# cleanup
if (checkMiss==FALSE){
rm(ADJ,covariates,y)
rm(statsvec)
	rm(directed)
	rm(canchange)
}
	#### begin estimation

ThetaCorrent <- matrix(0,Iterations,p)

theta <- tuneCons$theta
ThetaCorrent[1,] <- theta 

 if (!is.null(ALAAMobj$canchangeMiss) && checkMiss==TRUE){
		 		 	print('check canchange here at the middle')
		 		 	#browser()
		 		 }
#### allow for MPLE ###########
if (MPLE==TRUE)
{
	# 
	currentlike <- 0
	obsstats <- changestatsALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=100,theta=priorMu,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,directed=ALAAMobj$directed,interaction=interaction, thinning = 1000,burnin = 3000, contagion =contagion,
		RecEdgeList=ALAAMobj$RecEdgeList,
		Recdegree=ALAAMobj$Recdegree,
		RecRowIn=ALAAMobj$RecRowIn, 
		P2EdgeList=ALAAMobj$P2EdgeList,
		P2degree=ALAAMobj$P2degree, 
		P2RowIn=ALAAMobj$P2RowIn, 
		P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
		P2degreeIn=ALAAMobj$P2degreeIn, 
		P2RowInIn=ALAAMobj$P2RowInIn ,
		CP2EdgeList = ALAAMobj$CP2EdgeList,
		CP2degree = ALAAMobj$CP2degree,
		CP2RowIn = ALAAMobj$CP2RowIn,
		CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
	CP2degreeIn = ALAAMobj$CP2degreeIn,
	CP2RowInIn = ALAAMobj$CP2RowInIn,
	translist = ALAAMobj$translist,
		 		transFirst = ALAAMobj$transFirst,
		 		transMiddle = ALAAMobj$transMiddle,
		 		transEndnode = ALAAMobj$transEndnode)
		 		
		 		yobs <- ALAAMobj$y
		 		
if (!is.null(ALAAMobj$ALAAMobj$non.fix))
{
	obsstats <- obsstats[,ALAAMobj$non.fix]
	yobs <- yobs[ALAAMobj$non.fix]
	
}
		 		
		 	
obsstats <- t(obsstats)
n<- length(yobs)
# y,x,theta,n,p,variable=TRUE
# drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)
# if there are missing
# missingUpdate$y <- yobs
# missingUpdate$y[ALAAMobj$canchangeMiss] <- drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)[ALAAMobj$canchangeMiss]
# 
if (ignore.cont )
{
  theta[2] <- 0
}
if (!is.null(ALAAMobj$canchangeMiss) ){
  missingUpdate <- list()
  missingUpdate$y <- yobs
  missingUpdate$y[ALAAMobj$canchangeMiss] <- drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)[ALAAMobj$canchangeMiss]
  
    missFreq <- ceiling(Iterations /missFreq )
    save.miss.imp[,1 ] <- missingUpdate$y
  
}		 	

	currentlike <-straightlike(y=yobs,x=obsstats,theta=theta ,n=n)
		 		cat('Estimation using pseudo likelihood:\n')
		 		
		 		
		 for (iteration in c(2:Iterations))
	{
theta1 <- rmvnorm(1,ThetaCorrent[iteration-1,], sigma = sigma.epsilon)
if (ignore.cont )
{
  theta1[2] <- 0
}
likestar <- straightlike(y=yobs,x=obsstats,theta=theta1 ,n=n)
Hratio <- likestar-currentlike
if (Hratio >= log(runif(1))) {
                
                theta <- theta1
               currentlike <- likestar
            }
ThetaCorrent[iteration,] <- theta


### update missing if present
if (!is.null(ALAAMobj$canchangeMiss)){

  
  missingUpdate$y[ALAAMobj$canchangeMiss] <- drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)[ALAAMobj$canchangeMiss]
 
  
}

if (!is.null(ALAAMobj$canchangeMiss) & ((iteration %% missFreq)==0) )
{
  put.imp.here <- min(put.imp.here,dim(save.miss.imp)[2])
  save.miss.imp[,put.imp.here ] <- missingUpdate$y
  put.imp.here <- put.imp.here + 1
}

if ( (silent==FALSE) && ((iteration %% saveFreq)==0) )
{
  cat('\nyou have done ',iteration,' iterations out of ',Iterations,' \ntheta:', round(ThetaCorrent[iteration,] ,3) ) 	 
}

}		
	
}


################################
if (MPLE==FALSE)
{
  
  if (!is.null(ALAAMobj$canchangeMiss) )
  {
    missFreq <- ceiling(Iterations /missFreq )
  }
for (iteration in c(2:Iterations))
	{
theta1 <- rmvnorm(1,ThetaCorrent[iteration-1,], sigma = sigma.epsilon)

if (checkMiss==TRUE){
tempALAAMobj <- prepALAAMdata(y=ALAAMobj$y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =contagion)
		 		statsvec <- getALAAMstats(tempALAAMobj , contagion =contagion)
		 		#dicrepancy <- sum(statsvec != ALAAMobj$statsvec)
		 		dicrepancy <- sum( abs(statsvec - ALAAMobj$statsvec)>0.0001)

		 		if (dicrepancy!=0)
		 		{
		 			print('missmatch between observed and simulated stats in main loop')
		 			browser()
		 		}
		 	}

statStar <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=burnin,theta=theta1,statsvec=ALAAMobj$statsvec ,DoSave=FALSE,directed=ALAAMobj$directed, interaction=interaction,  contagion =contagion , 
RecEdgeList=ALAAMobj$RecEdgeList,
		Recdegree=ALAAMobj$Recdegree,
		RecRowIn=ALAAMobj$RecRowIn, 
		P2EdgeList=ALAAMobj$P2EdgeList,
		P2degree=ALAAMobj$P2degree, 
		P2RowIn=ALAAMobj$P2RowIn, 
		P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
		P2degreeIn=ALAAMobj$P2degreeIn, 
		P2RowInIn=ALAAMobj$P2RowInIn,
		CP2EdgeList = ALAAMobj$CP2EdgeList,
		CP2degree = ALAAMobj$CP2degree,
		CP2RowIn = ALAAMobj$CP2RowIn,
		CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
	CP2degreeIn = ALAAMobj$CP2degreeIn,
	CP2RowInIn = ALAAMobj$CP2RowInIn,
	translist = ALAAMobj$translist,
		 		transFirst = ALAAMobj$transFirst,
		 		transMiddle = ALAAMobj$transMiddle,
		 		transEndnode = ALAAMobj$transEndnode,
		 		canchange = ALAAMobj$non.fix)
		 		
		 		 if (!is.null(ALAAMobj$canchangeMiss) && checkMiss==TRUE){
		 		 	print('check canchange')
		 		 	browser()
		 		 }
if (!is.null(priorSigma))
{
	pr <- dmvnorm(rbind(theta1, ThetaCorrent[iteration-1,]), 
                          mean = priorMu, 
                          sigma = priorSigma,
                          log=TRUE)
    prirat <- pr[1] - pr[2]
}
if (is.null(priorSigma))
{
	prirat <- 0
}
Hratio <- (ThetaCorrent[iteration-1,] - theta1) %*% ( statStar-ALAAMobj$statsvec )+prirat
ThetaCorrent[iteration,] <- ThetaCorrent[iteration-1,]
if (Hratio >= log(runif(1))) {
                
                ThetaCorrent[iteration,] <- theta1
               
            }
            ### update missing if present
            if (!is.null(ALAAMobj$canchangeMiss)){
           # supplying values to canchange=NULL,returnNet=FALSE
            missingUpdate <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,
            EdgeListIn=ALAAMobj$EdgeListIn,
            RowInIn=ALAAMobj$RowInIn,
            degreein = ALAAMobj$degreein,
            covariates = ALAAMobj$covariates,
            NumIterations=burnin,
            theta=ThetaCorrent[iteration,],
            statsvec=ALAAMobj$statsvec ,
            DoSave=FALSE,
            directed=ALAAMobj$directed,
            canchange=ALAAMobj$canchangeMiss ,
            returnNet=TRUE,
            interaction=interaction,
            bias = bias,  
            contagion =contagion,
            RecEdgeList=ALAAMobj$RecEdgeList,
		Recdegree=ALAAMobj$Recdegree,
		RecRowIn=ALAAMobj$RecRowIn, 
		P2EdgeList=ALAAMobj$P2EdgeList,
		P2degree=ALAAMobj$P2degree, 
		P2RowIn=ALAAMobj$P2RowIn, 
		P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
		P2degreeIn=ALAAMobj$P2degreeIn, 
		P2RowInIn=ALAAMobj$P2RowInIn,
		CP2EdgeList = ALAAMobj$CP2EdgeList,
		CP2degree = ALAAMobj$CP2degree,
		CP2RowIn = ALAAMobj$CP2RowIn,
		CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
	CP2degreeIn = ALAAMobj$CP2degreeIn,
	CP2RowInIn = ALAAMobj$CP2RowInIn,
	translist = ALAAMobj$translist,
		 		transFirst = ALAAMobj$transFirst,
		 		transMiddle = ALAAMobj$transMiddle,
		 		transEndnode = ALAAMobj$transEndnode)
		 		if (checkMiss==TRUE){
		 		dicrepancy <- sum(missingUpdate$y[!(c(1:length(ALAAMobj$y)) %in% ALAAMobj$canchangeMiss)]!=ALAAMobj$y[!(c(1:length(ALAAMobj$y)) %in% ALAAMobj$canchangeMiss)])
		 		
		 		if (dicrepancy!=0)
		 		{
		 			print('missmatch between observed and simulated')
		 		}
		 		
		 		tempALAAMobj <- prepALAAMdata(y=missingUpdate$y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =contagion)
		 		statsvec <- getALAAMstats(tempALAAMobj , contagion =contagion)
		 		dicrepancy <- sum( abs(statsvec - missingUpdate$statsvec)>0.0001)
		 		if (dicrepancy!=0)
		 		{
		 			print('missmatch between observed and simulated stats')
		 			browser()
		 		}
            }
            ALAAMobj$y <- missingUpdate$y
            ALAAMobj$statsvec <- missingUpdate$statsvec
           
            }
            #### finished updating missing
            
  if ( (silent==FALSE) && ((iteration %% saveFreq)==0) )
{
	
	cat('\nyou have done ',iteration,' iterations out of ',Iterations,' \ntheta:', round(ThetaCorrent[iteration,] ,3) ) 	
	save(ThetaCorrent,file='BayesALAAMdump.RData')
	if (!is.null(ALAAMobj$canchangeMiss)){
		cat('\nobserved stats: ', round(ALAAMobj$statsvec ,3),'\n' ) 
		cat('\nimputed ones: ', sum(missingUpdate$y[ALAAMobj$canchangeMiss]),' out of ',length(ALAAMobj$canchangeMiss),'\n' ) 
	}
	
	
  }

# save imputed missing observations
if (!is.null(ALAAMobj$canchangeMiss) & ((iteration %% missFreq)==0) )
{
  put.imp.here <- min(put.imp.here,dim(save.miss.imp)[2])
  save.miss.imp[,put.imp.here ] <- ALAAMobj$y
  put.imp.here <- put.imp.here + 1
}


}
}

ResTab <- matrix(0,p,5)

ResTab[,1] <- colMeans(ThetaCorrent[floor(.1*Iterations):Iterations,])
ResTab[,2] <- apply(ThetaCorrent[floor(.1*Iterations):Iterations,],2,sd)
ResTab[,3] <- effectiveSize(ThetaCorrent[floor(.1*Iterations):Iterations,])
for (k in c(1:p)){
ResTab[k,4:5] <- as.numeric(acf(ThetaCorrent[,k],plot=FALSE)[c(10,30)][[1]])
}
ResTab <- as.table(ResTab)
colnames(ResTab) <- c('mean','sd','ESS','SACF 10','SACF 30')
rownames(ResTab)[1] <- c('intercept')
colnames(ThetaCorrent) <- rep('eff',dim(ThetaCorrent)[2])
colnames(ThetaCorrent)[1]<-c('intercept')
k <- 2
if ('simple' %in% contagion)
{
	rownames(ResTab)[k] <- c('contagion')
	colnames(ThetaCorrent)[k]<-c('contagion')
	k<- k+1
	
}

if ('recip' %in% contagion)
{
	rownames(ResTab)[k] <- c('recip cont')
	colnames(ThetaCorrent)[k]<-c('recip cont')
	k<- k+1
}

if ('indirect' %in% contagion)
{
	rownames(ResTab)[k] <- c('indirect cont')
	colnames(ThetaCorrent)[k]<-c('indirect cont')
	k<- k+1
}

if ('closedind' %in% contagion)
		{
	rownames(ResTab)[k] <- c('closed indirect cont')
	colnames(ThetaCorrent)[k]<-c('closed indirect cont')
	k<- k+1
}
	doTrans <- FALSE
if ('transitive' %in% contagion)
		 {
	rownames(ResTab)[k] <- c('transitiv cont')
	colnames(ThetaCorrent)[k]<-c('transitiv cont')
	k<- k+1
}

#if (useDegree){
#rownames(ResTab)[k] <- 'outdegree'
#k <- k+1
#} 

if (is.null(colnames(ALAAMobj$covariates)))
{
for (i in c(1:dim(ALAAMobj$covariates)[2] ) ) {
if (useDegree & i==1){
rownames(ResTab)[k] <- 'outdegree'
colnames(ThetaCorrent)[k]<-  'outdegree'

k <- k+1
} 
if (useDegree==FALSE & i>1){
rownames(ResTab)[k] <- paste('cov ',i)
colnames(ThetaCorrent)[k]<- paste('cov ',i)
k <- k+1
}
}

}
if (!is.null(colnames(ALAAMobj$covariates)) )
{
if (useDegree && is.null(colnames(ALAAMobj$covariates)[1] ))
{
colnames(ALAAMobj$covariates)[1] <- 'outdegree'

}
for (i in c(1:dim(ALAAMobj$covariates)[2] ) ){
rownames(ResTab)[k] <- colnames(ALAAMobj$covariates)[i]
colnames(ThetaCorrent)[k]<- colnames(ALAAMobj$covariates)[i]
k <- k +1
}

}

if (!is.null(interaction))
	{
		p2 <- length(interaction)
		for (i in c(1:p2))
		{
		  if (is.null(colnames(ALAAMobj$covariates) ) ) {
			rownames(ResTab)[k] <- paste('interaction ',interaction[i])
			colnames(ThetaCorrent)[k]<- paste('interaction ',interaction[i])
		  }
			if (!is.null(colnames(ALAAMobj$covariates))){
			  rownames(ResTab)[k] <- paste('interaction ',interaction[i],': ',colnames(ALAAMobj$covariates)[interaction[i]])
			  colnames(ThetaCorrent)[k]<- paste('interaction ',colnames(ALAAMobj$covariates)[interaction[i]])
			}
		  
			k <- k +1
		}
		
	}



cat('\nsummaries of the posterior draws:\n')
print(ResTab)
cat('\nif ESS is small you need to increase the number of Iterations')
cat('\nor adjust the proposals by improving the proposal variance-covariance:')
cat('\ne.g. set PropSigma equal to the covariance of the posterior draws, or')
cat('\nincrease the scaling of the porposals (defaul: scaling=1)')

if (is.null(ALAAMobj$canchangeMiss))
{
  # there is no missing
  imputed.obs <- NULL
}

if (!is.null(ALAAMobj$canchangeMiss))
{
  # there are missings
  imputed.obs <- save.miss.imp
}

Results <- list(Thetas=ThetaCorrent, 
                ALAAMobj = ALAAMobj , 
                sigma.epsilon = sigma.epsilon, 
                ResTab = ResTab , 
                priorSigma = priorSigma, 
                priorMu = priorMu,
                imputed.obs =imputed.obs )


	Results
}

sampleMuSigma <- function(p,Thetas,nGroup,priorDf,priorSigma,priorKappa,priorMu)
		{		
		##@protectedInverse internal sampleMuSigma inverse of p.s.d matrix
		protectedInverse <- function(x)
			{
				if (inherits(try(xinv <- chol2inv(chol(x)),
					silent=TRUE), "try-error"))
				{
			# Now make this x positive definite, if it is not.
			# See above for a more extensive treatment of the same.
			# Adapted from function make.positive.definite in package corpcor
			# which uses a method by Higham (Linear Algebra Appl 1988)
			# but changed to make the matrix positive definite (psd)
			# instead of nonnegative definite.
			# The idea is to left-truncate all eigenvalues to delta0.
			# The construction with tol, not used now,
			# is to ensure positive definiteness given numerical inaccuracy.
					es <- eigen(x)
					esv <- es$values
					delta0 <- 1e-6
					cat("Eigenvalues Sigma = ", sort(esv), "\n")
					if (min(esv) < delta0)
					{
						delta <- delta0
						tau <- pmax(delta, esv)
#						cat("Smallest eigenvalue of Sigma now is ",
#									min(esv),"; make posdef.\n")
						xinv <- es$vectors %*% 
									diag(1/tau, dim(x)[1]) %*% t(es$vectors)
					}
				}
				xinv
			}
		invWish <- function(v,S){
		# Draw from the inverse Wishart distribution with df v
		# and scale matrix S
		# inefficient if drawn multiply for the same S
		protectedInverse(rWishart(1,v,protectedInverse(S))[,,1])
		}

		muhat <- matrix( c(0), p, 1)
	
		muhat <- rowMeans( Thetas )# the average across groups
		matQ <- (nGroup - 1)*cov(t(Thetas))
# prior Lambda = z$priorDf*z$priorSigma
		SigmaTemp <- invWish(priorDf + nGroup,
				priorDf*priorSigma + matQ +
				(priorKappa*nGroup/(priorKappa + nGroup))*
				tcrossprod( muhat - priorMu, muhat - priorMu ) )
		muTemp <- t(chol( ( priorKappa + nGroup )^(-1)*SigmaTemp )) %*%
			rnorm( p , 0 , 1 ) +
			(nGroup/( priorKappa + nGroup ))*muhat +
			(priorKappa/( priorKappa + nGroup ))*priorMu
			
		MuAndSigma <- cbind(muTemp,SigmaTemp)
		MuAndSigma
}
	
multiALAAM <- function(y,ADJ,covariates,useDegree=FALSE,directed=FALSE, silent=FALSE,Iterations,saveFreq = 100,MHtuning=NULL,priorMu=NULL, priorSigma=NULL, priorDf=NULL, priorKappa=NULL,burnin=2000)
{
	require('mvtnorm')
	require('coda')
	# to interact contagion and network-level effect
	# set p2 to 1 and let eta be the parameter
	# this assumes theta*neigh = (eta*x_j)*neigh
	M <- length(y)
	networklist <- vector('list',M)
	if (useDegree)
	{
		p <- 3 + dim(covariates[[1]])[2]
	}
	if (useDegree==FALSE)
	{
		p <- 2 + dim(covariates[[1]])[2]
	}
	#### set up priors
	if (is.null(priorMu))
	{
		priorMu <- matrix( c(0), p, 1 )
		
	}
	
	if (is.null(priorKappa))
	{
		
			priorKappa <- 1
		
	}
	
	if (is.null(priorDf))
	{
		
			priorDf <- p+2
		
	}
	if (priorDf + M <= p)
	{
		priorDf <- p - M + 2
	}
	
	if (is.null(priorSigma))
	{
		priorSigma <- diag(p)
	}
	### DONE setting up priors
	Thetas <- array(0 , dim=c(p,M,Iterations))
	doTune <- FALSE
		if (is.null(MHtuning))
		{
			doTune <- TRUE
			MHtuning <- vector('list',M)
			cat('\ndefault tuning of MCMC initialised\n')
		}
	### go through groups and format
	for (group  in c(1:M))
		{
			ALAAMobj <- prepALAAMdata(y=y[[group]],
			ADJ=ADJ[[group]],
			covariates=covariates[[group]],
			directed=directed,
			useDegree = useDegree)
			ALAAMobj$directed <- directed
			#statsvec <- getALAAMstats(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,directed=directed)
		statsvec <- getALAAMstats(ALAAMobj)	
		ALAAMobj$statsvec <- statsvec
			if (any(is.na(statsvec)))
			{
				print('there are missing values in stufficient statistics: check your covariates and response variable')
				
			}
			p <- length(statsvec)
			networklist[[group]] <- ALAAMobj
			
			theta <- matrix(0,1,p)
		
	if (doTune)
	{
	tuneCons <- tuneALAAMprop(ALAAMobj,directed=directed,theta=theta)
			MHtuning[[ group ]] <- (1/sqrt(p))*tuneCons$PropSigma
			Thetas[,group,1] <- tuneCons$theta
			}
			}
	cat('\nTuning done\n')
	p <- length(statsvec)
	
	SigmaPost <- array(0 , dim=c(p,p,Iterations))
	MuPost <- matrix(0, p, Iterations )
	SigmaPost[,,1 ] <- diag(p)
	####
		cat('\nMain MCMC estimation phase initiated\n')
	for (iteration in c(2:Iterations))
	{
	# update group theta
		for (group  in c(1:M))
		{
		Thetas[,group,iteration] <- upDateALAAMTheta( ALAAMobj=networklist[[group]], theta = Thetas[,group,(iteration -1)] , m.prior = MuPost[,(iteration -1)], sigma.prior = SigmaPost[,,(iteration -1)], sigma.epsilon = MHtuning[[ group ]], burnin)
		}
		# update mu and sigma
		
		MuAndSigma <- sampleMuSigma(p,Thetas = Thetas[, ,iteration ],nGroup = M,priorDf,priorSigma,priorKappa,priorMu)
MuPost[,iteration ] <- MuAndSigma[,1]
SigmaPost[,,iteration ] <- MuAndSigma[,2:(p+1)]

if ( (silent==FALSE) && ((iteration %% saveFreq)==0) )
{
		
	cat('\nyou have done ',iteration,' iterations out of ',Iterations,' \nmu:', round(t(MuPost[,iteration ]) ,2) ) 
	mup <- list(MuPost=MuPost,SigmaPost=SigmaPost,Thetas=Thetas)
	save(mup,file='TempDumpMALAAM.RData')
	
	}

}

### some print functions
ResTab <- matrix(0,p,5)

ResTab[,1] <- colMeans(t(MuPost[,floor(.1*Iterations):Iterations]))
ResTab[,2] <- apply(t(MuPost[,floor(.1*Iterations):Iterations]),2,sd)
ResTab[,3] <- effectiveSize(t(MuPost[,floor(.1*Iterations):Iterations]) )
for (k in c(1:p)){
ResTab[k,4:5] <- as.numeric(acf(MuPost[k,floor(.1*Iterations):Iterations],plot=FALSE)[c(10,30)][[1]])
}
ResTab <- as.table(ResTab)
colnames(ResTab) <- c('mean','sd','ESS','SACF 10','SACF 30')
rownames(ResTab)[1:2] <- c('intercept','contagion')
k <- 3
if (useDegree){
rownames(ResTab)[k] <- 'outdegree'
k <- k+1
} 
for (i in c(k:p)){
rownames(ResTab)[i] <- paste('cov ',i-k+1)
}

cat('\nsummaries of the posterior draws for mu:\n')
print(ResTab)
cat('\nif ESS is small you need to increase the number of Iterations')
cat('\nor adjust the proposals by improving the proposal variance-covariance:')
cat('\ne.g. set PropSigma equal to the covariance of the posterior draws, or')
cat('\nincrease the scaling of the porposals (defaul: scaling=1)')


out.put <- list(MuPost=MuPost, SigmaPost=SigmaPost, Thetas=Thetas, MHtuning=MHtuning)
}


upDateALAAMTheta <- function(ALAAMobj,theta,m.prior,sigma.prior,sigma.epsilon,burnin)
{
	
	theta1 <- rmvnorm(1,theta, sigma = sigma.epsilon)
	
	pr <- dmvnorm(rbind(theta1, theta), 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)

statStar <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates
,NumIterations=burnin,theta=theta1,statsvec=ALAAMobj$statsvec ,DoSave=FALSE,directed=ALAAMobj$directed)

# if (any(is.na(c(theta,theta1))))
# {
	# print('corrupted theta values in main update of algortihm')
	
	# browser()
# }
# if (any(is.na(pr)))
# {
	# print('corrupted prior probabilities')
	
	# browser()
# }
# if (any(is.na(statStar)) )
# {
	# print('your statistics are weird')
	
	# browser()
# }


Hratio <- (theta - theta1) %*% ( statStar-ALAAMobj$statsvec )+ pr[1] - pr[2]
# if (is.na(Hratio))
# {
	# print('bonkers')
	# browser()
	
# }

if (Hratio >= log(runif(1))) {
                
                theta <- theta1
               
            }

	 theta
}


modelSelALAAM <- function(ALAAMresult, burnin=NULL,thinning=NULL,modalPoint=NULL,numYsamps=NULL,m.prior=NULL,sigma.prior=NULL, interaction=NULL,likelihoodOnly=FALSE,thetaRef=NULL,numbridges=20,saveGOF=FALSE,yburning=3000,iterationsDenominator=200)
{
	require('mvtnorm')
	saveFreq <- 10
	sigma.epsilon <- ALAAMresult$sigma.epsilon
	ALAAMobj <- ALAAMresult$ALAAMobj
	if (likelihoodOnly==FALSE)
	{
	if (is.null(m.prior) | is.null(sigma.prior))
	{
		error('you can only do modelselection with proper priors: supply prior mean and variance')
	}
	}
	# determine number of draws to use
	Iterations <- dim(ALAAMresult$Thetas)[1]
	
		if (is.null(thinning))
		{
			thinning <- 10
		
		}
		
		if (is.null(burnin))
		{
			burnin <- 0.1
			burnin <-floor(burnin*Iterations)
		
		}
		

	thinTheta <- ALAAMresult$Thetas[seq(burnin,Iterations,by=thinning),]
	newIterations <- dim(thinTheta)[1]
	if (is.null(numYsamps))
		{
			numYsamps <-newIterations
		
		}
	if (is.null(thetaRef))
	{
		cat('\nlogistic MLE will be used for reference parameters\n ')
	tuneCons <- tuneALAAMprop(ALAAMobj,directed=ALAAMobj$directed,theta=modalPoint,interaction=interaction,logitOnly=TRUE)

	modalPointNoDep <- tuneCons$theta
	modalPointNoDep[2] <- 0
	}
	if (!is.null(thetaRef))
	{
		cat('\nusersupplied parameter will be used for reference parameters\n ')
		modalPointNoDep <- thetaRef
		
	}
		
	patEstim <-pathsampALAAM(ALAAMobj=ALAAMobj,theta_tilde=modalPoint,theta_star=modalPointNoDep,burnin=10000,numbridges= numbridges,interaction=interaction)
	
	
if (likelihoodOnly==FALSE)
	{
		tempburning <- ceiling(yburning/3)
	STATS <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=numYsamps,theta=modalPoint,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,directed=ALAAMobj$directed,thinning=tempburning ,burnin=yburning,interaction=interaction)
	}
	if (!is.null(ALAAMobj$canchange))
	{
		
	if (likelihoodOnly==FALSE)
	{
		STATSobs <- matrix(0,length(ALAAMobj$statsvec),newIterations)
		
		tempStats <- ALAAMobj$statsvec
		tempyy <-ALAAMobj$y
		for (r in c(1:newIterations))
		{
		tempStats <- simulateALAAM(y=tempyy,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates, NumIterations=1,theta=thinTheta[r,],statsvec=tempStats ,DoSave=TRUE,returnNet=TRUE,directed=ALAAMobj$directed,  thinning=yburning, burnin=yburning,canchange=ALAAMobj$canchange,interaction=interaction)
		
		 tempyy <- tempStats$y[,1]
tempStats <- tempStats$statsvec[,1]

		STATSobs[,r] <- tempStats
		}
		}
		patEstimNumerator <-pathsampALAAM(ALAAMobj=ALAAMobj,theta_tilde=modalPoint,theta_star=modalPointNoDep,burnin=10000,numbridges=numbridges,interaction=interaction, canchange=ALAAMobj$canchange)
		
		
		logLike <- patEstimNumerator - patEstim 

		
	}
	if (is.null(ALAAMobj$canchange))
	{
		logLike <- (modalPoint-modalPointNoDep) %*% ALAAMobj$statsvec -patEstim
	}
	if (likelihoodOnly==FALSE)
	{
logLike <- logLike + dmvnorm(modalPoint, 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)

}

	
# make sure that sample is good



if (likelihoodOnly==FALSE)
	{
		
		p <- dim(STATS)[1]
if (p != length(modalPoint))
{
	warning('Dimansion missmatch of simulated statistics and modal point')
}
for (k in c(1:p))# lowlevel centering of simulated statistics
{
	STATSuse <- STATS
	if (is.null(ALAAMobj$canchange))
	{
	STATSuse[k,] <- STATS[k,]-ALAAMobj$statsvec[k]
	}
		
}

M <- iterationsDenominator

	accepts <- matrix(0,newIterations,1)
	denomaccepts <- matrix(0,newIterations,1)
	vectorzeros <- matrix(0,1,iterationsDenominator)
	# calculate numerator
	for (r in c(1:newIterations))
	{
		PropProb <- dmvnorm(modalPoint,thinTheta[r,], sigma = sigma.epsilon,log=TRUE)
	
	pr <- dmvnorm(rbind(modalPoint, thinTheta[r,]), 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)
        if (!is.null(ALAAMobj$canchange))
	{
		for (k in c(1:p))
		{
		STATSuse[k,] <- STATS[k,]-STATSobs[k,r]
		}
		}
	
		
		hratios <-  pmin(vectorzeros,((thinTheta[r,] - modalPoint) %*% STATSuse + pr[1] - pr[2] ) ) +PropProb
		accepts[r] <- mean(exp(hratios))
		
	}
	# calculate denominator, average of acceptance prob
	vectorzeros <- matrix(0,1,dim(STATS)[2])
	for (r in c(1:newIterations))
	{
		# theta ~ h(theta|modalPoint)
		theta1 <- rmvnorm(1,modalPoint, sigma = sigma.epsilon)
		# calculate acceptance probability
		pr <- dmvnorm(rbind(theta1, modalPoint), 
                          mean = m.prior, 
                          sigma = sigma.prior,
                          log=TRUE)
		# draw missing values
		 if (!is.null(ALAAMobj$canchange))
	{
		STATSobs <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=1,theta=thinTheta[r,],statsvec=ALAAMobj$statsvec ,DoSave=FALSE,directed=ALAAMobj$directed,burnin=yburning,canchange=ALAAMobj$canchange,interaction=interaction)
		}
		
		
			# need draws from p( | theta)
			STATS <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=M,theta=theta1,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,directed=ALAAMobj$directed,thinning=1000,burnin=yburning,interaction=interaction)
			
		
		 if (!is.null(ALAAMobj$canchange))
	{
		for (k in c(1:p))
		{
		STATS[k,] <- STATS[k,]-STATSobs[k]
		}
		}
		
		if (is.null(ALAAMobj$canchange))
	{
		for (k in c(1:p))
		{
		STATS[k,] <- STATS[k,]-ALAAMobj$statsvec[k]
		}
		}
		
		hratios <-  pmin(vectorzeros,((modalPoint-theta1 ) %*% STATS + pr[1] - pr[2] ) )
		denomaccepts[r] <- mean(exp(hratios))
		 if ( (r %% saveFreq)==0 )
{
	
	cat('\nyou have done ',r,' iterations out of ',newIterations) 	
	}

		}
	
	}
	
	
	if (likelihoodOnly==FALSE){
	modelres <- list(relloglike = logLike,accepts = accepts, denomaccepts= denomaccepts , modalPoint=modalPoint, modalPointNoDep=modalPointNoDep)
	}
	if (likelihoodOnly==TRUE)
	{
		modelres <- list(relloglike = logLike )
	}
	
	modelres
}



plotPost <- function(ALAAMresult,figname,
                     nameaxs=NULL,
                     sizearg=NULL,
                     doTable=FALSE,
                     thinning=NULL,
                     burnin = NULL,
                     showplot = FALSE)
{
	require('xtable')
	require('coda')
	Iterations <- dim(ALAAMresult$Thetas)[1]
	p <- dim(ALAAMresult$Thetas)[2]
	HPDlowhigh <- matrix(0,p,2)
	# return plot or save plot
	if (showplot==FALSE){
	if (is.null(sizearg))
	{
	pdf(paste(figname,'.pdf',sep=''))
	}
	if (!is.null(sizearg))
	{
	pdf(paste(figname,'.pdf',sep=''),width=sizearg[1],height=sizearg[2])
	}
	}
	####
	if (is.null(nameaxs))
	{
	  nameaxs <- colnames(ALAAMresult$Thetas)
	}
	if (p <=3)
	{
		par(mfrow=c(p,3))
	}
	if (p > 3 & p<10)
	{
		par(mfrow=c(ceiling(p/2),6),cex.lab=.75,
		oma = c(0,4,0,0) + 0.1,
          mar = c(5,0,1,1) + 0.1)
	}
	for (r in c(1:p))
	{
	  if (any(diff(ALAAMresult$Thetas[,r])!=0)){
				if (is.null(nameaxs))
	{
		plot(density(ALAAMresult$Thetas[,r]),bty='n',xlab=expression(theta[r]),main='',ylab='',yaxt='n', ann=FALSE)
		}
		if (!is.null(nameaxs))
	{
		plot(density(ALAAMresult$Thetas[,r]),bty='n',xlab=nameaxs[r],main='',ylab='',yaxt='n')
	}
	plot(ALAAMresult$Thetas[,r],type='l',bty='n',ylab='',xlab='iterations',cex.axis=0.5)
	muppet <- acf(ALAAMresult$Thetas[,r],lag.max=100,plot=FALSE)
	plot(muppet,bty='n',main='',ylab='Autocorrelation',cex.axis=0.5)
	  }
	  if (any(diff(ALAAMresult$Thetas[,r])!=0)==FALSE){
	    cat('Parameter ',r,' all zeros\n')
	  }
	}
	
	if (showplot==FALSE){
	dev.off()
	}
	
	if (doTable)
	{
		if (is.null(thinning))
		{
			thinning <- 10
		
		}
		
		if (is.null(burnin))
		{
			burnin <- 0.1
			burnin <-floor(burnin*Iterations)
		
		}
		pickupsamp <- seq(burnin,Iterations,by= thinning)

		for (r in c(1:p))
		{
			tempThet <- sort(ALAAMresult$Thetas[pickupsamp,r])
			bratwurst <- length(tempThet)
		HPDlowhigh[r,1] <- tempThet[floor(bratwurst*0.025)]
		HPDlowhigh[r,2] <- tempThet[ceiling(bratwurst*(1-0.025))]

		}
		
		print(paste('using every ',thinning,' draw and burnin is set to',  burnin))
		ResTab <- matrix(0,p,7)
ResTab[,1] <- colMeans(ALAAMresult$Thetas[pickupsamp,])
ResTab[,2] <- apply(ALAAMresult$Thetas[pickupsamp,],2,sd)
ResTab[,3] <- effectiveSize(ALAAMresult$Thetas[pickupsamp,])
for (k in c(1:p)){
ResTab[k,4:5] <- as.numeric(acf(ALAAMresult$Thetas[,k],plot=FALSE)[c(10,30)][[1]])
}
ResTab[,6:7] <- HPDlowhigh

ResTab <- as.table(ResTab)
colnames(ResTab) <- c('mean','sd','ESS','SACF 10','SACF 30','2.5 perc','97.5 perc')
rownames(ResTab)[1:2] <- c('intercept','contagion')
k <- 3
if (is.null(nameaxs)){

for (i in c(k:p)){
rownames(ResTab)[i] <- paste('cov ',i-k+1)
}
}
if (!is.null(nameaxs))
		{
			rownames(ResTab) <- nameaxs
		}
	
print(xtable(ResTab, type = "latex"), file = paste(figname,".tex",sep=''))

#toLatex(ResTab)
	}
	
}




plotHPDbayes <- function(theta,figname,nameaxs=NULL,sizearg=NULL,thinning=NULL,burnin = NULL)
{

if (!is.null(dim(theta)))
{
	p <- dim(theta)[2]
	Iterations <- dim(theta)[1]

}
if (is.null(dim(theta)))
{
	p <- 1
Iterations <- length(theta)[1]
}
	if (is.null(thinning))
		{
			thinning <- 10
		
		}
		
		if (is.null(burnin))
		{
			burnin <- 0.1
			burnin <-floor(burnin*Iterations)
		
		}
		
		pickupsamp <- seq(burnin,Iterations,by= thinning)
	polyColours <- gray.colors(8,alpha=.55)
	print(p)
	if (is.null(sizearg))
	{
	pdf(paste(figname,'.pdf',sep=''))
	}
	if (!is.null(sizearg))
	{
	pdf(paste(figname,'.pdf',sep=''),width=sizearg[1],height=sizearg[2])
	}
	
	if (p==1){
		useCol <-c(polyColours[3])
	tempVars <- sort(theta[pickupsamp])
	}
	if (p>1){
		useCol <-c(polyColours[8],polyColours[4],polyColours[2])
		varMin <- matrix(0,1,p)
		varMax <- matrix(0,1,p)
		varVari <- matrix(0,1,p)
		tempVars <- matrix(0,length(pickupsamp ),p)
		heightMax <- 0
		for (k in c(1:p))
		{
	tempVars[,k] <- sort(theta[pickupsamp,k])
	varMin[,k] <- min( tempVars[,k])
	varMax[,k] <- max( tempVars[,k])
	varVari[,k] <- var( tempVars[,k])
	heightMax <- max(heightMax,density(tempVars[,k])$y)
	}
	plotOrder <- order(varVari)
	plotOrder2 <- rev(plotOrder)
	
	
	}
	
	if (p==1){
		
	HK <- density(tempVars)
	if (is.null(nameaxs)){
	plot(HK$x,HK$y,type='l',main=expression(theta[1]),ylab='',xlab='',yaxt='n',bty='n',col= polyColours[1],xlim=range(theta))
	}
	if (!is.null(nameaxs)){
	plot(HK$x,HK$y,type='l',main=nameaxs[1],ylab='',xlab='',yaxt='n',bty='n',col= polyColours[1],xlim=range(theta))
	}
	numReps <- length(tempVars)
	lowerPoint <- tempVars[floor(0.025*numReps)]
	upperPoint <- tempVars[ceiling((1-0.025)*numReps)]
	densTies <- HK$y[HK$x>lowerPoint & HK$x<upperPoint]
	densXs <- HK$x[HK$x>lowerPoint & HK$x<upperPoint]
	polygon(c(densXs[1], densXs ,densXs[length(densXs)],rev(densXs)), c( 0,densTies,0, rep(0,length(densXs)) ), col = useCol[plotOrder[r]], border = NA)

	}
	if (p>1){
	for (j in c(1:p))
	{
		
		r <- which(plotOrder2==j)
		if (j==1){
		HK <- density(tempVars[,r])
	if (is.null(nameaxs)){
	plot(HK$x,HK$y,type='l',main=expression(theta[r]),ylab='',xlab='',yaxt='n',bty='n',col= polyColours[r],xlim=range(theta),ylim=c(0,heightMax))
	}
	if (!is.null(nameaxs)){
	plot(HK$x,HK$y,type='l',main='',ylab='',xlab='',yaxt='n',bty='n',col= polyColours[r],xlim=range(theta),ylim=c(0,heightMax))
	}
	

	
	}
	
	if (j>1){
		HK <- density(tempVars[,r])
	
	lines(HK$x,HK$y,type='l',col= useCol[plotOrder[j]])
	
	}
	
	numReps <- length(tempVars[,r])
	lowerPoint <- tempVars[floor(0.025*numReps),r]
	upperPoint <- tempVars[ceiling((1-0.025)*numReps),r]
	densTies <- HK$y[HK$x>lowerPoint & HK$x<upperPoint]
	densXs <- HK$x[HK$x>lowerPoint & HK$x<upperPoint]
	polygon(c(densXs[1], densXs ,densXs[length(densXs)],rev(densXs)), c( 0,densTies,0, rep(0,length(densXs)) ), col =  useCol[plotOrder[r]], border = NA)
				
	}
	legend('topright',nameaxs,fill=useCol[plotOrder])

	}
	
	dev.off()
	
}


pathsampALAAM <- function(ALAAMobj,
                          theta_tilde,
                          theta_star,
                          burnin=3000,
                          numbridges=20,
                          interaction=interaction,
                          sampleSize=100,
                          canchange=NULL,
                          silent=FALSE)
{
	
	theta <- theta_star 
	STATS <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=sampleSize, theta= theta,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,directed=ALAAMobj$directed,interaction=interaction,returnNet=TRUE,thinning = 1000, burnin = burnin,canchange=canchange)



ALAAMobj$y <- STATS$y[,sampleSize]
ALAAMobj$statsvec <- STATS$statsvec[,sampleSize]
rm(STATS)

	PathEstim <- 0
	
	for (RepSamp in  c(0:numbridges) )
	{
    theta <- theta_star * (1 - RepSamp/numbridges) + theta_tilde * (RepSamp/numbridges)
       
    STATS <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,NumIterations=sampleSize, theta= theta,statsvec=ALAAMobj$statsvec ,DoSave=TRUE,returnNet=TRUE,directed=ALAAMobj$directed,interaction=interaction,thinning = 1000, burnin = burnin, canchange=canchange)
  ALAAMobj$y <- STATS$y[,sampleSize]
ALAAMobj$statsvec <- STATS$statsvec[,sampleSize]
    #PathEstim = PathEstim + (1/K)*(theta_tilde - theta_star) *  A_obs;
    tempest <- (1/(numbridges+1))*(theta_tilde - theta_star) %*% rowMeans(STATS$statsvec)
    rm(STATS)
    PathEstim <-  PathEstim +  tempest
    if (silent==FALSE){
    print(paste('Path at bridge: ',RepSamp,' contributes: ', round(tempest,3)))
    }
}
	
 PathEstim
	
	
}


aitkinPostDev <- function(ALAAMresult, burnin=NULL,thinning=NULL,numYsamps=NULL, interaction=NULL,thetaRef=NULL,numbridges=20,Yburnin=NULL)
{
	ALAAMobj <- ALAAMresult$ALAAMobj
	saveFreq <- 10
	# determine number of draws to use
	Iterations <- dim(ALAAMresult$Thetas)[1]
	
		if (is.null(thinning))
		{
			thinning <- 10
		
		}
		
		if (is.null(burnin))
		{
			burnin <- 0.1
			burnin <-floor(burnin*Iterations)
		
		}
		if (is.null(Yburnin))
		{
		
			Yburnin <-length(ALAAMobj$y)^2
		
		}

	thinTheta <- ALAAMresult$Thetas[seq(burnin,Iterations,by=thinning),]
	newIterations <- dim(thinTheta)[1]
	if (is.null(numYsamps))
		{
			numYsamps <-newIterations
		
		}
	if (is.null(thetaRef))
	{
		cat('\nlogistic MLE will be used for reference parameters\n ')
		modalPoint<- colMeans(thinTheta)
	tuneCons <- tuneALAAMprop(ALAAMobj,directed=ALAAMobj$directed,theta=modalPoint,interaction=interaction,logitOnly=TRUE)
	modalPointNoDep <- tuneCons$theta
	modalPointNoDep[2] <- 0
	}
	if (!is.null(thetaRef))
	{
		cat('\nusersupplied parameter will be used for reference parameters\n ')
		modalPointNoDep <- thetaRef
		
	}
	
	RelLike <- matrix(0,newIterations,1)
	
	for (k in c(1:newIterations))
		{
			modalPoint <- thinTheta[k,] 
	patEstim <-pathsampALAAM(ALAAMobj=ALAAMobj,theta_tilde=modalPoint,theta_star=modalPointNoDep,burnin=Yburnin,numbridges=numbridges,interaction=interaction,sampleSize=numYsamps,silent=TRUE)

	if (!is.null(ALAAMobj$canchange))
	{
		
	
		patEstimNumerator <-pathsampALAAM(ALAAMobj=ALAAMobj,theta_tilde=modalPoint,theta_star=modalPointNoDep,burnin=Yburnin,numbridges=numbridges,interaction=interaction, canchange=ALAAMobj$canchange,sampleSize=numYsamps,silent=TRUE)
		
		
		logLike <- patEstimNumerator - patEstim 

		
	}
	if (is.null(ALAAMobj$canchange))
	{
		logLike <- (modalPoint-modalPointNoDep) %*% ALAAMobj$statsvec -patEstim
		}
		RelLike[k] <- logLike
	 if ( (k %% saveFreq)==0 )
{
	
	cat('\nyou have done ',k,' iterations out of ',newIterations ,'\n') 	
	}
	}

RelLike	
}

independLike <- function(ALAAMobj,theta)
{
	n <- length(ALAAMobj$y)
	doActors <- c(1:n)
	if (!is.null(ALAAMobj$canchange))
	 {
		 # there is missing - this only works for MAR
		 useRows <- !(doActors %in% ALAAMobj$canchange)
		 doActors <- doActors[useRows]
	 }
	loglike <- 0
	for (k in c(1:length(doActors)))
	{
		actor <- doActors[k]
		if (ALAAMobj$y[actor]==1)
		{
		loglike <- loglike - log(1+exp(  -theta %*%  c(1, ALAAMobj$covariates[actor,]) ) )
		}
		if (ALAAMobj$y[actor]==0)
		{
		loglike <- loglike - log(1+exp(  theta %*%  c(1, ALAAMobj$covariates[actor,]) ) )
		}

		
		
		
	}
	loglike 
}

nocontagionBayes <- function(y,x,iterations,priormu,priorsigma,thin=10,numYdenom=10,numYnum=100)
{
	p <- dim(x)[2]
	n <- dim(x)[1]
	BigTheta <- matrix(0,ceiling(iterations/thin),p)
	sigmaprop <- diag(p)*.1
	loglike <- matrix(0,1,iterations)
	printfreq <- ceiling(iterations/10)
#	browser()
	loglike[1] <- straightlike(y=y,x=x,theta=BigTheta[1,],n=n)
	ans.log.Int <- glm(y~ x-1, family = binomial(link = "logit"))
	mltheta <- ans.log.Int$coefficients
	mlthetalike <- straightlike(y=y,x=x,theta=mltheta,n=n)
	chibjelnum <-  matrix(0,1,ceiling(iterations/thin) )
	
	chibjelnumSAV <-  matrix(0,1,ceiling(iterations/thin) )
	StatsVec <- matrix(0,numYnum ,p)
	z <- matrix(0,1,p)

	for (actor in c(1:n) )
	{
		z <- z + y[actor]*t(x[actor,])
		
	}
	for (k in c(1:numYnum))
	{
		StatsVec[k,] <- drawstraight(y=y,x=x,theta=mltheta,n=n,p=p,variable=FALSE)-z
	}
	StatsVec <- t(StatsVec)
	vectorones <- matrix(1,1,numYnum )
	placehere <- 1
	theta <- BigTheta[1,]
	for (k in c(2:iterations) )
	{
		
		theta1 <- rmvnorm(1,theta, sigma = sigmaprop)
		pr <- dmvnorm(rbind(theta1, theta), 
                          mean = priormu, 
                          sigma = priorsigma,
                          log=TRUE)
     templike <- straightlike(y=y,x=x,theta=theta1,n=n)
		hratio <- templike-loglike[k-1] +pr[1] - pr[2]
		loglike[k]<- loglike[k-1]
		if (hratio >= log(runif( 1 ) ) )
		{
			theta <- theta1
			loglike[k]<- templike
		}
		
		
		if ( (k %% printfreq )==0 )
		{
			print(paste('iteration ',k))
		}
		
		if ( (k %% thin)==0 )
		{
			BigTheta[placehere,] <- theta
			pr <- dmvnorm(rbind(mltheta, theta), 
                          mean = priormu, 
                          sigma = priorsigma,
                          log=TRUE)
			alpha <- min(1,  exp(mlthetalike-loglike[k] + pr[1] - pr[2]))
			propprop <- dmvnorm(mltheta, 
                          mean = theta, 
                          sigma = sigmaprop)
			chibjelnum[placehere] <- alpha * propprop
			
			
			#### alternative:
			
			hratios <- exp( (theta - mltheta) %*% StatsVec + pr[1] - pr[2])
			alpha2 <-  pmin(vectorones,hratios)
		
			
		chibjelnumSAV[placehere] <-mean( alpha2  * propprop )
			placehere <- placehere+1
			
		}
		
		
	}
	M <- dim(BigTheta)[1]
	chibjeldennum <-  matrix(0,1,M )
	chibjeldennumSAV <-  matrix(0,1,M )
StatsVec <- matrix(0,numYdenom ,p)
vectorones <- matrix(1,1,numYdenom)
	for (k in c(1:M))
	{
		theta1 <- rmvnorm(1,mltheta, sigma = sigmaprop)
		pr <- dmvnorm(rbind(theta1, mltheta), 
                          mean = priormu, 
                          sigma = priorsigma,
                          log=TRUE)
        templike <- straightlike(y=y,x=x,theta=theta1,n=n)
		hratio <- templike-mlthetalike +pr[1] - pr[2]
		chibjeldennum[k] <- min(1, exp(hratio) )
		
		
		for (petter in c(1:numYdenom))
	{
		StatsVec[petter,] <- drawstraight(y=y,x=x,theta=theta1 ,n=n,p=p,variable=FALSE)-z
	}
	
		hratios <- exp( (mltheta-theta1 ) %*% t(StatsVec) + pr[1] - pr[2])
			alpha2 <-  pmin(vectorones,hratios)
		
		chibjeldennumSAV[k] <- mean( alpha2 )

		
	}
	
	
	Results <- list(BigTheta = BigTheta,chibjelnum=chibjelnum, chibjeldennum = chibjeldennum, chibjelnumSAV=chibjelnumSAV,chibjeldennumSAV=chibjeldennumSAV)
	Results
}

straightlike <- function(y,x,theta,n)
{
	loglike <- 0
	for (actor in c(1:n) )
	{
		
		if (y[actor]==1)
		{
		loglike <- loglike - log(1+exp(  -theta %*%  x[actor,] ) )
		}
		if (y[actor]==0)
		{
		loglike <- loglike - log(1+exp(  theta %*%   x[actor,] ) )
		}

		
	}
	loglike 
}

drawstraight <- function(y,x,theta,n,p,variable=TRUE)
{
	z <- matrix(0,1,p)

	for (actor in c(1:n) )
	{
		pProb <- 1/(1+exp(  -theta %*%  x[actor,] ) )
		y[actor] <- sample(c(0,1),size=1,prob=c(1-pProb,pProb))
		z <- z + y[actor]*t(x[actor,])
		
	}
	if (variable)
	{
	y
	}
	if (variable==FALSE)
	{
	y <- z
	}
y	
}

getlistTransbase <- function(ADJ,directed=FALSE)
{
	# compute an ege-list with edges i -> j if in ADJ i->j and i->k->j
	# these edges also have values #{k: i->k->j }
	# as for indirect paths we need the reverse j -> i 
	n <- dim(ADJ)[1]
	
	# list neighbours at distance 2
	if (directed==FALSE)
	{
		error('indirect contagion is not implemented for undirected networks')
	}
	# create matrix of paths i -> k -> j
	path2 <- matrix(0,n,n)
	for (i in c(1:n))
	{
		for (j in c(1:n))
		{
			if (i != j)
			{
			path2[i,j] <- ADJ[i,] %*% ADJ[,j] * ADJ[i,j]
			}
		}
	}
	# produce edgelist for neighbours at distance 2
	P2RowIn <- matrix(NA,n,1)
	P2degree <- matrix(0,n,1)
	P2EdgeList <- which(path2>=1,arr.ind=TRUE)
	
	P2EdgeList <- cbind(P2EdgeList[order(P2EdgeList[,1]),],matrix(0,dim(P2EdgeList)[1],1) )
		 for (actor in c(1:n))
		 {
		 	TheseRows <- which(P2EdgeList[,1]==actor)
		 	if (length( P2EdgeList[TheseRows,c(1:2)])==2)
		 	{
		 	#	print('pucko')
		 	#	browser()
		 		P2EdgeList[TheseRows,3] <-  path2[P2EdgeList[TheseRows,1],P2EdgeList[TheseRows,2]]
		 	}
		 	if (length( P2EdgeList[TheseRows,c(1:2)])>2)
		 	{
		 	P2EdgeList[TheseRows,3] <- path2[P2EdgeList[TheseRows,c(1:2)]]
		 	}
		 	if (length(TheseRows)>0)
		 	{
		 		P2RowIn[actor] <- min(TheseRows)
		 		P2degree[actor] <- length(TheseRows)
		 	}
		 	
		 }
	 if (directed){
	 	path2 <- matrix(0,n,n)
	for (i in c(1:n))
	{
		for (j in c(1:n))
		{
			if (i != j)
			{
			path2[i,j] <- ADJ[j,] %*% ADJ[,i] * ADJ[j,i]
			}
		}
	}
	 	
	#ALAAMobj$P2EdgeListIn <- P2EdgeListIn
	#ALAAMobj$P2degreeIn <- P2degreeIn
	#ALAAMobj$P2RowInIn <- P2RowInIn	 
	P2RowInIn <- matrix(NA,n,1)
	P2degreeIn <- matrix(0,n,1)
	P2EdgeListIn <- which(path2>=1,arr.ind=TRUE)
	P2EdgeListIn <- P2EdgeListIn[order(P2EdgeListIn[,1]),]
	P2EdgeListIn <- cbind(P2EdgeListIn[order(P2EdgeListIn[,1]),],matrix(0,dim(P2EdgeListIn)[1],1) )
	for (actor in c(1:n))
		 {
		 	TheseRows <- which(P2EdgeListIn[,1]==actor)
		 	if (length(TheseRows) >0){
		 		if (n==3 && sum(ADJ)<6)
		 		{
		 			if ( dim(P2EdgeListIn)[1] < max(TheseRows))
		 			{
		 			print('muppettery')
		 			browser()
		 			}
		 		}
		 	if (length( P2EdgeListIn[TheseRows,c(1:2)])==2)
		 	{
		 	#	print('pucko')
		 	#	browser()
		 		P2EdgeListIn[TheseRows,3] <-  path2[P2EdgeListIn[TheseRows,1],P2EdgeListIn[TheseRows,2]]
		 	}
		 	if (length( P2EdgeListIn[TheseRows,c(1:2)])>2)
		 	{


		 	P2EdgeListIn[TheseRows,3] <- path2[P2EdgeListIn[TheseRows,c(1:2)]]
		 	}
		 	if (length(TheseRows)>0)
		 	{
		 		P2RowInIn[actor] <- min(TheseRows)
		 		P2degreeIn[actor] <- length(TheseRows)
		 	}
		 	}
		 }

	}
	

	closedP <- list( CP2RowIn = P2RowIn, CP2degree = P2degree, CP2EdgeList = P2EdgeList , CP2RowInIn = P2RowInIn, CP2degreeIn = P2degreeIn, CP2EdgeListIn = P2EdgeListIn)
	closedP
}

getlistTrans <- function(ADJ,directed=FALSE)
{
n <- dim(ADJ)[1]
# create a k by 3 list translist of actors i, j, k such that i->j, i->k, k->j
# if y[i] is toggled, the change statistic is counts of:
# y[translist[indexF,2]]*y[translistF[index,3]]
# y[translist[indexM,1]]*y[translist[indexM,3]]
# y[translist[indexE,1]]*y[translist[indexE,2]]
# for indices indexF (index first in triple), indexM (index middle in triple), indexE (index end of triple)

translist <- matrix(c(9999,0,0),1,3)
First <- vector('list',n)
Middle <- vector('list',n)
Endnode <- vector('list',n)
rowIndex <-1
if (directed==FALSE)
{
for (i in c(1:(n-1) ))
	{
		
		for (j in c((i+1):(n-1)))
		{
			if (ADJ[i,j] )
			{
				for (k in c((j+1):n))
				{
					if (ADJ[i,k] & ADJ[k,j])
					{
						thisTrans <- matrix(c(i,k,j),1,3)
						translist <- rbind(translist,thisTrans )
						First[[i]] <- c(First[[i]] ,rowIndex)
						Middle[[k]] <- c(Middle[[k]],rowIndex)
						Endnode[[j]] <- c(Endnode[[j]] ,rowIndex)

						rowIndex <-rowIndex + 1
					}
					

				}
		 
		}
		}
	
	}
}

if (directed==TRUE)
{
for (i in c(1:n))
	{
		
		for (j in c(1:n))
		{
			if (ADJ[i,j] )
			{
				for (k in c(1:n))
				{
					if (ADJ[i,k] & ADJ[k,j])
					{
						thisTrans <- matrix(c(i,k,j),1,3)
						translist <- rbind(translist,thisTrans )
						First[[i]] <-c(First[[i]] ,rowIndex)
						Middle[[k]] <- c(Middle[[k]],rowIndex)
						Endnode[[j]] <-  c(Endnode[[j]] ,rowIndex)
						rowIndex <-rowIndex + 1
					}
					# if (ADJ[i,k] & ADJ[j,k])
					# {
						# thisTrans <- matrix(c(i,k,j),1,3)
						# translist <- rbind(translist,thisTrans )
						# First[[i]] <-c(First[[i]] ,rowIndex)
						# Middle[[k]] <-c(Middle[[k]],rowIndex)
						# Endnode[[j]] <- c(Endnode[[j]] ,rowIndex)
						# rowIndex <-rowIndex + 1
					# }
				}
		 
		}
		}
	
	}
}


# as we need both F, M, and E, we might need a list for each actor
translist <- translist[translist[,1]<9999,]
transLists <- list(translist= translist ,First = First, Middle = Middle, Endnode = Endnode )

transLists 
}

checkClosedP <- function(y,ADJ)
{
	n <- dim(ADJ)[1]
	countCP2 <- 0
	for (i in c(1:n))
	{
		if (y[i]==1){
		for (j in c(1:n))
		{
			if (ADJ[i,j] && y[j]==1 )
			{
		countCP2 <- countCP2 + ADJ[i,] %*% ADJ[,j] 
		}
		}
		}
	}
countCP2	
}


checkTrans <- function(y,ADJ,directed=FALSE)
{
	n <- dim(ADJ)[1]
	countT <- 0
	if (directed==TRUE)
	{
	for (i in c(1:n))
	{
		if (y[i]==1){
		for (j in c(1:n))
		{
			if (ADJ[i,j] && y[j]==1 )
			{
				for (k in c(1:n))
				{
					if (ADJ[i,k] & ADJ[k,j] & y[k]==1)
					{
					if (y[i]*y[j]*y[k]!=1)
					{
						print('maputo')
					}
						countT <- countT+1
					}
				}
		 
		}
		}
	}
	}
	}
	if (directed==FALSE)
	{
	for (i in c(1:n))
	{
		if (y[i]==1){
		for (j in c(1:n))
		{
			if (ADJ[i,j] && y[j]==1 )
			{
				for (k in c(1:n))
				{
					if (ADJ[i,k] & ADJ[k,j] & y[k]==1)
					{
						countT <- countT+1
					}
				}
		 
		}
		}
		}
	}
	countT <- countT/3
	}
countT	
}

checktrans2 <- function(y,ADJ)
{
	Tc <- 0
	n <- dim(ADJ)[1]
	for (i in c(1:(n-2)))
	{
		for (j in c((i+1):(n-1)))
		{
			for (k in c((j+1):n))
			{
				if (y[i]==1 && y[j]==1 && y[k]==1)
				{
					# i -> j ,i -> k -> j
					Tc <- Tc + ADJ[i,j]*ADJ[i,k]*ADJ[k,j]
					# j -> i ,j -> k -> i
					Tc <- Tc + ADJ[j,i]*ADJ[j,k]*ADJ[k,i]
					# k -> j ,k -> i -> j
					Tc <- Tc + ADJ[k,j]*ADJ[k,i]*ADJ[i,j]
					# i -> j, j -> k, i -> k
					Tc <- Tc + ADJ[i,j]*ADJ[j,k]*ADJ[i,k]
					# j -> i, i -> k, j -> k
					Tc <- Tc + ADJ[j,i]*ADJ[i,k]*ADJ[j,k]
					# k -> j, j - > i, k->i 
					Tc <- Tc + ADJ[k,j]*ADJ[j,i]*ADJ[k,i]

					
				}
				
			}
		}
	}
	Tc
}

unittestT <- function(n)
{
	y <- matrix(1,n,1)
	ADJ <- matrix(1,n,n)
	diag(ADJ) <- 0
	edges <- which(ADJ>0 ,arr.ind=TRUE)
	doThisorder <- sample(c(1:dim(edges)[1]),size=dim(edges)[1])
	edges <- edges[order(doThisorder),]
	compT <- matrix(0,dim(edges)[1],3)
	covs <- matrix(0,n,1)

	for (k in c(1:(dim(edges)[1]-1) ) )
	{
		ADJ[edges[k,1],edges[k,2]] <- 0
	mod.1 <- prepALAAMdata(y = y,ADJ = ADJ,covariates = covs,directed=TRUE, contagion =c('simple','recip','indirect','closedind','transitive'))
	
	if (k>=3)
	{
		print('muppetshow')
		browser()
		
	}
	stats2 <-  getALAAMstats(mod.1,contagion =c('simple','recip','indirect','closedind','transitive') )
	#stats3 <- checkClosedP(y = y, ADJ = ADJ)
	stats4 <- checkTrans(y = y, ADJ = ADJ, directed=TRUE)
	stats5 <- checktrans2(y = y, ADJ = ADJ)
	compT[k,1] <- stats2[6]
	compT[k,2] <- stats4
	compT[k,3] <- stats5
	}
	compT
}


prep.gof.cov <- function(ADJ,directed=TRUE){
if (directed)
	{# indegree:
# 	indegree <- colSums(ADJ)
# 	# outdegree
# 	outdegree <- rowSums(ADJ)
# 	# 2-paths
# 	twopaths <- colSums(ADJ)*rowSums(ADJ)
# 	# two-outstars
		n <- dim(ADJ)[1]
		
		### - Manual
# 	twoout <- matrix(0,n,1)
# 	twoin <- matrix(0,n,1)
# 	triangle <- matrix(0,n,3)
# 	indir <- ADJ %*% ADJ
# diag(indir) <- 0
# indirectTies <- rowSums(indir)
# 
# 
# 	for (k in c(1:n))
# 	{
# 		twoout[k] <- choose(outdegree[k],2)
# 		twoin[k] <- choose(indegree[k],2)
# 		triangle[k,1] <- sum(ADJ[ADJ[k,]==1,ADJ[k,]==1])
# 		triangle[k,2] <- sum(ADJ[ADJ[,k]==1,ADJ[,k]==1])
# 		triangle[k,3] <- 0
# 		triangle[k,3] <- 0
# 		for (j in c(1:n) )
# 		{
# 			if (ADJ[j,k]==1){
# 		#	triangle[k,3] <- triangle[k,3]+ sum(ADJ[j,]*ADJ[k,]) -ADJ[k,j]
# 			triangle[k,3] <- triangle[k,3]+ sum(ADJ[j,]*ADJ[k,]) -ADJ[k,j]
# 			}
# 		}
# 
# 	}
# 		
		#### - end manual
  out.degree <-matrix( rowSums(ADJ), n, 1) # number of ties sent
  in.degree <- matrix( colSums(ADJ) , n, 1 ) # number of ties received
  rec.ties <-  matrix( rowSums(ADJ * t(ADJ) ), n , 1) # number of ties that are mutual
  in.two.star <- matrix( choose(in.degree,2),n,1) #  in-stars refecting dispersion in popularity
  out.two.star <- matrix( choose(out.degree,2),n,1) #  out-stars refecting dispersion in activity
  mix.two.star <- in.degree*out.degree - rec.ties # correlation between indegree and outdegree
  in.three.star <- matrix( choose(in.degree,3),n,1) # furhter measure of in-degree heterogeneity
  out.three.star <- matrix( choose(out.degree,3),n,1) # furhter measure of out-degree heterogeneity
  triangles <- rowSums( ADJ* (ADJ %*% t(ADJ) )  ) # e
  cyclic.tri <- rowSums( t(ADJ)* (ADJ %*% ADJ )  )
  
  num.indirect.nonexcl <-  ADJ %*% out.degree - rec.ties
  
  num.indirect <- matrix(0,n,1)
  for (i in c(1:n))
  {
    if ( sum( ADJ[ i, ])==1 )
    {
      num.indirect[i] <- sum(  ADJ[ADJ[ i, ]==1, ADJ[ i, ]==0]  ) - (ADJ[i,] %*% ADJ[,i]>0)
    }
    
    if ( sum(ADJ[ i, ])>1 )
    {
      num.indirect[i] <- sum( colSums( ADJ[ADJ[ i, ]==1, ADJ[ i, ]==0] ) > 0  ) - (ADJ[i,] %*% ADJ[,i]>0)
    }
  }
  
  covariates <- cbind(out.degree, 
                in.degree,
                rec.ties,
                in.two.star,
                out.two.star,
                mix.two.star,
                in.three.star,
                out.three.star,
                triangles,
                cyclic.tri,
                num.indirect.nonexcl,
                num.indirect )
  colnames(covariates) <- c("outdegree",
                    "indegree",
                      "reciprochation" ,
                      "instar",
                      "outstar",
                      "twopath",
                      "in3star",
                      "out3star",
                      "transitive",
                    "cyclic",
                    "indirect",
                    "excl.indirect")
  
	#covariates <- cbind(indegree,outdegree,twopaths,twoout,twoin,triangle,indirectTies)
	#colnames(covariates) <- c('indegree','outdegree','twopaths','out2star','in2star','outtria','intria','transtri','indirect ties')
	}
	if (directed==FALSE)
	{
		degree <- colSums(ADJ)
		n <- length(y)
		outstar <- matrix(0,n,1)
		triangle <- matrix(0,n,1)
	for (k in c(1:n))
	{
		outstar[k] <- choose(degree[k],2)
		triangle[k,1] <- 0
		for (j in c(1:n) )
		{
			if (ADJ[j,k]==1){
			triangle[k,1] <- triangle[k,1]+ sum(ADJ[j,]*ADJ[k,]) -1
			}
		}


		}
		covariates <- cbind(degree,outstar,triangle)
		colnames(covariates) <- c('degree','star','triangle')
	}
	
	covariates
}

alaam.gof.statscalc <- function(y,ADJ,covariates,directed=TRUE)
{
		mod.1 <- prepALAAMdata(y = y,ADJ = ADJ,covariates = covariates,directed=directed, contagion =c('simple','recip','indirect','closedind','transitive'))
		stats <- getALAAMstats(ALAAMobj=mod.1,contagion =c('simple','recip','indirect','closedind','transitive'))
		stats
}

alaam.gof <- function(ALAAMobj,res.obj,thinning=10,burning=100,contagion='simple')
{
	# figure out where to put parameter values
	interaction <- ALAAMobj$interaction
	pcov <- dim(ALAAMobj$covariates)[2]
	pcont <- length(contagion)
	pinter <- length(interaction)
	pOther <- 5 # reciprocated ties, indegree, outegree, mixedpath, indirect ties
	non.zeros <- matrix(NA,1+5+pcov+pinter+pOther,1)
	non.zeros[1] <- 1 # intercept
	if ('simple' %in% contagion)
	{
		non.zeros[2] <- 1 # simplecont
	}
	if ('recip' %in% contagion)
	{
		non.zeros[3] <- 1 # simplecont
	}
	
	if ('indirect' %in% contagion)
	{
		non.zeros[4] <- 1 # recip
	}
	
	if ('closedind' %in% contagion)
	{
		non.zeros[5] <- 1 # closed
	}
	if ('transitive' %in% contagion)
	{
		non.zeros[6] <- 1 # transitive
	}
	k <- (7+pcont)
	non.zeros[7:k] <- 1
	k <- k+1
	if (pinter >0 )
	{
		non.zeros[k:(k+pinter)] <- 1
		k <- (k+pinter)+1
	}
	non.zeros[k:(k+pOther )] <- 0
	ADJ <- edglist.to.matrix(ALAAMobj = ALAAMobj,n = length(ALAAMobj$y))
	indegree <- colSums( ADJ )
	outdegree <- rowSums( ADJ )
	reciprocal <- ADJ*t(ADJ)
reciprocal[reciprocal>0] <- 1
diag(reciprocal) <- 0
recties <- rowSums(reciprocal)
	
	browser()
	sim.1 <- simulate.alaam(ALAAMobj=ALAAMobj,statsvec=stats,theta=theta1,contagion =c('simple','recip','indirect','closedind','transitive') , thinning = 100, NumIterations = 30, burnin = 3000, DoSave=TRUE, returnNet= TRUE)
	
}


edglist.to.matrix <- function(ALAAMobj,n)
{
	ADJ <- matrix(0,n,n)
	ADJ[ALAAMobj$EdgeList] <- 1
	
	ADJ
	
}

##### for pseufolikelihood
changestatsALAAM <- function(y,EdgeList,RowIn,degree,covariates,NumIterations,theta,statsvec,directed=FALSE,EdgeListIn=NULL, degreein=NULL,RowInIn=NULL,DoSave=FALSE,thinning = NULL,burnin=NULL,canchange=NULL,returnNet=FALSE,interaction=NULL,checkAlgorithm=FALSE, bias = NULL, contagion ='simple',RecEdgeList=NULL,Recdegree=NULL,RecRowIn=NULL,P2EdgeList=NULL,P2degree=NULL,P2RowIn=NULL, P2EdgeListIn=NULL,P2degreeIn=NULL,P2RowInIn	=NULL, CP2EdgeList=NULL, CP2degree = NULL, CP2RowIn = NULL , CP2EdgeListIn = NULL , CP2degreeIn = NULL, CP2RowInIn = NULL, translist = NULL , transFirst = NULL , transMiddle = NULL , transEndnode = NULL)
{
	####
		
	
	#####
	n <- length(y)
	p <- dim(covariates)[2]+1# add 1 for intercept
	p1 <- p
	if ('simple' %in% contagion){
		p <- p+1
		}
		
	if (('simple' %in% contagion)==FALSE){
		error('model without simple contagion not yet specified')
		}
		dorec <- FALSE
	if ('recip' %in% contagion){
		p <- p+1
		dorec <- TRUE
		if (is.null(RecEdgeList) | is.null(Recdegree) | is.null(RecRowIn))
		{
		error('edgelist, reciprocated degree, and index to reciprocated edgelist required if reciprocated contagion requested')	
		}
		
		}
		doindir <- FALSE
	if ('indirect' %in% contagion)
		{
			p <- p+1
			doindir <- TRUE
		
			if (is.null(P2EdgeList) | is.null(P2degree) | is.null(P2RowIn))
		{
		error('edgelist, 2-path degree, and index to 2-path edgelist required if indirect contagion requested')	
		}
		 }
	doCloseindir <- FALSE
	if ('closedind' %in% contagion)
		{
			p <- p+1
			doCloseindir <- TRUE
			}
	doTrans <- FALSE
	if ('transitive' %in% contagion)
		 	{
		 		p <- p+1
		 		doTrans <- TRUE
		 				 		}

	p2 <- 0
	if (!is.null(interaction))
	{
		p2 <- length(interaction)
		p <- p+p2
	}
	ystar <- y

	statsvecstar <- matrix(0,p,1)
	
	if (is.null(canchange))
	{
		# regular sampling
		canchange <- c(1:n)
		
	}
	
	
			BigTotal <- n
	
	SimStats <- matrix(0,p,n)
	
	NumIterations <- n
	saveHere <- 1
	

	
	
	for (actor in c(1:n))
	{
		
		statsvecstar[1] <- 1# - 2*ystar[actor]
		neig <- 0 # for contagion count
		statsvecstar[2] <- 0
		Recneig <- 0 # for reciprocated contagion count
		Indineig <- 0 # for indirect contagion count
		CIndineig <- 0 # for closed indirect contagion count
		TIndineig <- 0 # for transitive contagion count
		if (degree[actor] > 0)
		{
		for (alter in c(1:degree[actor]))
		{
			neig <- neig + ystar[ EdgeList[RowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		#neig <- neig + sum(ystar[ EdgeList[RowIn[actor]+c(1:degree[actor])-1,2] ])
		}
		
		if (directed)
		{
			if (degreein[actor]>0)
			{
			for (alter in c(1:degreein[actor]))
		{
			neig <- neig + ystar[ EdgeListIn[RowInIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row RowIn[actor]
		}
		#neig <- neig + sum(ystar[ EdgeListIn[RowInIn[actor]+c(1:degreein[actor])-1,2] ])
		}
		}
		
		### reciprocated contagion if requested
		if (dorec)
		{
				#ALAAMobj$RecEdgeList <- RecEdgeList
		 		#ALAAMobj$Recdegree <- Recdegree
		 		#ALAAMobj$RecRowIn <- RecRowIn

		if (Recdegree[actor]>0)
		{
			for (alter in c(1:Recdegree[actor]))
		{
			Recneig <- Recneig + ystar[ RecEdgeList[RecRowIn[actor]+alter-1,2] ]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
			
		}

		### indirect contagion if requested
		if (doindir)
		{
		 		#ALAAMobj$P2EdgeList <- P2EdgeList
		 		#ALAAMobj$P2degree <- P2degree
		 		#ALAAMobj$P2RowIn <- P2RowIn

		if (P2degree[actor]>0)
		{
			for (alter in c(1:P2degree[actor]))
		{
			Indineig <- Indineig + ystar[ P2EdgeList[P2RowIn[actor]+alter-1,2] ]*P2EdgeList[P2RowIn[actor]+alter-1,3]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		
		if (directed)
		{
			
			if (P2degreeIn[actor]>0)
		{
			for (alter in c(1:P2degreeIn[actor]))
		{
			Indineig <- Indineig + ystar[ P2EdgeListIn[P2RowInIn[actor]+alter-1,2] ]*P2EdgeListIn[P2RowInIn[actor]+alter-1,3]#assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
			
			}
			
		}
		
		### indirect + direct if requested
		if (doCloseindir)
		{
			if (CP2degree[actor]>0)
		{
			for (alter in c(1:CP2degree[actor]))
		{
			CIndineig <- CIndineig + ystar[ CP2EdgeList[CP2RowIn[actor]+alter-1,2] ]*CP2EdgeList[CP2RowIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		
		if (directed)
		{
			if (CP2degreeIn[actor]>0)
		{
			for (alter in c(1:CP2degreeIn[actor]))
		{
			CIndineig <- CIndineig + ystar[ CP2EdgeListIn[CP2RowInIn[actor]+alter-1,2] ]*CP2EdgeListIn[CP2RowInIn[actor]+alter-1,3] #assumes edgelist
			# where the ties of actor starts at row Row[actor]
		}
		}
		}

		}
		
		### transitive contagion
		if (doTrans)
		{
			
			copyTP <- transFirst[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],2]]==1 && ystar[translist[copyTP[rowind],3]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			copyTP <- transMiddle[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],3]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			copyTP <- transEndnode[[actor]]
			numTps <- length(copyTP)
			if (numTps>0)
			{
				for (rowind in c(1:numTps))
				{
					if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],2]]==1)
					{
						TIndineig <- TIndineig + 1
						}
				}
			}
			

			
		}
		
		
		statsvecstar[2] <- statsvecstar[1]*neig
			k <- 3
		if (dorec)
		{
			statsvecstar[k] <- statsvecstar[1]*Recneig
			k <- k+1
		}
		if (doindir)
		{
			statsvecstar[k] <- statsvecstar[1]*Indineig
			k <- k+1
		}
		
		if (doCloseindir)
		{
			statsvecstar[k] <- statsvecstar[1]*CIndineig
			k <- k+1
			
		}
		if (doTrans)
		{
			statsvecstar[k] <- statsvecstar[1]*TIndineig
			k <- k+1
		}
		
		for (t in c(1:(p1 - 1) ))
		{
			statsvecstar[k] <- statsvecstar[1]*covariates[actor,t]
			k <- k+1
		}
		
		if (!is.null(interaction))
	{
		#statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction]
		statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction[1]]
	}
	
	
		
		
		
		
		
		
	
	#print(paste('you have done ',iterations,' iterations'))
	SimStats[,actor] <- statsvecstar
	
	

	
	}
	
	
		
		statsvec <- SimStats
			
	
		
		
	statsvec 

}

returnCI <- function(paradraws)
{
xsort <- sort(paradraws)
N <- length(xsort)
lowervalue <- floor(N*0.025)
uppervalue <- ceiling(N*(1-0.025))
CI <- c(xsort[lowervalue],xsort[uppervalue])

CI
}

write.res.table <- function(ALAAMresult=NULL,burnin=1,thin=1,tabname=NULL,nameVec=NULL)
{
	# this function writes the posterior means and credibility intervals
	# to a csv file
	datamat <- ALAAMresult$Thetas
	N <- dim(datamat)[1]
	p <- dim(datamat)[2]
	# name, mean, sd, ci x 2
	if (is.null(nameVec))
	{
		nameVec <- colnames(datamat)
		
	}
CredibilityIntervals <- as.data.frame(matrix(0,p,5))
names(CredibilityIntervals) <- c('parameter','mean','sd','.025','0.975')
for (t in c(1:p))
{
 CredibilityIntervals[t,1] <- nameVec[t]
  CredibilityIntervals[t,2] <- round(mean(datamat[seq(burnin,N,by=thin),t]),3)
  CredibilityIntervals[t,3] <- round(sd(datamat[seq(burnin,N,by=thin),t]),3)
 CredibilityIntervals[t,4:5] <- returnCI(paradraws=datamat[seq(burnin,N,by=thin),t])
  CredibilityIntervals[t,4] <- round(CredibilityIntervals[t,4],3)
  CredibilityIntervals[t,5] <- round(CredibilityIntervals[t,5],3)
}
print(CredibilityIntervals)
write.csv(CredibilityIntervals,file=paste(tabname,'.csv',sep=''))

}

truncNormalIndep <- function(mu,variance,zlower,zupper)
{
	# zvector should be 1 times length(zvector)
	# the purpose is to draw a conditional normal variate
	# zcoord (the argument is not used)
	# conditional on
	# zvector[-coordindex]
	# truncated to the left in zlower and to the right in zupper
		sigma <- sqrt(variance)
	
	zcoord <- mu +sigma * qnorm( runif(1,pnorm((zlower - mu)/sigma),pnorm((zupper - mu)/sigma)) ,0,1)
	#	if (is.na(zcoord))
	#{print(zvector)
	#	print(zlower)
	#	print(zupper - zlower)
	#	print(pnorm((zlower - mu)/sigma)) }
	return(zcoord)

}

truncNormalpdf <- function(x,mu,variance,zlower,zupper)
{
	fx <- dnorm(x,mean = mu, sd = sqrt(variance))
	
	fx <- fx/(pnorm(zupper, mean = mu, sd = sqrt(variance) ) - pnorm(zlower, mean = mu, sd = sqrt(variance) ))
	
	fx
		
}

GenerateMultivariateNormLike <- function(model='neteffs',y,X,beta,u,rho,sigma,ADJ,alternative)
{
	if (model == 'netauto')
	{
		n <- length(y)
		
	
	A <- solve( diag( rep(1,n) ) - rho * ADJ )

	mu <- A %*% (X %*% beta + u)

	Sigma <- A %*% diag( rep(sigma,n) ) %*% t(A)
	y <- mu + sqrtm(Sigma)%*%rnorm(n,0,1)
	if (alternative==TRUE)
	{
		y <- A %*%(X %*% beta + u + rnorm(n,0,1))
	}
	}
	if (model == 'neteffs')
	{
		n <- length(y)
		
	A <- solve( diag( rep(1,n) ) - rho * ADJ )
	mu <- X %*% beta +  u
	
	Sigma <- A %*% diag( rep(sigma,n) ) %*% t(A)
	y <- mu + sqrtm(Sigma)%*%rnorm(n,0,1)
	}
	y
}


updateMu <- function(type = 'multi',rho,rhomuHyper,tau,c0,G,n0,rhoMu)
{
	
	#rhoMu <- ( c0*sum(rho ) + rhomuHyper )/(c0 *(G + 1) ) +rnorm(1,0,1) * sqrt((c0 * tau )/ (c0 *(G + 1) ) )
	#rhoMu
	rhoMuStar <- truncNormalIndep(mu=rhoMu,variance=tau/3,zlower=-1,zupper=1)
	LikeOld <- 0
	LikeNew <- 0
	PriorOld <- log(truncNormalpdf(x=rhoMu,mu=rhomuHyper,variance=tau/n0,zlower=-1,zupper=1))
	PriorNew <- log(truncNormalpdf(x=rhoMuStar,mu=rhomuHyper,variance=tau/n0,zlower=-1,zupper=1))
	 PropRatio <-log(truncNormalpdf(x=rhoMu,mu=rhoMuStar,variance=tau/3,zlower=-1,zupper=1))
	 PropRatio <- PropRatio-log(truncNormalpdf(x=rhoMuStar,mu=rhoMu,variance=tau/3,zlower=-1,zupper=1))
	if (type == 'multi')
	{
		for (g in c(1:G))
	{
	LikeOld <-LikeOld+log(truncNormalpdf(x=rho[g],mu=rhoMu,variance=tau,zlower=-1,zupper=1))
	LikeNew <- LikeNew+log(truncNormalpdf(x=rho[g],mu=rhoMuStar,variance=tau,zlower=-1,zupper=1))
	
	}
	 Hasting <- LikeNew-LikeOld+PriorNew-PriorOld+PropRatio
        if (log(runif(1,0,1)) < Hasting )
         {
         	rhoMu <- rhoMuStar
         }
         
        # browser()
	}
	if (type == 'regular')
	{
		
	LikeOld <-LikeOld+log(truncNormalpdf(x=rho,mu=rhoMu,variance=tau,zlower=-1,zupper=1))
	LikeNew <- LikeNew+log(truncNormalpdf(x=rho,mu=rhoMuStar,variance=tau,zlower=-1,zupper=1))
	
	
	 Hasting <- LikeNew-LikeOld+PriorNew-PriorOld+PropRatio
        if (log(runif(1,0,1)) < Hasting )
         {
         	rhoMu <- rhoMuStar
         }
         
        # browser()
	}
	
	
	
	rhoMu
}

boundrho <- function(ADJ,steplength)
{
	
		
		lambda <- Re(eigen(ADJ)$values)
    AllSmallestRho <- 0
    AllLargestRho <- 0
    RhoLamProd <- 1
    while (AllLargestRho < 1 && RhoLamProd>0)
       { # increment LargestRho
       	AllLargestRho <- AllLargestRho + steplength/10
        RhoLamProd <- prod( 1 - AllLargestRho * lambda)
}   

    RhoLamProd <- 1
    while ( AllSmallestRho  > -1 && RhoLamProd>0)
        {
       # increment SmallestRho
        AllSmallestRho   <- AllSmallestRho  - steplength/10
        RhoLamProd = prod( 1 - AllSmallestRho   * lambda )
        }
        
    TunigStuff <-  list(lambda = lambda,RhoMin=AllSmallestRho,RhoMax=AllLargestRho)
    TunigStuff 
}

updaterhoBrute <- function(model='neteffs',type='multi',y,X,beta,u,rho,sigma,ADJ,steplength,RhoMin,RhoMax,rhoMu,tau=NULL)
{
	G <- length(X)
	if (type == 'multi')
	{
		for (g in c(1:G))
		{
			
		lowerrhoStar <- max(c(rho[g] - steplength, RhoMin[g] ) )
        upperrhoStar <- min( c(rho[g] + steplength, RhoMax[g] ) )
        rhoStar <- runif(1,lowerrhoStar,upperrhoStar)
       # rhoStar <- runif(1,0,1)*(upperrhoStar - lowerrhoStar) + lowerrhoStar
        lowerrho <- max(c(rhoStar - steplength, RhoMin[g] ))
        upperrho <- min(c(rhoStar + steplength, RhoMax[g]) )
       # priorRate <- -(1/(2* tau))*( (rhoStar - rhoMu )^2 - ( rho[g]  - rhoMu )^2 )
   priorRate <- log(truncNormalpdf(x=rhoStar,mu=rhoMu,variance=tau,zlower= RhoMin[g] ,zupper=RhoMax[g] ))
   
   priorRate <- priorRate -log(truncNormalpdf(x=rho[g] ,mu=rhoMu,variance=tau,zlower= RhoMin[g] ,zupper=RhoMax[g] ))
        
        LikeOld <- MultivariateNormLike(model=model,y = y[[g]],X = X[[g]],beta,u=u[g],rho = rho[g],sigma,ADJ =ADJ[[g]])
   
        LikeNew <- MultivariateNormLike(model=model,y = y[[g]],X = X[[g]],beta,u=u[g],rho = rhoStar,sigma,ADJ =ADJ[[g]])
         #AccRatio <- sum( log( 1 - rhoStar * lambda[[g]])) - sum( log( 1 - rho[g]  * lambda[[g]]))
         #PropRatio <- log(upperrhoStar - lowerrhoStar) -  log(upperrho - lowerrho)
        PropRatio <- dunif(rho[g],lowerrho,upperrho,log=TRUE)-dunif(rhoStar,lowerrhoStar,upperrhoStar,log=TRUE)
        Hasting <- LikeNew-LikeOld+priorRate+PropRatio
        if (log(runif(1,0,1)) < Hasting )
         {
         	rho[g] <- rhoStar
         }
		}
	}
	if (type == 'notmulti')
	{
		g <- 1# all groups have the same rho
		lowerrhoStar <- max(c(rho[g] - steplength, max(RhoMin ) ))
        upperrhoStar <- min( c(rho[g] + steplength, min(RhoMax) ) )
        rhoStar <- runif(1,lowerrhoStar,upperrhoStar)
        lowerrho <- max(c(rhoStar - steplength, max(RhoMin )  ))
        upperrho <- min(c(rhoStar + steplength, min(RhoMax) ) )
        #priorRate <- -(1/(2* tau))*( (rhoStar - rhoMu )^2 - ( rho[g]  - rhoMu )^2 )
        priorRate <- log(truncNormalpdf(x=rhoStar,mu=rhoMu,variance=tau,zlower= max(RhoMin ) ,zupper=min(RhoMax)  ))
   
   priorRate <- priorRate -log(truncNormalpdf(x=rho[g] ,mu=rhoMu,variance=tau,zlower= max(RhoMin ),zupper=min(RhoMax)  ))
        Hasting <-0
       
        LikeOld <- combineLikeslihoods(model=model,y,X,beta,u,rho ,sigma,ADJ)
        LikeNew <- combineLikeslihoods(model=model,y,X,beta,u,rho = matrix(rep(rhoStar,G),G,1),sigma,ADJ)
        #PropRatio <- log(upperrhoStar - lowerrhoStar) -  log(upperrho - lowerrho)
        
        PropRatio <- dunif(rho[g],lowerrho,upperrho,log=TRUE)-dunif(rhoStar,lowerrhoStar,upperrhoStar,log=TRUE)
      
         Hasting <- LikeNew-LikeOld+priorRate+PropRatio
         if (log(runif(1,0,1)) < Hasting)
         {
         	rho <- matrix(rep(rhoStar,G),G,1)
         }

	}

if (type == 'regular')
	{
		g <- 1# all groups have the same rho
		lowerrhoStar <- max(c(rho - steplength, max(RhoMin ) ))
        upperrhoStar <- min( c(rho + steplength, min(RhoMax) ) )
        rhoStar <- runif(1,lowerrhoStar,upperrhoStar)
        lowerrho <- max(c(rhoStar - steplength, max(RhoMin )  ))
        upperrho <- min(c(rhoStar + steplength, min(RhoMax) ) )
        #priorRate <- -(1/(2* tau))*( (rhoStar - rhoMu )^2 - ( rho[g]  - rhoMu )^2 )
        priorRate <- log(truncNormalpdf(x=rhoStar,mu=rhoMu,variance=tau,zlower= max(RhoMin ) ,zupper=min(RhoMax)  ))
   
   priorRate <- priorRate -log(truncNormalpdf(x=rho ,mu=rhoMu,variance=tau,zlower= max(RhoMin ),zupper=min(RhoMax)  ))
        Hasting <-0
       
        LikeOld <- combineLikeslihoods(model=model,y,X,beta,u,rho ,sigma,ADJ)
        LikeNew <- combineLikeslihoods(model=model,y,X,beta,u,rho = matrix(rep(rhoStar,G),G,1),sigma,ADJ)
        #PropRatio <- log(upperrhoStar - lowerrhoStar) -  log(upperrho - lowerrho)
        
        PropRatio <- dunif(rho,lowerrho,upperrho,log=TRUE)-dunif(rhoStar,lowerrhoStar,upperrhoStar,log=TRUE)
      
         Hasting <- LikeNew-LikeOld+priorRate+PropRatio
         if (log(runif(1,0,1)) < Hasting)
         {
         	rho <- matrix(rep(rhoStar,G),G,1)
         }

	}
	
rho

}




BayesNamNec <- function(NumIter,y,x,W1=NULL,W2=NULL,alpha,saveFreq=10)
{
	require('expm')
	# estimation of beta for model
	# y =rho*ADJ*y+x*beta+epsilon
	yobs <-y # we will let y be the latent variable
	RESULTS <- vector("list", 4)
	tau <- 1
	rhoMu <-0 
	beta <- solve(t(as.matrix(x))%*%as.matrix(x))%*% t(as.matrix(x))%*%as.matrix(y)

	binary = FALSE
	y.vals <- range(y)
	if (y.vals[1]==0 & y.vals[2]==1)
	{
		cat('dependent variable binary. probit model will be used \n')
		binary = TRUE
		for (k in c(1:length(y))){
		if (yobs[k]==0){
		zlower <- -Inf
		zupper <- 0	
		}
		if (yobs[k]==1){
		zlower <- 0
		zupper <- Inf
		}
		y[k] <- truncNormalIndep(mu=x[k,]%*%beta,variance=1,zlower=zlower,zupper=zupper)
		}
	}
	
	saveNumIter <- ceiling(NumIter/saveFreq)
	printFreq <- 1000
	p <- dim(x)[2]
	RESULTS$beta <- matrix(0,p,saveNumIter)
	RESULTS$rho <- matrix(0,1,saveNumIter)
	RESULTS$alpha <- matrix(0,1,saveNumIter)
	RESULTS$sigma <- matrix(0,1,saveNumIter)
			sigma <- (1/(dim(x)[1] - p) )*t(y-x%*%beta)%*%(y-x%*%beta)

	PropVar <- as.numeric(sigma)*solve(t(as.matrix(x))%*%as.matrix(x))/(1+p)
	PropVarSqrt <- sqrtm(PropVar)
	PropVarInv <- solve(PropVar)
	PhiPriorInv <- diag(rep(1,p)) * 0.01
	muprior <- matrix(0,p,1)
	steplength <- .1
	
	if (!is.null(W1)){
	TunigStuffAlpha <- boundrho(W1,steplength =.1)
	#alpha <- 0
	}
	
	if (is.null(W1))
	{
		alpha <- NULL
	}
	
	if (!is.null(W2)){
	TunigStuffRho <- boundrho(W2,steplength =.1)
	rho<- 0 
	}
	
	if (is.null(W2)){
		rho=NULL
}

saveIndex <- 1
	
	for (inter in c(1:NumIter))
	{
		
		beta <- update.beta(y,x,beta,rho=rho,alpha=alpha,W1=W1,W2=W2,PsiInv=PhiPriorInv, muprior = muprior)
	
		if (is.null(alpha)==FALSE){
			
		alpha <- update.alpha(y,x,beta=beta,rho=rho,alpha=alpha,sigma,W1=W1,W2=W2,rho.min=TunigStuffAlpha$RhoMin,rho.max=TunigStuffAlpha$RhoMax,steplength,rhoMu,tau)
		}
		
		if (binary==TRUE)
		{
			y <- update.trunc.norm(y=y,yobs=yobs,x, beta=beta, rho=rho,alpha=alpha,W1=W1,W2=W2)
		}
		
	if (round(inter/saveFreq)==(inter/saveFreq))
				{
						RESULTS$beta[,saveIndex] <- beta
						RESULTS$alpha[,saveIndex] <- alpha
		#RESULTS$rho[,saveIndex] <- rho
	
	
	saveIndex <- saveIndex + 1

	}

if (round(inter/printFreq)==(inter/printFreq))
				{
						print(paste('iterations: ',inter,' out of ',NumIter))
						print('regression coeffs')
						print(beta)
						cat('correlation coeffs rho ',rho,' and alpha: ',alpha,'\n')
						
						
						
						}

}

RESULTS
	
}

net.reg.like <- function(y,x,beta,rho=NULL,alpha=NULL,sigma,W1=NULL,W2=NULL)
{
	
	if (is.null(W2))
	{
		if (is.null(W1)==FALSE)
		{
		N <- length(y)
	A <- solve( diag( N ) - alpha * W1 )
	mu <- A %*% (x %*% beta )

	Sigma <- A %*%  t(A)
	logLike <- -.5* log(det(Sigma))-0.5*t(y-mu) %*% solve(Sigma) %*% (y-mu)-(N/2)*log(2*pi)
	}
	if (is.null(W1))
		{
			# regular regression
			}
	
	}
	if (is.null(W2)==FALSE)
	{
		if (is.null(W1)==FALSE)
		{
			# both
			}
		if (is.null(W1))
		{
		N <- length(y)
	A <- solve( diag( N ) - rho * W1 )
	mu <- x %*% beta
	
	Sigma <- A %*% t(A)
	logLike <- -.5* log(det(Sigma))-0.5*t(y-mu) %*% solve(Sigma) %*% (y-mu)-(N/2)*log(2*pi)
	}
	}
	logLike
}

update.rho <- function(y,x,beta,rho=NULL,alpha=NULL,sigma,W1=NULL,W2=NULL,rho.min,rho.max,steplength)
{
	lowerrhoStar <- max(c(rho - steplength, rho.min ))
        upperrhoStar <- min( c(rho + steplength, rho.max ) )
        rhoStar <- runif(1,lowerrhoStar,upperrhoStar)
        lowerrho <- max(c(rhoStar - steplength, rho.min  ))
        upperrho <- min(c(rhoStar + steplength, rho.max ) )
        #priorRate <- -(1/(2* tau))*( (rhoStar - rhoMu )^2 - ( rho[g]  - rhoMu )^2 )
        priorRate <- log(truncNormalpdf(x=rhoStar,mu=rhoMu,variance=tau,zlower= rho.min ,zupper=rho.max ))
   
   priorRate <- priorRate -log(truncNormalpdf(x=rho ,mu=rhoMu,variance=tau,zlower= rho.min ,zupper=rho.max  ))
        Hasting <-0
       
         LikeOld <- net.reg.like(y=y,x=x,beta=beta,rho=rho,alpha=alpha,sigma,W1=W1,W2=W2)
           LikeNew <- net.reg.like(y=y,x=x,beta=beta,rho= rhoStar,alpha=alpha,sigma,W1=W1,W2=W2)
    
       
        #PropRatio <- log(upperrhoStar - lowerrhoStar) -  log(upperrho - lowerrho)
        
        PropRatio <- dunif(rho,lowerrho,upperrho,log=TRUE)-dunif(rhoStar,lowerrhoStar,upperrhoStar,log=TRUE)
      
         Hasting <- LikeNew-LikeOld+priorRate+PropRatio
         if (log(runif(1,0,1)) < Hasting)
         {
         	rho <- rhoStar
         }
	rho
}

update.alpha <- function(y,x,beta,rho=NULL,alpha=NULL,sigma,W1=NULL,W2=NULL,rho.min,rho.max,steplength,rhoMu,tau)
{
	loweralphaStar <- max(c(alpha - steplength, rho.min ))
        upperalphaStar <- min( c(alpha + steplength, rho.max ) )
        alphaStar <- runif(1,loweralphaStar,upperalphaStar)
        lowerrho <- max(c( alphaStar - steplength, rho.min  ))
        upperrho <- min(c( alphaStar + steplength, rho.max ) )
        #priorRate <- -(1/(2* tau))*( (rhoStar - rhoMu )^2 - ( rho[g]  - rhoMu )^2 )
        
    
        priorRate <- log(truncNormalpdf(x=alphaStar,mu=rhoMu,variance=tau,zlower= rho.min ,zupper=rho.max ))
   
   priorRate <- priorRate -log(truncNormalpdf(x=alpha ,mu=rhoMu,variance=tau,zlower= rho.min ,zupper=rho.max   ))
        Hasting <-0
       
         LikeOld <- net.reg.like(y=y,x=x,beta=beta,rho=rho,alpha=alpha,sigma,W1=W1,W2=W2)
           LikeNew <- net.reg.like(y=y,x=x,beta=beta,rho= rho,alpha= alphaStar,sigma,W1=W1,W2=W2)
    
       
        #PropRatio <- log(upperrhoStar - lowerrhoStar) -  log(upperrho - lowerrho)
        
        PropRatio <- dunif(alpha,lowerrho,upperrho,log=TRUE)-dunif(alphaStar,loweralphaStar,upperalphaStar,log=TRUE)
      
         Hasting <- LikeNew-LikeOld+priorRate+PropRatio
         if (log(runif(1,0,1)) < Hasting)
         {
         	alpha <-  alphaStar
         }
	alpha
}

update.beta <- function(y,x,beta,rho=NULL,alpha=NULL,W1=NULL,W2=NULL,PsiInv, muprior )
{
	p <- dim(x)[2]
	if (is.null(W2))
	{
		if (is.null(W1)==FALSE)
		{
			
			# network autocorrelation model
		n <- length(y)
	A <- solve( diag( n) - alpha * W1 )

	x <- A %*% x # to get regression of form y = X*beta+epsilon

	Sigma <- A %*%  t(A)
	C <- chol(Sigma ) 
	Xstar <- t(C) %*% x
	ystar <- t(C) %*% y
	PsiStar <- solve( t(Xstar) %*%  Xstar + PsiInv )
	mu <- PsiStar %*% ( t(Xstar) %*%  ystar + PsiInv %*% muprior )
	
	#SigmaStar <- A %*% t(A)
	
	
	
	
	}
	if (is.null(W1))
		{
			# regular regression
			}
	
	}
	if (is.null(W2)==FALSE)
	{
		if (is.null(W1)==FALSE)
		{
			# both
			}
		if (is.null(W1))
		{
		n <- length(y)
	A <- solve( diag( rep(1,n) ) - rho * W1 )
	Sigma <- A %*%  t(A)
	C <- chol(Sigma ) 
	Xstar <- t(C) %*% x
	ystar <- t(C) %*% y
	PsiStar <- solve( t(Xstar) %*%  Xstar + PsiInv )
	mu <- PsiStar %*% ( t(Xstar) %*%  ystar + PsiInv %*% muprior )
	

	}
	}
	
	 beta <- matrix(rmvnorm(1,mu, sigma = PsiStar),p,1)
	 beta
}

update.trunc.norm <- function(y,yobs,x, beta, rho=NULL,alpha=NULL,W1=NULL,W2=NULL)
{
	
	p <- dim(x)[2]
	if (is.null(W2))
	{
		if (is.null(W1)==FALSE)
		{
			
			# network autocorrelation model
			N <- length(y)
	A <- solve( diag( N ) - alpha * W1 )
		mu <- A %*% (x %*% beta )

	Sigma <- A %*%  t(A)	
	#SigmaStar <- A %*% t(A)
	
	
	
	
	}
	if (is.null(W1))
		{
			# regular regression
			}
	
	}
	if (is.null(W2)==FALSE)
	{
		if (is.null(W1)==FALSE)
		{
			# both
			}
		if (is.null(W1))
		{
	
		N <- length(y)
	A <- solve( diag( N ) - rho * W2 )
	mu <- x %*% beta
	
	Sigma <- A %*% t(A)
	

	}
	}
	
	actor.list <- c(1:N)
	for (k in c(1:N))
	{
		notAct <- actor.list!=k
		invSigma22 <- solve(Sigma[notAct,notAct])
		muhat <- mu[k]+Sigma[k,notAct] %*% invSigma22 %*% (y[notAct,1] -mu[notAct,1])
		Sigmahat <- Sigma[k,k]-Sigma[k,notAct] %*% invSigma22 %*% Sigma[notAct,k]
		if (yobs[k]==0){
		zlower <- -Inf
		zupper <- 0	
		}
		if (yobs[k]==1){
		zlower <- 0
		zupper <- Inf
		}
		y[k] <- truncNormalIndep(mu=muhat,variance=Sigmahat,zlower=zlower,zupper=zupper)
	}
	y
}



draw.nam.trunc.norm <- function(x, beta, alpha=NULL,W1=NULL)
{
	
	p <- dim(x)[2]
				
			# network autocorrelation model
			N <- dim(x)[1]
	A <- solve( diag( N ) - alpha * W1 )
y <- A %*% (x %*% beta + rnorm(N))
y[y<=0] <- 0
y[y>0] <- 1
y
}


gof.nam.trunc <- function(x, beta, alpha,W1)
{
	p <- dim(beta)[1]
	num.thetas <- dim(beta)[2]
	n <- dim(x)[1]
	gof.y <- matrix(0,n,num.thetas)
	for (i in c(1:num.thetas))
	{
		gof.y[,i] <- draw.nam.trunc.norm(x=x, beta=beta[,i], alpha=alpha[i],W1=W1)
	}
	gof.y
}


gof.table <- function(obs.stats,sim.stats,name.vec,tabname,compare=FALSE,sim.stats2=NULL,pvalues=FALSE,save.tab='latex',directed=TRUE, Imp.gof =NULL)
{
  
	p <- length(obs.stats)
	if (p != length(name.vec))
	{
		warning('length differs between observed statistics and name vector')
	}
	
	
	# if there are missing
	if (!is.null(Imp.gof)){
	  obs.stats <- matrix(0,p,1)
	  Imp.gof<- Imp.gof[,round(seq(1,dim(Imp.gof)[2], length.out = dim(sim.stats)[2] ) )]
	for ( i in c(1:p))
	{
	 
	  obs.stats[i] <- mean(Imp.gof[i,])
	}
	  
	}
	
	if (dim(sim.stats)[2] <=p )
	{
		sim.stats <- t(sim.stats)
		if (dim(sim.stats)[1] !=p)
		{
			warning('missmatch in dimensions')
			
		}
		
	}
	if (compare)
	{
		
	
	if (dim(sim.stats2)[2] <=p )
	{
		sim.stats2 <- t(sim.stats2)
		if (dim(sim.stats)[1] !=p)
		{
			warning('missmatch in dimensions')
			
		}
		
	}
	}
	# obs, mean, sd, z
	if (compare==FALSE)
	{
		if (pvalues==FALSE)
		{
	ResTab <- matrix(0,p,4)
	}
	if (pvalues==TRUE)
		{
		# obs, mean, pval
		ResTab <- matrix(0,p,3)
		}
	}
	if (compare)
	{
			if (pvalues==FALSE)
		{

	ResTab <- matrix(0,p,7)
	}
	if (pvalues==TRUE)
		{
ResTab <- matrix(0,p,5)
}
	}
	
#	browser()
	for (i in c(1:p))
	{
		if (pvalues==FALSE)
		{
		ResTab[i,1] <- obs.stats[i]
		ResTab[i,2] <- mean(sim.stats[i,])
		ResTab[i,3] <- sd(sim.stats[i,])
		ResTab[i,4] <- (ResTab[i,1]-ResTab[i,2])/ResTab[i,3] 
		if (compare)
	{
		ResTab[i,5] <- mean(sim.stats2[i,])
		ResTab[i,6] <- sd(sim.stats2[i,])
		ResTab[i,7] <- (ResTab[i,1] -ResTab[i,5])/ResTab[i,6]
		}
		}
		if (pvalues==TRUE)
		{
		ResTab[i,1] <- obs.stats[i]
		ResTab[i,2] <- mean(sim.stats[i,])
				if (ResTab[i,1]>ResTab[i,2] )
		{
		if (!is.null(Imp.gof))
		{
		  pobs <- mean(sim.stats[i,] > Imp.gof[i,])/2
		  
		}
				  if (is.null(Imp.gof))
				  {
			pobs <- mean(sim.stats[i,] > ResTab[i,1])/2
				  }
		}
		if (ResTab[i,1]<ResTab[i,2] )
		{
		  if (!is.null(Imp.gof))
		  {
		    pobs <- mean(sim.stats[i,] < Imp.gof[i,])/2
		    
		  }
		  if (is.null(Imp.gof))
		  {
			pobs <- mean(sim.stats[i,] < ResTab[i,1])/2
		  }
		}
		ResTab[i,3] <- pobs

		
		if (compare)
	{
		ResTab[i,4] <- mean(sim.stats2[i,])
		if (ResTab[i,1]>ResTab[i,4] )
		{
			pobs <- mean(sim.stats2[i,] > ResTab[i,1])/2
		}
		if (ResTab[i,1]<ResTab[i,4] )
		{
			pobs <- mean(sim.stats2[i,]< ResTab[i,1])/2
		}
		
		ResTab[i,5] <- pobs 
		}
		}
		
		
	}
	

	
	ResTab <- as.table(ResTab)
	if (pvalues==FALSE){
	if (compare)
	{
colnames(ResTab) <- c('obs','mean','sd','z-stat','mean','sd','z-stat')
}
if (compare==FALSE)
	{

colnames(ResTab) <- c('obs','mean','sd','z-stat')
}
}

if (pvalues==TRUE){
	if (compare)
	{
colnames(ResTab) <- c('obs','mean','p-val','mean','p-val')
}
if (compare==FALSE)
	{

colnames(ResTab) <- c('obs','mean','p-val')
}
}


rownames(ResTab) <- name.vec

if (directed==FALSE)
{
  # remove directed effects

  ResTab <- ResTab[!(rownames(ResTab) %in% c('recip','indegree','in2star','intria','transtri')),]
  rownames(ResTab)[rownames(ResTab)=='outdegree'] <- 'degree'
  rownames(ResTab)[rownames(ResTab)=='transitive'] <- 'triadic'
  rownames(ResTab)[rownames(ResTab)=='out2star'] <- '2star'
  rownames(ResTab)[rownames(ResTab)=='outtria'] <- 'triangle'
  
}

	print(ResTab)
	if (save.tab=='latex'){
print(xtable(ResTab, type = "latex"), file = paste(tabname,".tex",sep=''))
	}
	if (save.tab=='csv'){
	write.csv(ResTab,file= paste(tabname,".csv",sep=''))
	  }
}


test.missing <- function(n=NULL,propmiss=NULL,iterations=NULL,ave.ties=4)
{
	# create network data
	# source("/Users/johankoskinen/Desktop/frombackup/johanadmid/melbourne 2012/melbourne 2017/data hub position/Burnet/dean gof/MultivarALAAMalt.R")
	# RESULTS <- test.missing(n=150,propmiss=0.00265252,iterations=3000)
	# ,propmiss=.025
	require('network')
	require('sna')
	require('ergm')
	#ADJ <- rgraph(n,tprob = ave.ties/(n-1))
	#sim.net <- simulate(as.network(ADJ)~edges+mutual+twopath+dgwnsp(decay=log(2),fixed=TRUE),coef=c(-.45-log(dim(ADJ)[1]), 2,-.01,.25 ))
#	print( gden(	sim.net) )
#	ADJ <- as.matrix.network(sim.net)
	#data(coleman)
	#ADJ <- coleman[2,,]
	doscratch <- FALSE
	directed <- TRUE
	
	useNets <- c(101,103,114,154,156,212,213,248,328,464,475,476,516)
	out.put	 <- get.sbc.net.internal(useNets)
	ADJ <- out.put$BigAdj
	n <- dim(ADJ)[1]
	dim(out.put$BigCov)

	covars <- cbind(colSums(ADJ), out.put$BigCov[,c(1,8,12)] )
	covars[,4] <- (covars[,4]-mean(covars[,4]))/sd(covars[,4])
		colnames(covars) <- c('indegree','sex','attitude','marks')
	if (doscratch==TRUE){
	y <- matrix(runif(n)>.5,n,1)*1
#pdf('SBC_test.pdf')
#plot(as.network(ADJ))
#dev.off()
	ALAAMobj <- prepALAAMdata(y=y,ADJ=ADJ,covariates=covars,directed=directed, contagion ='simple')
	theta <- c(-5.5,.35,-0.1,0,.48,.9)
		
		# generate data
	sim.data <- simulate.alaam(ALAAMobj=ALAAMobj,statsvec=NULL,theta=theta,thinning = 100, NumIterations = 1000, burnin = 5000, returnNet= TRUE, DoSave=TRUE)
#	sim.data <- simulate.alaam(ALAAMobj=ALAAMobj,statsvec=NULL,theta=theta,thinning = 10, NumIterations = 10000, burnin = 1, returnNet= TRUE)
	PropSigma <- solve(cov(t(sim.data$statsvec)))
	# plot(ts(t(sim.data$statsvec)))
	 print( table(sim.data$y[,1000]) )
	y <- sim.data$y[,1000]
	#browser()
	#	print(table(sim.data$y))
	#	browser()
	#	browser()
	# fit model
	save(y , ADJ, covars, file='data_used_for_missing.RData')
	}
	if (doscratch==FALSE){
	load('data_used_for_missing.RData')
	n <- dim(ADJ)[1]
	}
		res.1 <- BayesALAAM(y=y,ADJ=ADJ ,directed=directed,covariates=covars ,Iterations=iterations,silent=FALSE, initcontagion = 0.5 , scaling = 1.3, burnin = 6000,saveFreq=100)#,PropSigma=PropSigma)
	
	# plot(ts(res.1$Thetas))
	
	# test standard MAR
	yMAR <- y
	
	set.missing.MAR <- sample(c(1:n),size=ceiling(propmiss*n), replace=FALSE)
		yMAR[set.missing.MAR] <- NA
	res.mar <- BayesALAAM(y=yMAR,ADJ=ADJ ,directed=directed,covariates=covars ,Iterations=iterations,silent=FALSE, initcontagion = 0.5 , scaling = 1.3, burnin = 6000,saveFreq=100)#,PropSigma=PropSigma)
	#plot(ts(res.mar$Thetas))
	# missingCovs = NULL,missingPhi =NULL
	# # 9. 'Safe in School'
# 10. 'interest in school'
# 12. 'Marks'

missingCovs <- cbind( y, matrix(1,n,1),out.put$BigCov[,c(9,10,12)] )
missingCovs[,5] <- (missingCovs [,5]-mean(missingCovs[,5]))/sd(missingCovs[,5])
	missingPhi <- matrix(c(-.15 ,-1.5, -.15, -0.15,0.2),5,1)
	pmiss <- 1/(1+exp( -missingCovs %*% missingPhi ))
#	hist(pmiss,main=paste('ave: ',round(mean(pmiss),3)))

	set.missing.NMAR <- runif(n)<pmiss
	#	hist(pmiss,main=paste('ave: ',round(mean(pmiss),3)))
print(table(set.missing.NMAR))
	yNMAR <- y
	yNMAR[set.missing.NMAR] <- NA
	res.mar.nmar <- BayesALAAM(y=yNMAR,ADJ=ADJ ,directed=directed,covariates=covars ,Iterations=iterations,silent=FALSE, initcontagion = 0.5 , scaling = 1.3, burnin = 6000)#,PropSigma=PropSigma)
res.nmar <- BayesALAAM(y=yNMAR,ADJ=ADJ ,directed=directed,covariates=covars ,Iterations=iterations,silent=FALSE, initcontagion = 0.5 , scaling = 1.3,missingPhi=missingPhi[1], burnin = 6000)#,PropSigma=PropSigma)

AllRes <- list(res.1=res.1, res.mar=res.mar ,res.mar.nmar = res.mar.nmar , res.nmar=res.nmar)	

AllRes
}

get.sbc.net.internal <- function(useNets)
{
	num_nets <- 619
sizeNet <- matrix(0,num_nets,1)
for (k in c(1:num_nets)){
filenameNet <- paste('/Users/johankoskinen/Desktop/frombackup/SBC networks/ SBCschoolNetwork',k,'.csv',sep='')
schoolNet <- as.matrix(read.csv(filenameNet,header=FALSE))
sizeNet[k] <- dim(schoolNet)[1]
}

# 101 103 114 154 156 212 213 248 328 464 475 476 516
num_nets <- length(useNets)
ADJ <- vector('list',length=num_nets)
y <- vector('list',length=num_nets)
covariates <- vector('list',length=num_nets)
for (k in c(1:num_nets))
{
chooseSchool <- useNets[k] # just an example if you wanted to use school number 605
filenameCov <- paste('/Users/johankoskinen/Desktop/frombackup/SBC networks/ SBCschoolCovars',chooseSchool,'.csv',sep='')

filenameNet <- paste('/Users/johankoskinen/Desktop/frombackup/SBC networks/ SBCschoolNetwork',chooseSchool,'.csv',sep='')
# filenameNet <- paste('/Users/johankoskinen/Documents/files from old macbook pro/johankoskinen/Documents/RSiena course in Zurich/SOFI/ERGM/SBCschoolNetworkBellaMiss',chooseSchool,'.csv',sep='')

schoolCov <- as.matrix(read.csv(filenameCov,header=FALSE))
schoolNet <- as.matrix(read.csv(filenameNet,header=FALSE))
# 1. intend to apply to higher secondary education (yes=1, no/don't know =0, missing=99)
# 2. 'Sex' (male=0, fem=1)
# 3. 'Welfare'
# 4.'Parent psych'
# 5.'Social class 1'
# 6. 'Social class 2'
# 7. 'Social class 3'
# 8. 'Social class 4'
# 9. 'Family attitude'
# 10. 'Safe in School'
# 11. 'interest in school'
# 12. test score on verb, spat & num test
# 13. 'Marks'
y[[k]] <-  schoolCov[,1]
#y[[k]][y[[k]]==99] <- 0
covariates[[k]] <- schoolCov[,c(2:13)]
covariates[[k]][covariates[[k]]==99] <- 0
ADJ[[k]] <- schoolNet
sizeNet[k] <- dim(schoolNet)[1]
if (k==1)
{
	BigAdj <- schoolNet
	BigCov <- covariates[[k]]
	
}
if (k>1)
{
	prevSize <- dim(BigAdj)[1]
	BigAdj <- cbind(BigAdj,matrix(0,prevSize, sizeNet[k] ) )
	BigAdj <- rbind(BigAdj,cbind(matrix(0,sizeNet[k], prevSize ), schoolNet) )
	BigCov <- rbind(BigCov,covariates[[k]])
	
}

}

out.put <- list(BigAdj = BigAdj, BigCov = BigCov)

out.put	
	
}


plot.post.stats <- function(n=NULL,cd = NULL, sim.out.4=NULL,N=NULL)
{
	xaxname <- c('L','S','T')
	lambdaT <- 2
lambdaS <- 1.1
tableT <- penging_ergm(lambdaT,n-1)
tableS <- penging_ergm(lambdaS,n-1)
take.out <- which(cd$csize==n)
if (length(take.out)>1)
{
	bbb <- length(take.out)
	Zobs <- matrix(0,3,bbb)
	for (r in c(1:length(take.out)))
	{
		verticesSamp <- which(cd$membership == take.out[r] )
smallNet <- BigNet[verticesSamp,verticesSamp]

lambdaT <- 2
lambdaS <- 1.1
tableT <- penging_ergm(lambdaT,n-1)
tableS <- penging_ergm(lambdaS,n-1)
Zobs[,r] <- InitializeLazegaERGM(X=smallNet, lambdaT=2,lambdaS=1.1, n = n ,p = 3 ,tableT=tableT, tableS =tableS)$Z_obs

}
}
if (length(take.out)==1)
{
	bbb <- 1
}
	for (k in c(1:3)){
		out.stats <- matrix(sim.out.4$simStats[,k+1,sim.out.3$simStats[1,1,]==n],N,bbb)
		for (m in c(1:bbb))
		{
			if (m==1)
			{
				plot( out.stats[,1], type='l',ylab='n',xlab=paste('n: ',n,' ',xaxname[k]))
				grid( nx = NULL)
	abline(h = Zobs[k,m], col='red')
			}
			if (m>1)
			{
				lines( out.stats[,1])
				
	abline(h = Zobs[k,m], col='red')
}
		}
		
		}
	

}


get.gof.distribution <- function(NumIterations=NULL,res=NULL,burnin=100,thinning= 1000, contagion ='simple',do.loop=FALSE,user.covars=NULL)
{
  
  
  
  if (contagion == 'none')
  {
    contagion <- 'simple' 
  }
  num.thetas <- dim(res$Thetas)[2]

  start_time <- Sys.time()
  Big.stats <- simulate.alaam(ALAAMobj=res$ALAAMobj,theta=res$Thetas[floor(seq(burnin,num.thetas, length.out =NumIterations)),],contagion =contagion, thinning = thinning, NumIterations = NumIterations, burnin = burnin, DoSave=TRUE, returnNet= TRUE, doGOF=TRUE)
  
  end_time <- Sys.time()
  cat('Simulating GOF took ',end_time - start_time,'\n')
  
  
  ADJ <- edglist.to.matrix(ALAAMobj = res$ALAAMobj,n = length(res$ALAAMobj$y))
  gof.covs <- prep.gof.cov(ADJ=ADJ,directed=TRUE)
  if (!is.null(user.covars))
  {
    if (class(user.covars)=="data.frame")
    {
      cat('user supplied covaraites, user.covars, must be class matrix \n')
    }
    gof.covs <- cbind(gof.covs,user.covars)
    
  
  }
  
  #browser()
  #browser()
  
  #Sav.gof <- matrix(0,length(stats),dim(Big.stats$y)[2] )
  
  # start_time <- Sys.time()
  # for (i in c(1:dim(Big.stats$y)[2] ))
  # {
  #   Sav.gof[,i] <-  alaam.gof.statscalc(y=Big.stats$y[,i],ADJ=ADJ,covariates=gof.covs,directed=TRUE)
  # }
  # end_time <- Sys.time()
  # cat('Calculating statistics took ',end_time - start_time,'\n')


  mod.1 <- prepALAAMdata(y = Big.stats$y[,1],ADJ = ADJ,covariates = gof.covs,directed=TRUE, contagion =c('simple','recip','indirect','closedind','transitive'))
  stats.1 <- getALAAMstats(ALAAMobj=mod.1,contagion =c('simple','recip','indirect','closedind','transitive'))
  start_time <- Sys.time()
  
  Sav.gof.2 <- simulate.alaam.NULL(ALAAMobj=mod.1,
                                   contagion =c('simple','recip','indirect','closedind','transitive'), 
                                   DoSave=TRUE, returnNet= TRUE, doGOF=TRUE, sim.data= Big.stats$y)
  
  Sav.gof <- Sav.gof.2$statsvec
  end_time <- Sys.time()
  cat('Calculating statistics took ',end_time - start_time,'\n')


  #pairs(~Sav.gof[1,]+Sav.gof[2,]+Sav.gof.2$statsvec[1,]+Sav.gof.2$statsvec[2,],
   #     main="Simple Scatterplot Matrix")
  
 # pairs(~Sav.gof[3,]+Sav.gof[4,]+Sav.gof.2$statsvec[3,]+Sav.gof.2$statsvec[4,],
  #      main="Simple Scatterplot Matrix")
  stats <- alaam.gof.statscalc(y=res$ALAAMobj$y,ADJ=ADJ,covariates=gof.covs,directed=TRUE)
  #browser()

  if (!is.null(res$ALAAMobj$canchangeMiss))
  {
   # 
    if (do.loop==FALSE){
      mod.1 <- prepALAAMdata(y = res$imputed.obs[,1],
                             ADJ = ADJ,covariates = gof.covs,directed=TRUE, 
                             contagion =c('simple','recip','indirect','closedind','transitive'))
      
    Imp.gof <- simulate.alaam.NULL(ALAAMobj=mod.1,
                           contagion =c('simple','recip','indirect','closedind','transitive'), 
                           DoSave=TRUE, returnNet= TRUE, doGOF=TRUE, 
                           sim.data= res$imputed.obs)
    Imp.gof <- Imp.gof$statsvec
    
    
    }
    if (do.loop==TRUE){
      Imp.gof <- matrix( 0,length(stats),dim(res$imputed.obs )[2] )
    for (i in c(1:dim(res$imputed.obs )[2] ))
    {
      Imp.gof[,i] <-  alaam.gof.statscalc(y=res$imputed.obs[,i],ADJ=ADJ,covariates=gof.covs,directed=TRUE)
    }
    }
  }
  
  if (is.null(res$ALAAMobj$canchangeMiss))
  {
    Imp.gof<- NULL
    #
  }
  

  gof.stats.names<- c('intercept' ,'simple cont.','recip cont.','indirect cont.','closedind cont.','transitive cont.',
                      colnames(gof.covs))
  gof.out <- list(Sav.gof=Sav.gof,ADJ=ADJ,stats=stats, y = res$ALAAMobj$y, gof.stats.names=gof.stats.names, Imp.gof=Imp.gof,gof.y =Sav.gof.2$y)
  
  gof.out 
}



simulate.alaam.NULL <- function(ALAAMobj,statsvec=NULL,contagion ='simple', DoSave=FALSE, returnNet= FALSE,doGOF=FALSE,sim.data=NULL)
{
  y <- ALAAMobj$y
  n <- length(y)
  directed <- ALAAMobj$directed
  interaction <- ALAAMobj$interaction
  covariates <- ALAAMobj$covariates
  p <- dim(covariates)[2]+1# add 1 for intercept
  p1 <- p
  
  if (is.null(statsvec))
  {
    statsvec <- getALAAMstats(ALAAMobj,contagion =contagion)
  }
  if ('simple' %in% contagion){
    p <- p+1
    
    EdgeList <- ALAAMobj$EdgeList
    RowIn <- ALAAMobj$RowIn
    degree <- ALAAMobj$degree
    covariates <- ALAAMobj$covariates
    if (directed==TRUE){
      EdgeListIn <- ALAAMobj$EdgeListIn 
      degreein <- ALAAMobj$degreein
      RowInIn <- ALAAMobj$RowInIn
    }
    
  }
  
  if (('simple' %in% contagion)==FALSE){
    error('model without simple contagion not yet specified')
  }
  dorec <- FALSE
  if ('recip' %in% contagion){
    p <- p+1
    dorec <- TRUE
    
    
    RecEdgeList <- ALAAMobj$RecEdgeList
    Recdegree <- ALAAMobj$Recdegree
    RecRowIn <- ALAAMobj$RecRowIn
  }
  doindir <- FALSE
  if ('indirect' %in% contagion)
  {
    p <- p+1
    doindir <- TRUE
    
    
    P2EdgeList <- ALAAMobj$P2EdgeList
    P2degree <- ALAAMobj$P2degree
    P2RowIn <- ALAAMobj$P2RowIn
    if (directed==TRUE){
      P2EdgeListIn <- ALAAMobj$P2EdgeListIn
      P2degreeIn <- ALAAMobj$P2degreeIn
      P2RowInIn <- ALAAMobj$P2RowInIn
    }
    if (is.null(P2EdgeList) | is.null(P2degree) | is.null(P2RowIn))
    {
      error('edgelist, 2-path degree, and index to 2-path edgelist required if indirect contagion requested')	
    }
  }
  doCloseindir <- FALSE
  if ('closedind' %in% contagion)
  {
    p <- p+1
    doCloseindir <- TRUE
    CP2EdgeList <- ALAAMobj$CP2EdgeList
    CP2degree <- ALAAMobj$CP2degree
    CP2RowIn <- ALAAMobj$CP2RowIn
    
    if (directed==TRUE){
      CP2EdgeListIn <- ALAAMobj$CP2EdgeListIn
      CP2degreeIn <- ALAAMobj$CP2degreeIn
      CP2RowInIn <- ALAAMobj$CP2RowInIn
    }
  }
  doTrans <- FALSE
  if ('transitive' %in% contagion)
  {
    p <- p+1
    doTrans <- TRUE
    translist <- ALAAMobj$translist
    transFirst <- ALAAMobj$transFirst
    transMiddle <- ALAAMobj$transMiddle
    transEndnode <- ALAAMobj$transEndnode
    
  }
  
  
  p2 <- 0
  if (!is.null(interaction))
  {
    p2 <- length(interaction)
    p <- p+p2
  }
  
  out <- simulateALAAM.NULL(y=y,EdgeList=EdgeList,
                            RowIn=RowIn,
                            degree = degree,
                            contagion = contagion,
                            EdgeListIn=EdgeListIn,
                            RowInIn=RowInIn,
                            degreein = degreein,
                            covariates = covariates,
                            statsvec=statsvec ,
                       DoSave=DoSave,
                       returnNet=returnNet, 
                       directed=directed,
                       interaction=interaction, 
                       RecEdgeList=RecEdgeList,
                       Recdegree=Recdegree,
                       RecRowIn=RecRowIn, 
                       P2EdgeList=P2EdgeList,
                       P2degree=P2degree, 
                       P2RowIn=P2RowIn, 
                       P2EdgeListIn=P2EdgeListIn, 
                       P2degreeIn=P2degreeIn, 
                       P2RowInIn=P2RowInIn,
                       CP2EdgeList = CP2EdgeList,
                       CP2degree = CP2degree,
                       CP2RowIn =CP2RowIn,
                       CP2EdgeListIn = CP2EdgeListIn,
                       CP2degreeIn = CP2degreeIn,
                       CP2RowInIn = CP2RowInIn,
                       translist = translist,
                       transFirst = transFirst,
                       transMiddle = transMiddle,
                       transEndnode = transEndnode,
                       doGOF=doGOF,
                       sim.data=sim.data)
  
  
  out
}

simulateALAAM.NULL <- function(y,EdgeList,RowIn,degree,covariates,statsvec,directed=FALSE,EdgeListIn=NULL, degreein=NULL,RowInIn=NULL,DoSave=FALSE,canchange=NULL,returnNet=FALSE,interaction=NULL,checkAlgorithm=FALSE, bias = NULL, contagion ='simple',RecEdgeList=NULL,Recdegree=NULL,RecRowIn=NULL,P2EdgeList=NULL,P2degree=NULL,P2RowIn=NULL, P2EdgeListIn=NULL,P2degreeIn=NULL,P2RowInIn	=NULL, CP2EdgeList=NULL, CP2degree = NULL, CP2RowIn = NULL , CP2EdgeListIn = NULL , CP2degreeIn = NULL, CP2RowInIn = NULL, translist = NULL , transFirst = NULL , transMiddle = NULL , transEndnode = NULL,doGOF=FALSE, sim.data = NULL)
{
  #browser()
  ####
  
  if (doGOF==TRUE){
   
    num.thetas <- dim(sim.data)[2]
  }
  
  #####
  n <- length(y)
  p <- dim(covariates)[2]+1# add 1 for intercept
  p1 <- p
  if ('simple' %in% contagion){
    p <- p+1
  }
  
  if (('simple' %in% contagion)==FALSE){
    error('model without simple contagion not yet specified')
  }
  dorec <- FALSE
  if ('recip' %in% contagion){
    p <- p+1
    dorec <- TRUE
    if (is.null(RecEdgeList) | is.null(Recdegree) | is.null(RecRowIn))
    {
      error('edgelist, reciprocated degree, and index to reciprocated edgelist required if reciprocated contagion requested')	
    }
    
  }
  doindir <- FALSE
  if ('indirect' %in% contagion)
  {
    p <- p+1
    doindir <- TRUE
    
    if (is.null(P2EdgeList) | is.null(P2degree) | is.null(P2RowIn))
    {
      error('edgelist, 2-path degree, and index to 2-path edgelist required if indirect contagion requested')	
    }
  }
  doCloseindir <- FALSE
  if ('closedind' %in% contagion)
  {
    p <- p+1
    doCloseindir <- TRUE
  }
  doTrans <- FALSE
  if ('transitive' %in% contagion)
  {
    p <- p+1
    doTrans <- TRUE
  }
  
  p2 <- 0
  if (!is.null(interaction))
  {
    
    p2 <- length(interaction)
    p <- p+p2
  }
  ystar <- y
  
  statsvecstar <- matrix(0,p,1)
  
  
  
  if (is.null(canchange))
  {
    # regular sampling
    canchange <- c(1:n)
    
  }
  
  if ( DoSave )# we output statistics
  {
    
    SimStats <- matrix(0,p,num.thetas)
    SampleNet<- matrix(0,n,num.thetas)
    saveHere <- 1
    SimStats[,1] <- statsvec
    
    
    SampleNet[,1] <-  ystar
    
  }
  
  for (iterations in c(2:num.thetas))
  {
    actors <- which( sim.data[ ,(iterations-1) ] != sim.data[,iterations] )
    if (length(actors)>0){
    for (do.this.dude in c(1:length(actors)))
    {
    actor <- actors[do.this.dude]
    statsvecstar[1] <- 1 - 2*ystar[actor]
    neig <- 0 # for contagion count
    statsvecstar[2] <- 0
    Recneig <- 0 # for reciprocated contagion count
    Indineig <- 0 # for indirect contagion count
    CIndineig <- 0 # for closed indirect contagion count
    TIndineig <- 0 # for transitive contagion count
    if (degree[actor] > 0)
    {
      for (alter in c(1:degree[actor]))
      {
        neig <- neig + ystar[ EdgeList[RowIn[actor]+alter-1,2] ]#assumes edgelist
        # where the ties of actor starts at row RowIn[actor]
      }
      #neig <- neig + sum(ystar[ EdgeList[RowIn[actor]+c(1:degree[actor])-1,2] ])
    }
    
    if (directed)
    {
      if (degreein[actor]>0)
      {
        for (alter in c(1:degreein[actor]))
        {
          neig <- neig + ystar[ EdgeListIn[RowInIn[actor]+alter-1,2] ]#assumes edgelist
          # where the ties of actor starts at row RowIn[actor]
        }
        #neig <- neig + sum(ystar[ EdgeListIn[RowInIn[actor]+c(1:degreein[actor])-1,2] ])
      }
    }
    
    ### reciprocated contagion if requested
    if (dorec)
    {
      #ALAAMobj$RecEdgeList <- RecEdgeList
      #ALAAMobj$Recdegree <- Recdegree
      #ALAAMobj$RecRowIn <- RecRowIn
      
      if (Recdegree[actor]>0)
      {
        for (alter in c(1:Recdegree[actor]))
        {
          Recneig <- Recneig + ystar[ RecEdgeList[RecRowIn[actor]+alter-1,2] ]#assumes edgelist
          # where the ties of actor starts at row Row[actor]
        }
      }
      
    }
    
    ### indirect contagion if requested
    if (doindir)
    {
      #ALAAMobj$P2EdgeList <- P2EdgeList
      #ALAAMobj$P2degree <- P2degree
      #ALAAMobj$P2RowIn <- P2RowIn
      if (is.na(P2degree[actor]))
      {
        print('degree distribution is undefined')
        #browser()
        
      }
      if (P2degree[actor]>0)
      {
        for (alter in c(1:P2degree[actor]))
        {
          Indineig <- Indineig + ystar[ P2EdgeList[P2RowIn[actor]+alter-1,2] ]*P2EdgeList[P2RowIn[actor]+alter-1,3]#assumes edgelist
          # where the ties of actor starts at row Row[actor]
        }
      }
      
      if (directed)
      {
        
        if (P2degreeIn[actor]>0)
        {
          for (alter in c(1:P2degreeIn[actor]))
          {
            Indineig <- Indineig + ystar[ P2EdgeListIn[P2RowInIn[actor]+alter-1,2] ]*P2EdgeListIn[P2RowInIn[actor]+alter-1,3]#assumes edgelist
            # where the ties of actor starts at row Row[actor]
          }
        }
        
      }
      
    }
    
    ### indirect + direct if requested
    if (doCloseindir)
    {
      if (CP2degree[actor]>0)
      {
        for (alter in c(1:CP2degree[actor]))
        {
          CIndineig <- CIndineig + ystar[ CP2EdgeList[CP2RowIn[actor]+alter-1,2] ]*CP2EdgeList[CP2RowIn[actor]+alter-1,3] #assumes edgelist
          # where the ties of actor starts at row Row[actor]
        }
      }
      
      if (directed)
      {
        if (CP2degreeIn[actor]>0)
        {
          for (alter in c(1:CP2degreeIn[actor]))
          {
            CIndineig <- CIndineig + ystar[ CP2EdgeListIn[CP2RowInIn[actor]+alter-1,2] ]*CP2EdgeListIn[CP2RowInIn[actor]+alter-1,3] #assumes edgelist
            # where the ties of actor starts at row Row[actor]
          }
        }
      }
      
    }
    
    ### transitive contagion
    if (doTrans)
    {
      
      copyTP <- transFirst[[actor]]
      numTps <- length(copyTP)
      if (numTps>0)
      {
        for (rowind in c(1:numTps))
        {
          if (ystar[translist[copyTP[rowind],2]]==1 && ystar[translist[copyTP[rowind],3]]==1)
          {
            TIndineig <- TIndineig + 1
          }
        }
      }
      copyTP <- transMiddle[[actor]]
      numTps <- length(copyTP)
      if (numTps>0)
      {
        for (rowind in c(1:numTps))
        {
          if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],3]]==1)
          {
            TIndineig <- TIndineig + 1
          }
        }
      }
      copyTP <- transEndnode[[actor]]
      numTps <- length(copyTP)
      if (numTps>0)
      {
        for (rowind in c(1:numTps))
        {
          if (ystar[translist[copyTP[rowind],1]]==1 && ystar[translist[copyTP[rowind],2]]==1)
          {
            TIndineig <- TIndineig + 1
          }
        }
      }
      
      
      
    }
    
    
    statsvecstar[2] <- statsvecstar[1]*neig
    k <- 3
    if (dorec)
    {
      statsvecstar[k] <- statsvecstar[1]*Recneig
      k <- k+1
    }
    if (doindir)
    {
      statsvecstar[k] <- statsvecstar[1]*Indineig
      k <- k+1
    }
    
    if (doCloseindir)
    {
      statsvecstar[k] <- statsvecstar[1]*CIndineig
      k <- k+1
      
    }
    if (doTrans)
    {
      statsvecstar[k] <- statsvecstar[1]*TIndineig
      k <- k+1
    }
    
    for (t in c(1:(p1 - 1) ))
    {
      statsvecstar[k] <- statsvecstar[1]*covariates[actor,t]
      k <- k+1
    }
    
    if (!is.null(interaction))
    {
      #statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction]
      statsvecstar[p] <- statsvecstar[1]*neig*covariates[actor,interaction[1]]
    }
    
   
      
      statsvec <- statsvec+statsvecstar
      ystar[actor] <- 1 - ystar[actor]
    }
    }
    
    
    
    
      SimStats[,iterations] <- statsvec
      
      
        SampleNet[,iterations] <-  ystar
      
    
  }
  
  if ( DoSave )# we output statistics
  {
    
    statsvec <- SimStats
  }
  if ( DoSave==FALSE )# we output statistics
  {
    
  }
  
  if (returnNet==TRUE)
  {
    if (DoSave)
    {
      
      ystar <- SampleNet
    }
    statsvec <-list(y=ystar,statsvec=statsvec)
  }
  
  
  statsvec 
  
}

get_indep_cov <- function(ALAAMobj,directed=FALSE,theta=NULL,interaction=NULL,logitOnly=FALSE,contagion ='simple')
{
  n <- length(ALAAMobj$y)
  p <- dim(ALAAMobj$covariates)[2]+1# add 1 for intercept
  p1 <- p
  pstart <- 2
  if ('simple' %in% contagion){
    p <- p+1
    pstart <- pstart +1
  }
  
  
  dorec <- FALSE
  if ('recip' %in% contagion){
    p <- p+1
    pstart <- pstart +1
    dorec <- TRUE
    
  }
  doindir <- FALSE
  if ('indirect' %in% contagion)
  {
    p <- p+1
    doindir <- TRUE
    pstart <- pstart +1
  }
  
  doCloseindir <- FALSE
  if ('closedind' %in% contagion)
  {
    p <- p+1
    doCloseindir <- TRUE
    pstart <- pstart +1
    
  }
  doTrans <- FALSE
  if ('transitive' %in% contagion)
  {
    p <- p+1
    doTrans <- TRUE
    pstart <- pstart +1
    
  }		 
  
  p2 <- 0
  if (!is.null(interaction))
  {
    p2 <- length(interaction)
    p <- p+p2
  }
  
  ind.covs <- matrix(0,n,length(theta))
  ind.covs[,pstart:(p-p2)] <- ALAAMobj$covariates
  ind.covs[,1] <- 1
  
  
  if (!is.null(interaction))
  {
    ind.covs[,p] <- 0
    
  }
  
  ind.covs
  
}

BayesALAAM.formula <- function(ALAAMobj,Iterations=1000,silent=FALSE,
                               PropSigma=NULL,scaling = 1,
                               do.scaling = TRUE,
                               initcontagion = NULL,
                               burnin = NULL,
                               missingCovs = NULL,missingPhi =NULL,
                               priorSigma=NULL,
                               priorMu=NULL,
                               scalePrior=NULL,
                               MPLE=FALSE, 
                               saveFreq=NULL,
                               missFreq=100,
                               initTheta = NULL,
                               thinning=1,
                               par.burnin=1)
{
  #  ALAAMobj <- prepALAAMdata(y=y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =contagion, interaction= interaction)
  
  checkMiss <- FALSE
  require(MASS)
  require('mvtnorm')
  require('coda')
  
  ignore.cont <- FALSE
  
  
  if ('none' %in% ALAAMobj$contagion)
  {
    # run independent 
    MPLE <- TRUE
    ignore.cont <- TRUE
    #contagion <- 'simple'
    ALAAMobj$contagion <- 'simple'
    ## but we want to reset this later on
    
  }
  
  
  if (is.null(saveFreq))
  {	
    saveFreq <-  ceiling(Iterations/10)
  }
  if (is.null(burnin))
  {
    burnin <- 3000
  }
  if (any(is.na(ALAAMobj$y)))
  {
    # first: save the original data, with missings
    y.with.missings <- ALAAMobj$y
    # set missing to 0 and add loop
    canchangeMiss <- which(is.na(ALAAMobj$y))
    ALAAMobj$y[is.na(ALAAMobj$y)] <- 0
    print(paste('you have ',length(canchangeMiss),' missing entries in response'))
    print(ALAAMobj$canchange)
    bias <- matrix(0,length(ALAAMobj$y),1)
    
    missFreq <- missFreq
    save.miss.imp <- matrix(0,length(ALAAMobj$y),missFreq)
    
    put.imp.here <- 1
    save.miss.imp[,put.imp.here ] <- ALAAMobj$y
    put.imp.here <- put.imp.here + 1
  }
  else {canchangeMiss <- NULL}
  
 # ALAAMobj <- prepALAAMdata(y=y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =contagion, interaction= interaction)
  
  # statsvec <- getALAAMstats(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,covariates = ALAAMobj$covariates,directed=directed,interaction=interaction, contagion =contagion,
  # RecEdgeList=ALAAMobj$RecEdgeList,
  # Recdegree=ALAAMobj$Recdegree,
  # RecRowIn=ALAAMobj$RecRowIn, 
  # P2EdgeList=ALAAMobj$P2EdgeList,
  # P2degree=ALAAMobj$P2degree,
  # P2RowIn=ALAAMobj$P2RowIn, 
  # P2EdgeListIn=ALAAMobj$P2EdgeListIn,
  # P2degreeIn=ALAAMobj$P2degreeIn,
  # P2RowInIn=ALAAMobj$P2RowInIn)
  
  statsvec <- getALAAMstats(ALAAMobj , contagion =ALAAMobj$contagion)
  
  ALAAMobj$statsvec <- statsvec
  #ALAAMobj$directed <- directed
  #ALAAMobj$canchangeMiss <- NULL
  #ALAAMobj$canchange <- NULL
  if (!is.null(canchangeMiss))
  {
    ALAAMobj$canchangeMiss <- canchangeMiss
    
  }
  
  if (!is.null( ALAAMobj$canchange))
  {
    ALAAMobj$non.fix <-  ALAAMobj$canchange
  }
  #if (is.null(canchange))
  #{
  #	tem.obj <- NULL
  #ALAAMobj$canchange <- tem.obj
  #}
  
  
  
  # clean up
  if (!is.null(ALAAMobj$canchangeMiss) && checkMiss==TRUE){
    print('check canchange here at the beginning')
    browser()
  }
  ## === below reported in estimate.alaam instead === ####
  #	print(paste('you have chosen directed equals ',ALAAMobj$directed))
  #	print(paste('you have chosen degree equals ',useDegree))
  
  p <- length(ALAAMobj$statsvec)
  print(paste('you have p: ', p))
  cat('\nobserved stats: ', round(ALAAMobj$statsvec ,3) ) 
  theta <- matrix(0,1,p)
  cat('\nnumber of covariates: ', dim(ALAAMobj$covariates)[2] ) 
  
  
  tuneCons <- tuneALAAMprop(ALAAMobj,directed=ALAAMobj$directed,theta=theta,
                            interaction=ALAAMobj$interaction, contagion =ALAAMobj$contagion)
  if (is.null(PropSigma))
  {
    sigma.epsilon <-(scaling/sqrt(p))*tuneCons$PropSigma
  }
  if (is.null(PropSigma)==FALSE)
  {
    if (do.scaling){
    sigma.epsilon <-(scaling/sqrt(p))*PropSigma
    }
    if (do.scaling==FALSE){
      sigma.epsilon <- PropSigma
    }
  }
  if (!is.null(initcontagion))
  {
    tuneCons$theta[2] <- initcontagion
  }
  if (!is.null(missingPhi))
  {
    # this used to be detemined by  (!is.null(missingCovs))
    # but note that only missingPhi is needed
    print(paste('MNAR mechanism specified based on response and ', dim(missingCovs)[2],' covariates'))
    # define logit(p_miss_1) = missingPhi[1]*y[i]+missingPhi[2]*missingCovs[i,1]
    # define logit(p_miss_0) = missingPhi[1]*y[i]+missingPhi[2]*missingCovs[i,1]
    #bias <- missingPhi[2]+log(1+exp(missingPhi[1]+missingPhi[3]* colSums( ADJ )))-log(1+exp(missingPhi[1]+missingPhi[2]+missingPhi[3]* colSums( ADJ )))
    
    for (i in c(1:length(ALAAMobj$y)))
    {
      # missingPhi[1]*y[i]
      # missingPhi[2] is intercept
      # missingPhi[3:dim(missingCovs)[2]] are other covariates
      if (i %in% canchangeMiss)
      {
        bias[i] <-  missingPhi[1]
      }
    }
    
  }
  #if (is.null(missingCovs))
  #{
  #	bias <- NULL
  #	}
  if (!is.null(scalePrior))
  {
    cat('\npseudo-conjugate prior will be used with scaling factor: ',scalePrior)
    if (is.null(priorMu))
    {
      # use 0 prior mean
      priorMu <- matrix(0,1,p)
    }
    zStats <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,
                            degree = ALAAMobj$degree,
                            EdgeListIn=ALAAMobj$EdgeListIn,
                            RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,
                            covariates = ALAAMobj$covariates,NumIterations=100,
                            theta=priorMu,statsvec=ALAAMobj$statsvec ,
                            DoSave=TRUE,directed=ALAAMobj$directed,
                            interaction=ALAAMobj$interaction, thinning = 1000,
                            burnin = 3000,
                            contagion = ALAAMobj$contagion,
                            RecEdgeList=ALAAMobj$RecEdgeList,
                            Recdegree=ALAAMobj$Recdegree,
                            RecRowIn=ALAAMobj$RecRowIn, 
                            P2EdgeList=ALAAMobj$P2EdgeList,
                            P2degree=ALAAMobj$P2degree, 
                            P2RowIn=ALAAMobj$P2RowIn, 
                            P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
                            P2degreeIn=ALAAMobj$P2degreeIn, 
                            P2RowInIn=ALAAMobj$P2RowInIn ,
                            CP2EdgeList = ALAAMobj$CP2EdgeList,
                            CP2degree = ALAAMobj$CP2degree,
                            CP2RowIn = ALAAMobj$CP2RowIn,
                            CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
                            CP2degreeIn = ALAAMobj$CP2degreeIn,
                            CP2RowInIn = ALAAMobj$CP2RowInIn,
                            translist = ALAAMobj$translist,
                            transFirst = ALAAMobj$transFirst,
                            transMiddle = ALAAMobj$transMiddle,
                            transEndnode = ALAAMobj$transEndnode)
    priorSigma <- solve(cov(t(zStats)))*scalePrior
  }
  
  if (!is.null(priorSigma))
  {
    require('mvtnorm')
  }
  # cleanup
  if (checkMiss==FALSE){
    #rm(ADJ,covariates,y)
    rm(statsvec)
    #rm(directed)
    #rm(canchange)
  }
  #### begin estimation
  # calculate thin.Iterations
  if (par.burnin>Iterations)
  {
    stop('Burnin (par.burnin) cannot be greater than the number of iterations')
  }

  store.inds <- seq(from=par.burnin,to=Iterations,by=thinning)
  thin.Iterations <- length(store.inds)
  cat('\nA thinning of ',thinning,', ')
  cat('\n(parameter) burn-in of ',par.burnin,', and ')
  cat('\niterations of ',Iterations,',')
  cat('\nyields a sample size of ',thin.Iterations,'.\n')
  ThetaCorrent <- matrix(0,thin.Iterations,p)
  
  if (is.null(initTheta))
  {
  theta <- tuneCons$theta
  }
  if (is.null(initTheta)==FALSE)
  {
    theta <- initTheta
  }
  #ThetaCorrent[1,] <- theta 
  theta.current  <- theta 
  if (!is.null(ALAAMobj$canchangeMiss) && checkMiss==TRUE){
    print('check canchange here at the middle')
    #browser()
  }
  
  #### Count accepts
  acc.ratio <- 0
  
  # Initialise theta
  put.theta.here <- 1
  ThetaCorrent[put.theta.here,] <- theta.current
  #### allow for MPLE ###########
  if (MPLE==TRUE)
  {
    # 
    currentlike <- 0
    obsstats <- changestatsALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,
                                 RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,EdgeListIn=ALAAMobj$EdgeListIn,
                                 RowInIn=ALAAMobj$RowInIn,degreein = ALAAMobj$degreein,
                                 covariates = ALAAMobj$covariates,NumIterations=100,theta=priorMu,
                                 statsvec=ALAAMobj$statsvec ,DoSave=TRUE,directed=ALAAMobj$directed,
                                 interaction=ALAAMobj$interaction, thinning = 1000,burnin = 3000,
                                 contagion = ALAAMobj$contagion,
                                 RecEdgeList=ALAAMobj$RecEdgeList,
                                 Recdegree=ALAAMobj$Recdegree,
                                 RecRowIn=ALAAMobj$RecRowIn, 
                                 P2EdgeList=ALAAMobj$P2EdgeList,
                                 P2degree=ALAAMobj$P2degree, 
                                 P2RowIn=ALAAMobj$P2RowIn, 
                                 P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
                                 P2degreeIn=ALAAMobj$P2degreeIn, 
                                 P2RowInIn=ALAAMobj$P2RowInIn ,
                                 CP2EdgeList = ALAAMobj$CP2EdgeList,
                                 CP2degree = ALAAMobj$CP2degree,
                                 CP2RowIn = ALAAMobj$CP2RowIn,
                                 CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
                                 CP2degreeIn = ALAAMobj$CP2degreeIn,
                                 CP2RowInIn = ALAAMobj$CP2RowInIn,
                                 translist = ALAAMobj$translist,
                                 transFirst = ALAAMobj$transFirst,
                                 transMiddle = ALAAMobj$transMiddle,
                                 transEndnode = ALAAMobj$transEndnode)
    
    yobs <- ALAAMobj$y
    
    if (!is.null(ALAAMobj$ALAAMobj$non.fix))
    {
      obsstats <- obsstats[,ALAAMobj$non.fix]
      yobs <- yobs[ALAAMobj$non.fix]
      
    }
    
    
    obsstats <- t(obsstats)
    n<- length(yobs)
    # y,x,theta,n,p,variable=TRUE
    # drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)
    # if there are missing
    # missingUpdate$y <- yobs
    # missingUpdate$y[ALAAMobj$canchangeMiss] <- drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)[ALAAMobj$canchangeMiss]
    # 
    if (ignore.cont )
    {
      theta[2] <- 0
    }
    if (!is.null(ALAAMobj$canchangeMiss) ){
      missingUpdate <- list()
      missingUpdate$y <- yobs
      missingUpdate$y[ALAAMobj$canchangeMiss] <- drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)[ALAAMobj$canchangeMiss]
      
      missFreq <- ceiling(Iterations /missFreq )
      save.miss.imp[,1 ] <- missingUpdate$y
      
    }		 	
    
    currentlike <-straightlike(y=yobs,x=obsstats,theta=theta ,n=n)
    cat('Estimation using pseudo likelihood:\n')
    
    #### change: from storing every 'theta' to storing theta.current
    
    theta.current <- theta 
    
    cat('\n(Printing every ',saveFreq,' iterations)\n')
    
    for (iteration in c(1:Iterations))
    {
     
      theta1 <- rmvnorm(1,theta.current, sigma = sigma.epsilon)
     # theta1 <- rmvnorm(1,ThetaCorrent[iteration-1,], sigma = sigma.epsilon)
      if (ignore.cont )
      {
        theta1[2] <- 0
      }
      likestar <- straightlike(y=yobs,x=obsstats,theta=theta1 ,n=n)
      Hratio <- likestar-currentlike
      if (Hratio >= log(runif(1))) {
        
        theta.current  <- theta1
        currentlike <- likestar
        acc.ratio <- acc.ratio +(1/Iterations)
      }
      #ThetaCorrent[iteration,] <- theta
      
      
      ### update missing if present
      if (!is.null(ALAAMobj$canchangeMiss)){
        
        
        missingUpdate$y[ALAAMobj$canchangeMiss] <- drawstraight(y=yobs,x=obsstats,theta=theta ,n=n,p=p)[ALAAMobj$canchangeMiss]
        
        
      }
      
      if (!is.null(ALAAMobj$canchangeMiss) & ((iteration %% missFreq)==0) )
      {
        put.imp.here <- min(put.imp.here,dim(save.miss.imp)[2])
        save.miss.imp[,put.imp.here ] <- missingUpdate$y
        put.imp.here <- put.imp.here + 1
      }
      
      if ( (silent==FALSE) && ((iteration %% saveFreq)==0) )
      {
        cat('\nyou have done ',iteration,' iterations out of ',Iterations,' \ntheta:', round(theta.current  ,3) ) 	 
      }
      ## Store  <- theta : but 
      if ((iteration %in% store.inds))
      {
        ThetaCorrent[put.theta.here,] <- theta.current
        put.theta.here <- put.theta.here + 1
      } 
    }		
   
    

  }
  
  
  ################################
  
  # check ((iteration %in% store.inds))
  
  if (MPLE==FALSE)
  {
    # theta.current <- theta 
    if (!is.null(ALAAMobj$canchangeMiss) )
    {
      missFreq <- ceiling(Iterations /missFreq )
    }
    for (iteration in c(1:Iterations))
    {
      theta1 <- rmvnorm(1,theta.current, sigma = sigma.epsilon)
      
      if (checkMiss==TRUE){
        tempALAAMobj <- prepALAAMdata(y=ALAAMobj$y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =contagion)
        statsvec <- getALAAMstats(tempALAAMobj , contagion =contagion)
        #dicrepancy <- sum(statsvec != ALAAMobj$statsvec)
        dicrepancy <- sum( abs(statsvec - ALAAMobj$statsvec)>0.0001)
        
        if (dicrepancy!=0)
        {
          print('missmatch between observed and simulated stats in main loop')
          browser()
        }
      }
      
      statStar <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,
                                RowIn=ALAAMobj$RowIn,
                                degree = ALAAMobj$degree,
                                EdgeListIn=ALAAMobj$EdgeListIn,
                                RowInIn=ALAAMobj$RowInIn,
                                degreein = ALAAMobj$degreein,
                                covariates = ALAAMobj$covariates,
                                NumIterations=burnin,
                                theta=theta1,
                                statsvec=ALAAMobj$statsvec,
                                DoSave=FALSE,directed=ALAAMobj$directed, 
                                interaction=ALAAMobj$interaction, 
                                contagion =ALAAMobj$contagion , 
                                RecEdgeList=ALAAMobj$RecEdgeList,
                                Recdegree=ALAAMobj$Recdegree,
                                RecRowIn=ALAAMobj$RecRowIn, 
                                P2EdgeList=ALAAMobj$P2EdgeList,
                                P2degree=ALAAMobj$P2degree, 
                                P2RowIn=ALAAMobj$P2RowIn, 
                                P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
                                P2degreeIn=ALAAMobj$P2degreeIn, 
                                P2RowInIn=ALAAMobj$P2RowInIn,
                                CP2EdgeList = ALAAMobj$CP2EdgeList,
                                CP2degree = ALAAMobj$CP2degree,
                                CP2RowIn = ALAAMobj$CP2RowIn,
                                CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
                                CP2degreeIn = ALAAMobj$CP2degreeIn,
                                CP2RowInIn = ALAAMobj$CP2RowInIn,
                                translist = ALAAMobj$translist,
                                transFirst = ALAAMobj$transFirst,
                                transMiddle = ALAAMobj$transMiddle,
                                transEndnode = ALAAMobj$transEndnode,
                                canchange = ALAAMobj$non.fix)
      
      if (!is.null(ALAAMobj$canchangeMiss) && checkMiss==TRUE){
        print('check canchange')
        browser()
      }
      if (!is.null(priorSigma))
      {
        pr <- dmvnorm(rbind(theta1, theta.current ), 
                      mean = priorMu, 
                      sigma = priorSigma,
                      log=TRUE)
        prirat <- pr[1] - pr[2]
      }
      if (is.null(priorSigma))
      {
        prirat <- 0
      }
      Hratio <- (theta.current - theta1) %*% ( statStar-ALAAMobj$statsvec )+prirat
     # ThetaCorrent[iteration,] <- ThetaCorrent[iteration-1,]
      if (Hratio >= log(runif(1))) {
        
        theta.current <- theta1
        acc.ratio <- acc.ratio +(1/Iterations)
      }
      ### update missing if present
      if (!is.null(ALAAMobj$canchangeMiss)){
        # supplying values to canchange=NULL,returnNet=FALSE
        missingUpdate <- simulateALAAM(y=ALAAMobj$y,EdgeList=ALAAMobj$EdgeList,RowIn=ALAAMobj$RowIn,degree = ALAAMobj$degree,
                                       EdgeListIn=ALAAMobj$EdgeListIn,
                                       RowInIn=ALAAMobj$RowInIn,
                                       degreein = ALAAMobj$degreein,
                                       covariates = ALAAMobj$covariates,
                                       NumIterations=burnin,
                                       theta=theta.current,
                                       statsvec=ALAAMobj$statsvec ,
                                       DoSave=FALSE,
                                       directed=ALAAMobj$directed,
                                       canchange=ALAAMobj$canchangeMiss ,
                                       returnNet=TRUE,
                                       interaction=ALAAMobj$interaction,
                                       bias = bias,  
                                       contagion =ALAAMobj$contagion,
                                       RecEdgeList=ALAAMobj$RecEdgeList,
                                       Recdegree=ALAAMobj$Recdegree,
                                       RecRowIn=ALAAMobj$RecRowIn, 
                                       P2EdgeList=ALAAMobj$P2EdgeList,
                                       P2degree=ALAAMobj$P2degree, 
                                       P2RowIn=ALAAMobj$P2RowIn, 
                                       P2EdgeListIn=ALAAMobj$P2EdgeListIn, 
                                       P2degreeIn=ALAAMobj$P2degreeIn, 
                                       P2RowInIn=ALAAMobj$P2RowInIn,
                                       CP2EdgeList = ALAAMobj$CP2EdgeList,
                                       CP2degree = ALAAMobj$CP2degree,
                                       CP2RowIn = ALAAMobj$CP2RowIn,
                                       CP2EdgeListIn = ALAAMobj$CP2EdgeListIn,
                                       CP2degreeIn = ALAAMobj$CP2degreeIn,
                                       CP2RowInIn = ALAAMobj$CP2RowInIn,
                                       translist = ALAAMobj$translist,
                                       transFirst = ALAAMobj$transFirst,
                                       transMiddle = ALAAMobj$transMiddle,
                                       transEndnode = ALAAMobj$transEndnode)
        if (checkMiss==TRUE){
          dicrepancy <- sum(missingUpdate$y[!(c(1:length(ALAAMobj$y)) %in% ALAAMobj$canchangeMiss)]!=ALAAMobj$y[!(c(1:length(ALAAMobj$y)) %in% ALAAMobj$canchangeMiss)])
          
          if (dicrepancy!=0)
          {
            print('missmatch between observed and simulated')
          }
          
          tempALAAMobj <- prepALAAMdata(y=missingUpdate$y,ADJ=ADJ,covariates=covariates,directed=directed,useDegree = useDegree, contagion =ALAAMobj$contagion)
          statsvec <- getALAAMstats(tempALAAMobj , contagion =contagion)
          dicrepancy <- sum( abs(statsvec - missingUpdate$statsvec)>0.0001)
          if (dicrepancy!=0)
          {
            print('missmatch between observed and simulated stats')
            browser()
          }
        }
        ALAAMobj$y <- missingUpdate$y
        ALAAMobj$statsvec <- missingUpdate$statsvec
        
      }
      #### finished updating missing
      
      if ( (silent==FALSE) && ((iteration %% saveFreq)==0) )
      {
        
        cat('\nyou have done ',iteration,' iterations out of ',Iterations,' \ntheta:', round(theta.current ,3) ) 	
        save(ThetaCorrent,file='BayesALAAMdump.RData')
        if (!is.null(ALAAMobj$canchangeMiss)){
          cat('\nobserved stats: ', round(ALAAMobj$statsvec ,3),'\n' ) 
          cat('\nimputed ones: ', sum(missingUpdate$y[ALAAMobj$canchangeMiss]),' out of ',length(ALAAMobj$canchangeMiss),'\n' ) 
        }
        
        
      }
      
      # save imputed missing observations
      if (!is.null(ALAAMobj$canchangeMiss) & ((iteration %% missFreq)==0) )
      {
        put.imp.here <- min(put.imp.here,dim(save.miss.imp)[2])
        save.miss.imp[,put.imp.here ] <- ALAAMobj$y
        put.imp.here <- put.imp.here + 1
      }
      
      if ((iteration %in% store.inds))
      {
        ThetaCorrent[put.theta.here,] <- theta.current
        put.theta.here <- put.theta.here + 1
      }
    }
  }
  
  ResTab <- matrix(0,p,5)
 
  # no longer needed to thin or subtract burnin post hoc
 # ResTab[,1] <- colMeans(ThetaCorrent[floor(.1*Iterations):Iterations,])
  ResTab[,1] <- colMeans(ThetaCorrent)
 # ResTab[,2] <- apply(ThetaCorrent[floor(.1*Iterations):Iterations,],2,sd)
  ResTab[,2] <- apply(ThetaCorrent,2,sd)
  # ResTab[,3] <- effectiveSize(ThetaCorrent[floor(.1*Iterations):Iterations,])
  ResTab[,3] <- effectiveSize(ThetaCorrent)
  for (k in c(1:p)){
    ResTab[k,4:5] <- as.numeric(acf(ThetaCorrent[,k],plot=FALSE)[c(10,30)][[1]])
  }
  ResTab <- as.table(ResTab)
  colnames(ResTab) <- c('mean','sd','ESS','SACF 10','SACF 30')
  rownames(ResTab)[1] <- c('intercept')
  colnames(ThetaCorrent) <- rep('eff',dim(ThetaCorrent)[2])
  colnames(ThetaCorrent)[1]<-c('intercept')
  k <- 2
  if ('simple' %in% ALAAMobj$contagion)
  {
    rownames(ResTab)[k] <- c('contagion')
    colnames(ThetaCorrent)[k]<-c('contagion')
    k<- k+1
    
  }
  
  if ('recip' %in% ALAAMobj$contagion)
  {
    rownames(ResTab)[k] <- c('recip cont')
    colnames(ThetaCorrent)[k]<-c('recip cont')
    k<- k+1
  }
  
  if ('indirect' %in% ALAAMobj$contagion)
  {
    rownames(ResTab)[k] <- c('indirect cont')
    colnames(ThetaCorrent)[k]<-c('indirect cont')
    k<- k+1
  }
  
  if ('closedind' %in% ALAAMobj$contagion)
  {
    rownames(ResTab)[k] <- c('closed indirect cont')
    colnames(ThetaCorrent)[k]<-c('closed indirect cont')
    k<- k+1
  }
  doTrans <- FALSE
  if ('transitive' %in% ALAAMobj$contagion)
  {
    rownames(ResTab)[k] <- c('transitiv cont')
    colnames(ThetaCorrent)[k]<-c('transitiv cont')
    k<- k+1
  }
  
  #if (useDegree){
  #rownames(ResTab)[k] <- 'outdegree'
  #k <- k+1
  #} 
  useDegree <- FALSE
  
  if (is.null(colnames(ALAAMobj$covariates)))
  {
    for (i in c(1:dim(ALAAMobj$covariates)[2] ) ) {
      if (useDegree & i==1){
        rownames(ResTab)[k] <- 'outdegree'
        colnames(ThetaCorrent)[k]<-  'outdegree'
        
        k <- k+1
      } 
      if (useDegree==FALSE & i>1){
        rownames(ResTab)[k] <- paste('cov ',i)
        colnames(ThetaCorrent)[k]<- paste('cov ',i)
        k <- k+1
      }
    }
    
  }
  if (!is.null(colnames(ALAAMobj$covariates)) )
  {
    if (useDegree && is.null(colnames(ALAAMobj$covariates)[1] ))
    {
      colnames(ALAAMobj$covariates)[1] <- 'outdegree'
      
    }
    for (i in c(1:dim(ALAAMobj$covariates)[2] ) ){
      rownames(ResTab)[k] <- colnames(ALAAMobj$covariates)[i]
      colnames(ThetaCorrent)[k]<- colnames(ALAAMobj$covariates)[i]
      k <- k +1
    }
    
  }
  
  if (!is.null(ALAAMobj$interaction))
  {
    p2 <- length(ALAAMobj$interaction)
    for (i in c(1:p2))
    {
      if (is.null(colnames(ALAAMobj$covariates) ) ) {
        rownames(ResTab)[k] <- paste('interaction ',ALAAMobj$interaction[i])
        colnames(ThetaCorrent)[k]<- paste('interaction ',ALAAMobj$interaction[i])
      }
      if (!is.null(colnames(ALAAMobj$covariates))){
        rownames(ResTab)[k] <- paste('interaction ',ALAAMobj$interaction[i],': ',colnames(ALAAMobj$covariates)[ALAAMobj$interaction[i]])
        colnames(ThetaCorrent)[k]<- paste('interaction ',colnames(ALAAMobj$covariates)[ALAAMobj$interaction[i]])
      }
      
      k <- k +1
    }
    
  }
  
  
  
  cat('\nsummaries of the posterior draws:\n')
  
  print(round(ResTab,4))
  cat('\nif ESS is small you need to increase the number of Iterations')
  cat('\nif and/or increase the parameter thinning which is currently ',thinning,'.')
  cat('\nAlternatively, adjust the proposals by improving the proposal variance-covariance:')
  cat('\ne.g. set PropSigma equal to the covariance of the posterior draws, or')
  cat('\nincrease/decrease the scaling of the proposals (default: scaling=1/sqrt(',p,'); current: ',scaling,'/sqrt(',p,')')
  cat('\nAcceptance ratio: ',acc.ratio,' (idealy around 0.25)')
 
  
  if (is.null(ALAAMobj$canchangeMiss))
  {
    # there is no missing
    imputed.obs <- NULL
  }
  
  if (!is.null(ALAAMobj$canchangeMiss))
  {
    # there are missings
    imputed.obs <- save.miss.imp
    # restore the missings
    ALAAMobj$y <- y.with.missings
    # if this is not done, then missings will magically dissapear
  }
  
  if (ignore.cont)
  {
    # reset contagion to 'none'
    ALAAMobj$contagion <- 'none'
  }
  

  
  Results <- list(Thetas=ThetaCorrent, 
                  ALAAMobj = ALAAMobj , 
                  ResTab = ResTab , 
                  priorSigma = priorSigma, 
                  priorMu = priorMu,
                  imputed.obs =imputed.obs,
                  silent=silent,
                  do.scaling = FALSE,
                  PropSigma = sigma.epsilon,  # note that this will only work with do.scaling==FALSE
                  MPLE=MPLE,
                  saveFreq=saveFreq,
                  missFreq=missFreq,
                  initTheta=  theta.current,
                  acc.ratio = acc.ratio)
  
  
  
  # BayesALAAM.formula <- function(ALAAMobj,Iterations=1000,silent=FALSE,
  #                                PropSigma=NULL,scaling = 1,
  #                                initcontagion = NULL,
  #                                burnin = NULL,
  #                                missingCovs = NULL,missingPhi =NULL,
  #                                priorSigma=NULL,
  #                                priorMu=NULL,
  #                                scalePrior=NULL,
  #                                canchange=NULL,MPLE=FALSE, saveFreq=NULL,missFreq=100)
  # 
  
  Results
}

get.alaam.obj.from.form <- function(temp.alaam.obj,adjacency=NULL,canchange=NULL)
{
  # in getALAAMstats interactions defined as:
  # statsvecstar[1]*neig*covariates[actor,interaction[1]]
  
  ALAAMobj <- prepALAAMdata(y=temp.alaam.obj$y,
                            ADJ=adjacency,
                            covariates=as.matrix(temp.alaam.obj$covariates),
                            directed=temp.alaam.obj$directed,
                            useDegree = FALSE,
                            contagion =temp.alaam.obj$contagion,
                            interaction= temp.alaam.obj$interaction.by.col)
  
  ALAAMobj$formula <- temp.alaam.obj$formula
  ALAAMobj$canchange <- canchange
  ALAAMobj
}


# calling BayesALAAM by function
# BayesALAAM.formula <- function(ALAAMobj,Iterations=1000,silent=FALSE,
#                                PropSigma=NULL,scaling = 1,
#                                initcontagion = NULL,
#                                burnin = NULL,
#                                missingCovs = NULL,missingPhi =NULL,
#                                priorSigma=NULL,
#                                priorMu=NULL,
#                                scalePrior=NULL,
#                                canchange=NULL,MPLE=FALSE, saveFreq=NULL,missFreq=100)
estimate.alaam <- function(formula, data, adjacency=NULL,
                           Iterations=1000,
                           prevBayes=NULL,
                           recalibrate =FALSE,
                           silent=FALSE,
                           PropSigma=NULL,
                           scaling = 1,
                           do.scaling = TRUE,
                           initcontagion = NULL,
                           burnin = NULL,
                           missingCovs = NULL,
                           missingPhi =NULL,
                           priorSigma=NULL,
                           priorMu=NULL,
                           scalePrior=NULL,
                           canchange=NULL,
                           MPLE=FALSE,
                           saveFreq=NULL,
                           missFreq=100,
                           initTheta=NULL,
                           thinning=1,
                           par.burnin=1)
{
  form <- as.formula(formula)
  ## for future robustness:
  ## check if 'intercept' is specified or a covariate is constant
  ## then remove from formula
  if (is.null(prevBayes)){
    if (class(data)=="alaam.obj")
    {
      # Data formated and we do not need to create edgelist, or anything else for that matter
      
    }
    if (class(data)=="data.frame")
    {
      # Data in data frame and we need to format stuff
      if (is.null(adjacency)==TRUE)
      {
        stop('An adjacency matrix is required')
      }
      directed <- sum(adjacency != t(adjacency))
      if (directed)
      {
        cat('Network is directed with ',dim(adjacency)[1],'nodes \n')
      }
      else
      {
        cat('Network is undirected  ',dim(adjacency)[1],'nodes \n')
      }
      temp.alaam.obj <- match.alaam.formula(form, data, adjacency)
      ALAAMobj <- get.alaam.obj.from.form(temp.alaam.obj,adjacency)
      ALAAMobj$canchange <- canchange
      
    }
    
    
  }
  ### if prevBayes - continue in previous estimation
  
  if (is.null(prevBayes)==FALSE )
  {
    
    if (recalibrate)
    {
      ## Use posterior to recalibrate proposal 
      # Very weird: contagion = 'none' implies theta1[2] <- 0 in algorithm
      # consequently dim(prevBayes$Thetas)[1] and prevBayes$PropSigma includes contagoin but set to 0
      if (prevBayes$ALAAMobj$contagion=='none')
      {
        non.zero.pars <- c(1,3:dim(prevBayes$Thetas)[2])
      }
      if (prevBayes$ALAAMobj$contagion!='none')
      {
        non.zero.pars <- c(1:dim(prevBayes$Thetas)[2])
      }
      ## first check that acf not too bad
      my.max.lag <- determ.lag.num(prevBayes$Thetas[,non.zero.pars],10,4)[1]
      ROmeg <- robust.cov(prevBayes$Thetas[,non.zero.pars],max.lag = my.max.lag)
      PropSigma <- matrix(0,dim(prevBayes$Thetas)[2],dim(prevBayes$Thetas)[2])
      PropSigma[non.zero.pars,non.zero.pars] <- ROmeg
      ## or simply restart with new theta, burnin, 
      if (is.null(scaling))
      {
        # keep at null
        do.scaling <- FALSE
      }
      if (!is.null(scaling))
      {
        # keep at null
        do.scaling <- TRUE
      }
    }
    # propvar
    if (recalibrate==FALSE)
    {
      if (is.null(PropSigma))
      {
        PropSigma<- prevBayes$PropSigma
        if (is.null(scaling))
        {
          do.scaling <- FALSE
          scaling<- prevBayes$scaling
        }
      }
    }
    
    if (is.null(burnin))
    {
      burnin <- prevBayes$burnin
    }
    
    ALAAMobj <- prevBayes$ALAAMobj
    if (form != ALAAMobj$formula)
    {
      stop('prevBayes estimated with model:\n ',
           as.character(ALAAMobj$formula)[2],as.character(ALAAMobj$formula)[1],as.character(ALAAMobj$formula)[3] )
    }
    silent<- prevBayes$silent
    
    
    initcontagion <- prevBayes$initcontagion
    
    missingCovs <- prevBayes$missingCovs
    missingPhi <- prevBayes$missingPhi
    priorSigma<- prevBayes$priorSigma
    priorMu<- prevBayes$priorMu
    
    scalePrior<- prevBayes$scalePrior
    # might not be necessary:
    # canchange<- prevBayes$canchange
    # cause it shoud be in ALAAMobj$canchange
    # so remove it
    MPLE<- prevBayes$MPLE
    # saveFreq<- prevBayes$saveFreq
    missFreq<- prevBayes$missFreq
    initTheta <- prevBayes$initTheta
    
  }
  
  ### CALL the estimation function
  
  
  if (is.binary(ALAAMobj$y)==FALSE)
  {
    stop('Outcome variable must be dichotomous')
  }
  
  
  
  Results <- BayesALAAM.formula(ALAAMobj,
                                Iterations=Iterations,
                                silent=silent,
                                PropSigma=PropSigma,
                                scaling = scaling,
                                do.scaling = do.scaling,
                                initcontagion = initcontagion,
                                burnin = burnin,
                                missingCovs = missingCovs,
                                missingPhi =missingPhi,
                                priorSigma=priorSigma,
                                priorMu=priorMu,
                                scalePrior=scalePrior,
                                MPLE=MPLE,
                                saveFreq=saveFreq,
                                missFreq=missFreq,
                                initTheta= initTheta,
                                thinning=thinning,
                                par.burnin= par.burnin)
  Results
}

match.alaam.formula <- function(formula, data, adjacency=NULL)
{
  # NOTE: if (useDegree==TRUE) in prepALAAMdata needs to be disabled
  ### - Example usage - ########
  # A mock dataset
  # my.dat <- data.frame(agree=runif(10)>.5,peter=rnorm(10),john=runif(10),sex=runif(10)>.5)
  # Calling the function without an adjacency matrix:
  #
  # match.alaam.formula(as.formula(agree ~ john+peter+simple+recip), my.dat)
  #
  # Calling the function with interactions of covariates:
  #
  # match.alaam.formula(as.formula(agree ~ john+peter+simple+recip+john*peter+peter*sex), my.dat)
  # Calling the function with a network adjacency matrix:
  # 
  # match.alaam.formula(as.formula(agree ~ idegree+john+peter+simple+recip+john*peter+peter*sex), my.dat, rgraph(10,1))
  # Calling the function with a network adjacency matrix AND interaction of contagion:
  # 
  # match.alaam.formula(as.formula(agree ~ idegree+john+peter+simple+recip+john*peter+peter*sex+sex*simple), my.dat, rgraph(10,1))
  
  n <- dim(data)[1]
  mod.terms <- terms(formula)
  if (is.null(adjacency)==FALSE)
  {
    require(sna)
    d.cens <- dyad.census(adjacency)
    if (n!=dim(adjacency)[1])
    {
      stop('adjacency matrix has different dimensions from data')
    }
    directed <- d.cens[2]!=0
    if (directed)
    {
      idegree <- matrix(colSums(adjacency),nrow=n,ncol=1)
      odegree <- matrix(rowSums(adjacency),nrow=n,ncol=1)
      recipties <- matrix( rowSums(adjacency * t(adjacency) ), n , 1)
    }
    if (directed==FALSE)
    {
      degree <- matrix(rowSums(adjacency),nrow=n,ncol=1)
    }
  }
  if (is.null(adjacency)){
    directed <- NULL
  }
  
  # if attr(mod.terms,'factors')
  # has a row that is greater than 1
  is.interaction <- any(rowSums(attr(mod.terms,'factors'))>1)
  # print(is.interaction)
  # there is an interaction
  vars.w.interact <- all.vars(formula)[rowSums(attr(mod.terms,'factors'))>1]
  all.interacts <- attr(mod.terms,'factors')[,colSums(attr(mod.terms,'factors'))>1]
  # Check what contagion type (-s) are requested
  contagion.list <- c('simple','recip','indirect','closedind','transitive')
  contagion <- contagion.list[contagion.list %in% all.vars(formula)] 
  
  
  ## Check if there are any contation interactions
  is.contagion.interaction <- any(contagion %in% vars.w.interact)
  
  if (is.contagion.interaction)
  {
    cont.interaction.terms <- c()
    interaction.by.col <- c()
  }
  if (is.contagion.interaction==FALSE)
  {
    cont.interaction.terms <- NULL # actually the same as c()
    interaction.by.col <- NULL
  }
  ###
  indep.var.names <- all.vars(formula)
  indep.var.names <- indep.var.names[!(indep.var.names %in% c(contagion,as.character(formula[[2]]))) ]
  
  if ('degree' %in% indep.var.names)
  {
    if (is.null(adjacency))
    {
      stop('An adjacency matrix is required')
    }
    if (directed)
    {
      stop('The degree term can only be used with symmetric networks (use idegree or odegree instead)')
    }
    
    data <- cbind(data,degree)
  }
  if ('idegree' %in% indep.var.names)
  {
    if (is.null(adjacency))
    {
      stop('An adjacency matrix is required')
    }
    if (directed==FALSE)
    {
      stop('The idegree term can only be used with directed networks (use degree  instead)')
    }
    
    data <- cbind(idegree,data)
  }
  if ('odegree' %in% indep.var.names)
  {
    if (is.null(adjacency))
    {
      stop('An adjacency matrix is required')
    }
    if (directed==FALSE)
    {
      stop('The odegree term can only be used with directed networks (use degree  instead)')
    }
    
    data <- cbind(odegree,data)
  }
  if ('recipties' %in% indep.var.names)
  {
    if (is.null(adjacency))
    {
      stop('An adjacency matrix is required')
    }
    if (directed==FALSE)
    {
      stop('The recipties term can only be used with directed networks')
    }
    
    data <- cbind(recipties,data)
  }
  
  if ('twostar' %in% indep.var.names)
  {
    
    twostar <- matrix( choose(idegree,2),n,1) #  in-stars refecting dispersion in popularity
    data <- cbind(twostar,data)
  }
  if ('intwostar' %in% indep.var.names)
  {
    
    intwostar <- matrix( choose(idegree,2),n,1) #  in-stars refecting dispersion in popularity
    data <- cbind(intwostar,data)
  }
  
  if ('outtwostar' %in% indep.var.names)
  {
    
    outtwostar <- matrix( choose(odegree,2),n,1) #  out-stars refecting dispersion in activity
    data <- cbind(outtwostar,data)
  }
  
  
  if ('twopath' %in% indep.var.names)
  {
    
    twopath <- idegree*odegree - recipties # correlation between indegree and outdegree
    data <- cbind(twopath,data)
  }
  
  if ('inthreestar' %in% indep.var.names)
  {
    
    inthreestar <- matrix( choose(idegree,3),n,1) # furhter measure of in-degree heterogeneity
    data <- cbind(inthreestar,data)
  }
  
  if ('outthreestar' %in% indep.var.names)
  {
    outthreestar <- matrix( choose(odegree,3),n,1) # furhter measure of out-degree heterogeneity
    data <- cbind(outthreestar,data)
  }
  
  if ('transties'  %in% indep.var.names)
  {
    transties <- rowSums( adjacency* (adjacency %*% t(adjacency) )  ) # embedded in transitive triads
    data <- cbind(transties,data)
  }
  
  if ('indirties'  %in% indep.var.names)
  {
    indirties <- matrix(0,n,1)
    for (i in c(1:n))
    {
      if ( sum( adjacency[ i, ])==1 )
      {
        indirties[i] <- sum(  adjacency[adjacency[ i, ]==1, adjacency[ i, ]==0]  ) - (adjacency[i,] %*% adjacency[,i]>0)
      }
      
      if ( sum(adjacency[ i, ])>1 )
      {
        indirties[i] <- sum( colSums( adjacency[adjacency[ i, ]==1, adjacency[ i, ]==0] ) > 0  ) - (adjacency[i,] %*% adjacency[,i]>0)
      }
    }
    
    data <- cbind(indirties,data)# number of indirect ties
  }
  
  
  y <- getResponse(formula, data)
  covariates <- data[,indep.var.names]
  
  if (is.interaction)
  {
    if (length(all.vars(formula))==length(all.interacts))
    {
      num.interacts <- 1
    }
    if (length(all.vars(formula))<length(all.interacts))
    {
      num.interacts <- dim(all.interacts)[2]
    }
    for (k in c(1:num.interacts))
    {
      if (num.interacts>1){
        # check if any of them a contagion effect
        if (any(contagion %in% names(which(all.interacts[,k]>0)))==FALSE)
        {
          these.interaction <- matrix(apply(data[,names(which(all.interacts[,k]>0))],1,prod),nrow=n,ncol=1)
          
          covariates <- cbind(covariates ,these.interaction)
          #covariates <- cbind(covariates ,assign( paste(names(which(all.interacts[,k]>0)),collapse ='*'), these.interaction))
          names(covariates)[dim(covariates)[2]] <- paste(names(which(all.interacts[,k]>0)),collapse ='*')
        }
        if (any(contagion %in% names(which(all.interacts[,k]>0)))){
          cont.interaction.terms <- c(cont.interaction.terms,names(which(all.interacts[,k]>0))[!(names(which(all.interacts[,k]>0))%in% contagion)])
        }
      }
      if (num.interacts==1){
        if (any(contagion %in% names(which(all.interacts>0)))==FALSE){
          these.interaction <- matrix(apply(data[,names(which(all.interacts>0))],1,prod),nrow=n,ncol=1)
          covariates <- cbind(covariates ,these.interaction)
          #covariates <- cbind(covariates ,assign( paste(names(which(all.interacts[,k]>0)),collapse ='*'), these.interaction))
          names(covariates)[dim(covariates)[2]] <- paste(names(which(all.interacts>0)),collapse ='*')
        }
        if (any(contagion %in% names(which(all.interacts>0))))
        {
          cont.interaction.terms <- c(cont.interaction.terms,names(which(all.interacts>0))[!(names(which(all.interacts>0)) %in% contagion)] )
        }
      }
      
      #assign( paste(names(which(all.interacts[,k]>0)),collapse ='*'), these.interaction)
      
    }
  }
  
  
  
  if (length(cont.interaction.terms)>0)
  {
    # warning('only one covariate can be interacted with one contagion term')
    interaction.by.col <- which(names(covariates) %in% cont.interaction.terms)
    #browser()
  }
  if (length(contagion)==0)
  {
    contagion <- 'none'
  }
  
  # we should also save the independent formula, which would be the covariates... but in formula form
  temp.alaam.obj <- list(y=y,covariates=covariates,
                         contagion=contagion,
                         formula=formula,
                         directed=directed, 
                         cont.interaction.terms=cont.interaction.terms,
                         interaction.by.col=interaction.by.col)
  temp.alaam.obj
}


getResponse <- function(form, data) {
  
  # mf <- match.call(expand.dots = FALSE)
  # m <- match(c("formula", "data"), names(mf), 0L)
  # mf <- mf[c(1L, m)]
  # mf$drop.unused.levels <- TRUE
  # mf[[1L]] <- as.name("model.frame")
  # mf <- eval(mf, parent.frame())
  # y <- model.response(mf, "numeric")
  y <- data[,as.character(form[[2]])]
  y
} 

is.binary <- function(y){
  # check if binary
  y.copy <- y[!is.na(y)]
  out.of.range <- any(range( y.copy)!=c(0,1))
  non.discrete <- length(unique( y.copy))!=2
  (out.of.range==FALSE) | ( non.discrete==FALSE)
}


robust.cov <- function(x,max.lag)
{
  n <- dim(x)[1]
  m <- dim(x)[2]
  if (m>n)
  {
    x <- t(x)
    n <- dim(x)[1]
    m <- dim(x)[2]
  }
  x <- scale(x, center =TRUE,scale=FALSE)
  initial.cov <- get.omega.j(x,0)
  if (max.lag>0){
    for (j in c(1:max.lag))
    {
      Ome <- get.omega.j(x,j)
      initial.cov <- initial.cov + (1-j/(max.lag+1))*( Ome +t( Ome ))
    }
  }
  #initial.cov <- initial.cov -matrix(colMeans(x),m,1) %*% matrix(colMeans(x),1,m)
  initial.cov
}


get.omega.j <- function(x,j)
{
  n <- dim(x)[1]
  m <- dim(x)[2]
  if (m>n)
  {
    x <- t(x)
    n <- dim(x)[1]
    m <- dim(x)[2]
  }
  Omega.j <- matrix(0,m,m)
  for (t in c((j+1):n))
  {
    Omega.j <- Omega.j +matrix(x[t,],m,1) %*% matrix(x[t-j,],1,m)/n
    
  }
  Omega.j
}

determ.lag <- function(x,lag.max)
{
  n <- dim(x)[1]
  m <- dim(x)[2]
  if (m>n)
  {
    x <- t(x)
    n <- dim(x)[1]
    m <- dim(x)[2]
  }
  Lags <- matrix(0,m,lag.max)
  for (i in c(1:m))
  {
    
    Lags[i,]<- acf(x[,i],lag.max=lag.max,plot=FALSE)$acf[1:lag.max]
  }
  Lags
}

determ.lag.num <- function(x,lag.max,tolerance)
{
  
  LAGS <- determ.lag(x,lag.max)
  accep1<- 1
  accep2 <-1
  if (any( colSums(LAGS>.3)>tolerance ))
  {
    accep1 <- max(which(colSums(LAGS>.3)>tolerance))
  }
  if (any(colSums(LAGS>.4)>tolerance ))
  {
    accep2 <-       max(which(colSums(LAGS>.4)>tolerance))
  }
  accep  <- c(accep1,accep2)
  accep 
}

post.deviance.alaam <- function(ALAAMresult,
                                numBridges=20,
                                thinning.like = 5000,
                                sample.size = 200,
                                cov.sample.burnin = NULL,
                                printFreq=10,
                                mult.fact = 30,
                                num.outs=100)
{
  post.devs <- matrix(0,num.outs,1)# we would like to have a measure of uncertainty in there (next year)
  # We assume ALAAMresult$Thetas 
  # nice and thinned
  Thetas <- ALAAMresult$Thetas[seq(from=1,to=dim(ALAAMresult$Thetas)[1],length.out = num.outs) ,]
  p <- dim(Thetas)[2]
  
  indep.index <- get.indep.vars(ALAAMresult$ALAAMobj)
  # determine if independent
  if (ALAAMresult$ALAAMobj$contagion=='none')
  {
    
    # independent model: no need for reference
    for (i in c(1:num.outs))
    {
      # calculate -2log.like for i=1,...,K
      post.devs[i] <- -2*independLike(ALAAMresult$ALAAMobj, Thetas[i,indep.index] )
      #currentlike <-straightlike(y=yobs,x=obsstats,theta=theta ,n=n)
    }
    
  }
  
  if (ALAAMresult$ALAAMobj$contagion!='none')
  {
    # non-independent model
    theta.logit.mle <- matrix(0,1,p)
    # get MLE
    theta.logit.mle[1,indep.index] <- glm(ALAAMresult$ALAAMobj$y~ ALAAMresult$ALAAMobj$covariates, family = binomial(link = "logit"))$coefficients
    # get log-likelihood
    logit.mle.loglike <- independLike(ALAAMresult$ALAAMobj, theta.logit.mle[1,indep.index] )
    # call path sampler
    for (i in c(1:num.outs))
    {
      temp.path.like <- eval.like.path.alaam(ALAAMresult$ALAAMobj,
                                             theta_tilde=Thetas[i,],
                                             theta_star=theta.logit.mle,
                                             burnin=cov.sample.burnin,
                                             numbridges=numBridges,
                                             sampleSize=sample.size,
                                             thinning.like=thinning.like)
      # cat('log-likelihood to be evaluated with reference to constrained MLE:\n')
      # cat(eta.0,'\n')
      # log.like <- t( eta.0 - eta.tilde ) %*% Z.obs+PathEstim-ML.res$Like.vals[15]
      # eta.0: MLE
      temp.log.like <- -( (theta.logit.mle-Thetas[i,])%*% ALAAMresult$ALAAMobj$statsvec+temp.path.like- logit.mle.loglike)
      post.devs[i] <- -2*temp.log.like
      if ((i %% ceiling(num.outs/10))==0)
      {
        cat('\nLikelihood evaluated for draw ',i,' (out of ',num.outs,')')
      }
    }
    
  }
  
  post.devs
}


get.indep.vars <- function(ALAAMobj)
{
  p <- dim(ALAAMobj$covariates)[2]+1# add 1 for intercept
  p1 <- p
  pstart <- 2
  if ('none' %in% ALAAMobj$contagion){
    p <- p+1
    pstart <- pstart +1
    # Why? Well, beucase if 'none' second column set to 0
  }
  if ('simple' %in% ALAAMobj$contagion){
    p <- p+1
    pstart <- pstart +1
  }
  
  
  dorec <- FALSE
  if ('recip' %in% ALAAMobj$contagion){
    p <- p+1
    pstart <- pstart +1
    dorec <- TRUE
    
  }
  doindir <- FALSE
  if ('indirect' %in% ALAAMobj$contagion)
  {
    p <- p+1
    doindir <- TRUE
    pstart <- pstart +1
  }
  
  doCloseindir <- FALSE
  if ('closedind' %in% ALAAMobj$contagion)
  {
    p <- p+1
    doCloseindir <- TRUE
    pstart <- pstart +1
    
  }
  doTrans <- FALSE
  if ('transitive' %in% ALAAMobj$contagion)
  {
    p <- p+1
    doTrans <- TRUE
    pstart <- pstart +1
    
  }		 
  
  p2 <- 0
  if (!is.null(ALAAMobj$interaction))
  {
    p2 <- length(ALAAMobj$interaction)
    p <- p+p2
  }
  
  # theta[pstart:(p-p2)] <- ans.log.Int$coefficients[2:length(ans.log.Int$coefficients)]
  # theta[1] <- ans.log.Int$coefficients[1]
  idep.index <- c(1,pstart:(p-p2))
  idep.index
}

eval.like.path.alaam <- function(ALAAMobj,
                                 theta_tilde,
                                 theta_star,
                                 burnin=3000,
                                 numbridges=20,
                                 sampleSize=100,
                                 thinning.like=3000)
{
  
  theta <- theta_star 
  STATS <- simulateALAAM(y=ALAAMobj$y,
                         EdgeList=ALAAMobj$EdgeList,
                         RowIn=ALAAMobj$RowIn,
                         degree = ALAAMobj$degree,
                         EdgeListIn=ALAAMobj$EdgeListIn,
                         RowInIn=ALAAMobj$RowInIn,
                         degreein = ALAAMobj$degreein,
                         covariates = ALAAMobj$covariates,
                         NumIterations=sampleSize,
                         theta= theta,
                         statsvec=ALAAMobj$statsvec ,
                         DoSave=TRUE,directed=ALAAMobj$directed,
                         interaction=ALAAMobj$interaction,
                         returnNet=TRUE,
                         thinning = thinning.like,
                         burnin = burnin,
                         canchange=ALAAMobj$canchange)
  
  
  
  ALAAMobj$y <- STATS$y[,sampleSize]
  ALAAMobj$statsvec <- STATS$statsvec[,sampleSize]
  delta.theta <-theta_tilde - theta_star
  PathEstim <-  (1/(numbridges+1))*delta.theta  %*% rowMeans(STATS$statsvec)
  rm(STATS)
  
  #PathEstim <- 0
  
  for (RepSamp in  c(1:numbridges) )
  {
    theta <- theta_star * (1 - RepSamp/numbridges) + theta_tilde * (RepSamp/numbridges)
    # update parameter on path from eta.0 to eta.tilde
    #  eta.current <- eta.0 * (1 - t/numBridges) + eta.tilde * (t/numBridges)
    
    STATS <- simulateALAAM(y=ALAAMobj$y,
                           EdgeList=ALAAMobj$EdgeList,
                           RowIn=ALAAMobj$RowIn,
                           degree = ALAAMobj$degree,
                           EdgeListIn=ALAAMobj$EdgeListIn,
                           RowInIn=ALAAMobj$RowInIn,
                           degreein = ALAAMobj$degreein,
                           covariates = ALAAMobj$covariates,
                           NumIterations=sampleSize,
                           theta= theta,
                           statsvec=ALAAMobj$statsvec,
                           DoSave=TRUE,
                           returnNet=TRUE,
                           directed=ALAAMobj$directed,
                           interaction=ALAAMobj$interaction,
                           thinning = thinning.like,
                           burnin = burnin,
                           canchange=ALAAMobj$canchange)
    ALAAMobj$y <- STATS$y[,sampleSize]
    ALAAMobj$statsvec <- STATS$statsvec[,sampleSize]
    #PathEstim = PathEstim + (1/K)*(theta_tilde - theta_star) *  A_obs;
    tempest <- (1/(numbridges+1))*delta.theta %*% rowMeans(STATS$statsvec)
    # delta.path <- (1/(numBridges+1))*delta.eta %*% Zave
    rm(STATS)
    PathEstim <-  PathEstim +  tempest
    # if (silent==FALSE){
    #    print(paste('Path at bridge: ',RepSamp,' contributes: ', round(tempest,3)))
    # }
  }
  
  PathEstim
  
  
  
}

plot.deviance.alaam <- function(dev.1,
                                dev.2=NULL,
                                dev.3=NULL,
                                dev.4=NULL,
                                colPal = c('black','red','darkorange','darkblue'),
                                line.tps = c(1,2,3,4),
                                mod.names = c('M1','M2','M3','M4'))
{
  dev <- dev.1
  dev.sort <- dev[order(dev.1)]
  muppet <- ecdf(dev.sort)
  #xa <- seq(min(RESULTS$deviances),max(RESULTS$deviances),length.out = 200)
  if (is.null(dev.2)==FALSE)
  {
    hor.range <- range(dev.1,dev.2)
    if (is.null(dev.3)==FALSE)
    {
      hor.range <- range(hor.range,dev.3)
      if (is.null(dev.4)==FALSE)
      {
        hor.range <- range(hor.range,dev.4)
        
      }
    }
    
  }
  if (is.null(dev.2))
  {
    hor.range <- range(dev.1)
  }
  xticks.labs <- seq(min(hor.range),max(hor.range),length.out=5)
  muppet.kermit <- c(muppet( dev.sort  ))
  plot(dev.sort ,muppet.kermit , type='l',bty='n' ,
       xlim=hor.range,xaxt='n',
       xlab='deviance',ylab='CDF',cex.lab=.7,cex.axis=0.7,
       col=colPal[1],
       lty=line.tps[1])
  axis(1, at = xticks.labs, cex.axis=0.9, srt=45, col.ticks = "grey",labels = FALSE, las=2)
  text(x = xticks.labs,
       ## Move labels to just below bottom of chart.
       y = par("usr")[3] - 0.085,
       ## Use names from the data list.
       labels = round(xticks.labs,3),
       ## Change the clipping region.
       xpd = NA,
       ## Rotate the labels by 35 degrees.
       srt = 35,
       ## Adjust the labels to almost 100% right-justified.
       adj = 0.965,
       ## Increase label size.
       cex = .5)
  #dev.sort <- dev[order(dev)]
  #se.dev.sort <- se.dev[order(dev)]
  #lines(dev.sort-conf.fact*se.dev.sort,muppet(dev.sort),col='blue')
  #lines(dev.sort+conf.fact*se.dev.sort,muppet(dev.sort),col='blue')
  if (is.null(dev.2)==FALSE)
  {
    dev <- dev.2
    dev.sort <- dev[order(dev)]
    muppet.show <- ecdf(dev.sort)
    #xa <- seq(min(RESULTS.2$deviances),max(RESULTS.2$deviances),length.out = 200)
    lines(dev.sort ,muppet.show( dev.sort  ), col=colPal[2],
          lty=line.tps[2])
    
    if (is.null(dev.3)==FALSE)
    {
      dev <- dev.3
      dev.sort <- dev[order(dev)]
      muppet.show <- ecdf(dev.sort)
      #xa <- seq(min(RESULTS.2$deviances),max(RESULTS.2$deviances),length.out = 200)
      lines(dev.sort ,muppet.show( dev.sort  ), col=colPal[3],
            lty=line.tps[3])
      if (is.null(dev.4)==FALSE)
      {
        dev <- dev.4
        dev.sort <- dev[order(dev)]
        muppet.show <- ecdf(dev.sort)
        #xa <- seq(min(RESULTS.2$deviances),max(RESULTS.2$deviances),length.out = 200)
        lines(dev.sort ,muppet.show( dev.sort  ), col=colPal[4],
              lty=line.tps[4])
      }
    }
    #lines(dev.sort+conf.fact*se.dev.sort,muppet(dev.sort),col='darkorange')
  }
  if (is.null(dev.4)==FALSE)
  {
    # legend for all four
    legend('bottomright',col=colPal,lty=line.tps,mod.names)
  }
  if (is.null(dev.4) &  (is.null(dev.3)==FALSE))
  {
    # legend for all four
    legend('bottomright',col=colPal[1:3],lty=line.tps[1:3],mod.names[1:3])
  }
  if (is.null(dev.3) &  (is.null(dev.2)==FALSE))
  {
    # legend for all four
    legend('bottomright',col=colPal[1:2],lty=line.tps[1:2],mod.names[1:2])
  }
  
}

alaam.dic <- function(Post.dev)
{
  # The idea is that models with smaller DIC should be preferred to models with larger DIC
  dev.bar <- mean(Post.dev)
  pV1 <- var(Post.dev)/2
  DIC <- dev.bar + pV1
  DIC
}
