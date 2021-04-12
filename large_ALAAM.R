# A script for testing 1000+ network
require('sna')
n <- 1050
g <- rgraph(n,m=1, tprob=0.025)# create a large network
y <- matrix(runif(n)<0.5,n,1)# create a binary outcome
x <- matrix(rnorm(n),n.1)# create a covariate
out.degree <-matrix( rowSums(g), n, 1) # number of ties sent
covs <- matrix(NA,n,2)
covs[,1] <- x
covs[,2] <- out.degree
colnames(covs) <- c("x",
                    "outdegree")
##== read source code
source('MultivarALAAMalt.R')
res.0 <- BayesALAAM(y = y,           # dependent variable
                    ADJ = g,           # network
                    covariates = covs,   # covariates
                    directed = TRUE,    # directed / undirecred network
                    Iterations = 500,   # number of iterations
                    saveFreq = 100,      # print and save frequency
                    contagion = 'none')  # type of contagion
##== create a covariate that is the average 'x' of alters
user.covars <- as.data.frame(g %*% matrix(x,n,1 )/(out.degree +1 ))
names(user.covars) <- 'prop.alc.alter'

##== Run GOF
sim.0 <- get.gof.distribution(NumIterations=500, # number of vectors to draw
                              res=res.0, # the ALAAM estimation object that contains model and results
                              burnin=100, # no. iterations discarded from GOF distribution
                              thinning = 1000, # no. iterations between sample points
                              contagion ='none',# should be the same as for model fitted
                              user.covars = user.covars) # the gof will now evaluate the fit for this covariate
##== Create a goodness-of-fit table
gof.table(obs.stats=	sim.0$stats, # observed statistics included  not fitted statistics
          sim.stats=	sim.0$Sav.gof, # simulated goodness-of-fit statistics
          name.vec=	sim.0$gof.stats.names, # names of statistics calculate, not all will be used if undirected
          tabname='ALAAMGofalt', # name of file saved
          pvalues=TRUE, # posterior predictive p-values
          save.tab ='csv', # save a csv file or a LaTex file
          directed=TRUE,
          Imp.gof = sim.0$Imp.gof)