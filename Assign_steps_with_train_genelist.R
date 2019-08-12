library(ASSIGN)

dir.create("tempdir")
tempdir <- "tempdir"

data(trainingData1)
data(testData1)
data(geneList1)
trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10,
                                      myc=1:10, ras=1:10, src=1:10),
                       bcat = 11:19, e2f3 = 20:28, myc= 29:38,
                       ras = 39:48, src = 49:55)
testLabel1 <- rep(c("Adeno", "Squamous"), c(53,58))

#dir.create(file.path(tempdir,"wrapper_example1"))
#assign.wrapper(trainingData=trainingData1, testData=testData1,
#               trainingLabel=trainingLabel1, testLabel=testLabel1,
#               geneList=NULL, n_sigGene=rep(200,5), adaptive_B=TRUE,
#               adaptive_S=FALSE, mixture_beta=TRUE,
#               outputDir=file.path(tempdir,"wrapper_example1"),
#               iter=2000, burn_in=1000)

#dir.create(file.path(tempdir,"wrapper_example2"))
#assign.wrapper(trainingData=trainingData1, testData=testData1,
#               trainingLabel=trainingLabel1, testLabel=NULL,
#               geneList=geneList1, n_sigGene=NULL, adaptive_B=TRUE,
#               adaptive_S=FALSE, mixture_beta=TRUE,
#               outputDir=file.path(tempdir,"wrapper_example2"),
#               iter=2000, burn_in=1000)

#dir.create(file.path(tempdir,"wrapper_example3"))
#assign.wrapper(trainingData=NULL, testData=testData1,
#               trainingLabel=NULL, testLabel=NULL,
#               geneList=geneList1, n_sigGene=NULL, adaptive_B=TRUE,
#               adaptive_S=TRUE, mixture_beta=TRUE,
#               outputDir=file.path(tempdir,"wrapper_example3"),
#               iter=2000, burn_in=1000)

# training dataset is available;
# the gene list of pathway signature is available
processed.data <- assign.preprocess(trainingData=trainingData1,
                                    testData=testData1,
                                    trainingLabel=trainingLabel1,
                                    geneList=geneList1)

mcmc.chain <- assign.mcmc(Y=processed.data$testData_sub,
                          Bg = processed.data$B_vector,
                          X=processed.data$S_matrix,
                          Delta_prior_p = processed.data$Pi_matrix,
                          iter = 2000, adaptive_B=TRUE,
                          adaptive_S=FALSE, mixture_beta=TRUE)

trace.plot <- assign.convergence(test=mcmc.chain, burn_in=0, iter=2000,
                                 parameter="B", whichGene=1,
                                 whichSample=NA, whichPath=NA)

mcmc.pos.mean <- assign.summary(test=mcmc.chain, burn_in=1000,
                                iter=2000, adaptive_B=TRUE,
                                adaptive_S=FALSE, mixture_beta=TRUE)

assign.output(processed.data=processed.data,
              mcmc.pos.mean.testData=mcmc.pos.mean,
              trainingData=trainingData1, testData=testData1,
              trainingLabel=trainingLabel1,
              testLabel=testLabel1, geneList=NULL,
              adaptive_B=TRUE, adaptive_S=FALSE,
              mixture_beta=TRUE, outputDir=tempdir)



