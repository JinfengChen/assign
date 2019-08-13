library(ASSIGN)
args <- commandArgs(trailingOnly = TRUE)
pathway=as.numeric(args[1])

pathwayList <- c("AKT", "BAD", "IGF1R", "ERK", "HER2", "EGFR", "RAF", "KRASQH", "KRASGV", "KRASWT")

pathway

dir.create(pathwayList[pathway])
tempdir <- pathwayList[pathway]

#data(trainingData1)
#data(testData1)
#data(geneList1)
#trainingLabel1 <- list(control = list(bcat=1:10, e2f3=1:10,
#                                      myc=1:10, ras=1:10, src=1:10),
#                       bcat = 11:19, e2f3 = 20:28, myc= 29:38,
#                       ras = 39:48, src = 49:55)
#testLabel1 <- rep(c("Adeno", "Squamous"), c(53,58))

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

######################
#getting the genelist#
######################
gene_list_dir <- "/home/jichen/Projects/Cellines/ABCB1/ASSIGN/GeneLists" 
load(paste(gene_list_dir,"akt_75_gene_list/adapB_single/output.rda", sep="/"))
akt_75_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"bad_200_gene_list/adapB_single/output.rda", sep="/"))
bad_200_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"igf1r_75_gene_list/adapB_single/output.rda", sep="/"))
igf1r_75_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"erk_250_gene_list/adapB_single/output.rda", sep="/"))
erk_250_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"her2_15_gene_list/adapB_single/output.rda", sep="/"))
her2_15_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"egfr_25_gene_list/adapB_single/output.rda", sep="/"))
egfr_25_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"raf_100_gene_list/adapB_single/output.rda", sep="/"))
raf_100_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"krasqh_300_gene_list/adapB_single/output.rda", sep="/"))
krasqh_300_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"krasgv_300_gene_list/adapB_single/output.rda", sep="/"))
krasgv_300_genelist<-output.data$processed.data$diffGeneList
load(paste(gene_list_dir,"kraswt_300_gene_list/adapB_single/output.rda", sep="/"))
kraswt_300_genelist<-output.data$processed.data$diffGeneList

##geting training data
source("Key_ASSIGN_functions.Rmd")
expr<-as.matrix(read.table("./Trainning_data/GFP18_AKT_BAD_HER2_IGF1R_RAF_ERK.tpmlog",sep='\t',row.names=1, header=1))
control<-subset(expr, select=GFP.1:GFP.12)
her2<-subset(expr, select=HER2.1:HER2.6)
akt<-subset(expr,select=AKT.1:AKT.6)
bad<-subset(expr,select=BAD.1:BAD.6)
igf1r<-subset(expr,select=IGF1R.1:IGF1R.6)
raf<-subset(expr,select=RAF.1:RAF.6)
erk<-subset(expr,select=ERK.1:ERK.6)
expr_all<-cbind(control,akt,bad,her2,igf1r,raf, erk)
dim(expr_all)
#expr_all_f <-expr_all[apply(expr_all[,1:47]==0,1,mean) < 0.85,]
control_egfr_l<-read.table("./Trainning_data/18_GFP_EGFR_TPMlog2.txt", sep='\t', header=1, row.names=1)
gfp_egfr_multi_f <- merge_drop(expr_all, control_egfr_l)
dim(gfp_egfr_multi_f)
gfp_kras<-read.table("./Trainning_data/36_GFP_KRAS_TPMlog2.txt", sep='\t', header=1, row.names=1)
#head(gfp_kras)
gfp_egfr_kras_multi_f<-merge_drop(gfp_egfr_multi_f, gfp_kras)
dim(gfp_egfr_kras_multi_f)
#test<-data.frame(fread(testFile), check.names=F,row.names=1)
#dim(test)
#expr_all_test_f<-merge_drop(gfp_egfr_multi_f,test,by=0)
trainningData1 <- gfp_egfr_kras_multi_f
trainningLabel1 <- list(control = list(AKT = 1:12, BAD = 1:12, HER2 = 1:12, IGF1R = 1:12, RAF = 1:12, ERK = 1:12, EGFR = 48:53, KRASWT = 60:68, KRASQH = 60:68, KRASGV = 60:68), AKT = 13:18, BAD = 19:24, HER2 = 25:29, IGF1R = 30:35, RAF = 36:41, ERK = 42:47, EGFR1 = 54:59, KRASWT = 69:77, KRASQH = 78:86, KRASGV = 87:95)

##read ABCB1 data
x <- read.table("expression.gene.tpm", sep="\t", header=T)
rownames(x) <- x[,1]
x <- x[,4:ncol(x)]
x[1:3,1:4]
x_label <- c("uninf_ctrl", "GFP_50", "GFP_50", "GFP_50", "ABCB1_5", "ABCB1_5", "ABCB1_5", "ABCB1_5", "ABCB1_50", "ABCB1_50", "ABCB1_50", "uninf_ctrl", "ABCB1_50", "uninf_ctrl", "uninf_ctrl", "GFP_5", "GFP_5", "GFP_5", "GFP_5", "GFP_50")

testData1  <- x
testLabel1 <- x_label
geneList1  <- list(AKT = akt_75_genelist[[1]], BAD = bad_200_genelist[[1]], IGF1R = igf1r_75_genelist[[1]], ERK = erk_250_genelist[[1]], HER2 = her2_15_genelist[[1]], EGFR = egfr_25_genelist[[1]], RAF = raf_100_genelist[[1]], KRASQH = krasqh_300_genelist[[1]], KRASGV = krasgv_300_genelist[[1]], KRASWT = kraswt_300_genelist[[1]])

#pathwayList <- c("AKT", "BAD", "IGF1R", "ERK", "HER2", "EGFR", "RAF", "KRASQH", "KRASGV", "KRASWT")

if ( TRUE ){
# training dataset is available;
# the gene list of pathway signature is available
processed.data <- assign.preprocess(trainingData=trainningData1,
                                    testData=testData1,
                                    trainingLabel=trainningLabel1,
                                    geneList=geneList1, n_sigGene=rep(200, 10))

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
              trainingData=trainningData1, testData=testData1,
              trainingLabel=trainningLabel1,
              testLabel=testLabel1, geneList=geneList1,
              adaptive_B=TRUE, adaptive_S=FALSE,
              mixture_beta=TRUE, outputDir=tempdir)

}

