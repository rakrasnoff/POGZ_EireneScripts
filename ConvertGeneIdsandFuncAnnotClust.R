######### Install required packages ###########

##Check if required packages are installed; if not, install them. If more packages become necessary, add them to the package list
packagelist <- c('gplots', 'ggplot2', 'gridExtra', 'WGCNA', 
                 'impute', 'GO.db', 'AnnotationDbi', 'org.Hs.eg.db', 
                 'annotate', 'compare')

for (i in 1:length(packagelist)) { 
  if(packagelist[i] %in% rownames(installed.packages()) == FALSE) 
  {install.packages(packagelist[i])
    cat("
        installing", packagelist[i])}
  if (packagelist[i] %in% rownames(installed.packages()) == TRUE) 
  {cat(packagelist[i], "is already installed
       ")}
  }

#####
source("http://bioconductor.org/biocLite.R")
biocLitepackagelist <- c("AnnotationDbi", 'org.Hs.eg.db', 'org.Mm.eg.db', 'annotate', 'cummeRbund')

for (i in 1:length(biocLitepackagelist)) { 
  if(biocLitepackagelist[i] %in% rownames(installed.packages()) == FALSE) 
  {biocLite(biocLitepackagelist[i])
    cat("
        installing", biocLitepackagelist[i]) }
  if (biocLitepackagelist[i] %in% rownames(installed.packages()) == TRUE) 
  {cat(biocLitepackagelist[i], "is already installed
       ")}
  }

########## Load required packages ##########

library(gplots)
library(ggplot2)
library(gridExtra)
library(WGCNA)
library(impute)
library(GO.db)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(annotate)
library(compare)
library(cummeRbund)


####### Set working directory
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cuffdiff_out") ##write path to file here


###### Load dataset
gene.exp.full <- read.table("gene_exp.diff", header=TRUE)
head(gene.exp.full)

###### Create gene list to use: remove genes without names, expression value of less than one for both sample 1 and 2, with status !OK
gene.exp <- gene.exp.full[
  (gene.exp.full$value_1 >= 1 | gene.exp.full$value_2 >= 1)
  & gene.exp.full$status == "OK"
  & gene.exp.full$gene != "-", ]
nrow(gene.exp)
nrow(gene.exp.full)

##### Convert symbols to entrez ids
keytypes(org.Mm.eg.db)
columns(org.Mm.eg.db)
head(keys(org.Mm.eg.db, keytype="SYMBOL"))
head(keys(org.Mm.eg.db, keytype="ENTREZID"))
head(select(org.Mm.eg.db, as.vector(gene.exp$gene), "ENTREZID", "SYMBOL"),300)
## Append column with entrez id to dataset
convertedids<- select(org.Mm.eg.db, as.vector(gene.exp$gene), "ENTREZID", "SYMBOL")
gene.exp$entrez_id <- convertedids$ENTREZID[match(gene.exp$gene, convertedids$SYMBOL)]

##### Define foreground ("gene list") and background genes
#Parameters for FG genes: difference in expression was >+2x, q value was less than .05,
#Note: if not done above, also filter out genes without names, with expression values less than one both sample 1 and 2, and with status !OK
gene.FG <- gene.exp[
  (gene.exp$log2.fold_change. >= 1 | gene.exp$log2.fold_change. <= -1) 
  & gene.exp$q_value < .05, ]
head(gene.FG)
#BG genes: using whole list of genes with names, expression values >1, status OK
gene.BG <- gene.exp

#### Gene lists are ready for functional annotation
### FUNCTIONAL ANNOTATION: UP-REGULATED AND DOWN-REGULATED GENES 

dyn.load('/Library/Java/JavaVirtualMachines/jdk1.8.0_91.jdk/Contents/Home/jre/lib/server/libjvm.dylib')
require(rJava)
.jinit()
.jcall("java/lang/System", "S", "getProperty", "java.runtime.version")
library("RDAVIDWebService")

###Connect to webservice
# Create a DAVIDWebService object connected to David, using your registration email.
# To register, go to: http://david.abcc.ncifcrf.gov/content.jsp?file=WS.html.
david <- DAVIDWebService(email="willseylab@ucsf.edu", url="https://david.ncifcrf.gov/webservice/services/DAVIDWebService.DAVIDWebServiceHttpSoap12Endpoint/")


####RUN UPREGULATED AND DOWN REGULATED GENES THROUGH DAVID SIMULTANEOUSLY
# Define foreground and background gene lists.
# The foreground list should be contained within the background list.
# note: for possible inputs for idType, use getIdTypes(david)
FG <- addList(david, gene.FG$entrez_id, idType="ENTREZ_GENE_ID", listName="foreground", listType="Gene")
BG <- addList(david, gene.BG$entrez_id, idType="ENTREZ_GENE_ID", listName="background", listType="Background")

# Inspect FG and BG to see the proportion of genes recognized by David, and those that are unmapped.
FG 
BG

# Inspect "david" object to see the gene lists selected as foreground and background.
david

# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart <- getFunctionalAnnotationChart(david)
View(FuncAnnotChart)
# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
head(summary(FuncAnnotClust))
dev.new()
plot2D(FuncAnnotClust, 1)

# Print functional annotation clustering to file (limited to 3000 genes).
getClusterReportFile(david, "FuncAnnotClust.tsv")

###############
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[1]], pvalueCutoff=0.05, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")
dev.new()
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[2]], pvalueCutoff=0.1, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")


####RUN ONLY UP REGULATED GENES THROUGH DAVID 
gene.u <- gene.FG[
  (gene.FG$log2.fold_change. >= 1) , ]
nrow(gene.u)
FGu <- addList(david, gene.u$entrez_id, idType="ENTREZ_GENE_ID", listName="upregulated", listType="Gene")
BG <- addList(david, gene.BG$entrez_id, idType="ENTREZ_GENE_ID", listName="background", listType="Background") #this line is only here because it wouldn't let me set the same backgroun as above
david

# Specifiy annotation categories.
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart.up <- getFunctionalAnnotationChart(david)
View(FuncAnnotChart.up)
nterms <- .1*(nrow(FuncAnnotChart.up))
FuncAnnotGraph.up <- data.frame(FuncAnnotChart.up$Term[1:nterms])
FuncAnnotGraph.up[,2] <- data.frame(FuncAnnotChart.up$PValue[1:nterms])
#PLOTTTTTTT

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart_up.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
head(summary(FuncAnnotClust))
dev.new()
plot2D(FuncAnnotClust, 1)
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[1]], pvalueCutoff=0.05, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")


####RUN ONLY DOWN REGULATED GENES THROUGH DAVID 
gene.d <- gene.FG[
  (gene.FG$log2.fold_change. <= (-1)) , ]
nrow(gene.d)
FGd <- addList(david, gene.d$entrez_id, idType="ENTREZ_GENE_ID", listName="downregulated", listType="Gene")
BG <- addList(david, gene.BG$entrez_id, idType="ENTREZ_GENE_ID", listName="backgroundu", listType="Background")
david


#set categories
setAnnotationCategories(david, c("GOTERM_BP_ALL", "GOTERM_MF_ALL", "GOTERM_CC_ALL"))

# Get functional annotation chart as R object.
FuncAnnotChart.d <- getFunctionalAnnotationChart(david)
View(FuncAnnotChart.d)
nterms <- .1*(nrow(FuncAnnotChart.d))
FuncAnnotGraph.d <- data.frame(FuncAnnotChart.d$Term[1:nterms])
FuncAnnotGraph.d[,2] <- data.frame(FuncAnnotChart.d$PValue[1:nterms])
##pick up here

# Print functional annotation chart to file.
getFunctionalAnnotationChartFile(david, "FuncAnnotChart_d.tsv")

# Get functional annotation clustering (limited to 3000 genes).
FuncAnnotClust <- getClusterReport(david)
head(summary(FuncAnnotClust))
dev.new()
plot2D(FuncAnnotClust, 1)
davidGODag<-DAVIDGODag(members(FuncAnnotClust)[[1]], pvalueCutoff=0.05, "CC")
plotGOTermGraph(g=goDag(davidGODag),r=davidGODag, max.nchar=40, node.shape="ellipse")


###Create file for GSEA analysis
##this is set up for 3 samples of from each group
#doesn't give me the name of the genes....
#HETfpkm1 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cufflinks_out/pogzP2-het1_ACTTGA_L001_001/genes.fpkm_tracking", header=TRUE)
#HETfpkm2 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cufflinks_out/pogzP2-het2_GATCAG_L001_001/genes.fpkm_tracking", header=TRUE)
#HETfpkm3 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cufflinks_out/pogzP2-het3_GGCTAC_L001_001/genes.fpkm_tracking", header=TRUE)
#WTfpkm1 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cufflinks_out/pogzP2-wt1_ATCACG_L001_001/genes.fpkm_tracking", header=TRUE)
#WTfpkm2 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cufflinks_out/pogzP2-wt2_TTAGGC_L001_001/genes.fpkm_tracking", header=TRUE)
#WTfpkm3 <- read.table("/Users/rebeccakrasnoff/Documents/Current/Willsey/POGZProject/data/cufflinks_out/pogzP2-wt3_TAGCTT_L001_001/genes.fpkm_tracking", header=TRUE)
###Was seeing what happened to error when I removed NAs
#head(gene.exp)
#gene.exp.nona <- gene.exp[gene.exp$entrez_id != "NA",]
#NAME = gene.exp.nona$entrez_id #give column with gene name
#Description <- gene.exp.nona$gene #give column with description of gene or gene symbol
#WT <- gene.exp.nona$value_1 #name for sample 1
#HET <- gene.exp.nona$value_2 #name for sample 2
#Create table
#df.glite <- data.frame(NAME, WTfpkm1, WTfpkm2, WTfpkm3, HETfpkm1, HETfpkm2, HETfpkm3)


##CURRENTLY SET UP LIKE THIS JUST TO EXPLORE ERROR
NAME = gene.exp$gene #give column with gene name
Description <- gene.exp$gene #give column with description of gene or gene symbol
WT <- gene.exp$value_1 #name for sample 1
HET <- gene.exp$value_2 #name for sample 2
df.gsea.fake <- data.frame(NAME, Description, 
                      WTfpkm1, WTfpkm2, WTfpkm3, HETfpkm1, HETfpkm2, HETfpkm3)
n <- nrow(df.gsea)
cl <- length(df.gsea)
n.samples <- 6 #number of sample groups
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/POGZ_dataset")


######create fake gct withougt 0s
gene.exp.nonzero <- gene.exp[gene.exp$value_1 != 0 & gene.exp$value_2 != 0,]
NAME = gene.exp.nonzero$gene #give column with gene name
Description <- gene.exp.nonzero$entrez_id #give column with description of gene or gene symbol
WT <- gene.exp.nonzero$value_1 #name for sample 1
HET <- gene.exp.nonzero$value_2 #name for sample 2
df.gsea <- data.frame(NAME, Description, 
                      WT, WT, WT, HET, HET, HET)
n <- nrow(df.gsea)
cl <- length(df.gsea)
n.samples <- 6 #number of sample groups
setwd("/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/POGZ_dataset")

write.table("#1.2", file="POGZ_GSEA_fake.gct", sep=",", row.names=FALSE, col.names=FALSE)
write.table(t((c(n, n.samples))), file="POGZ_GSEA_fake.gct", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(df.gsea, file="POGZ_GSEA_fake.gct", sep=",", row.names=FALSE, col.names=TRUE, append=TRUE)



############ original/good code
write.table("#1.2", file="POGZ_GSEA.gct", sep=",", row.names=FALSE, col.names=FALSE)
write.table(t((c(n, n.samples))), file="POGZ_GSEA.gct", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(df.gsea, file="POGZ_GSEA.gct", sep=",", row.names=FALSE, col.names=TRUE, append=TRUE)

#write as csv so I can look at it
write.table("#1.2", file="POGZ_GSEA.csv", sep=",", row.names=FALSE, col.names=FALSE)
write.table(t((c(n, n.samples))), file="POGZ_GSEA.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
write.table(df.gsea, file="POGZ_GSEA.csv", sep=",", row.names=FALSE, col.names=TRUE, append=TRUE)


#####create cls file
write.table(t((c(n.samples, 2, 1))), file="POGZ_GSEA.cls", sep=",", row.names=FALSE, col.names=FALSE)
sample1.name <- "WT"
sample2.name <- "HET"
write.table(t((c("#", sample1.name, sample2.name))), file="POGZ_GSEA.cls", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
sample1.list <- rep("sample_1", (n.samples/2))
sample2.list <- rep("sample_2", (n.samples/2))
write.table(t((c(sample1.list, sample2.list))), file="POGZ_GSEA.cls", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)

#write as csv so I can look at it
write.table(t((c(n.samples, 2, 1))), file="POGZ_GSEA.csv", sep=",", row.names=FALSE, col.names=FALSE)
sample1.name <- "WT"
sample2.name <- "HET"
write.table(t((c("#", sample1.name, sample2.name))), file="POGZ_GSEA.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)
sample1.list <- rep("sample_1", (n.samples/2))
sample2.list <- rep("sample_2", (n.samples/2))
write.table(t((c(sample1.list, sample2.list))), file="POGZ_GSEA.csv", sep=",", row.names=FALSE, col.names=FALSE, append=TRUE)


###################################GSEA
# GSEA 1.0 -- Gene Set Enrichment Analysis / Broad Institute 
#
# R script to run GSEA Analysis of the Leukemia ALL/AML vs C1 example (cut and paste into R console)

GSEA.program.location <- "/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/GSEA.1.0.R"   #  R source program (change pathname to the rigth location in local machine)
source(GSEA.program.location, echo=T, verbose=T, max.deparse.length=9999)


GSEA(# Input/Output Files :-------------------------------------------
     input.ds =  "/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/POGZ_dataset/GSEA_dummy.gct",        # Input gene expression Affy dataset file in RES or GCT format
     input.cls = "/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/Datasets/GSEA_dummy.cls",          # Input class vector (phenotype) file in CLS format
     gs.db =     "/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/GeneSetDatabases/C1.gmt",        # Gene set database in GMT format
     output.directory      = "/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/POGZ_dataset",        # Directory where to store output and results (default: "")
     
     #  Program parameters :-------------------------------------------------------------------------------------------------------------------------
     doc.string            = "SAMPLE",   # Documentation string used as a prefix to name result files (default: "GSEA.analysis")
     non.interactive.run   = F,               # Run in interactive (i.e. R GUI) or batch (R command line) mode (default: F)
     reshuffling.type      = "gene.labels", # Type of permutation reshuffling: "sample.labels" or "gene.labels" (default: "sample.labels" 
     nperm                 = 1000,            # Number of random permutations (default: 1000)
     weighted.score.type   =  1,              # Enrichment correlation-based weighting: 0=no weight (KS), 1= weigthed, 2 = over-weigthed (default: 1)
     nom.p.val.threshold   = -1,              # Significance threshold for nominal p-vals for gene sets (default: -1, no thres)
     fwer.p.val.threshold  = -1,              # Significance threshold for FWER p-vals for gene sets (default: -1, no thres)
     fdr.q.val.threshold   = .25,            # Significance threshold for FDR q-vals for gene sets (default: 0.25)
     topgs                 = 10,              # Besides those passing test, number of top scoring gene sets used for detailed reports (default: 10)
     adjust.FDR.q.val      = F,               # Adjust the FDR q-vals (default: F)
     gs.size.threshold.min = 15,              # Minimum size (in genes) for database gene sets to be considered (default: 25)
     gs.size.threshold.max = 30000,             # Maximum size (in genes) for database gene sets to be considered (default: 500)
     reverse.sign          = F,               # Reverse direction of gene list (pos. enrichment becomes negative, etc.) (default: F)
     preproc.type          = 0,               # Preproc.normalization: 0=none, 1=col(z-score)., 2=col(rank) and row(z-score)., 3=col(rank). (def: 0)
     random.seed           = 3338,            # Random number generator seed. (default: 123456)
     perm.type             = 0,               # For experts only. Permutation type: 0 = unbalanced, 1 = balanced (default: 0)
     fraction              = 1.0,             # For experts only. Subsampling fraction. Set to 1.0 (no resampling) (default: 1.0)
     replace               = F,               # For experts only, Resampling mode (replacement or not replacement) (default: F)
     save.intermediate.results = F,           # For experts only, save intermediate results (e.g. matrix of random perm. scores) (default: F)
     OLD.GSEA              = F,               # Use original (old) version of GSEA (default: F)
     use.fast.enrichment.routine = T          # Use faster routine to compute enrichment for random permutations (default: T)
)
#-----------------------------------------------------------------------------------------------------------------------------------------------

# Overlap and leading gene subset assignment analysis of the GSEA results

GSEA.Analyze.Sets(
  directory           = "/Users/rebeccakrasnoff/Documents/Current/Willsey/Scripts/GSEAPR/POGZ_dataset",        # Directory where to store output and results (default: "")
  topgs = 20,                                                           # number of top scoring gene sets used for analysis
  height = 16,
  width = 16
)



































