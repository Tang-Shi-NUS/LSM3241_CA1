if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("GEOquery", version = "3.8")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("affy", version = "3.8")
libif (!requireNamespace("BiocManager", quietly = TRgetUE))
  install.packages("BiocManager")
BiocManager::install("limma", version = "3.8")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("hgu133plus2.db", version = "3.8")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("org.Hs.eg.db", verinsion = "3.8")
gse <- GEOquery::getGEO(filename='C:/Users/tangshi/Documents/R/LSM3241_CA1/Data/GSE50697_family.soft.gz')
names(GEOquery::GSMList(gse))
list.celfiles('data') #generates a vector containing the names of all the CEL files in the named directory
data <- read.affybatch(paste0('data/',list.celfiles('data')))
image(data[,1])

gsm <- GSMList(gse)[[1]]
gsm50697 <- GSMList(gse50697)[[1]]
names(Meta(gsm))[!(names(Meta(gsm)) %in% names(Meta(gse)))]
for (gsm in GSMList(gse)) { 
  print(Meta(gsm)[['characteristics_ch1']])
}pd
treatment_type <- function(gsm) { 
  Meta(gsm)[['characteristics_ch1']][2]
}
sapply(GSMList(gse),treatment_type)
pd_treatment <- data.frame(treatment=as.factor(sapply(GSMList(gse),treatment_type)))
levels(pd_treatment$treatment)
pd_treatment$treatment <- as.factor(pd_treatment$treatment)
levels(pd_treatment$treatment) <- c("control", "miR203")
trcelfiles <- paste0('data/',rownames(pd_treatment),'.CEL.gz')
affydata_treat <- read.affybatch(celfiles,phenoData = new("AnnotatedDataFrame",pd_treatment))
phenoData(affydata_treat)
eset_treat <- rma(affydata_treat)
head(exprs(eset_treat))
plotDensity(exprs(eset_treat),xlab='log intensity',main=" GSE50697 Feature Level Densities After RMA",lwd=2)
pData(eset_treat)

eset_treat_no_bg <- rma(affydata_treat, subset=NULL, verbose=TRUE, destructive=TRUE, normalize=TRUE, background=FALSE, bgversion=2)
plotDensity(exprs(eset_treat_no_bg),xlab='log intensity',main=" GSE50697 Feature Level Densities w/o Background Correction",lwd=2)
eset_treat_no_norm <- rma(affydata_treat, subset=NULL, verbose=TRUE, destructive=TRUE, normalize=FALSE, background=TRUE, bgversion=2)
plotDensity(exprs(eset_treat_no_norm),xlab='log intensity',main=" GSE50697 Feature Level Densities w/o Normalization",lwd=2)
plotDensity(exprs(affydata_treat),xlab='log intensity',main=" GSE50697 Feature Level Densities w/o RMA",lwd=2)

featureData(eset_treat)
annotation(eset_treat)

treatment_model <- model.matrix( ~ 0 + eset_treat$treatment)
treatment_model
colnames(treatment_model) <- levels(eset_treat$treatment)
treatment_contrasts <- makeContrasts(miR203 - control, levels=treatment_model)
treatment_contrasts
fit <- lmFit(eset_treat, treatment_model)
fit
fitted.contrast <- contrasts.fit(fit,treatment_contrasts)
fitted.contrast 
fitted.ebayes <- eBayes(fitted.contrast)
topTable(fitted.ebayes) #gives a table of genes that BEST FIT the linear model; in this case, upregulated
ps_upreg <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1) #retrieves all upregulated genes with adjusted p < 0.05
ps_upreg
ps_upreg_log <-rownames(ps_upreg[ps_upreg$logFC > 0,]) #log
ps_upreg_log
library(hgu133plus2.db)
columns(hgu133plus2.db)
head(keys(hgu133plus2.db,keytype="PROBEID"))
dysregulated_genes_csv <- AnnotationDbi::select(hgu133plus2.db,rownames(upregulated_genes),c("SYMBOL","ENTREZID","GENENAME","PATH"),keytype="PROBEID")
# upregulated_genes_csv <- AnnotationDbi::select(hgu133plus2.db,ps_upreg_log,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
# Generates a dataframe of all the upregulated genes that fit our prior criteria 
write.csv(dysregulated_genes_csv, file = "dysregulated_genes.csv") #for easy viewing

upregulated_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=2) #set a cutoff and identify the genes that pass them;
# select for genes that will likely exhibit STATISTICALLY SIGNIFICANT upregulation -> stringent
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(upregulated_genes)))
points(upregulated_genes[['logFC']],-log10(upregulated_genes[['P.Value']]),col='blue')
abline(v =c(-2, 2), col="red", lwd = 2, lty = 2)
abline(h = -log10(0.05), col = "red", lwd=2, lty = 2)
#statistically significant (true positive) up&downregulation is highlighted; (p > 0.05?)
eset_upregulated_genes <-eset_treat[rownames(upregulated_genes),]
heatmap(exprs(eset_upregulated_genes),
        # labCol=eset_treat$treatment,labRow=NA, #PROBEID will not be shown; comment this line for PROBEID
        col       = rev(brewer.pal(10, "RdBu")), #better for red-green colourblindness
        distfun   = function(x) as.dist(1-cor(t(x)))) #arrange heatmap by correlation

dysregulated_genes_lfc1.5 <- topTable(fitted.ebayes,number=Inf,p.value = 0.01,lfc=1.5)
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(dysregulated_genes_lfc1.5)))
points(dysregulated_genes_lfc1.5[['logFC']], -log10(dysregulated_genes_lfc1.5[['P.Value']]), col = 'blue')
abline(v =c(-1.5, 1.5), col="red", lwd = 2, lty = 2)
abline(h = -log10(0.05), col = "red", lwd=2, lty = 2)
dysregulated_genes_lfc1.5_csv <- AnnotationDbi::select(hgu133plus2.db,rownames(dysregulated_genes_lfc1.5),c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(dysregulated_genes_lfc1.5_csv, file = "dysregulated_genes_lfc1.5.csv")
write.csv(dysregulated_genes_lfc1.5, file = "dysregulated_genes_lfc1.5_probs.csv")

is_dkk1 <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=1)
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(is_dkk1)))
points(is_dkk1[['logFC']], -log10(is_dkk1[['P.Value']]), col = 'blue')
abline(v =c(-1, 1), col="red", lwd = 2, lty = 2)
abline(h = -log10(0.05), col = "red", lwd=2, lty = 2)

####IGNORE EVERYTHING AFTER THIS LINE ######
#for downregulated genes after miR203 treatment: (not sure how to get)
ps_downreg <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc=-4)
ps_downreg_log <- rownames(ps_downreg[ps_downreg$logFC < 0,]) 
# log fold change is positive, expression is greater in miR203 than control, and VICE VERSA
# log fold change (LFC) is NEGATIVE (LFC < 0), expr is greater in control than miR203; miR203 downregulates these genes!
head(keys(hgu133plus2.db,keytype="PROBEID"))
downregulated_genes_csv <- AnnotationDbi::select(hgu133plus2.db,ps_downreg_log,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
downregulated_genes <- topTable(fitted.ebayes,number=Inf,p.value = 0.05,lfc = -2)
eset_downregulated_genes <- eset_treat[rownames(downregulated_genes),]
write.csv(downregulated_genes_csv, file = "downregulated_genes.csv")

volcanoplot(fitted.ebayes, main=sprintf("%d features pass our cutoffs",nrow(downregulated_genes)))
points(downregulated_genes[['logFC']],-log10(downregulated_genes[['P.Value']]),col='red')
heatmap(exprs(eset_downregulated_genes),
        labCol=eset_treat$treatment,labRow=NA, #PROBEID will not be shown; comment this line for PROBEID
        col       = rev(brewer.pal(10, "RdBu")), #better for red-green colourblindness
        distfun   = function(x) as.dist(1-cor(t(x)))) #arrange heatmap by correlation
# very stringent isolation; p = 0.0001; 0.01% chance that it's a false positive
ps_downreg_strin <- topTable(fitted.ebayes,number=Inf,p.value = 0.0001,lfc=-2)
ps_downreg_strin_log <- rownames(ps_downreg_strin[ps_downreg_strin$logFC < 0,])
downregulated_genes_strin_csv <- AnnotationDbi::select(hgu133plus2.db,ps_downreg_strin_log,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
write.csv(downregulated_genes_strin_csv, file = "downregulated_genes_strin.csv")

downregulated_genes_strin <- topTable(fitted.ebayes,number=Inf,p.value = 0.0001,lfc= -4)
upregulated_genes_strin <- topTable(fitted.ebayes, number=Inf, p.value = 0.0001, lfc = 2)
# ps_upreg_strin_log <- rownames(ps_upreg_strin[ps_upreg_strin$logFC < 0,])
# upregulated_genes_strin_csv <- AnnotationDbi::select(hgu133plus2.db,ps_upreg_strin_log,c("SYMBOL","ENTREZID","GENENAME"),keytype="PROBEID")
eset_downregulated_genes_strin <- eset_treat[rownames(downregulated_genes_strin),]
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our stringent cutoffs",nrow(downregulated_genes_strin)))
points(downregulated_genes_strin[['logFC']],-log10(downregulated_genes_strin[['P.Value']]),col='red')
eset_upregulated_genes_strin <- eset_treat[rownames(upregulated_genes_strin),]
volcanoplot(fitted.ebayes, main=sprintf("%d features pass our stringent cutoffs",nrow(upregulated_genes_strin)))
points(upregulated_genes_strin[['logFC']],-log10(upregulated_genes_strin[['P.Value']]),col='blue')
heatmap(exprs(eset_downregulated_genes_strin),
        # labCol=eset_treat$treatment,labRow=NA, #PROBEID will not be shown; comment this line for PROBEID
        col       = rev(brewer.pal(10, "RdBu")), #better for red-green colourblindness
        distfun   = function(x) as.dist(1-cor(t(x)))) #arrange heatmap by correlation