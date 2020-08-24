library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(PharmacoGx)
#data("CCLEsmall")
#data("GDSCsmall")

data.name = "CCLE"
data.dir = "D:/rws/DATA/Data/"

data.path = paste(data.dir,data.name,".RData",sep="")


#load(file = "D:/rws/CTRP/CTRPv2.RData")
#load(file = "D:/rws/CTRP/CCLE.RData")
DATA<-load(file = data.path)
DATA <- get(DATA)

# commonGenes <- intersect(fNames(GDSCsmall, "rna"),
#                          fNames(CCLEsmall,"rna"))
# common <- intersectPSet(list('CCLE'=CCLEsmall,
#                              'GDSC'=GDSCsmall),
#                         intersectOn=c("cell.lines", "drugs"),
#                         strictIntersect=TRUE)

# The cellid - drug matrix (values inside are sensitivities)
DATA.auc <- summarizeSensitivityProfiles(
  DATA,
  sensitivity.measure='auc_published',
  summary.stat="median",
  verbose=FALSE)

## Example for other datasets
# CCLE.auc <- summarizeSensitivityProfiles(
#   CCLE,
#   sensitivity.measure='auc_published',
#   summary.stat="median",
#   verbose=FALSE)
# CTRP.auc <- summarizeSensitivityProfiles(
#   CTRPv2,
#   sensitivity.measure='auc_published',
#   summary.stat="median",
#   verbose=FALSE)
# GDSC.ic50 <- summarizeSensitivityProfiles(
#   common$GDSC,
#   sensitivity.measure='ic50_published',
#   summary.stat="median",
#   verbose=FALSE)
# CCLE.ic50 <- summarizeSensitivityProfiles(
#   common$CCLE,
#   sensitivity.measure='ic50_published',
#   summary.stat="median",
#   verbose=FALSE)

# GDSCsummarized <- summarizeMolecularProfiles(GDSCsmall,
#                                         mDataType = "rna", cell.lines=cellNames(GDSCsmall),
#                                         summary.stat = 'median', fill.missing = TRUE, verbose=TRUE)


# Print dimensions of the metadata and the expression data
print(dim(DATA@molecularProfiles$rna@phenoData@data))
print(dim(DATA@molecularProfiles$rna@assayData$exprs))


# Extract sample name and cell id relation
DATA.samp2ID <- DATA@molecularProfiles$rna@phenoData@data[,c("samplename","cellid")]
print(head.DataTable(DATA.samp2ID))

# The cell - gene matrix, values inside are expression values
DATA.expression.transpose <- t(DATA@molecularProfiles$rna@assayData$exprs)
DATA.samp2ID$samplename<-toupper(DATA.samp2ID$samplename)


# Check if the sample name order in the expression match the order of the metadata
print(head.DataTable(rownames(DATA.expression.transpose)== DATA.samp2ID))
print(sum(rownames(DATA.expression.transpose)== DATA.samp2ID)==dim(DATA.expression.transpose)[1])

# Because I know their order are the same, I just replace the rowname with cell id
if(data.name %in% c("GDSC","CCLE")){
  rownames(DATA.expression.transpose) <- DATA.samp2ID$cellid
  print(head.DataTable(DATA.expression.transpose))
  
}

DATA.auc.T<- t(DATA.auc)
# Todo: please check if cell id in  and GDSC.expression.transpose match each other 
# If yes, megrge two according to their id 
# If no, find and method to calcute string simiarities to match them
DATA.label <- DATA.auc.T[rownames(DATA.expression.transpose),]


# Todo: please check if cell id in GDSC.auc and GDSC.expression.transpose match each other 
# If yes, megrge two according to their id 
# If no, find and method to calcute string simiarities to match them

write.csv(DATA.expression.transpose,file = "DATA/CCLEexpression.csv")

write.csv(DATA.label,file = "DATA/CCLElabel.csv")





