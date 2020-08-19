library(Biobase)
library(SummarizedExperiment)
library(S4Vectors)
library(PharmacoGx)
#data("CCLEsmall")
#data("GDSCsmall")

## TODO cross check the CCLE and L1000 CTRPV2s

data.name = "L1000_compounds"
data.dir = "Data/"

data.path = paste(data.dir,data.name,".RData",sep="")


#load(file = "D:/rws/CTRP/CTRPv2.RData")
#load(file = "D:/rws/CTRP/CCLE.RData")
DATA<-load(file = data.path)
DATA <- get(DATA)

CTRPV2<-load(file = "Data/CTRPv2.RData")
CTRPV2<- get(CTRPV2)

CCLE<-load(file = "Data/CCLE.RData")
CCLE<- get(CCLE)


# commonGenes <- intersect(fNames(GDSCsmall, "rna"),
#                          fNames(CCLEsmall,"rna"))
# common <- intersectPSet(list('CCLE'=CCLEsmall,
#                              'GDSC'=GDSCsmall),
#                         intersectOn=c("cell.lines", "drugs"),
#                         strictIntersect=TRUE)

# The cellid - drug matrix (values inside are sensitivities)
CTRPV2.auc <- summarizeSensitivityProfiles(
  CTRPV2,
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


DATA.meta<-DATA@molecularProfiles$rna@phenoData@data


# Extract sample name and cell id relation
DATA.samp2ID <- DATA@molecularProfiles$rna@phenoData@data[,c("distil_id","cellid","xptype")]
print(head.DataTable(DATA.samp2ID))

DATA.samp2ID.control <- DATA.samp2ID[DATA.samp2ID$xptype=="control",]
DATA.control.index <- DATA.samp2ID$xptype=="control"
print(head.DataTable(DATA.samp2ID.control))

# The cell - gene matrix, values inside are expression values
DATA.control.expression.transpose <- t(DATA@molecularProfiles$rna@assayData$exprs)[DATA.control.index,]

# Check if the sample name order in the expression match the order of the metadata
print(head.DataTable(rownames(DATA.control.expression.transpose)== DATA.samp2ID.control))
print(sum(rownames(DATA.control.expression.transpose)== DATA.samp2ID.control)==dim(DATA.control.expression.transpose)[1])

# Because I know their order are the same, I just replace the rowname with cell id
if(data.name == "GDSC"){
  rownames(DATA.expression.transpose) <- DATA.samp2ID$cellid
  print(head.DataTable(DATA.expression.transpose))
  
}else if(data.name = "L1000_compounds"){
  rownames(DATA.control.expression.transpose) <- DATA.samp2ID.control$cellid
  print(head.DataTable(DATA.control.expression.transpose))
  
}



# Todo: please check if cell id in GDSC.auc and GDSC.expression.transpose match each other 
# If yes, megrge two according to their id 
# If no, find and method to calcute string simiarities to match them


rownames(DATA.control.expression.transpose)

CTRPV2.auc.T <- t(CTRPV2.auc)
rownames(CTRPV2.auc.T)<-toupper(rownames(CTRPV2.auc.T))

badchars <- "[\xb5]|[\n]|[,]|[;]|[:]|[-]|[+]|[*]|[%]|[$]|[#]|[{]|[}]|[[]|[]]|[|]|[\\^]|[/]|[\\]|[.]|[(]|[)]|[_]|[ ]"
rownames(CTRPV2.auc.T)<-gsub(badchars, "", rownames(CTRPV2.auc.T))

save(DATA.control.expression.transpose,file="DATA/L1000train.RData")


DF.temp <- DATA.control.expression.transpose[,1]

DATA.result <- merge(DF.temp,CTRPV2.auc.T,all.x=TRUE)


save(DATA.result,file="DATA/L1000target.RData")