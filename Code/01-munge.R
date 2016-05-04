library(Biobase)
library(GEOquery)
library(rio)
library(sva)
# load series and platform data from GEO

gset <- getGEO("GSE14814", GSEMatrix =TRUE)
if (length(gset) > 1) idx <- grep("GPL96", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]


gdat <- as(gset, "data.frame")
gdat$ID <- sampleNames(phenoData(gset))

gexpr <- exprs(gset)
fdata <- fData(gset)
#fdex <- fdata$ID[fdata$ENTREZ_GENE_ID != "" & fdata$`Gene Ontology Cellular Component` != ""]
#gexpr <- gexpr[fdex, ]

anyna <- rowSums(is.na(gexpr)) > 0
gexpr <- gexpr[!anyna, ]
## 19,619 genes, close enough

parse_pheno <- function(x){

  pair <- strsplit(paste(x), ": ", fixed = TRUE)
  if(length(pair[[1]]) == 1){
    x
  } else {
    name <- pair[[1]]
    name <- name[!is.na(name)][1]

    val <- sapply(pair, "[", 2)

  ## try converting to numeric

    nval <- as.numeric(val)
    if(all(is.na(nval))) nval <- val

    out <- data.frame(nval)
    colnames(out) <- name
    out
}

}

pheno0 <- gdat[, c("ID", grep("characteristics", colnames((gdat[, (ncol(gdat) - 55):ncol(gdat)])[1:4, ]), value = TRUE))]
pheno0$characteristics_ch1.12 <- paste(pheno0$characteristics_ch1.12)
pheno0$characteristics_ch1.12[pheno0$characteristics_ch1.12 == ""] <- paste(pheno0$characteristics_ch1.11[pheno0$characteristics_ch1.12 == ""])
pheno0$characteristics_ch1.11 <- NULL

pheno <- as.data.frame(lapply(pheno0, parse_pheno))

#gdat <- merge(gdat, pheno, by = "ID")
gdat$batch <- ifelse(gdat$data_processing == "The data were preprocessed using RMA v0.5.", -1, 1)

## remove batch effect

modcom <- model.matrix( ~ 1, data = gdat)
batchadj <- ComBat(gexpr, batch = gdat$batch, mod = modcom, prior.plots = TRUE)

expradj <- as.data.frame(t(batchadj))

## pca
pca.fit <- prcomp(t(batchadj), center = TRUE, scale. = TRUE)
pca.dat <- as.data.frame(pca.fit$x[, 1:25])
pca.dat$ID <- rownames(pca.dat)
##

expradj$ID <- rownames(expradj)

gdat <- merge(expradj, pheno, by = "ID")

pca.dat2 <- merge(pca.dat, pheno, by = "ID")
## make sure no merging errors in data

pheno.xls <- import("../Data/Zhu_2010_JCO_microarray_pts_clin_info.xls", sheet = 1)

check <- merge(pheno, pheno.xls, by.x = "age", by.y = "Age")
check2 <- subset(check, OS.time == `OS time`)
# nrow(check2) == nrow(pheno)
# plot(OS.time ~ `OS time`, data = check2)
# we good

save(gdat, file = paste0("../Data/jbl10-data-", Sys.Date(), ".RData"))
save(pca.dat2, file = paste0("../Data/jbl10-pca-", Sys.Date(), ".RData"))


