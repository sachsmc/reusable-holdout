library(survival)

load("../Data/jbl10-data-2016-02-04.RData")

## develop signature in 62 OBS patient subset

psdex <- 2:(grep("tissue", colnames(gdat)) - 1)
colnames(gdat)[psdex] <- paste0("X", colnames(gdat)[psdex])

table(gdat$Post.Surgical.Treatment)
obsgrp <- subset(gdat, Post.Surgical.Treatment == "OBS")

## identify subset of 172 probes with P < .005 in cox model with outcome
# results from Zhu, supplementary table 6S
if(FALSE){
zhu <- readLines("../Data/sig172-probes-zhu.txt")
header <- zhu[1:8]
data0 <- zhu[-c(1:8)]
data1 <- unlist(strsplit(data0, split = " ", fixed = TRUE))
data2 <- as.data.frame(matrix(data1, ncol = 8, byrow = TRUE), stringsAsFactors = FALSE)
data2[, 5:8] <- lapply(data2[, 5:8], as.numeric)
colnames(data2) <- header
}
## can i reproduce it?

obsgrp[, psdex] <- lapply(obsgrp[, psdex], function(x) (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE))

ps172dex <- NULL
for(i in psdex){

  form <- as.formula(paste("Surv(OS.time, OS.status == 'Dead') ~ ", colnames(obsgrp)[i]))
  fit <- tryCatch(coxph(form, data = obsgrp, y = FALSE, ties = "breslow", robust = FALSE), error = function(e) NA)
  if(is.na(fit[[1]])) next
  fs <- summary(fit)$coefficients
  ps172dex <- rbind(ps172dex, data.frame(probe = colnames(obsgrp)[i], index = i,
                                p.value = fs[5], hr = fs[2], lhr = exp(fs[1] - 1.96*fs[3]),
                                uhr = exp(fs[1] + 1.96*fs[3]), stringsAsFactors = FALSE))

}


ps172dex$probe2 <- ifelse(substr(ps172dex$probe, 1, 1) == "X",
                         substr(ps172dex$probe, 2, nchar(ps172dex$probe)), ps172dex$probe)


###

table(ps172dex$p.value < .005)

matches <- subset(ps172dex, p.value < .005)
rownames(matches) <- matches$probe
## now do the forwards and backward selection
subgrp <- obsgrp[, c("ID", paste(matches$probe), "OS.time", "OS.status")]




