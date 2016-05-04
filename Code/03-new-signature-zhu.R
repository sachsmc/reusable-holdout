library(survival)
library(tidyr)
library(dplyr)

load("../Data/jbl10-data-2016-02-09.RData")

## supervised pca fit for heldout
##    - select x most correlated genes
##    - subset those, form 3 pcas
##    - Fit multivariate model
## cutpoint at median
## predict for held in dataset

scl.cen <- function(x){

  (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)

}

## dichotomize survival
gdat$DSS.status[gdat$DSS.status == "dead"] <- "Dead"
gdat$surv2.5 <- with(gdat, ifelse(DSS.time >= 5, 1, -1))
psdex <- 2:(grep("tissue", colnames(gdat)) - 1)
gdat[, psdex] <- lapply(gdat[, psdex], scl.cen)
ldat <- gather_(gdat, "gene", "expression", colnames(gdat[, psdex]))
ldat <- ldat %>% group_by(gene) %>% mutate(zexpression = (expression - mean(expression, na.rm = TRUE))/sd(expression, na.rm = TRUE))
cldat <- gdat[, -psdex]

fit.superpc <- function(dat0, cldat){

  cors <- dat0 %>% group_by(gene) %>%
    do(data.frame(
      t(summary(coxph(Surv(DSS.time, DSS.status == "Dead") ~ zexpression,
                           data = .,
                           eps = 1e-4))$coefficients[c(1,5)]))
      )

  ## select number of genes

  cors <- cors[cors$X2 < .005, ]
  seldat0 <- dat0 %>% select(ID, gene, zexpression) %>%
    filter(gene %in% cors$gene) %>% spread(gene, zexpression)

  seldat0[, 2:ncol(seldat0)] <- lapply(2:ncol(seldat0), function(it){

    unlist(cors %>% filter(gene == colnames(seldat0[, it])) %>% `$`("X1") * seldat0[, it])

  })

  seldat0 <- merge(seldat0, cldat, by = "ID", all.y = FALSE)

  ngen.max <- 25

  ## start with most significant gene and add 1 at a time

  gen.cand <- gen.remain <- cors$gene[order(cors$X2)]
  dateval <- datfin <- seldat0
  dateval$cand.score <- datfin$score <- 0
  ngen.0 <- ngen.cur <- 0
  cors.num <- gen.in <- NULL
  while(ngen.cur < ngen.max){

    cors.try <- rep(NA, length(gen.remain))
    names(cors.try) <- gen.remain
    for(gen.try in gen.remain){
      dateval$cand.score <- unlist(dateval$cand.score + dateval[, gen.try])
      cors.try[gen.try] <- survConcordance(Surv(DSS.time, DSS.status == "Dead") ~ cand.score, data = dateval)$concordance
    }
    # select score with highest concordance
    gen.sel <- names(sort(cors.try, decreasing = TRUE)[1])
    datfin$score <- unlist(datfin$score + datfin[, gen.sel])

    gen.remain <- setdiff(gen.remain, gen.sel)
    gen.in <- append(gen.in, gen.sel)
    dateval$cand.score <- datfin$score

    cors.num <- append(cors.num, survConcordance(Surv(DSS.time, DSS.status == "Dead") ~ score, data = datfin)$concordance)
    ngen.cur <- ngen.cur + 1

  }

  fin.sig <- gen.in[1:which(cors.num == max(cors.num))]

  seldat1 <- dat0 %>% select(ID, gene, zexpression) %>%
    filter(gene %in% fin.sig) %>% spread(gene, zexpression)
  seldat1 <- merge(seldat1, cldat, by = "ID", all.y = FALSE)

  form.make <- as.formula(paste0("Surv(DSS.time, DSS.status == 'Dead') ~ ", paste(paste0("`", fin.sig, "`"), collapse = "+")))
  fit.cox <- coxph(form.make, data = seldat1)

  datfin$score <- lps <- predict(fit.cox, type = "lp")
  ## middle 75%
  c.cands <- sort(lps)[floor(.25 * length(lps)):floor(.75 * length(lps))]
  p.cands <- sapply(c.cands, function(c.cand){
    ltest <- survdiff(Surv(DSS.time, DSS.status == "Dead") ~ I(score > c.cand), data = datfin)
    pchisq(ltest$chisq, df = 1, lower.tail = FALSE)
  })
  cutoff <- c.cands[which(order(p.cands) == 1)]

  list(cox.fit = fit.cox, cutoff = cutoff)

}


predict.superpc <- function(dat0, las.fit){

  coxests <- predict(las.fit$cox.fit, newdata = dat0, type = 'lp')
  list(lps = coxests, riskgrp = coxests > las.fit$cutoff)


}

eval.predict <- function(dat0, las.fit){

  lpreds <- predict.superpc(dat0, las.fit)
  sc <- survConcordance(Surv(dat0$DSS.time, dat0$DSS.status == "Dead") ~ lpreds$lps)

  disc <- coxph(Surv(dat0$DSS.time, dat0$DSS.status == "Dead") ~ lpreds$riskgrp)$coefficients

  c(calib = sc$concordance, discrim = disc)

}








