
args <- commandArgs(TRUE)
print(args)
ntasks <- as.numeric(args[1])

library(snow)
cl <- makeCluster(ntasks, type = "MPI")

clusterEvalQ(cl, {

auc <- function(lll, ddd){

  ROCR::performance(ROCR::prediction(lll, ddd), "auc")@y.values[[1]]

}

or <- function(sc, resp){

  sc01 <- sc > median(sc, na.rm = TRUE)
  tab <- table(sc01, resp)
  tab[1 ,1] * tab[2, 2] / (tab[1, 2] * tab[2, 1])

}

createnosignaldata <- function(n, d){

  data.frame(Y = rbinom(n, 1, p = .3), X = matrix(rnorm(n * d), nrow = n))

}

selectvars <- function(data0){

  pees <- sapply(data0[, -1], function(v1){

    fit0 <- glm(data0$Y ~ v1, family = "binomial", model = FALSE, y = FALSE)
    summary(fit0)$coefficients[2, 4]

  })

  order(pees)[1:25]

}

runclassifier <- function(data0, selected){

  dat2 <- data0[, c(1, selected + 1)]
  fit1 <- glm(Y ~ ., data = dat2, family = "binomial", model = FALSE, y = FALSE)

  fit1

}

fitclassifier <- function(fit1, data1){

  scores <- predict(fit1, newdata = data1, type = "response")
  response <- as.factor(data1$Y)

  c(AUC = auc(scores, response), OR = or(scores, response))

}

naiveest <- function(n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)

  train <- dat0

  selecttr <- selectvars(train)
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, train)

}

zhutest.hold <- function(trratio = .5, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  npart.tr <- floor(trratio * n)
  tr.dex <- sample(1:n, npart.tr)
  ho.dex <- setdiff(1:n, tr.dex)

  train <- dat0[tr.dex, ]
  hold <- dat0[ho.dex, ]

  selecttr <- selectvars(dat0) ## selection performed on entire dataset
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, hold)

}


holdoutest <- function(trratio = .5, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  npart.tr <- floor(trratio * n)
  tr.dex <- sample(1:n, npart.tr)
  ho.dex <- setdiff(1:n, tr.dex)

  train <- dat0[tr.dex, ]
  hold <- dat0[ho.dex, ]

  selecttr <- selectvars(train)
  fit <- runclassifier(train, selecttr)
  fitclassifier(fit, hold)

}

zhu.cv <- function(k = 10, K = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  selectin <- selectvars(dat0)

  cvests <- replicate(K, {

    oot.dex <- sample((1:n)[which(dat0$Y == 1)], floor(k * mean(dat0$Y)))
    oot.dex <- c(oot.dex, sample((1:n)[which(dat0$Y == 0)], ceiling(k * mean(dat0$Y == 0))))
    oot <- dat0[oot.dex, ]
    estin <- dat0[setdiff(1:n, oot.dex), ]

    fit <- runclassifier(estin, selectin)
    fitclassifier(fit, oot)

  })
  rowMeans(cvests, na.rm = TRUE)

}


cvest <- function(k = 10, K = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  cvests <- replicate(K, {

    oot.dex <- sample((1:n)[which(dat0$Y == 1)], floor(k * mean(dat0$Y)))
    oot.dex <- c(oot.dex, sample((1:n)[which(dat0$Y == 0)], ceiling(k * mean(dat0$Y == 0))))
    oot <- dat0[oot.dex, ]
    estin <- dat0[setdiff(1:n, oot.dex), ]

    selectin <- selectvars(estin)
    fit <- runclassifier(estin, selectin)
    fitclassifier(fit, oot)

  })
  rowMeans(cvests, na.rm = TRUE)

}


bootest <- function(B = 50, n = 2000, d = 500){

  dat0 <- createnosignaldata(n, d)
  bootests <- replicate(B, {

    boot.dex <- sample(1:n, n, replace = TRUE)
    bin <- dat0[boot.dex, ]
    notbin <- setdiff(1:n, unique(boot.dex))

    selectin <- selectvars(bin)
    fit <- runclassifier(bin, selectin)
    fitclassifier(fit, dat0[notbin, ])

  })

  rowMeans(bootests)

}

runsim <- function(n = 1000, d = 500, K = 200, B = 200){

  c(holdout.5 = holdoutest(trratio = .5, n = n, d = d),
    holdout.3 = holdoutest(trratio = .666, n = n, d = d),
    cv.10 = cvest(k = 50, K = K, n = n, d = d),
    cv.100 = cvest(k = 100, K = K, n = n, d = d),
    boot = bootest(B = B, n = n, d = d),
    zhu.hold = zhutest.hold(trratio = .5, n = n, d = d),
    zhu.cv = zhu.cv(k = 50, K = K, n = n, d = d),
    naive = naiveest(n = n, d = d))

}

})

runsim <- function(n = 1000, d = 500, K = 200, B = 200){

  c(holdout.5 = holdoutest(trratio = .5, n = n, d = d),
    holdout.3 = holdoutest(trratio = .666, n = n, d = d),
    cv.10 = cvest(k = 10, K = K, n = n, d = d),
    cv.100 = cvest(k = 100, K = K, n = n, d = d),
    boot = bootest(B = B, n = n, d = d),
    zhu.hold = zhutest.hold(trratio = .5, n = n, d = d),
    zhu.cv = zhu.cv(k = 10, K = K, n = n, d = d),
    naive = naiveest(n = n, d = d))

}


cvsims <- clusterApplyLB(cl, 1:1000, function(i) runsim())
stopCluster(cl)
save(cvsims, file = "cvsim-result.RData")

