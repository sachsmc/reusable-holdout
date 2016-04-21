createnosignaldata <- function(n, d){

  train <- data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n))
  holdout <- data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n))
  test <- data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n))

  list(train = train, holdout = holdout, test = test)

}

add_bias <- function(set, nbiased = 50, size = .1){

  set[, 2:(2 + nbiased)] <- set[, 2:(2 + nbiased)] + set[, 1] * size
  set

}

createhighsignaldata <- function(n, d){

  train <- add_bias(data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n)))
  holdout <- add_bias(data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n)))
  test <- add_bias(data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n)))

  list(train = train, holdout = holdout, test = test)

}

runClassifier <- function(n, d, krange, dataList, threshold = .1, tolerance = .05){

  # select based on correlation with combined sample

  cor.train <- cor(dataList$train)[1, -1]
  cor.hold <- cor(dataList$holdout)[1, -1]

  check <- cor.train[sign(cor.train) == sign(cor.hold)]

  ## thresholdout selection

  diffs <- abs(cor.train - cor.hold)
  abovethresh <- diffs > threshold + rnorm(d, 0, tolerance)
  cor.threshold <- cor.hold + rnorm(d, 0, tolerance)
  cor.threshold[!abovethresh] <- cor.train[!abovethresh]

  check.thresh <- cor.threshold[sign(cor.train) == sign(cor.threshold)]

  vals <- data.frame(krange = krange, ftrain = NA, fholdout = NA, ftest = NA,
                     fthresholdout = NA, fstrict = NA, ftest.thresh = NA,
                     ftrain.lasso = NA, ftest.lasso = NA)
  for(k in krange){

    vars <- names(sort(abs(check), decreasing = TRUE)[1:k])
    form <- as.formula(paste("Y ~ ", paste(vars, collapse = " + ")))
    fit.train <- lm(form, data = dataList$train)

    ftrain <- mean((dataList$train$Y - fit.train$fitted)^2)
    fholdout <- mean((dataList$holdout$Y - predict(fit.train, newdata = dataList$holdout))^2)
    ftest <- mean((dataList$test$Y - predict(fit.train, newdata = dataList$test))^2)

    vals[vals$krange == k, ]$ftrain <- ftrain
    vals[vals$krange == k, ]$fholdout <- fholdout
    vals[vals$krange == k, ]$ftest <- ftest

    ## implement thresholdout
    vars2 <- names(sort(abs(check.thresh), decreasing = TRUE)[1:k])
    form2 <- as.formula(paste("Y ~ ", paste(vars2, collapse = " + ")))
    fit.train2 <- lm(form2, data = dataList$train)

    ftrain <- mean((dataList$train$Y - fit.train2$fitted)^2)
    fholdout <- mean((dataList$holdout$Y - predict(fit.train2, newdata = dataList$holdout))^2)
    ftest <- mean((dataList$test$Y - predict(fit.train2, newdata = dataList$test))^2)

    threshtest <- abs(ftrain - fholdout) < threshold + rnorm(1, 0, tolerance)
    fthresholdout <- ifelse(threshtest, ftrain, fholdout + rnorm(1, 0, tolerance))

    vals[vals$krange == k, ]$fthresholdout <- fthresholdout
    vals[vals$krange == k, ]$ftest.thresh <- ftest

    ## strict variable selection
    vars3 <- names(sort(abs(cor.train), decreasing = TRUE)[1:k])
    form3 <- as.formula(paste("Y ~ ", paste(vars3, collapse = " + ")))
    fit.train3 <- lm(form3, data = dataList$train)

    fholdout <- mean((dataList$holdout$Y - predict(fit.train3, newdata = dataList$holdout))^2)

    vals[vals$krange == k, ]$fstrict <- fholdout

     ## lasso

    fit.l <- glmnet(model.matrix(form3, dataList$train)[, -1], dataList$train$Y)
    vals[vals$krange == k, ]$ftrain.lasso <- mean((dataList$train$Y -
                                                     predict(fit.l, newx = model.matrix(form3, dataList$train)[, -1]))^2)
    vals[vals$krange == k, ]$ftest.lasso <- mean((dataList$test$Y -
                                                    predict(fit.l, newx = model.matrix(form3, dataList$test)[, -1]))^2)



  }

  vals

}


listMean <- function(lmats, keepcol = 1){

  res <- matrix(0, nrow = nrow(lmats[[1]]), ncol = ncol(lmats[[1]]))
  for(i in 1:length(lmats)){
    res <- res + lmats[[i]]
  }
  ret <- res/length(lmats)
  ret[, keepcol] <- lmats[[1]][, keepcol]
  ret

}

listSd <- function(lmats, keepcol = 1){

  res <- matrix(0, nrow = nrow(lmats[[1]]), ncol = ncol(lmats[[1]]))
  means <- listMean(lmats)
  for(i in 1:length(lmats)){
    res <- res + (lmats[[i]] - means)^2
  }

  ret <- sqrt(res/length(res))
  ret[, keepcol] <- lmats[[1]][, keepcol]
  ret

}

repeatexp <- function(cl, n, d, krange, reps, datafn, threshold = .1, tolerance = .04){

  res <- clusterApplyLB(cl, 1:reps, function(i){

    dataList <- datafn(n, d)
    runClassifier(n, d, krange, dataList, threshold, tolerance)

  })

  list(mean = listMean(res), sd = listSd(res))

}


args <- commandArgs(TRUE)
print(args)
ntasks <- as.numeric(args[1])

library(glmnet)
library(snow)
cl <- makeCluster(ntasks, type = "MPI")

clusterEvalQ(cl, {
library(glmnet)
createnosignaldata <- function(n, d){

  train <- data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n))
  holdout <- data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n))
  test <- data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n))

  list(train = train, holdout = holdout, test = test)

}

add_bias <- function(set, nbiased = 50, size = .1){

  set[, 2:(2 + nbiased)] <- set[, 2:(2 + nbiased)] + set[, 1] * size
  set

}

createhighsignaldata <- function(n, d){

  train <- add_bias(data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n)))
  holdout <- add_bias(data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n)))
  test <- add_bias(data.frame(Y = rnorm(n), X = matrix(rnorm(n * d), nrow = n)))

  list(train = train, holdout = holdout, test = test)

}

runClassifier <- function(n, d, krange, dataList, threshold = .1, tolerance = .05){

  # select based on correlation with combined sample

  cor.train <- cor(dataList$train)[1, -1]
  cor.hold <- cor(dataList$holdout)[1, -1]

  check <- cor.train[sign(cor.train) == sign(cor.hold)]

  ## thresholdout selection

  diffs <- abs(cor.train - cor.hold)
  abovethresh <- diffs > threshold + rnorm(d, 0, tolerance)
  cor.threshold <- cor.hold + rnorm(d, 0, tolerance)
  cor.threshold[!abovethresh] <- cor.train[!abovethresh]

  check.thresh <- cor.threshold[sign(cor.train) == sign(cor.threshold)]

  vals <- data.frame(krange = krange, ftrain = NA, fholdout = NA, ftest = NA,
                     fthresholdout = NA, fstrict = NA, ftest.thresh = NA,
                     ftrain.lasso = NA, ftest.lasso = NA)
  for(k in krange){

    vars <- names(sort(abs(check), decreasing = TRUE)[1:k])
    form <- as.formula(paste("Y ~ ", paste(vars, collapse = " + ")))
    fit.train <- lm(form, data = dataList$train)

    ftrain <- mean((dataList$train$Y - fit.train$fitted)^2)
    fholdout <- mean((dataList$holdout$Y - predict(fit.train, newdata = dataList$holdout))^2)
    ftest <- mean((dataList$test$Y - predict(fit.train, newdata = dataList$test))^2)

    vals[vals$krange == k, ]$ftrain <- ftrain
    vals[vals$krange == k, ]$fholdout <- fholdout
    vals[vals$krange == k, ]$ftest <- ftest

    ## implement thresholdout
    vars2 <- names(sort(abs(check.thresh), decreasing = TRUE)[1:k])
    form2 <- as.formula(paste("Y ~ ", paste(vars2, collapse = " + ")))
    fit.train2 <- lm(form2, data = dataList$train)

    ftrain <- mean((dataList$train$Y - fit.train2$fitted)^2)
    fholdout <- mean((dataList$holdout$Y - predict(fit.train2, newdata = dataList$holdout))^2)
    ftest <- mean((dataList$test$Y - predict(fit.train2, newdata = dataList$test))^2)

    threshtest <- abs(ftrain - fholdout) < threshold + rnorm(1, 0, tolerance)
    fthresholdout <- ifelse(threshtest, ftrain, fholdout + rnorm(1, 0, tolerance))

    vals[vals$krange == k, ]$fthresholdout <- fthresholdout
    vals[vals$krange == k, ]$ftest.thresh <- ftest

    ## strict variable selection
    vars3 <- names(sort(abs(cor.train), decreasing = TRUE)[1:k])
    form3 <- as.formula(paste("Y ~ ", paste(vars3, collapse = " + ")))
    fit.train3 <- lm(form3, data = dataList$train)

    fholdout <- mean((dataList$holdout$Y - predict(fit.train3, newdata = dataList$holdout))^2)

    vals[vals$krange == k, ]$fstrict <- fholdout

     ## lasso

    fit.l <- glmnet(model.matrix(form3, dataList$train)[, -1], dataList$train$Y)
    vals[vals$krange == k, ]$ftrain.lasso <- mean((dataList$train$Y -
                                                     predict(fit.l, newx = model.matrix(form3, dataList$train)[, -1]))^2)
    vals[vals$krange == k, ]$ftest.lasso <- mean((dataList$test$Y -
                                                    predict(fit.l, newx = model.matrix(form3, dataList$test)[, -1]))^2)



  }

  vals

}


})
nosignal <- repeatexp(cl, 2000, 2000, krange = c(5, 10, 20, 50, 100, 150, 200, 250, 300, 400, 500), reps = 1000, datafn = createnosignaldata)
highsignal <- repeatexp(cl, 2000, 2000, krange = c(5, 10, 20, 50, 100, 150, 200, 250, 300, 400, 500), reps = 1000, datafn = createhighsignaldata)
stopCluster(cl)
save(nosignal, file = "nosignal-result.RData")
save(highsignal, file = "highsignal-result.RData")

