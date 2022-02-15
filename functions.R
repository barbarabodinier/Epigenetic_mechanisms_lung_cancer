PlotROCCurves <- function(mymodels, niter = 100) {
  mod <- mymodels[1]
  print(mod)
  model <- eval(parse(text = mod))
  plot(1 - apply(model$specificity, 1, mean), apply(model$sensitivity, 1, mean),
    type = "l", lwd = 2, col = colors[1],
    xlab = "False Positive Rate", ylab = "True Positive Rate", las = 1, cex.lab = 1.5,
    panel.first = abline(0, 1, lty = 3)
  )

  for (i in 1:length(mymodels)) { # area between most extreme points (quantiles)
    mod <- mymodels[i]
    print(mod)
    model <- eval(parse(text = mod))
    polygon(c(
      1 - apply(model$specificity, 1, FUN = function(x) {
        sort(x)[quantiles[1] * niter]
      }),
      rev(1 - apply(model$specificity, 1, FUN = function(x) {
        sort(x)[quantiles[2] * niter]
      }))
    ),
    c(
      apply(model$sensitivity, 1, FUN = function(x) {
        sort(x)[quantiles[1] * niter]
      }),
      rev(apply(model$sensitivity, 1, FUN = function(x) {
        sort(x)[quantiles[2] * niter]
      }))
    ),
    col = tcolors[i], border = NA
    )
  }

  for (i in 1:length(mymodels)) { # curve with most extreme points
    mod <- mymodels[i]
    model <- eval(parse(text = mod))
    qs <- sort.list(model$AUC)[c(quantiles[1] * niter, quantiles[2] * niter)]
    lines(1 - apply(model$specificity, 1, FUN = function(x) {
      sort(x)[quantiles[1] * niter]
    }),
    apply(model$sensitivity, 1, FUN = function(x) {
      sort(x)[quantiles[1] * niter]
    }),
    type = "l", lwd = 0.3, col = colors[i]
    )
    lines(1 - apply(model$specificity, 1, FUN = function(x) {
      sort(x)[quantiles[2] * niter]
    }),
    apply(model$sensitivity, 1, FUN = function(x) {
      sort(x)[quantiles[2] * niter]
    }),
    type = "l", lwd = 0.3, col = colors[i]
    )
  }

  for (i in 1:length(mymodels)) {
    mod <- mymodels[i]
    model <- eval(parse(text = mod))
    lines(1 - apply(model$specificity, 1, mean), apply(model$sensitivity, 1, mean), type = "l", lwd = 2, col = colors[i])
  }
}


rocfoldslogistic <- function(xdata, ydata, M = 5, niter = 1) {
  t0 <- Sys.time()
  sensitivity <- specificity <- AUC <- NULL
  for (i in 1:niter) {
    # scv=sample(rep(1:M,length.out=nrow(xdata)), size=nrow(xdata))
    tmp <- split(1:length(ydata), f = as.factor(ydata))
    cacofolds <- lapply(tmp, FUN = function(x) {
      sample(rep(1:M, length.out = length(x)), size = length(x))
    })
    scv <- rep(NA, length(ydata))
    for (k in 1:2) {
      scv[tmp[[k]]] <- cacofolds[[k]]
    } # balanced folds
    beta <- NULL
    fitted <- NULL
    for (k in 1:M) { # CV
      xdata_train <- xdata[scv != k, , drop = FALSE]
      y_train <- ydata[scv != k]
      xdata_test <- xdata[scv == k, , drop = FALSE]
      colnames(xdata_train) <- colnames(xdata_test) <- colnames(xdata)
      y_test <- ydata[scv == k]
      f <- paste("y_train ~", paste(colnames(xdata), collapse = " + "))
      # print(f)
      model <- glm(as.formula(f), data = as.data.frame(xdata_train), family = binomial)
      beta <- c(beta, coefficients(model)[-1])
      fitted <- c(fitted, predict(model, as.data.frame(xdata_test)))
    }
    roc_iter <- roc(response = ydata, predictor = fitted[names(ydata)])
    sensitivity <- cbind(sensitivity, roc_iter$sensitivities)
    specificity <- cbind(specificity, roc_iter$specificities)
    AUC <- c(AUC, roc_iter$auc[1])
  }
  t1 <- Sys.time()
  print(t1 - t0)
  return(list(
    sensitivity = sensitivity, specificity = specificity, AUC = AUC,
    nobs = nrow(xdata), ncases = sum(as.character(ydata) == "1")
  ))
}
