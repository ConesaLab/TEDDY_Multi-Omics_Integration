###########################################################################
#######################     NPLSDAfunctions.R     #########################
###########################################################################
# Author: Leandro Balzano-Nogueira
# Genetics Institute, University of Florida (Gainesville)

# These are the functions created for performing all Multiway and multiblock strategies used for the project

###########################################################################
# Functions: 

Imputemethod<-function (X, fac=c(2, 2, 2), conver = 1e-07, max.iter = 1000) 
{
  # Function to impute missing data from data array based of the Tucker3 method. It performs iterations until reach 
  # a convergence limit
  # You can control the convergence (conver) value and the number of iterations (max.iter)
  # Also you can control the number of components per mode to be retain in the Tucker3 model (fac) 
  # Return: An Array with no missing values
  if (any(is.na(X))){
    print("Imputing missing values with Tucker3 method")
    xdims<-dim(X)
  NA.position <- which(is.na(X))
  NA.values <- rnorm(length(NA.position), mean(X, na.rm = TRUE),sd(X, na.rm = TRUE))
  
  X[NA.position] <- NA.values
  SSTc <- 0
  for (it in 1:max.iter) {
    #print(paste("Iteration #", it, sep = ""))
    SSTc.old <- SSTc
    Xe <- tucker3mod2(X, fac, verbose = F)
    X[NA.position] <- Xe$Xhat[NA.position]
    SSTc <- sum(Xe$G^2)
    delta.SSTc <- abs(SSTc - SSTc.old)/SSTc
    if (delta.SSTc < conver) 
      break
  }
  if (it == max.iter) {
    print(paste("max. iterations (", max.iter, ")reached"))
    print(paste("SST increment:", formatC(delta.SSTc)))
  }
  } else{ 
    print("No NA values present in the 3D array, no imputation performed")
  }
  return(X)
}
###########################################################################
tucker3mod2<-function (X, COMP, conver = 1e-07, max.iter = 10000, verbose = FALSE) 
{
  # Function to perform a Tucker3 analysis from arrays nxpxq where n are individuals, p are variables and
  # q are time points
  # You can control the convergence (conver) value and the number of iterations (max.iter)
  # Also you can control the number of components per mode to be retain in the Tucker3 model (COMP) 
  # Return: Ten useful parameters
  # 1.- FactorsX: The Scores of the individuals (Mode1), the loadings of the Variables (Mode2)
  #               and the loadings of the Time (Mode3) that you decided to retain
  # 2.- G: The Core array of the Tucker3 which retains the correlations among the components of the modes
  # 3.- Xhat: The array of the best fitted model
  # 4.- SST: Total Sum of squares
  # 5.- expl.var: Explained variance of the model
  # 6.- SSF: Sum of squares of the Fitted Model
  # 7.- SSE: Sum of squares of the Error
  # 8.- GCV: The Generalized cross-validation. Is a method for evaluating the smoothing parameter 
  #          the more small the value, the better the model.
  # 9.- edf: Effective degrees of Freedom. is a measure of the degrees of freedom in each mode of the Tucker3 analysis and in the core
  #          array. For values obtained, for Mode1, Mode2, Mode3 and Core array
  # 10.- tdf: Total degrees of Freedom of the N-way model
  
  
  
  if (any(is.na(X))) {
    print("Missing values are taken care of")
    X <- Imputemethod(X, COMP)
  }
  xdims<-dim(X)
  XA <- matrix(X, xdims[1], xdims[2] * xdims[3])
  XB <- matrix(aperm(X, perm = c(2, 1, 3)), xdims[2], xdims[1] * xdims[3])
  XC <- matrix(aperm(X, perm = c(3, 1, 2)), xdims[3],xdims[1] * xdims[2])
  A <- svd(matrix(rnorm(xdims[1] * COMP[1]), xdims[1],COMP[1]), nu = COMP[1], nv = 0)$u
  B <- svd(matrix(rnorm(xdims[2] * COMP[2]), xdims[2],COMP[2]), nu = COMP[2], nv = 0)$u
  C <- svd(matrix(rnorm(xdims[3] * COMP[3]), xdims[3],COMP[3]), nu = COMP[3], nv = 0)$u
  SST <- 0
  SSF<-0
  if (verbose) 
    print(paste("iterations","               SST","                   SSF","                    SSE", 
                "                          expl.var", sep = ""))
  for (it in 1:max.iter) {
    count <- it
    SST.old <- SST
    SSF.old <- SSF
    A1 <- svd(XA %*% kronecker(C, B),nu = COMP[1], nv=0)$u
    B1 <- svd(XB %*% kronecker(C, A),nu = COMP[2], nv=0)$u
    C1 <- svd(XC %*% kronecker(B, A),nu = COMP[3], nv=0)$u
    Gu <- t(A1) %*% XA %*% kronecker(C1, B1)
    Xhat.u <- A1 %*% Gu %*% kronecker(t(C1), t(B1))
    #options(digits = 3)
    SST<-sum(X^2)
    SSE <- SST - sum(rowSums(Gu^2))
    SSF <- sum(Xhat.u^2)
    expl.var <- SSF/sum(XA^2)
    
    Adf <- xdims[1] * COMP[1] - COMP[1] * (COMP[1] + 1)/2
    Bdf <- xdims[2] * COMP[2] - COMP[2] * (COMP[2] + 1)/2
    Cdf <- xdims[3] * COMP[3] - COMP[3] * (COMP[3] + 1)/2
    edf <- c(Adf, Bdf, Cdf, prod(COMP))  # edf= Effective degrees of freedom
    tdf<-(xdims[1]*COMP[1]) + (xdims[2]*COMP[2]) + (xdims[3]*COMP[3]) + COMP[1] + COMP[2] + COMP[3] - COMP[1]^2 - COMP[2]^2 - COMP[3]^2  
    pxdim <- prod(xdims)
    GCV <- (SSE/pxdim)/(1 - sum(edf)/pxdim)^2  # Generalized Cross-validation
    
    if (verbose) 
      print(paste(it, "            ", SST, "      ",SSF,"        ", SSE, "               ",expl.var, 
                  sep = ""))
    if (abs(SSF - SSF.old)/SST < conver || (sum((A1 - A)^2, 
                                                na.rm = T) < conver & sum((B1 - B)^2, na.rm = T) < 
                                            conver & sum((C1 - C)^2, na.rm = T) < conver)) 
      break
    A <- A1
    B <- B1
    C <- C1
  }
  G <- array(as.vector(Gu), COMP)
  if (count == 10000) 
    print("Maximum number of iterations reached")
  Xhat <- array(as.vector(Xhat.u), dim = dim(X), dimnames = dimnames(X))
  SSF <- sum(Xhat.u^2)
  expl.var <- SSF/sum(XA^2)
  
  rownames(A) <- dimnames(X)[[1]]
  rownames(B) <- dimnames(X)[[2]]
  rownames(C) <- dimnames(X)[[3]]
  Factors <- list(Mode1 = A, Mode2 = B, Mode3 = C)
  results <- list(FactorsX = Factors, G = G, Xhat = Xhat,
                  SST = SST, expl.var = expl.var, SSF=SSF, SSE=SSE,GCV=GCV,edf=edf,tdf=tdf)
  return(results)
}

###########################################################################
# Function to decompose the model column to PQR components separatelly
listTOdata.frame <- function(x, row.names=NULL, optional=FALSE, ...) {
  # Function to transform a list to a data.frame so we can be able to sum some parameters
  
  if(!all(unlist(lapply(x, class)) %in% 
          c('raw','character','complex','numeric','integer','logical'))) {
    warning('All elements of the list must be a vector.')
    NextMethod(x, row.names=row.names, optional=optional, ...)
  }
  allequal <- all(unlist(lapply(x, length)) == length(x[[1]]))
  havenames <- all(unlist(lapply(x, FUN=function(x) !is.null(names(x)))))
  if(havenames) { #All the vectors in the list have names we can use
    colnames <- unique(unlist(lapply(x, names)))
    df <- data.frame(matrix(
      unlist(lapply(x, FUN=function(x) { x[colnames] })),
      nrow=length(x), byrow=TRUE))
    names(df) <- colnames
  } else if(allequal) { #No names, but are of the same length
    df <- data.frame(matrix(unlist(x), nrow=length(x), byrow=TRUE), ...)
    hasnames <- which(unlist(lapply(x, FUN=function(x) !is.null(names(x)))))
    if(length(hasnames) > 0) { #We'll use the first element that has names
      names(df) <- names(x[[ hasnames[1] ]])
    }
  } else { 
    warning(paste("The length of vectors are not the same and do not ",
                  "are not named, the results may not be correct.", sep=''))
    #Find the largest
    lsizes <- unlist(lapply(x, length))
    start <- which(lsizes == max(lsizes))[1]
    df <- x[[start]]
    for(i in (1:length(x))[-start]) {
      y <- x[[i]]
      if(length(y) < length(x[[start]])) {
        y <- c(y, rep(NA, length(x[[start]]) - length(y)))
      }
      if(i < start) {
        df <- rbind(y, df)
      } else {
        df <- rbind(df, y)
      }
    }
    df <- as.data.frame(df, row.names=1:length(x))
    names(df) <- paste('Col', 1:ncol(df), sep='')
  }
  if(missing(row.names)) {
    row.names(df) <- names(x)
  } else {
    row.names(df) <- row.names
  }
  return(df)
}

###########################################################################
plsreg<-function (predictors, responses, comps = 2, crosval = TRUE) {
  # Function to perform a Partial least squares regression
  # Input: 
  # It needs a data.frame of predictor variables
  # Response(s) variable(s)
  # You can control the number of need components (comps)
  # You can perform or not a cross-validation type K-fold setted at one tenth.
  
  # Return: Nineteen useful parameters
  # 1.- x.scores: The Scores of the individuals in terms of the predictor variables
  # 2.- x.loads: The Loadings of the predictor variables
  # 3.- y.scores: The Scores of the individuals in terms of the response(s) variable(s)
  # 4.- y.loads: The Loadings of the response(s) variable(s)
  # 5.- cor.xt: Correlation amongst x scores and t (decomposed predictor matrix)
  # 6.- cor.yt: Correlation amongst y scores and t (decomposed predictor matrix)
  # 7.- cor.xu: Correlation amongst x scores and u (decomposed response matrix)
  # 8.- cor.yu: Correlation amongst y scores and u (decomposed response matrix)
  # 9.- cor.tu: Correlation amongst t (decomposed predictor matrix) and u (decomposed response matrix)
  # 10.- raw.wgs: Raw weights of the variables
  # 11.- mod.wgs: Weights of the variables in the model
  # 12.- std.coefs: Standard coefficients
  # 13.- reg.coefs: Regularized coefficients
  # 14.- y.pred: Predicted Y (Yhat)
  # 15.- resid: Residuals
  # 16.- expvar: Explained variance of the model
  # 17.- VIP: Variable importance for the projection
  # 18.- Q2: Q squared of the model
  # 19.- Q2cum: accumulated Q squared of the model
  
  
    X = as.matrix(predictors)
  if (any(is.na(X))) 
    stop("\nSorry, no missing values are allowed")
  n = nrow(X)
  p = ncol(X)
  if (is.null(colnames(X))) 
    colnames(X) = paste(rep("X", p), 1:p, sep = "")
  if (is.null(rownames(X))) 
    rownames(X) = 1:n
  Y = as.matrix(responses)
  if (any(is.na(Y))) 
    stop("\nSorry, no missing values are allowed")
  if (nrow(X) != nrow(Y)) 
    stop("\ndifferent number of rows in predictors and responses")
  q = ncol(Y)
  if (is.null(colnames(Y))) 
    colnames(Y) = paste(rep("Y", q), 1:q, sep = "")
  if (is.null(rownames(Y))) 
    rownames(Y) = 1:n
  if (p < 2 || q < 2) 
    stop("predictors and responses must have more than one column")
  if (!is.null(comps)) {
    nc = comps
    if (mode(nc) != "numeric" || length(nc) != 1 || nc <= 
        1 || (nc%%1) != 0 || nc > min(n, p)) 
      nc = min(n, p)
    if (nc == n) 
      nc = n - 1
  }
  else {
    if (n >= 10) {
      crosval = TRUE
      nc = min(n, p)
    }
    else {
      crosval = FALSE
      nc = 2
      message("\nSorry, no cross-validation with less than 10 observations")
    }
  }
  if (!is.logical(crosval)) 
    crosval = FALSE
  X.old = scale(X)
  Y.old = scale(Y)
  Wh = matrix(0, p, nc)
  Uh = matrix(0, n, nc)
  Th = matrix(0, n, nc)
  Ch = matrix(0, q, nc)
  Ph = matrix(0, p, nc)
  bh = rep(0, nc)
  if (crosval) {
    RSS = rbind(rep(n - 1, q), matrix(NA, nc, q))
    PRESS = matrix(NA, nc, q)
    Q2 = matrix(NA, nc, q)
    sets_size = c(rep(n%/%10, 9), n - 9 * (n%/%10))
    obs = sample(1:n, size = n)
    segments = vector("list", length = 10)
    ini = cumsum(sets_size) - sets_size + 1
    fin = cumsum(sets_size)
    for (k in 1:10) segments[[k]] = obs[ini[k]:fin[k]]
  }
  h = 1
  repeat {
    u.new = Y.old[, 1]
    w.old = rep(1, p)
    iter = 1
    repeat {
      w.new = t(X.old) %*% u.new/sum(u.new^2)
      w.new = w.new/sqrt(sum(w.new^2))
      t.new = X.old %*% w.new
      c.new = t(Y.old) %*% t.new/sum(t.new^2)
      u.new = Y.old %*% c.new/sum(c.new^2)
      w.dif = w.new - w.old
      w.old = w.new
      if (sum(w.dif^2) < 1e-06 || iter == 100) 
        break
      iter = iter + 1
    }
    p.new = t(X.old) %*% t.new/sum(t.new^2)
    if (crosval) {
      RSS[h + 1, ] = colSums((Y.old - t.new %*% t(c.new))^2)
      press = matrix(0, 10, q)
      for (i in 1:10) {
        aux = segments[[i]]
        uh.si = Y.old[-aux, 1]
        wh.siold = rep(1, p)
        itcv = 1
        repeat {
          wh.si = t(X.old[-aux, ]) %*% uh.si/sum(uh.si^2)
          wh.si = wh.si/sqrt(sum(wh.si^2))
          th.si = X.old[-aux, ] %*% wh.si
          ch.si = t(Y.old[-aux, ]) %*% th.si/sum(th.si^2)
          uh.si = Y.old[-aux, ] %*% ch.si/sum(ch.si^2)
          wsi.dif = wh.si - wh.siold
          wh.siold = wh.si
          if (sum(wsi.dif^2) < 1e-16 || itcv == 100) 
            break
          itcv = itcv + 1
        }
        Yhat.si = (X.old[aux, ] %*% wh.si) %*% t(ch.si)
        press[i, ] = colSums((Y.old[aux, ] - Yhat.si)^2)
      }
      PRESS[h, ] = colSums(press)
      Q2[h, ] = 1 - (PRESS[h, ]/RSS[h, ])
    }
    X.old = X.old - (t.new %*% t(p.new))
    Y.old = Y.old - (t.new %*% t(c.new))
    Wh[, h] = w.new
    Uh[, h] = u.new
    Th[, h] = t.new
    Ch[, h] = c.new
    Ph[, h] = p.new
    bh[h] = t(u.new) %*% t.new
    if (is.null(comps) && crosval) {
      if (sum(Q2[h, ] < 0.0975) == q || h == nc) 
        break
    }
    else {
      if (h == nc) 
        break
    }
    h = h + 1
  }
  
  ##################################################
  Th = Th[, 1:h]
  Ph = Ph[, 1:h]
  Wh = Wh[, 1:h]
  Uh = Uh[, 1:h]
  Ch = Ch[, 1:h]
  Ph = Ph[, 1:h]
  Ws = Wh %*% solve(t(Ph) %*% Wh,tol=1e-20)
  Bs = Ws %*% t(Ch)
  ricky<-apply(X, 2, sd)
  head(ricky)
  invricky<-1/ricky
  class(invricky)
  head(invricky)
  
  
  Br = diag(invricky) %*% Bs %*% diag(apply(Y, 2, sd))
  cte = as.vector(apply(Y, 2, mean) - apply(X, 2, mean) %*% Br)
  Y.hat = X %*% Br + matrix(rep(cte, each = n), n, q)
  resids = Y - Y.hat
  cor.xt = cor(X, Th)
  cor.yt = cor(Y, Th)
  cor.tu = cor(Th, Uh)
  cor.xu = cor(X, Uh)
  cor.yu = cor(Y, Uh)
  R2x = cor(X, Th)^2
  R2y = cor(Y, Th)^2
  Rdx = colMeans(R2x)
  Rdy = colMeans(R2y)
  EV = cbind(Rdx, cumsum(Rdx), Rdy, cumsum(Rdy))
  Rd.mat = matrix(0, h, h)
  for (j in 1:h) Rd.mat[1:j, j] = Rdy[1:j]
  VIP = sqrt((Wh^2) %*% Rd.mat %*% diag(p/cumsum(Rdy), h, h))
  if (crosval) {
    PRESS = PRESS[1:h, ]
    RSS = RSS[1:(h + 1), ]
    Q2 = Q2[1:h, ]
    Q2G = 1 - (rowSums(PRESS)/rowSums(RSS[1:h, ]))
    Q2cum = Q2
    Q2cum[1, ] = 1 - (PRESS[1, ]/RSS[1, ])
    for (i in 2:h) Q2cum[i, ] = 1 - apply(PRESS[1:i, ]/RSS[1:i, 
                                                           ], 2, prod)
    Q2Gcum = Q2G
    for (i in 1:h) Q2Gcum[i] = 1 - prod((rowSums(PRESS)/rowSums(RSS[-h, 
                                                                    ]))[1:i])
    Q2T = cbind(Q2, Q2G)
    Q2TC = cbind(Q2cum, Q2Gcum)
    q2 = c(paste(rep("Q2", q), colnames(Y), sep = "."), "Q2")
    q2c = c(paste(rep("Q2cum", q), colnames(Y), sep = "."), 
            "Q2cum")
    dimnames(Q2T) = list(paste(rep("t", h), 1:h, sep = ""), 
                         q2)
    dimnames(Q2TC) = list(paste(rep("t", h), 1:h, sep = ""), 
                          q2c)
  }
  else {
    Q2T = NULL
    Q2TC = NULL
  }
  dimnames(Wh) = list(colnames(X), paste(rep("w", h), 1:h, 
                                         sep = ""))
  dimnames(Ws) = list(colnames(X), paste(rep("w*", h), 1:h, 
                                         sep = ""))
  dimnames(Uh) = list(rownames(Y), paste(rep("u", h), 1:h, 
                                         sep = ""))
  dimnames(Th) = list(rownames(X), paste(rep("t", h), 1:h, 
                                         sep = ""))
  dimnames(Ch) = list(colnames(Y), paste(rep("c", h), 1:h, 
                                         sep = ""))
  dimnames(Ph) = list(colnames(X), paste(rep("p", h), 1:h, 
                                         sep = ""))
  dimnames(Bs) = list(colnames(X), colnames(Y))
  dimnames(Br) = list(colnames(X), colnames(Y))
  dimnames(cor.xt) = list(colnames(X), paste(rep("t", h), 1:h, 
                                             sep = ""))
  dimnames(cor.yt) = list(colnames(Y), paste(rep("t", h), 1:h, 
                                             sep = ""))
  dimnames(cor.xu) = list(colnames(X), paste(rep("u", h), 1:h, 
                                             sep = ""))
  dimnames(cor.yu) = list(colnames(Y), paste(rep("u", h), 1:h, 
                                             sep = ""))
  dimnames(cor.tu) = list(paste(rep("t", h), 1:h, sep = ""), 
                          paste(rep("u", h), 1:h, sep = ""))
  dimnames(EV) = list(paste(rep("t", h), 1:h, sep = ""), c("R2X", 
                                                           "R2Xcum", "R2Y", "R2Ycum"))
  dimnames(Y.hat) = list(rownames(Y), colnames(Y))
  dimnames(resids) = list(rownames(Y), colnames(Y))
  dimnames(VIP) = list(colnames(X), paste(rep("t", h), 1:h, 
                                          sep = ""))
  coeffs = rbind(Br, INTERCEPT = cte)
  structure(list(x.scores = Th, x.loads = Ph, y.scores = Uh, 
                 y.loads = Ch, cor.xt = cor.xt, cor.yt = cor.yt, cor.xu = cor.xu, 
                 cor.yu = cor.yu, cor.tu = cor.tu, raw.wgs = Wh, mod.wgs = Ws, 
                 std.coefs = Bs, reg.coefs = coeffs, y.pred = Y.hat, resid = resids, 
                 expvar = EV, VIP = VIP, Q2 = Q2T, Q2cum = Q2TC))
}




##################################################################################################################
# Determine the best fitted model

bestfittedmodel<- function (X, centering=2)
  
  # Function to determine the best fitted model for a Tucker3 analysis and its posterior use in NPLS-DA
  # This function was based on the paper:
  
  # Timmerman, M. E. & Kiers, H. A. L. (2000). Three-mode Principal Components Analysis: Choosing the numbers of
  # components and sensitivity to local optima. British Journal of Mathematical and Statistical Psychology. 53: 1-16.
  
  # In this paper they describe a very interesting strategy that i implemented here.
  
  # The input is an array with or without missing values.
  # You can center the data by individuals (1), by variables (2), by time (3) or not center it at all (0)
# The default would be centering by variables which is the more commonly used centering method.
# Return: a table with the combinations of number of components per mode that gathers most of the variability.
# The table has the following information
# 1.- sumnumberofcomponents: Is the sum of the number of components in all modes. The result just delivers the best
#                            representative of each subset.
# 2.- Model: The models with more explained variance
# 3.- P: Number of components in the mode1
# 4.- Q: Number of components in the mode2
# 5.- R: Number of components in the mode3
# 6.- EXPL.VAR: Explained Variance of each model
# 7.- DF: Degrees of Freedom of each model
# 8.- Deviance (SSE): Sum of squares of the Error of each model
# 9.- SS(FIT): Sum of squares of each Fitted Model
# 10.- Fitper: Explained Variance of each model in percentage
# 11.- diftm: The difference of Explained Variance between a model and the next model with more components
# 12.- btm: The ratio of the diftm between a model and the next model with more components
# 13.- CriVal: Critical value based on the euclidean norm of the array and the number of components.
#              If btm>CriVal, that is the best minimal fitted model. Nevertheless, it is important to analyze your own results


{
  if (centering==0) {
    xdim <- dim(X)
    Xcentrado<-X
    A3D<-Imputemethod(Xcentrado)
    xdimA3D<-dim(A3D)
    A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
  }
  
  if (centering==1) {
    xdim <- dim(X)
    #length(xdim) == 3L
    individuos<-rownames(X)
    variables<-colnames(X)
    X <- matrix(aperm(X, perm = c(1, 2, 3)), xdim[1], 
                xdim[2] * xdim[3])
    #dim(X) #2D por Col
    #X[1:10,1:10]
    X<-apply(X, 2, as.numeric)
    #X[1:10,1:10]
    xmns <- rowMeans(X, na.rm = T)
    xdim2<-dim(X)
    Y<- matrix(xmns, xdim2[1], xdim2[2])
    #Y[1:10,1:10]
    
    X <- X - Y
    #dim(X)
    X <- aperm(array(X, dim = xdim[c(1, 2, 3)]), 
               perm = c(1, 2, 3))
    dim(X)
    xdim3 <- dim(X)
    Xcporind <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
    rownames(Xcporind)<-individuos
    colnames(Xcporind)<-variables
    Xcentrado<-Xcporind
    A3D<-Imputemethod(Xcentrado)
    xdimA3D<-dim(A3D)
    A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
  }
  
  if (centering==2) {
    xdim <- dim(X)
    #length(xdim) == 3L
    individuos<-rownames(X)
    variables<-colnames(X)
    X <- matrix(aperm(X, perm = c(2, 1, 3)), xdim[2], 
                xdim[1] * xdim[3])
    #dim(X) #2D por Col
    #X[1:10,1:10]
    X<-apply(X, 2, as.numeric)
    #X[1:10,1:10]
    xmns <- rowMeans(X, na.rm = T)
    xdim2<-dim(X)
    Y<- matrix(xmns, xdim2[1], xdim2[2])
    #Y[1:10,1:10]
    
    X <- X - Y
    #dim(X)
    X <- aperm(array(X, dim = xdim[c(2, 1, 3)]), 
               perm = c(2, 1, 3))
    dim(X)
    xdim3 <- dim(X)
    Xcporvar <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
    rownames(Xcporvar)<-individuos
    colnames(Xcporvar)<-variables
    Xcentrado<-Xcporvar
    A3D<-Imputemethod(Xcentrado)
    xdimA3D<-dim(A3D)
    A<-matrix(A3D, xdimA3D[2], xdimA3D[1] * xdimA3D[3])
    
  }
  if (centering==3) {
    xdim <- dim(X)
    #length(xdim) == 3L
    individuos<-rownames(X)
    variables<-colnames(X)
    X <- matrix(aperm(X, perm = c(3, 1, 2)), xdim[3], 
                xdim[1] * xdim[2])
    #dim(X) #2D por Col
    #X[1:10,1:10]
    X<-apply(X, 2, as.numeric)
    #X[1:10,1:10]
    xmns <- rowMeans(X, na.rm = T)
    xdim2<-dim(X)
    Y<- matrix(xmns, xdim2[1], xdim2[2])
    #Y[1:10,1:10]
    
    X <- X - Y
    #dim(X)
    X <- aperm(array(X, dim = xdim[c(3, 2, 1)]), 
               perm = c(3, 2, 1))
    dim(X)
    xdim3 <- dim(X)
    Xcportime <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
    rownames(Xcportime)<-individuos
    colnames(Xcportime)<-variables
    Xcentrado<-Xcportime
    A3D<-Imputemethod(Xcentrado)
    xdimA3D<-dim(A3D)
    A<-matrix(A3D, xdimA3D[3], xdimA3D[1] * xdimA3D[2])
  }
  
  nf<-list(c(2,2,2),c(2,2,3),c(2,2,4),c(2,3,2),c(2,3,3),c(2,3,4),
           c(2,4,2),c(2,4,3),c(2,4,4),c(3,2,2),c(3,2,3),c(3,2,4),
           c(3,3,2),c(3,3,3),c(3,3,4),c(3,4,2),c(3,4,3),c(3,4,4),
           c(4,2,2),c(4,2,3),c(4,2,4),c(4,3,2),c(4,3,3),c(4,3,4),
           c(4,4,2),c(4,4,3),c(4,4,4)
  )
  vectorderows<-c(222:224,232:234,242:244,
                  322:324,332:334,342:344,
                  422:424,432:434,442:444
  )
  
  sumnumberofcomponents<-c(6,7,8,7,8,9,8,9,10,
                           7,8,9,8,9,10,9,10,11,
                           8,9,10,9,10,11,10,11,12
                           
  )
  
  vectitoRsqs<-NULL
  vectitoSSE<-NULL
  vectitoTDF<-NULL
  vectitoSSF<-NULL
  vectitopropSSF<-NULL
  valorrsq<-NULL
  valorSSE<-NULL
  valorTDF<-NULL
  valorSSF<-NULL
  valorpropSSF<-NULL
  
  for (a in 1:length(nf)) {
    tucker<-tucker3mod2(A3D, COMP = nf[[a]])
    valorrsq<-tucker$expl.var
    vectitoRsqs<-rbind(vectitoRsqs,valorrsq)
    valorSSE<-tucker$SSE
    vectitoSSE<-rbind(vectitoSSE,valorSSE)
    valorTDF<-tucker$tdf
    vectitoTDF<-rbind(vectitoTDF,valorTDF)
    valorSSF<-tucker$SSF
    vectitoSSF<-rbind(vectitoSSF,valorSSF)
  }
  rownames(vectitoRsqs)<-vectorderows;rownames(vectitoSSE)<-vectorderows;rownames(vectitoTDF)<-vectorderows;
  rownames(vectitoSSF)<-vectorderows
  vectitoRsqs<-data.frame(rownames(vectitoRsqs),vectitoRsqs)
  colnames(vectitoRsqs)<-c("Model","EXPL.VAR")
  vectitoSSE<-data.frame(rownames(vectitoSSE),vectitoSSE)
  colnames(vectitoSSE)<-c("Model","Deviance(SSE)")
  vectitoTDF<-data.frame(rownames(vectitoTDF),vectitoTDF)
  colnames(vectitoTDF)<-c("Model","DF")
  vectitoSSF<-data.frame(rownames(vectitoSSF),vectitoSSF)
  colnames(vectitoSSF)<-c("Model","SS(FIT)")
  
  #dim(vectitoRsqs);dim(vectitoSSE);dim(vectitoTDF);dim(vectitoSSF)
  tabla1<-merge(vectitoRsqs,vectitoTDF, by="Model"); dim(tabla1)
  tabla2<-merge(tabla1,vectitoSSE, by="Model"); dim(tabla2)
  tabla<-merge(tabla2,vectitoSSF, by="Model"); dim(tabla)
  tabla<-cbind(sumnumberofcomponents,tabla)
  tabla<-tabla[with(tabla,order(tabla[,1],-tabla[,6])),]
  #tabla
  #tabla<-tabla[with(tabla,order(-tabla[,3])),]
  #summarytable<-tabla
  tabla2<-tabla[!duplicated(tabla$sumnumberofcomponents),]
  #tabla2
  
  PQR<-strsplit(as.character(tabla2$Model), "")
  #PQR
  PQR2<-listTOdata.frame(PQR)
  colnames(PQR2)<-c("P","Q","R")
  tabla3<-cbind(tabla2,PQR2)
  tabla3<-tabla3[,c(1,2,7,8,9,3:6)]
  
  #class(tabla3)
  #sapply(tabla3, class)
  tabla3$Fitper<-tabla3$`EXPL.VAR`*100
  
  diftm<-c(0,tabla3[2:nrow(tabla3),10] - tabla3[1:nrow(tabla3)-1,10])
  tabla4<-cbind(tabla3,diftm)
  btm<-c(tabla4[1:nrow(tabla4),11] / tabla4[1:nrow(tabla4)+1,11])
  tabla5<-cbind(tabla4,btm)
  #sapply(tabla5, class)
  
  as.numeric.factor <- function(x) {as.numeric(levels(x))[x]}
  tabla5$Model<-as.numeric.factor(tabla5$Model)
  tabla5$P<-as.numeric.factor(tabla5$P)
  tabla5$Q<-as.numeric.factor(tabla5$Q)
  tabla5$R<-as.numeric.factor(tabla5$R)
  
  Vectitodetotales<-colSums(tabla5)
  totals<-as.data.frame(Vectitodetotales)
  tabla6<-rbind(tabla5,Totals=Vectitodetotales)
  
  
  # Calculating the Critical Value
  xdimsA<-dim(A3D)
  norma<-norm(A, type = c("O"))
  Smax<-totals[1,1]
  CriVal<-norma/(Smax-3)
  CriVal
  # If btm>CriVal, that is the best minimal fitted model
  
  return(list(summarytable=tabla6,CriVal=CriVal, Note="If btm>CriVal, that is the best minimal fitted model"))
}

###############################################################################################################
# NPLSmodel

NPLSmod<- function (XN, YN, factors=3, COMP= c(2,2,2), centering=2, conver = 1e-16, max.iteration = 10000) 
{
  # Function to perform the NPLS analysis from arrays nxpxq where n are individuals, p are variables and
  # q are time points
  # You can input the array with missing values, or also the Tucker3 array already processed
  # You can control the convergence (conver) value and the number of iterations (max.iter)
  # You can center by any of the three modes or not center at all
  # Also you can control the number of components per mode to be retain in the Tucker3 model to get the most quantity
  # of information for the NPLS (COMP). This has to be based on the results of the bestfittedmodel function 
  # You can control the number of factors to be retained per mode in the NPLS (this is different to 
  # the number of components per mode that best fit the Tucker3 model (COMP) 
  # Finally you can perform Comparative NPLS between two datasets or classic NPLS to explain response variables (Y)
  # with predictors variables (X)
  
  # Return: Seventeen useful parameters
  # 1.- FactorsX: The Scores of the individuals (Mode1) (wsuprai), the loadings of the Variables (Mode2) (wsupraj)
  #               and the loadings of the Time (Mode3) (wsuprak) that you decided to retain in the predictors matrix (X)
  # 2.- FactorsY: The Scores of the individuals (Mode1) (qsuprai), the loadings of the Variables (Mode2) (qsupraj)
  #               and the loadings of the Time (Mode3) (qsuprak) that you decided to retain in the response matrix (Y)
  # 3.- T: Scores matrix of the predictors matrix (X)
  # 4.- WsupraJ: Loadings matrix of the variables components of the predictors matrix (X)
  # 5.- WsupraK: Loadings matrix of the Time component of the predictors matrix (X)
  # 6.- U: Scores matrix of the response matrix (Y)
  # 7.- QsupraJ: Loadings matrix of the variables components of the response matrix (Y)
  # 8.- QsupraK: Loadings matrix of the Time component of the response matrix (Y)
  # 9.- B: Matrix of rfactors. These are the regression coefficients of the model that correlates the factors of X with Y
  # 10.- Gu: Deconvoluted Core array
  # 11.- G: Core matrix retaining the the correlations among the components of the modes 
  # 12.- Ypred: Matrix of predicted values of response (Y) matrix by predictors (X) matrix. 
  # 13.- explvar: Table of explained variance of X and explained variance of Y by X
  # 14.- VIP2D: Variable importance for the projection as a measure of how much a variable explains the bidimensional model
  #             in terms of variable in every timepoint studied
  # 15.- VIP3D: Variable importance for the projection as a measure of how much a variable explains the bidimensional model
  #             in terrms of variables, no matter the timeoint
  # 16.- tucorrelation: t-u correlations
  # 17.- residuals: Table of residuals
  
  require (plsdepot)
  require (mixOmics)
  require (MASS)
  require(abind)
  
    Ydim<-dim(YN)
  if (Ydim[2]>1){   # This part performs the centering in matrix X and Y
    print ("Performing NPLS regression class 2")
    if (any(is.na(XN))) {
      print("Missing values are taken care of")
      XN <- Imputemethod(XN, COMP)
      #matrizXnoNA<-XN
      #save(matrizXnoNA, file ="matrizXnoNA.RData")
    }
    if (any(is.na(YN))) {
      print("Missing values are taken care of")
      YN<- Imputemethod(YN, COMP)
      #matrizYnoNA<-YN
      #save(matrizYnoNA, file ="matrizYnoNA.RData")
    }
    if (centering==0) {
      xdim <- dim(XN)
      A3D<-XN
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      # MatrixY
      ydim <- dim(YN)
      B3D<-YN
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
    }
    
    if (centering==1) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(1, 2, 3)), xdim[1], 
                  xdim[2] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(1, 2, 3)]), 
                 perm = c(1, 2, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporind <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporind)<-individuos
      colnames(Xcporind)<-variables
      Xcentrado<-Xcporind
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      
      # MatrixY
      ydim <- dim(YN)
      #length(xdim) == 3L
      individuosY<-rownames(YN)
      variablesY<-colnames(YN)
      Y <- matrix(aperm(YN, perm = c(1, 2, 3)), ydim[1], 
                  ydim[2] * ydim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      Y<-apply(Y, 2, as.numeric)
      #X[1:10,1:10]
      ymns <- rowMeans(Y, na.rm = T)
      ydim2<-dim(Y)
      YY<- matrix(ymns, ydim2[1], ydim2[2])
      #Y[1:10,1:10]
      
      Y <- Y - YY
      #dim(X)
      Y <- aperm(array(Y, dim = ydim[c(1, 2, 3)]), 
                 perm = c(1, 2, 3))
      dim(Y)
      ydim3 <- dim(Y)
      Ycporind <- array(data = Y, dim = c(ydim3[1],ydim3[2],ydim3[3]),dimnames = list(NULL, NULL, colnames(Y[1,,])))
      rownames(Ycporind)<-individuosY
      colnames(Ycporind)<-variablesY
      Ycentrado<-Ycporind
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
      
    }
    
    if (centering==2) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(2, 1, 3)), xdim[2], 
                  xdim[1] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(2, 1, 3)]), 
                 perm = c(2, 1, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporvar <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporvar)<-individuos
      colnames(Xcporvar)<-variables
      Xcentrado<-Xcporvar
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[2], xdimA3D[1] * xdimA3D[3])
      
      
      #matrixY
      ydim <- dim(YN)
      #length(xdim) == 3L
      individuosY<-rownames(YN)
      variablesY<-colnames(YN)
      Y <- matrix(aperm(YN, perm = c(2, 1, 3)), ydim[2], 
                  ydim[1] * ydim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      Y<-apply(Y, 2, as.numeric)
      #X[1:10,1:10]
      ymns <- rowMeans(Y, na.rm = T)
      ydim2<-dim(Y)
      YY<- matrix(ymns, ydim2[1], ydim2[2])
      #Y[1:10,1:10]
      
      Y <- Y - YY
      #dim(X)
      Y <- aperm(array(Y, dim = ydim[c(2, 1, 3)]), 
                 perm = c(2, 1, 3))
      #dim(Y)
      ydim3 <- dim(Y)
      Ycporvar <- array(data = Y, dim = c(ydim3[1],ydim3[2],ydim3[3]),dimnames = list(NULL, NULL, colnames(Y[1,,])))
      rownames(Ycporvar)<-individuosY
      colnames(Ycporvar)<-variablesY
      Ycentrado<-Ycporvar
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[2], ydimB3D[1] * ydimB3D[3])
      
    }
    if (centering==3) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(3, 1, 2)), xdim[3], 
                  xdim[1] * xdim[2])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(3, 2, 1)]), 
                 perm = c(3, 2, 1))
      dim(X)
      xdim3 <- dim(X)
      Xcportime <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcportime)<-individuos
      colnames(Xcportime)<-variables
      Xcentrado<-Xcportime
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[3], xdimA3D[1] * xdimA3D[2])
      
      
      #MatrixY
      ydim <- dim(YN)
      #length(xdim) == 3L
      individuosY<-rownames(YN)
      variablesY<-colnames(YN)
      Y <- matrix(aperm(YN, perm = c(3, 1, 2)), ydim[3], 
                  ydim[1] * ydim[2])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      Y<-apply(Y, 2, as.numeric)
      #X[1:10,1:10]
      ymns <- rowMeans(Y, na.rm = T)
      ydim2<-dim(Y)
      YY<- matrix(ymns, ydim2[1], ydim2[2])
      #Y[1:10,1:10]
      
      Y <- Y - YY
      #dim(X)
      Y <- aperm(array(Y, dim = ydim[c(3, 2, 1)]), 
                 perm = c(3, 2, 1))
      dim(Y)
      ydim3 <- dim(Y)
      Ycportime <- array(data = Y, dim = c(ydim3[1],ydim3[2],ydim3[3]),dimnames = list(NULL, NULL, colnames(Y[1,,])))
      rownames(Ycportime)<-individuosY
      colnames(Ycportime)<-variablesY
      Ycentrado<-Ycportime
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[3], ydimB3D[1] * ydimB3D[2])
   }
    
    xdims<-dim(A3D)
    ydims<-dim(B3D)
    Tt <- U <- WsupraJ <- WsupraK <- QsupraJ <- QsupraK <- NULL
    B <- G <- matrix(0, ncol = factors, nrow = factors)
    Gu <- vector("list", factors)
    X <- matrix(A3D, xdims[1], xdims[2] * xdims[3])
    ordX <- length(dim(A3D))
    Y <- matrix(B3D, ydims[1], ydims[2] * ydims[3])
    ordY <- length(dim(B3D))
    Uf <- svd(Y)$u
    u <- Uf[, 1]
    f = 1
   #
    for (f in 1:factors) {
      it = 1
      while (it < max.iteration) {
        tX<-t(X)
        #tX<-tX[1:313]
        Zrow <- tX %*% u
        Z <- matrix(Zrow, nrow = dim(XN)[2], ncol = dim(XN)[3])
        svd.z <- svd(Z)
        wsupraj <- svd.z$u[, 1]
        wsuprak <- svd.z$v[, 1]
        tf <- X %*% kronecker(wsuprak, wsupraj)
        Vrow <- t(Y) %*% tf
        V <- matrix(Vrow, nrow = dim(YN)[2], ncol = dim(YN)[3])
        svd.v <- svd(V)
        qsupraj <- svd.v$u[, 1]
        qsuprak <- svd.v$v[, 1]
        uf <- Y %*% kronecker(qsuprak, qsupraj)
        if (sum((uf - u)^2) < conver) {
          print(paste("component number ", f))
          print(paste("number of iterations: ", it))
          it <- max.iteration
          Tt <- cbind(Tt, tf)
          WsupraJ <- cbind(WsupraJ, wsupraj)
          WsupraK <- cbind(WsupraK, wsuprak)
          QsupraJ <- cbind(QsupraJ, qsupraj)
          QsupraK <- cbind(QsupraK, qsuprak)
          bf <- ginv(t(Tt) %*% Tt) %*% t(Tt) %*% uf
          B[1:length(bf), f] <- bf
          U <- cbind(U, uf)
          TM <- ginv(t(Tt) %*% Tt) %*% t(Tt)
          WkM <- ginv(t(WsupraK) %*% WsupraK) %*% t(WsupraK)
          WjM <- ginv(t(WsupraJ) %*% WsupraJ) %*% t(WsupraJ)
          Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
          Y <- Y - Tt %*% bf %*% t(kronecker(qsuprak, qsupraj))
          f <- f + 1
          Uf <- svd(Y)$u
          u <- Uf[, 1]
        }
        else {
          u <- uf
          it <- it + 1
        }
      }
    }  
    rownames(Tt) <- dimnames(XN)[[1]]
    rownames(WsupraJ) <- dimnames(XN)[[2]]
    rownames(WsupraK) <- dimnames(XN)[[3]]
    rownames(U) <- dimnames(YN)[[1]]
    rownames(QsupraJ) <- dimnames(YN)[[2]]
    rownames(QsupraK) <- dimnames(YN)[[3]]
    G <- array(as.vector(Gu[[2]]), c(factors, factors,factors))
    FactorsX = list(Mode1 = Tt, Mode2 = WsupraJ, Mode3 = WsupraK)
    FactorsY = list(Mode1 = U, Mode2 = QsupraJ, Mode3 = QsupraK)
    #####
    # Evaluation of the model Vars
    PLSmodel<-plsreg(predictors = X, responses=Y, comps = factors, crosval = TRUE)
    Ypred<-PLSmodel$y.pred
    explvar<-PLSmodel$expvar
    VIP2D<-PLSmodel$VIP 
    tucorr<-PLSmodel$cor.tu  
    residuals<-PLSmodel$resid
    Q2<-PLSmodel$Q2
    Q2cum<-PLSmodel$Q2cum
    ####
    modo1 <- rowSums(FactorsX$Mode1[,c(1:factors)])
    modo2 <- rowSums(FactorsX$Mode2[,c(1:factors)])
    
    newx<-as.matrix(modo1) %*% as.matrix(t(modo2))
    #newx<-as.matrix(FactorsX$Mode1[,1]) %*% as.matrix(t(FactorsX$Mode2[,1]))
    PLSmodelNway<-plsreg(predictors = newx, responses=Y, comps = factors, crosval = FALSE)
    VIP3D<-PLSmodelNway$VIP 
    # head(VIP3D)
    
    
    
    ######
    result <- list(FactorsX = FactorsX, FactorsY = FactorsY, 
                   T = Tt, WsupraJ = WsupraJ, WsupraK = WsupraK, U = U, QsupraJ = QsupraJ, 
                   QsupraK = QsupraK, B = B, Gu = Gu, G = G,Ypred =Ypred,
                   explvar=explvar, VIP2D=VIP2D, VIP3D=VIP3D,tucorrelation=tucorr,residuals=residuals,
                   Q2=Q2,Q2cum=Q2cum
    )
    return(result)
    
}
    
    
    
    
  ###
  if (Ydim[2]==1){ # Classic NPLS
    print ("Performing Classic NPLS")
    if (any(is.na(XN))) {
      print("Missing values are taken care of")
      XN <- Imputemethod(XN, COMP)
      #matrizXnoNA<-XN
      #save(matrizXnoNA, file ="matrizXnoNA.RData")
    }
    if (any(is.na(YN))) {
      print("Missing values are taken care of")
      YN<- Imputemethod(YN, COMP)
      #matrizYnoNA<-YN
      #save(matrizYnoNA, file ="matrizYnoNA.RData")
    }
    if (centering==0) {
      xdim <- dim(XN)
      Xcentrado<-XN
      A3D<-Imputemethod(Xcentrado)
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      # MatrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Imputemethod(Ycentrado)
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
    }
    
    if (centering==1) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(1, 2, 3)), xdim[1], 
                  xdim[2] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(1, 2, 3)]), 
                 perm = c(1, 2, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporind <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporind)<-individuos
      colnames(Xcporind)<-variables
      Xcentrado<-Xcporind
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      
      # MatrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      
    }
    
    if (centering==2) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(2, 1, 3)), xdim[2], 
                  xdim[1] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(2, 1, 3)]), 
                 perm = c(2, 1, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporvar <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporvar)<-individuos
      colnames(Xcporvar)<-variables
      Xcentrado<-Xcporvar
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[2], xdimA3D[1] * xdimA3D[3])
      
      
      #matrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
    }
    if (centering==3) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(3, 1, 2)), xdim[3], 
                  xdim[1] * xdim[2])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(3, 2, 1)]), 
                 perm = c(3, 2, 1))
      dim(X)
      xdim3 <- dim(X)
      Xcportime <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcportime)<-individuos
      colnames(Xcportime)<-variables
      Xcentrado<-Xcportime
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[3], xdimA3D[1] * xdimA3D[2])
      
      
      #MatrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      
    }
    
    
    xdims<-dim(A3D)
    ydims<-dim(B3D)
    B3D2 <- abind(B3D, B3D, along=2) 
    ydims2<-dim(B3D2)
    Tt <- U <- WsupraJ <- WsupraK <- QsupraJ <- QsupraK <- NULL
    B <- G <- matrix(0, ncol = factors, nrow = factors)
    Gu <- vector("list", factors)
    X <- matrix(A3D, xdims[1], xdims[2] * xdims[3])
    ordX <- length(dim(A3D))
    Y <- matrix(B3D, ydims[1], ydims[2] * ydims[3])
    Y2 <- matrix(B3D2, ydims[1], ydims[2] * ydims[3])
    ordY <- length(dim(B3D))
    Uf <- svd(Y)$u
    u <- Uf[, 1]
    f = 1
    #
    for (f in 1:factors) {
      it = 1
      while (it < max.iteration) {
        tX<-t(X)
        #tX<-tX[1:313]
        Zrow <- tX %*% u
        Z <- matrix(Zrow, nrow = dim(XN)[2], ncol = dim(XN)[3])
        svd.z <- svd(Z)
        wsupraj <- svd.z$u[, 1]
        wsuprak <- svd.z$v[, 1]
        tf <- X %*% kronecker(wsuprak, wsupraj)
        Vrow <- t(Y) %*% tf
        V <- matrix(Vrow, nrow = dim(YN)[2], ncol = dim(YN)[3])
        svd.v <- svd(V)
        qsupraj <- svd.v$u[, 1]
        qsuprak <- svd.v$v[, 1]
        uf <- Y %*% kronecker(qsuprak, qsupraj)
        if (sum((uf - u)^2) < conver) {
          print(paste("component number ", f))
          print(paste("number of iterations: ", it))
          it <- max.iteration
          Tt <- cbind(Tt, tf)
          WsupraJ <- cbind(WsupraJ, wsupraj)
          WsupraK <- cbind(WsupraK, wsuprak)
          QsupraJ <- cbind(QsupraJ, qsupraj)
          QsupraK <- cbind(QsupraK, qsuprak)
          bf <- ginv(t(Tt) %*% Tt) %*% t(Tt) %*% uf
          B[1:length(bf), f] <- bf
          U <- cbind(U, uf)
          TM <- ginv(t(Tt) %*% Tt) %*% t(Tt)
          WkM <- ginv(t(WsupraK) %*% WsupraK) %*% t(WsupraK)
          WjM <- ginv(t(WsupraJ) %*% WsupraJ) %*% t(WsupraJ)
          Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
          Y <- Y - Tt %*% bf %*% t(kronecker(qsuprak, qsupraj))
          f <- f + 1
          Uf <- svd(Y)$u
          u <- Uf[, 1]
        }
        else {
          u <- uf
          it <- it + 1
        }
      }
    }  
    rownames(Tt) <- dimnames(XN)[[1]]
    rownames(WsupraJ) <- dimnames(XN)[[2]]
    rownames(WsupraK) <- dimnames(XN)[[3]]
    rownames(U) <- dimnames(YN)[[1]]
    rownames(QsupraJ) <- dimnames(YN)[[2]]
    rownames(QsupraK) <- dimnames(YN)[[3]]
    G <- array(as.vector(Gu[[2]]), c(factors, factors,factors))
    FactorsX = list(Mode1 = Tt, Mode2 = WsupraJ, Mode3 = WsupraK)
    FactorsY = list(Mode1 = U, Mode2 = QsupraJ, Mode3 = QsupraK)
    #####
    # Evaluation of the model Vars
    PLSmodel<-plsreg(predictors = X, responses=Y2, comps = factors, crosval = TRUE)
    Ypred<-PLSmodel$y.pred
    explvar<-PLSmodel$expvar
    VIP2D<-PLSmodel$VIP 
    tucorr<-PLSmodel$cor.tu  
    residuals<-PLSmodel$resid
    Q2<-PLSmodel$Q2
    Q2cum<-PLSmodel$Q2cum
    ####
    modo1 <- rowSums(FactorsX$Mode1[,c(1:factors)])
    modo2 <- rowSums(FactorsX$Mode2[,c(1:factors)])
    
    newx<-as.matrix(modo1) %*% as.matrix(t(modo2))
    #newx<-as.matrix(FactorsX$Mode1[,1]) %*% as.matrix(t(FactorsX$Mode2[,1]))
    PLSmodelNway<-plsreg(predictors = newx, responses=Y2, comps = factors, crosval = FALSE)
    VIP3D<-PLSmodelNway$VIP 
    
    ######
    result <- list(FactorsX = FactorsX, FactorsY = FactorsY, 
                   T = Tt, WsupraJ = WsupraJ, WsupraK = WsupraK, U = U, QsupraJ = QsupraJ, 
                   QsupraK = QsupraK, B = B, Gu = Gu, G = G,Ypred =Ypred,
                   explvar=explvar, VIP2D=VIP2D, VIP3D=VIP3D,tucorrelation=tucorr,residuals=residuals,
                   Q2=Q2,Q2cum=Q2cum
    )
    return(result)
 }
} 
##################################################################################################################
#Explcore

explcore<-function (G, factors, n = prod(dim(G))) 
  # Function to analize the core array of the Tucker3 model
  # You can control the number of factors to be evaluated per mode, (3) means three per mode, (2) Means two per mode
  
{
  if (factors == 3) {
    
    one.one.one<-G[1,1,1];one.one.two<-G[1,1,2];one.one.three<-G[1,1,3];one.two.one<-G[1,2,1];one.two.two<-G[1,2,2];
    one.two.three<-G[1,2,3];one.three.one<-G[1,3,1];one.three.two<-G[1,3,2];one.three.three<-G[1,3,3];two.one.one<-G[2,1,1];
    two.one.two<-G[2,1,2];two.one.three<-G[2,1,3];two.two.one<-G[2,2,1];two.two.two<-G[2,2,2];two.two.three<-G[2,2,3];
    two.three.one<-G[2,3,1];two.three.two<-G[2,3,2];two.three.three<-G[2,3,3];three.one.one<-G[3,1,1];three.one.two<-G[3,1,2];
    three.one.three<-G[3,1,3];three.two.one<-G[3,2,1];three.two.two<-G[3,2,2];three.two.three<-G[3,2,3];three.three.one<-G[3,3,1];
    three.three.two<-G[3,3,2];three.three.three<-G[3,3,3]
    
    preweight<-rbind(one.one.one,one.one.two, one.one.three, one.two.one,one.two.two,one.two.three,one.three.one,
                     one.three.two,one.three.three,two.one.one,two.one.two,two.one.three,two.two.one,two.two.two,
                     two.two.three,two.three.one,two.three.two,two.three.three,three.one.one,three.one.two,
                     three.one.three,three.two.one,three.two.two,three.two.three,three.three.one,
                     three.three.two,three.three.three)
    
    preweight<-data.frame(rownames(preweight),preweight)
    colnames(preweight)[1:2]<-c("Element", "Weight")
    weight<-preweight[with(preweight, order(-abs(preweight[,2]))),]
    weight$SS <- weight$Weight^2
    SST <- sum(weight$SS)
    weight$Expl.Var <- formatC((weight$SS/SST)*100)
    rownames(weight)<-NULL
    H <- weight
    H <- H[1:n, ]
    return(H)
  }
  if (factors == 2) {
    one.one.one<-G[1,1,1];one.one.two<-G[1,1,2];one.two.one<-G[1,2,1];one.two.two<-G[1,2,2];
    two.one.one<-G[2,1,1];two.one.two<-G[2,1,2];two.two.one<-G[2,2,1];two.two.two<-G[2,2,2]
    
    preweight<-rbind(one.one.one,one.one.two,one.two.one,one.two.two,two.one.one,
                     two.one.two,two.two.one,two.two.two
    )
    
    preweight<-data.frame(rownames(preweight),preweight)
    colnames(preweight)[1:2]<-c("Element", "Weight")
    weight<-preweight[with(preweight, order(-abs(preweight[,2]))),]
    weight$SS <- weight$Weight^2
    SST <- sum(weight$SS)
    weight$Expl.Var <- formatC((weight$SS/SST)*100)
    rownames(weight)<-NULL
    H <- weight
    H <- H[1:n, ]
    return(H)
  }
}



####################################################################################################################
# Plotting things
plotNPLSmod<-function (X, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                       cutoff = length(X), factors=2) 
  # Function to plot all necesary figures to interpret NPLS results
  # Input: Array of the best model
  # Output:
  # Plot of the Individuals (Mode 1) of block X in the selected components
  # Heatmap of the behavior of all variables in all individuals (it is turned off, but you can turn it on)
  # Plot of the Individuals (Mode 1) of block Y in the selected components
  # Joint plot of Components 1 of both blocks (X and Y) to see how related both blocks of variables are (the more linear, the better)
  # Plot of the Variables (Mode 2) of block X in the selected components (Here you can play with the cuttoff to know which are the
  #            variables with mor variance)
  # Plot of the Variables (Mode 2) of block Y in the selected components (Here you can play with the cuttoff to know which are the
  #            variables with mor variance)
  # Plot of the time (Mode 3) of block X in the selected components 
  # Plot of the time (Mode 3) of block Y in the selected components, When NPLS class 1 it makes no sense
  # Biplot of variables (Mode 2) and time (Mode 3) in block X, to see how X variables describe the behavior of time points in the subspace
  # Biplot of variables (Mode 2) and time (Mode 3) in block Y, to see how Y variables describe the behavior of time points in the subspace
  # Joint plot of Components 1 of both blocks (X and Y) to see how related both blocks of time (Mode 3) are
  #             (the more linear, the better explained)
  # Analysis of the Core Array of the model to determine explained variances and the contributions of the components in each mode
  # Here you can play with the script to obtain what you desire.

{
  library(sfsmisc)
  par(mar=c(5,4,4,2))
  
  #
  pc1 <- PCs[1]
  pc2 <- PCs[2]
  a = paste("Factors", "X", sep = "")
  label = NULL
  Factors <- X[[a]]
  x1 <- Factors[[1]][, pc1]
  y1 <- Factors[[1]][, pc2]
  x2 <- Factors[[2]][, pc1]
  y2 <- Factors[[2]][, pc2]
  x3 <- Factors[[3]][, pc1]
  y3 <- Factors[[3]][, pc2]
  
  b = paste("Factors", "Y", sep = "")
  label = NULL
  FactorsY <- X[[b]]
  x1Y <- FactorsY[[1]][, pc1]
  y1Y <- FactorsY[[1]][, pc2]
  x2Y <- FactorsY[[2]][, pc1]
  y2Y <- FactorsY[[2]][, pc2]
  x3Y <- FactorsY[[3]][, pc1]
  y3Y <- FactorsY[[3]][, pc2]
  par(mfrow=c(1,1))
  
  ######################################## Mode 1 Block X
  plot(x1, y1, type = "p", col = "black", pch=21,
       bg="dodgerblue2", 
       xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), las=1,
       main = paste("Block X Mode 1: NPLS", sep= " "), cex.main=2,
       cex.lab= 1.2,
       xlim = c(min(x1), max(x1)),
       ylim = c(min(y1), max(y1)))
       #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
       #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
  abline(v=0,h=0)
  text (x1,y1, labels=rownames(Factors$Mode1), col="dodgerblue2",cex=1, pos=3)
  ######################################## Mode 1 Block Y
  plot(x1Y, y1Y, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), las=1,
       main = paste("Block Y Mode 1: NPLS", sep= " "), cex.main=2,
       cex.lab= 1.2,
       xlim = c(min(x1Y), max(x1Y)),
       ylim = c(min(y1Y), max(y1Y)))
       #xlim = c((min(x1Y)+ (0.2*min(x1Y))), (max(x1Y)+ (0.2*max(x1Y)))),
       #ylim = c((min(y1Y)+ (0.2*min(y1Y))), (max(y1Y)+ (0.2*max(y1Y)))))
  abline(v=0,h=0)
  text (x1Y,y1Y, labels=rownames(FactorsY$Mode1), col="dodgerblue2",cex=1, pos=3)
  
  ######################################## Mode 1 Joint Biplot 
  plot(x1, x1Y, type = "p", col = "black", pch=21,bg="dodgerblue2", 
       xlab = paste("Component", pc1, "Block X", sep = " "),
       ylab = paste("Component", pc1, "Block Y", sep = " "), las=1,
       main = paste("Biplot Mode1 Blocks X & Y: NPLS", sep= " "), cex.main=2,
       cex.lab= 1.2,
       xlim = c((min(x1)), (max(x1)+ abs(0.2*max(x1)))),
       ylim = c((min(x1Y)), (max(x1Y)+ abs(0.2*max(x1Y)))))
  abline(v=0,h=0)
  text (x1,x1Y, labels=rownames(Factors$Mode1), col="dodgerblue2",cex=1, pos=3)
  ######################################## 
  plot(x2, y2, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block X Mode 2: NPLS", sep= " "), 
       xlim = c((min(x2)), (max(x2))),
       ylim = c((min(y2)), (max(y2))))
       #xlim = c((min(x2)), (max(x2)+ abs(0.2*max(x2)))),
       #ylim = c((min(y2)), (max(y2)+ abs(0.2*max(y2)))))
  
  abline(v=0,h=0)
  labels=NULL
  x<-NULL
  y<-NULL
  if (is.null(labels)) {
    labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
                   rownames(Factors$Mode3))
  }
  ABS <- abs(x2 * y2)
  r <- rank(-ABS)
  value = match(c(1:cutoff), r)
  for (n in 1:cutoff) {
    x[c(n)] <- Factors[[2]][value[c(n)], 1]
    y[c(n)] <- Factors[[2]][value[c(n)], 2]
    label[n] <- labels[[2]][value[c(n)]]
  }
  x <- x[-c((n + 1):(length(x) + 1))]
  y <- y[-c((n + 1):(length(y) + 1))]
  text(x, y, labels = label, adj = c(0, -1),col="dodgerblue2", cex=0.7)
  
  #text(x2,y2, labels= ifelse (x2 <= -0.1, names(x2), ifelse (x2>= 0.1,names(x2), "")),
  #     cex= 1, pos=1, col="dodgerblue2" )
  #text(x2,y2, labels= ifelse (y2 <= -0.1, names(y2), ifelse (y2>= 0.1, names(y2), "")),
  #     cex= 1, pos=1, col="dodgerblue2")
  ######################################## 
  plot(x2Y, y2Y, type = "p",col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block Y Mode 2: NPLS", sep= " "), 
       xlim = c((min(x2Y)), (max(x2Y)+ abs(0.2*max(x2Y)))),
       ylim = c((min(y2Y)), (max(y2Y)+ abs(0.2*max(y2Y)))))
  abline(v=0,h=0)
  labels=NULL
  x<-NULL
  y<-NULL
  if (is.null(labels)) {
    labels <- list(rownames(FactorsY$Mode1), rownames(FactorsY$Mode2),
                   rownames(FactorsY$Mode3))
  }
  ABS <- abs(x2Y * y2Y)
  r <- rank(-ABS)
  value = match(c(1:cutoff), r)
  for (n in 1:cutoff) {
    x[c(n)] <- FactorsY[[2]][value[c(n)], 1]
    y[c(n)] <- FactorsY[[2]][value[c(n)], 2]
    label[n] <- labels[[2]][value[c(n)]]
  }
  x <- x[-c((n + 1):(length(x) + 1))]
  y <- y[-c((n + 1):(length(y) + 1))]
  text(x, y, labels = label, adj = c(0, -1),col="dodgerblue2")
  
  #text(x2Y,y2Y, labels= ifelse (x2Y <= -0.1, names(x2Y), ifelse (x2Y>= 0.1,names(x2Y), "")),
  #     cex= 1, pos=1, col="dodgerblue2" )
  #text(x2Y,y2Y, labels= ifelse (y2Y <= -0.1, names(y2Y), ifelse (y2Y>= 0.1, names(y2Y), "")),
  #     cex= 1, pos=1, col="dodgerblue2")
  
  
  ######################################## 
  ###########      Mode 3     ############
  ######################################## 
  # Without variable vectors
  plot(x3, y3, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block X Mode 3: NPLS"), 
       xlim = c((min(x3)), (max(x3))),
       ylim = c((min(y3)), (max(y3))))
       #xlim = c((min(x3)), (max(x3)+ abs(0.2*max(x3)))),
       #ylim = c((min(y3)), (max(y3)+ abs(0.2*max(y3)))))
  abline(v=0,h=0)
  text (x3,y3, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1, pos=3)
  ######################################## 
  plot(x3Y, y3Y, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block Y Mode 3: NPLS"), 
       xlim = c((min(x3Y)), (max(x3Y))),
       ylim = c((min(y3Y)), (max(y3Y))))
       #xlim = c((min(x3Y)), (max(x3Y)+ abs(0.2*max(x3Y)))),
       #ylim = c((min(y3Y)), (max(y3Y)+ abs(0.2*max(y3Y)))))
  abline(v=0,h=0)
  text (x3Y,y3Y, labels=rownames(FactorsY$Mode3), col="dodgerblue2",cex=1, pos=3)
  
  # With variable vectors
  plot(x3, y3, type = "p", col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1, sep = " "), las=1,
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Biplot Block X Modes 2 and 3: NPLS"), 
       xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
  abline(v=0,h=0)
  text (x3,y3, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1, pos=3)
  arrows(0,0,x2,y2, pch=16, lwd=2,col = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"orange2", "transparent") )
  labels<-NULL
  x<-NULL
  y<-NULL
  if (is.null(labels)) {
    labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
                   rownames(Factors$Mode3))
  }
  ABS <- abs(x2 * y2)
  r <- rank(-ABS)
  value = match(c(1:cutoff), r)
  for (n in 1:cutoff) {
    x[c(n)] <- Factors[[2]][value[c(n)], 1]
    y[c(n)] <- Factors[[2]][value[c(n)], 2]
    label[n] <- labels[[2]][value[c(n)]]
  }
  x <- x[-c((n + 1):(length(x) + 1))]
  y <- y[-c((n + 1):(length(y) + 1))]
  text(x, y, labels = label, adj = c(0, -1),col="orange2", cex=0.7)
  
  #text (x2,y2, labels= ifelse ((abs(x2)>0.05 | abs(y2)>0.05),names(x2), ""),
  #      col = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"orange2", "transparent"),    
  #      cex=1.5, pos=1)
  
  ######################################## 
  plot(x3Y, y3Y, type = "p",  col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1, sep = " "), las=1,
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Biplot Block Y Modes 2 and 3: NPLS"), 
       xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
  abline(v=0,h=0)
  text (x3Y,y3Y, labels=rownames(FactorsY$Mode3), col="dodgerblue2",cex=1, pos=3)
  arrows(0,0,x2Y,y2Y, pch=16, lwd=2,col = ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),"orange2", "transparent") )
  labels=NULL
  x<-NULL
  y<-NULL
  if (is.null(labels)) {
    labels <- list(rownames(FactorsY$Mode1), rownames(FactorsY$Mode2),
                   rownames(FactorsY$Mode3))
  }
  ABS <- abs(x2Y * y2Y)
  r <- rank(-ABS)
  value = match(c(1:cutoff), r)
  for (n in 1:cutoff) {
    x[c(n)] <- FactorsY[[2]][value[c(n)], 1]
    y[c(n)] <- FactorsY[[2]][value[c(n)], 2]
    label[n] <- labels[[2]][value[c(n)]]
  }
  x <- x[-c((n + 1):(length(x) + 1))]
  y <- y[-c((n + 1):(length(y) + 1))]
  text(x, y, labels = label, adj = c(0, -1),col="orange2", cex=0.7)
  
  #text (x2Y,y2Y, labels= ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),names(x2Y), ""),
  #      col = ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),"orange2", "transparent"),    
  #      cex=1, pos=1)

  #    Joint Mode 3 ####################################### 
  plot(x3, x3Y, type = "p", col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1,"Block X", sep = " "), 
       ylab = paste("Component", pc1, "Block Y",sep = " "), las =1,
       main = paste("Joint plot: NPLS Mode 3"), 
       xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
  abline(v=0,h=0)
  abline(a=0,b=1,col= "gray47", lwd=1)
  text (x3,x3Y, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1, pos=3)
  
  ############################################################
  ex <- as.matrix(explcore(X$G, factors, n=5))
  plot.new()
  mtext(paste("Expl.var of X: ", round(sum(X$explvar[c(1:factors),1])*100,digits=2),"%"))
  mtext(paste("Expl.var of Y by X: ", round(sum(X$explvar[c(1:factors),3])*100, digits=2),"%"),line= -2)
  #mtext(paste("model expl.var: ", sum(as.numeric(ex[c(1:2),4]))))
  #mtext(paste("model expl.var: ", formatC(X$expl.var)))
  mtext(paste(colnames(ex), collapse = "                "), line = -4)
  par(cex = 0.8)
  for (i in 1:nrow(ex)) {
    mtext(paste(as.vector(ex[i, ]), collapse = "       "), line = -(i + 
                                                                      3) * 2, adj = 0)
  }
  
  }
  #######################################

# NPLSDAmodel

NPLSDAmod<- function (XN, YN, outcome.Y=NULL, factors=3, COMP= c(2,2,2), conver = 1e-16, max.iteration = 10000, centering=2) 
{
  # Function to perform the NPLS-DA analysis from arrays nxpxq where n are individuals, p are variables and
  # q are time points
  # You can input the array with missing values, or also the Tucker3 array already processed
  # You can control the convergence (conver) value and the number of iterations (max.iter)
  # You can center the data by any of the three modes or not center at all
  # Also you can control the number of components per mode to be retain in the Tucker3 model to get the most quantity
  # of information for the NPLS-DA (COMP). This has to be based on the results of the bestfittedmodel function 
  # You can control the number of factors to be retained per mode in the NPLS-DA (this is different to 
  # the number of components per mode that best fit the Tucker3 model (COMP) 
  # Finally you can perform Comparative NPLS-DA between two datasets to explain response variables (Y)
  # with predictors variables (X)
  # Return: Twenty one useful parameters
  # 1.- FactorsX: The Scores of the individuals (Mode1) (wsuprai), the loadings of the Variables (Mode2) (wsupraj)
  #               and the loadings of the Time (Mode3) (wsuprak) that you decided to retain in the predictors matrix (X)
  # 2.- FactorsY: The Scores of the individuals (Mode1) (qsuprai), the loadings of the Variables (Mode2) (qsupraj)
  #               and the loadings of the Time (Mode3) (qsuprak) that you decided to retain in the response matrix (Y)
  # 3.- T: Scores matrix of the predictors matrix (X)
  # 4.- WsupraJ: Loadings matrix of the variables components of the predictors matrix (X)
  # 5.- WsupraK: Loadings matrix of the Time component of the predictors matrix (X)
  # 6.- U: Scores matrix of the response matrix (Y)
  # 7.- QsupraJ: Loadings matrix of the variables components of the response matrix (Y)
  # 8.- QsupraK: Loadings matrix of the Time component of the response matrix (Y)
  # 9.- B: Matrix of rfactors. These are the regression coefficients of the model that correlates the factors of X with Y
  # 10.- Gu: Deconvoluted Core array
  # 11.- G: Core matrix retaining the the correlations among the components of the modes 
  # 12.- Ypred: Matrix of predicted values of response (Y) matrix by predictors (X) matrix. 
  # 13.- explvar: Table of explained variance of X and explained variance of Y by X
  # 14.- VIP2D: Variable importance for the projection as a measure of how much a variable explains the bidimensional model
  #             in terms of variable in every timepoint studied
  # 15.- VIP3D: Variable importance for the projection as a measure of how much a variable explains the bidimensional model
  #             in terrms of variables, no matter the timeoint
  # 16.- tucorrelation: t-u correlations
  # 17.- residuals: Table of residuals
  # 18.- NPLSDAvariates: NPLSDA scores for all studied blocks and all components evaluated
  # 19.- NPLSDAloadings: NPLSDA loadings for all studied blocks and all components evaluated
  # 20.- NPLSDAexplVar: NPLSDA explained variances for all studied blocks and all components evaluated
  # 21.- Design: The design of the analysis
  
  require (plsdepot)
  require (mixOmics)
  require (MASS)
  require(abind)
  
  Ydim<-dim(YN)
  if (Ydim[2]>1){   # This part performs the centering in matrix X and Y
    if (is.null(outcome.Y)){   
      stop("\nFor NPLS-DA regression class 2 you have to define the outcome design")
    }
    print ("Performing NPLS-DA regression class 2")
    if (any(is.na(XN))) {
      print("Missing values are taken care of")
      XN <- Imputemethod(XN, COMP)
      #matrizXnoNA<-XN
      #save(matrizXnoNA, file ="matrizXnoNA.RData")
    }
    if (any(is.na(YN))) {
      print("Missing values are taken care of")
      YN<- Imputemethod(YN, COMP)
      #matrizYnoNA<-YN
      #save(matrizYnoNA, file ="matrizYnoNA.RData")
    }
    if (centering==0) {
      xdim <- dim(XN)
      A3D<-XN
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      # MatrixY
      ydim <- dim(YN)
      B3D<-YN
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
    }
    
    if (centering==1) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(1, 2, 3)), xdim[1], 
                  xdim[2] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(1, 2, 3)]), 
                 perm = c(1, 2, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporind <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporind)<-individuos
      colnames(Xcporind)<-variables
      Xcentrado<-Xcporind
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      
      # MatrixY
      ydim <- dim(YN)
      #length(xdim) == 3L
      individuosY<-rownames(YN)
      variablesY<-colnames(YN)
      Y <- matrix(aperm(YN, perm = c(1, 2, 3)), ydim[1], 
                  ydim[2] * ydim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      Y<-apply(Y, 2, as.numeric)
      #X[1:10,1:10]
      ymns <- rowMeans(Y, na.rm = T)
      ydim2<-dim(Y)
      YY<- matrix(ymns, ydim2[1], ydim2[2])
      #Y[1:10,1:10]
      
      Y <- Y - YY
      #dim(X)
      Y <- aperm(array(Y, dim = ydim[c(1, 2, 3)]), 
                 perm = c(1, 2, 3))
      dim(Y)
      ydim3 <- dim(Y)
      Ycporind <- array(data = Y, dim = c(ydim3[1],ydim3[2],ydim3[3]),dimnames = list(NULL, NULL, colnames(Y[1,,])))
      rownames(Ycporind)<-individuosY
      colnames(Ycporind)<-variablesY
      Ycentrado<-Ycporind
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
      
    }
    
    if (centering==2) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(2, 1, 3)), xdim[2], 
                  xdim[1] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(2, 1, 3)]), 
                 perm = c(2, 1, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporvar <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporvar)<-individuos
      colnames(Xcporvar)<-variables
      Xcentrado<-Xcporvar
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[2], xdimA3D[1] * xdimA3D[3])
      
      
      #matrixY
      ydim <- dim(YN)
      #length(xdim) == 3L
      individuosY<-rownames(YN)
      variablesY<-colnames(YN)
      Y <- matrix(aperm(YN, perm = c(2, 1, 3)), ydim[2], 
                  ydim[1] * ydim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      Y<-apply(Y, 2, as.numeric)
      #X[1:10,1:10]
      ymns <- rowMeans(Y, na.rm = T)
      ydim2<-dim(Y)
      YY<- matrix(ymns, ydim2[1], ydim2[2])
      #Y[1:10,1:10]
      
      Y <- Y - YY
      #dim(X)
      Y <- aperm(array(Y, dim = ydim[c(2, 1, 3)]), 
                 perm = c(2, 1, 3))
      #dim(Y)
      ydim3 <- dim(Y)
      Ycporvar <- array(data = Y, dim = c(ydim3[1],ydim3[2],ydim3[3]),dimnames = list(NULL, NULL, colnames(Y[1,,])))
      rownames(Ycporvar)<-individuosY
      colnames(Ycporvar)<-variablesY
      Ycentrado<-Ycporvar
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[2], ydimB3D[1] * ydimB3D[3])
      
    }
    if (centering==3) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(3, 1, 2)), xdim[3], 
                  xdim[1] * xdim[2])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(3, 2, 1)]), 
                 perm = c(3, 2, 1))
      dim(X)
      xdim3 <- dim(X)
      Xcportime <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcportime)<-individuos
      colnames(Xcportime)<-variables
      Xcentrado<-Xcportime
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[3], xdimA3D[1] * xdimA3D[2])
      
      
      #MatrixY
      ydim <- dim(YN)
      #length(xdim) == 3L
      individuosY<-rownames(YN)
      variablesY<-colnames(YN)
      Y <- matrix(aperm(YN, perm = c(3, 1, 2)), ydim[3], 
                  ydim[1] * ydim[2])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      Y<-apply(Y, 2, as.numeric)
      #X[1:10,1:10]
      ymns <- rowMeans(Y, na.rm = T)
      ydim2<-dim(Y)
      YY<- matrix(ymns, ydim2[1], ydim2[2])
      #Y[1:10,1:10]
      
      Y <- Y - YY
      #dim(X)
      Y <- aperm(array(Y, dim = ydim[c(3, 2, 1)]), 
                 perm = c(3, 2, 1))
      dim(Y)
      ydim3 <- dim(Y)
      Ycportime <- array(data = Y, dim = c(ydim3[1],ydim3[2],ydim3[3]),dimnames = list(NULL, NULL, colnames(Y[1,,])))
      rownames(Ycportime)<-individuosY
      colnames(Ycportime)<-variablesY
      Ycentrado<-Ycportime
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[3], ydimB3D[1] * ydimB3D[2])
    }
    
    xdims<-dim(A3D)
    ydims<-dim(B3D)
    Tt <- U <- WsupraJ <- WsupraK <- QsupraJ <- QsupraK <- NULL
    B <- G <- matrix(0, ncol = factors, nrow = factors)
    Gu <- vector("list", factors)
    X <- matrix(A3D, xdims[1], xdims[2] * xdims[3])
    ordX <- length(dim(A3D))
    Y <- matrix(B3D, ydims[1], ydims[2] * ydims[3])
    ordY <- length(dim(B3D))
    Uf <- svd(Y)$u
    u <- Uf[, 1]
    f = 1
    #
    for (f in 1:factors) {
      it = 1
      while (it < max.iteration) {
        tX<-t(X)
        #tX<-tX[1:313]
        Zrow <- tX %*% u
        Z <- matrix(Zrow, nrow = dim(XN)[2], ncol = dim(XN)[3])
        svd.z <- svd(Z)
        wsupraj <- svd.z$u[, 1]
        wsuprak <- svd.z$v[, 1]
        tf <- X %*% kronecker(wsuprak, wsupraj)
        Vrow <- t(Y) %*% tf
        V <- matrix(Vrow, nrow = dim(YN)[2], ncol = dim(YN)[3])
        svd.v <- svd(V)
        qsupraj <- svd.v$u[, 1]
        qsuprak <- svd.v$v[, 1]
        uf <- Y %*% kronecker(qsuprak, qsupraj)
        if (sum((uf - u)^2) < conver) {
          print(paste("component number ", f))
          print(paste("number of iterations: ", it))
          it <- max.iteration
          Tt <- cbind(Tt, tf)
          WsupraJ <- cbind(WsupraJ, wsupraj)
          WsupraK <- cbind(WsupraK, wsuprak)
          QsupraJ <- cbind(QsupraJ, qsupraj)
          QsupraK <- cbind(QsupraK, qsuprak)
          bf <- ginv(t(Tt) %*% Tt) %*% t(Tt) %*% uf
          B[1:length(bf), f] <- bf
          U <- cbind(U, uf)
          TM <- ginv(t(Tt) %*% Tt) %*% t(Tt)
          WkM <- ginv(t(WsupraK) %*% WsupraK) %*% t(WsupraK)
          WjM <- ginv(t(WsupraJ) %*% WsupraJ) %*% t(WsupraJ)
          Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
          Y <- Y - Tt %*% bf %*% t(kronecker(qsuprak, qsupraj))
          f <- f + 1
          Uf <- svd(Y)$u
          u <- Uf[, 1]
        }
        else {
          u <- uf
          it <- it + 1
        }
      }
    }  
    rownames(Tt) <- dimnames(XN)[[1]]
    rownames(WsupraJ) <- dimnames(XN)[[2]]
    rownames(WsupraK) <- dimnames(XN)[[3]]
    rownames(U) <- dimnames(YN)[[1]]
    rownames(QsupraJ) <- dimnames(YN)[[2]]
    rownames(QsupraK) <- dimnames(YN)[[3]]
    G <- array(as.vector(Gu[[2]]), c(factors, factors,factors))
    FactorsX = list(Mode1 = Tt, Mode2 = WsupraJ, Mode3 = WsupraK)
    FactorsY = list(Mode1 = U, Mode2 = QsupraJ, Mode3 = QsupraK)
    #####
    # Evaluation of the model Vars
    PLSmodel<-plsreg(predictors = X, responses=Y, comps = factors, crosval = TRUE)
    Ypred<-PLSmodel$y.pred
    Explvar2D<-PLSmodel$expvar
    VIP2D<-PLSmodel$VIP 
    tucorr<-PLSmodel$cor.tu  
    residuals<-PLSmodel$resid
    Q22D<-PLSmodel$Q2
    Q2cum2D<-PLSmodel$Q2cum
    
    ################################################################################################################
    ####VIP3D model 1: Calculating the mean of all times per each t
    
    explvar3d<-list()
    explvar3dtita<-NULL
    #
    VIP3Dsas<-NULL
    vip3d<-list()
    #
    Q23Dsas<-NULL
    Q23D<-list()
    #
    Q23Dcumsas<-NULL
    Q23Dcum<-list()
    
    for (i in 1:dim(A3D)[3]) {
    PLSmodel<-plsreg(predictors = A3D[,,i], responses=B3D[,,i], comps = factors, crosval = TRUE)
    explvar3dtita<-PLSmodel$expvar
    nameexplvar <- paste('item:',i,sep='')
    explvar3dsas<-list(explvar3dtita)
    explvar3d[[nameexplvar]]<-explvar3dsas
    #
    VIP3Dtita<-PLSmodel$VIP 
    namevip <- paste('item:',i,sep='')
    vip3dsas<-list(VIP3Dtita)
    vip3d[[namevip]]<-vip3dsas
    #
    #
    Q23Dtita<-PLSmodel$Q2
    nameQ23D <- paste('item:',i,sep='')
    Q23Dsas<-list(Q23Dtita)
    Q23D[[nameQ23D]]<-Q23Dsas
    #
    Q23Dcumtita<-PLSmodel$Q2cum
    nameQ23Dcum <- paste('item:',i,sep='')
    Q23Dcumsas<-list(Q23Dcumtita)
    Q23Dcum[[nameQ23Dcum]]<-Q23Dcumsas
    
    
    }
    R2Xmean<-NULL
    R2Ymean<-NULL
    for (i in 1:length(explvar3d)) {
      a=as.matrix(explvar3d[i][[1]][[1]][,1])
      R2Xmean<-cbind(R2Xmean,a)
      b=as.matrix(explvar3d[i][[1]][[1]][,3])
      R2Ymean<-cbind(R2Ymean,b)
    }
    R2Xmean3D<-as.matrix(rowMeans(R2Xmean)); colnames(R2Xmean3D)<-"R2Xmean3D"
    R2Ymean3D<-as.matrix(rowMeans(R2Ymean)); colnames(R2Ymean3D)<-"R2Ymean3D"
    Explvar3D<-cbind(R2Xmean3D,R2Ymean3D)
    ###############################################
    # Calculating Q2mean3D
    Q2mean<-NULL
    for (i in 1:length(Q23D)) {
      a=as.matrix(Q23D[i][[1]][[1]][,1])
      Q2mean<-cbind(Q2mean,a)
    }
    Q2mean3D<-as.matrix(rowMeans(Q2mean)); colnames(Q2mean3D)<-"Q2mean3D"
    ###############################################
    # Calculating Q2cummean3D
    Q2cummean<-NULL
    for (i in 1:length(Q23Dcum)) {
      a=as.matrix(Q23Dcum[i][[1]][[1]][,1])
      Q2cummean<-cbind(Q2cummean,a)
    }
    Q2cummean3D<-as.matrix(rowMeans(Q2cummean)); colnames(Q2cummean3D)<-"Q2cummean3D"
    
    ####################################################################################
    
    t1mean<-t2mean<-t3mean<-t4mean<-NULL
   if(factors==2) {
    for (i in 1:factors) {
      t1=as.matrix(vip3d[i][[1]][[1]][,1])
      t1mean<-cbind(t1mean,t1)
      t2=as.matrix(vip3d[i][[1]][[1]][,2])
      t2mean<-cbind(t2mean,t2)
      #t3=as.matrix(vip3d[i][[1]][[1]][,3])
      #t3mean<-cbind(t3mean,t3)
    }
    t1mean3D<-as.matrix(rowMeans(t1mean)); colnames(t1mean3D)<-"t1mean3D"
    t2mean3D<-as.matrix(rowMeans(t2mean)); colnames(t2mean3D)<-"t2mean3D"
    #t3mean3D<-as.matrix(rowMeans(t3mean)); colnames(t3mean3D)<-"t3mean3D"
    VIP3Dmodel1<-cbind(t1mean3D,t2mean3D)#,t3mean3D)
  }
    if(factors==3) {
      for (i in 1:factors) {
        t1=as.matrix(vip3d[i][[1]][[1]][,1])
        t1mean<-cbind(t1mean,t1)
        t2=as.matrix(vip3d[i][[1]][[1]][,2])
        t2mean<-cbind(t2mean,t2)
        t3=as.matrix(vip3d[i][[1]][[1]][,3])
        t3mean<-cbind(t3mean,t3)
      }
      t1mean3D<-as.matrix(rowMeans(t1mean)); colnames(t1mean3D)<-"t1mean3D"
      t2mean3D<-as.matrix(rowMeans(t2mean)); colnames(t2mean3D)<-"t2mean3D"
      t3mean3D<-as.matrix(rowMeans(t3mean)); colnames(t3mean3D)<-"t3mean3D"
      VIP3Dmodel1<-cbind(t1mean3D,t2mean3D,t3mean3D)
    } 
    if(factors==4) {
      for (i in 1:factors) {
        t1=as.matrix(vip3d[i][[1]][[1]][,1])
        t1mean<-cbind(t1mean,t1)
        t2=as.matrix(vip3d[i][[1]][[1]][,2])
        t2mean<-cbind(t2mean,t2)
        t3=as.matrix(vip3d[i][[1]][[1]][,3])
        t3mean<-cbind(t3mean,t3)
        t4=as.matrix(vip3d[i][[1]][[1]][,4])
        t4mean<-cbind(t4mean,t4)
      }
      t1mean3D<-as.matrix(rowMeans(t1mean)); colnames(t1mean3D)<-"t1mean3D"
      t2mean3D<-as.matrix(rowMeans(t2mean)); colnames(t2mean3D)<-"t2mean3D"
      t3mean3D<-as.matrix(rowMeans(t3mean)); colnames(t3mean3D)<-"t3mean3D"
      t4mean3D<-as.matrix(rowMeans(t4mean)); colnames(t4mean3D)<-"t4mean3D"
      VIP3Dmodel1<-cbind(t1mean3D,t2mean3D,t3mean3D,t4mean3D)
    } 
    ##########################################################################################################
    ####VIP3D model 2: Calculating VIP per time, sum all t's generating a VIPtable per timepoint
    
    time<-NULL
    timetable<-NULL
    #t1mean<-t2mean<-t3mean<-NULL
    for (i in 1:length(vip3d)) {
      time<- as.matrix(rowSums(vip3d[i][[1]][[1]]))
      timetable<-cbind(timetable,time)
    }
    #timetable
    VIP3Dmodel2<-timetable
    #########
    xdimA3D<-dim(A3D)
    A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
    xdimB3D<-dim(B3D)
    B<-matrix(B3D, xdimB3D[1], xdimB3D[2] * xdimB3D[3])
    plsdaforplotting<-block.plsda (X=list(Block.X=A, Block.Y=B), Y = outcome.Y, ncomp=factors)
    NPLSDAvariates<-plsdaforplotting$variates       # *Para graficar plotIndiv
    NPLSDAloadings<-plsdaforplotting$loadings       # * Para graficar plotloadings
    ##########
    tabla<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,2))
    tabla4<-tabla[,c(4,1,2)]
    transformedNPLSDAloadings<-tabla4
    
    ##########
    NPvariates<-list()
    NPvariatestita<-NULL
    NPloadingstita<-NULL
    NPloadings<-list()
    
    for (i in 1:dim(A3D)[3]) {
      plsdaforplotting<-block.plsda (X=list(Block.X=A3D[,,i], Block.Y=B3D[,,i]), Y = outcome.Y, ncomp=factors)
      NPvariatestita<-plsdaforplotting$variates       # *Para graficar plotIndiv
      namevariatestita<-paste("item:",i,sep = "")
      variatestitatemp<-list(NPvariatestita)
      NPvariates[[namevariatestita]]<-variatestitatemp
      
      NPloadingstita<-plsdaforplotting$loadings       # * Para graficar plotloadings
      nameloadingdtita<-paste("item:",i,sep = "")
      loadingstitatemp<-list(NPloadingstita)
      NPloadings[[nameloadingdtita]]<-loadingstitatemp
    }
    NPLSDAvariatesperMode3<-NPvariates
    NPLSDAloadingsperMode3<-NPloadings
    ###### Consensus configuration
    
    
    Xcomp1variatesmean<-Xcomp2variatesmean<-Xcomp3variatesmean<-Xcomp4variatesmean<-NULL
    Ycomp1variatesmean<-Ycomp2variatesmean<-Ycomp3variatesmean<-Ycomp4variatesmean<-NULL
    if(factors==2){
    for (i in 1:factors) {
      a=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,1])
      Xcomp1variatesmean<-cbind(Xcomp1variatesmean,a)
      b=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,2])
      Xcomp2variatesmean<-cbind(Xcomp2variatesmean,b)
      #c=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,3])
      #Xcomp3variatesmean<-cbind(Xcomp3variatesmean,c)
      
      d=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,1])
      Ycomp1variatesmean<-cbind(Ycomp1variatesmean,d)
      e=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,2])
      Ycomp2variatesmean<-cbind(Ycomp2variatesmean,e)
      #f=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,3])
      #Ycomp3variatesmean<-cbind(Ycomp3variatesmean,f)
      
      }
    comp1mean3D<-as.matrix(rowMeans(Xcomp1variatesmean)); colnames(comp1mean3D)<-"Comp1mean3D"
    comp2mean3D<-as.matrix(rowMeans(Xcomp2variatesmean)); colnames(comp2mean3D)<-"Comp2mean3D"
    #comp3mean3D<-as.matrix(rowMeans(Xcomp3variatesmean)); colnames(comp3mean3D)<-"Comp3mean3D"
    Block.XConsensus<-cbind(comp1mean3D,comp2mean3D)#,comp3mean3D)
    
    Ycomp1mean3D<-as.matrix(rowMeans(Ycomp1variatesmean)); colnames(Ycomp1mean3D)<-"Comp1mean3D"
    Ycomp2mean3D<-as.matrix(rowMeans(Ycomp2variatesmean)); colnames(Ycomp2mean3D)<-"Comp2mean3D"
    #Ycomp3mean3D<-as.matrix(rowMeans(Ycomp3variatesmean)); colnames(Ycomp3mean3D)<-"Comp3mean3D"
    Block.YConsensus<-cbind(Ycomp1mean3D,Ycomp2mean3D)#,Ycomp3mean3D)
    
    NPLSDAConsensusvariates<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    if(factors==3){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1variatesmean<-cbind(Xcomp1variatesmean,a)
        b=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2variatesmean<-cbind(Xcomp2variatesmean,b)
        c=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3variatesmean<-cbind(Xcomp3variatesmean,c)
        
        d=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1variatesmean<-cbind(Ycomp1variatesmean,d)
        e=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2variatesmean<-cbind(Ycomp2variatesmean,e)
        f=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3variatesmean<-cbind(Ycomp3variatesmean,f)
        
      }
      comp1mean3D<-as.matrix(rowMeans(Xcomp1variatesmean)); colnames(comp1mean3D)<-"Comp1mean3D"
      comp2mean3D<-as.matrix(rowMeans(Xcomp2variatesmean)); colnames(comp2mean3D)<-"Comp2mean3D"
      comp3mean3D<-as.matrix(rowMeans(Xcomp3variatesmean)); colnames(comp3mean3D)<-"Comp3mean3D"
      Block.XConsensus<-cbind(comp1mean3D,comp2mean3D,comp3mean3D)
      
      Ycomp1mean3D<-as.matrix(rowMeans(Ycomp1variatesmean)); colnames(Ycomp1mean3D)<-"Comp1mean3D"
      Ycomp2mean3D<-as.matrix(rowMeans(Ycomp2variatesmean)); colnames(Ycomp2mean3D)<-"Comp2mean3D"
      Ycomp3mean3D<-as.matrix(rowMeans(Ycomp3variatesmean)); colnames(Ycomp3mean3D)<-"Comp3mean3D"
      Block.YConsensus<-cbind(Ycomp1mean3D,Ycomp2mean3D,Ycomp3mean3D)
      
      NPLSDAConsensusvariates<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    if(factors==4){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1variatesmean<-cbind(Xcomp1variatesmean,a)
        b=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2variatesmean<-cbind(Xcomp2variatesmean,b)
        c=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3variatesmean<-cbind(Xcomp3variatesmean,c)
        d=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,4])
        Xcomp4variatesmean<-cbind(Xcomp4variatesmean,d)
        
        e=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1variatesmean<-cbind(Ycomp1variatesmean,e)
        f=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2variatesmean<-cbind(Ycomp2variatesmean,f)
        g=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3variatesmean<-cbind(Ycomp3variatesmean,g)
        h=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,4])
        Ycomp4variatesmean<-cbind(Ycomp4variatesmean,h)
      }
      comp1mean3D<-as.matrix(rowMeans(Xcomp1variatesmean)); colnames(comp1mean3D)<-"Comp1mean3D"
      comp2mean3D<-as.matrix(rowMeans(Xcomp2variatesmean)); colnames(comp2mean3D)<-"Comp2mean3D"
      comp3mean3D<-as.matrix(rowMeans(Xcomp3variatesmean)); colnames(comp3mean3D)<-"Comp3mean3D"
      comp4mean3D<-as.matrix(rowMeans(Xcomp4variatesmean)); colnames(comp4mean3D)<-"Comp4mean3D"
      Block.XConsensus<-cbind(comp1mean3D,comp2mean3D,comp3mean3D,comp4mean3D)
      
      Ycomp1mean3D<-as.matrix(rowMeans(Ycomp1variatesmean)); colnames(Ycomp1mean3D)<-"Comp1mean3D"
      Ycomp2mean3D<-as.matrix(rowMeans(Ycomp2variatesmean)); colnames(Ycomp2mean3D)<-"Comp2mean3D"
      Ycomp3mean3D<-as.matrix(rowMeans(Ycomp3variatesmean)); colnames(Ycomp3mean3D)<-"Comp3mean3D"
      Ycomp4mean3D<-as.matrix(rowMeans(Ycomp4variatesmean)); colnames(Ycomp4mean3D)<-"Comp4mean3D"
      Block.YConsensus<-cbind(Ycomp1mean3D,Ycomp2mean3D,Ycomp3mean3D,Ycomp4mean3D)
      
      NPLSDAConsensusvariates<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    
    
    
    
    
    
    
    
    ### CONSENSUS LOADINGS
    Xcomp1loadingsmean<-Xcomp2loadingsmean<-Xcomp3loadingsmean<-Xcomp4loadingsmean<-NULL
    Ycomp1loadingsmean<-Ycomp2loadingsmean<-Ycomp3loadingsmean<-Ycomp4loadingsmean<-NULL
    
    if (factors==2){
    for (i in 1:factors) {
      a=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,1])
      Xcomp1loadingsmean<-cbind(Xcomp1loadingsmean,a)
      b=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,2])
      Xcomp2loadingsmean<-cbind(Xcomp2loadingsmean,b)
      #c=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,3])
      #Xcomp3loadingsmean<-cbind(Xcomp3loadingsmean,c)
      
      d=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,1])
      Ycomp1loadingsmean<-cbind(Ycomp1loadingsmean,d)
      e=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,2])
      Ycomp2loadingsmean<-cbind(Ycomp2loadingsmean,e)
      #f=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,3])
      #Ycomp3loadingsmean<-cbind(Ycomp3loadingsmean,f)
      
    }
    comp1meanloadings3D<-as.matrix(rowMeans(Xcomp1loadingsmean)); colnames(comp1meanloadings3D)<-"Comp1loadingsmean3D"
    comp2meanloadings3D<-as.matrix(rowMeans(Xcomp2loadingsmean)); colnames(comp2meanloadings3D)<-"Comp2loadingsmean3D"
    #comp3meanloadings3D<-as.matrix(rowMeans(Xcomp3loadingsmean)); colnames(comp3meanloadings3D)<-"Comp3loadingsmean3D"
    Block.XConsensus<-cbind(comp1meanloadings3D,comp2meanloadings3D)#,comp3meanloadings3D)
    
    Ycomp1meanloadings3D<-as.matrix(rowMeans(Ycomp1loadingsmean)); colnames(Ycomp1meanloadings3D)<-"YComp1loadingsmean3D"
    Ycomp2meanloadings3D<-as.matrix(rowMeans(Ycomp2loadingsmean)); colnames(Ycomp2meanloadings3D)<-"YComp2loadingsmean3D"
    #Ycomp3meanloadings3D<-as.matrix(rowMeans(Ycomp3loadingsmean)); colnames(Ycomp3meanloadings3D)<-"YComp3loadingsmean3D"
    Block.YConsensus<-cbind(Ycomp1meanloadings3D,Ycomp2meanloadings3D)#,Ycomp3meanloadings3D)
    
    NPLSDAConsensusloadings<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    if (factors==3){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1loadingsmean<-cbind(Xcomp1loadingsmean,a)
        b=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2loadingsmean<-cbind(Xcomp2loadingsmean,b)
        c=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3loadingsmean<-cbind(Xcomp3loadingsmean,c)
        
        d=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1loadingsmean<-cbind(Ycomp1loadingsmean,d)
        e=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2loadingsmean<-cbind(Ycomp2loadingsmean,e)
        f=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3loadingsmean<-cbind(Ycomp3loadingsmean,f)
        
      }
      comp1meanloadings3D<-as.matrix(rowMeans(Xcomp1loadingsmean)); colnames(comp1meanloadings3D)<-"Comp1loadingsmean3D"
      comp2meanloadings3D<-as.matrix(rowMeans(Xcomp2loadingsmean)); colnames(comp2meanloadings3D)<-"Comp2loadingsmean3D"
      comp3meanloadings3D<-as.matrix(rowMeans(Xcomp3loadingsmean)); colnames(comp3meanloadings3D)<-"Comp3loadingsmean3D"
      Block.XConsensus<-cbind(comp1meanloadings3D,comp2meanloadings3D,comp3meanloadings3D)
      
      Ycomp1meanloadings3D<-as.matrix(rowMeans(Ycomp1loadingsmean)); colnames(Ycomp1meanloadings3D)<-"YComp1loadingsmean3D"
      Ycomp2meanloadings3D<-as.matrix(rowMeans(Ycomp2loadingsmean)); colnames(Ycomp2meanloadings3D)<-"YComp2loadingsmean3D"
      Ycomp3meanloadings3D<-as.matrix(rowMeans(Ycomp3loadingsmean)); colnames(Ycomp3meanloadings3D)<-"YComp3loadingsmean3D"
      Block.YConsensus<-cbind(Ycomp1meanloadings3D,Ycomp2meanloadings3D,Ycomp3meanloadings3D)
      
      NPLSDAConsensusloadings<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    if (factors==4){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1loadingsmean<-cbind(Xcomp1loadingsmean,a)
        b=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2loadingsmean<-cbind(Xcomp2loadingsmean,b)
        c=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3loadingsmean<-cbind(Xcomp3loadingsmean,c)
        d=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,4])
        Xcomp4loadingsmean<-cbind(Xcomp4loadingsmean,d)
        
        e=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1loadingsmean<-cbind(Ycomp1loadingsmean,e)
        f=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2loadingsmean<-cbind(Ycomp2loadingsmean,f)
        g=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3loadingsmean<-cbind(Ycomp3loadingsmean,g)
        h=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,4])
        Ycomp4loadingsmean<-cbind(Ycomp4loadingsmean,h)
        
      }
      comp1meanloadings3D<-as.matrix(rowMeans(Xcomp1loadingsmean)); colnames(comp1meanloadings3D)<-"Comp1loadingsmean3D"
      comp2meanloadings3D<-as.matrix(rowMeans(Xcomp2loadingsmean)); colnames(comp2meanloadings3D)<-"Comp2loadingsmean3D"
      comp3meanloadings3D<-as.matrix(rowMeans(Xcomp3loadingsmean)); colnames(comp3meanloadings3D)<-"Comp3loadingsmean3D"
      comp4meanloadings3D<-as.matrix(rowMeans(Xcomp4loadingsmean)); colnames(comp4meanloadings3D)<-"Comp4loadingsmean3D"
      Block.XConsensus<-cbind(comp1meanloadings3D,comp2meanloadings3D,comp3meanloadings3D,comp4meanloadings3D)
      
      Ycomp1meanloadings3D<-as.matrix(rowMeans(Ycomp1loadingsmean)); colnames(Ycomp1meanloadings3D)<-"YComp1loadingsmean3D"
      Ycomp2meanloadings3D<-as.matrix(rowMeans(Ycomp2loadingsmean)); colnames(Ycomp2meanloadings3D)<-"YComp2loadingsmean3D"
      Ycomp3meanloadings3D<-as.matrix(rowMeans(Ycomp3loadingsmean)); colnames(Ycomp3meanloadings3D)<-"YComp3loadingsmean3D"
      Ycomp4meanloadings3D<-as.matrix(rowMeans(Ycomp4loadingsmean)); colnames(Ycomp4meanloadings3D)<-"YComp4loadingsmean3D"
      Block.YConsensus<-cbind(Ycomp1meanloadings3D,Ycomp2meanloadings3D,Ycomp3meanloadings3D,Ycomp4meanloadings3D)
      
      NPLSDAConsensusloadings<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    
    
    #if (factors>2){
    #  tabla<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,2))
    #  tabla2<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,3))
    #  tabla3<-cbind(tabla,tabla2)
    #  tabla4<-tabla3[,c(4,1,2,10)]
    #  colnames(tabla4)[4]<-"z"
    #  transformedNPLSDAloadings<-tabla4
      
    #}
    #if (factors==2){
    #  tabla<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,2))
    #  tabla4<-tabla[,c(4,1,2)]
    #  transformedNPLSDAloadings<-tabla4
      
    #}
    
    DAR2blockX<-plsdaforplotting$explained_variance$Block.X
    DAR2blockY<-plsdaforplotting$explained_variance$Block.Y
    blox<-sum(Explvar2D[c(1:factors),1]) * DAR2blockX[1:factors]
    bloy<-sum(Explvar2D[c(1:factors),3]) * DAR2blockX[1:factors] #BlockX tambien porque es que tanto X explica a Y
    NPLSDAexplVar<-t(rbind(R2Block.X=blox,R2Block.Y=bloy))
    
    #pch=c(0,1)
    #cex=c(1,2)
    #plotArrow(plsdaforplotting, ind.names=TRUE, abline = TRUE, pch=pch,cex=cex) # Esta si
    #plotDiablo(plsdaforplotting)
    # plotVar(plsdaforplotting) Mejorar esta
    # plotLoadings(plsdaforplotting)  Esta si
    #network(plsdaforplotting) # Muy pesada
    #selectVar ((plsdaforplotting))
    
    ######
    result <- list(FactorsX = FactorsX, FactorsY = FactorsY, 
                   T = Tt, WsupraJ = WsupraJ, WsupraK = WsupraK, U = U, QsupraJ = QsupraJ, 
                   QsupraK = QsupraK, B = B, Gu = Gu, G = G,Ypred =Ypred,
                   Explvar2D=Explvar2D, VIP2D=VIP2D, Explvar3D=Explvar3D,VIP3Dmodel1=VIP3Dmodel1,VIP3Dmodel2=VIP3Dmodel2,
                   tucorrelation=tucorr,residuals=residuals,
                   NPLSDAvariates=NPLSDAvariates, NPLSDAloadings=NPLSDAloadings, transformedNPLSDAloadings=transformedNPLSDAloadings,
                   NPLSDAexplVar=NPLSDAexplVar,
                   NPLSDAQ22D=Q22D,
                   NPLSDAQ2cum2D=Q2cum2D,
                   NPLSDAQ2mean3D=Q2mean3D,
                   NPLSDAQ2cummean3D=Q2cummean3D,
                   NPLSDAvariatesperMode3=NPLSDAvariatesperMode3,NPLSDAloadingsperMode3=NPLSDAloadingsperMode3,
                   NPLSDAConsensusvariates=NPLSDAConsensusvariates,NPLSDAConsensusloadings=NPLSDAConsensusloadings,
                   Design=outcome.Y
                   
                   
    )
    return(result)
    
  }
    
  ###################################################
  if (Ydim[2]==1){   # Classic NPLSDA
    print ("Performing Classic NPLS-DA")
    if (any(is.na(XN))) {
      print("Missing values are taken care of")
      XN <- Imputemethod(XN, COMP)
      #matrizXnoNA<-XN
      #save(matrizXnoNA, file ="matrizXnoNA.RData")
    }
    if (any(is.na(YN))) {
      print("Missing values are taken care of")
      YN<- Imputemethod(YN, COMP)
      #matrizYnoNA<-YN
      #save(matrizYnoNA, file ="matrizYnoNA.RData")
    }
    if (centering==0) {
      xdim <- dim(XN)
      Xcentrado<-XN
      A3D<-Imputemethod(Xcentrado)
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      # MatrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Imputemethod(Ycentrado)
      ydimB3D<-dim(B3D)
      #B<-matrix(B3D, ydimB3D[1], ydimB3D[2] * ydimB3D[3])
    }
    
    if (centering==1) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(1, 2, 3)), xdim[1], 
                  xdim[2] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(1, 2, 3)]), 
                 perm = c(1, 2, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporind <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporind)<-individuos
      colnames(Xcporind)<-variables
      Xcentrado<-Xcporind
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
      
      # MatrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      
    }
    
    if (centering==2) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(2, 1, 3)), xdim[2], 
                  xdim[1] * xdim[3])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(2, 1, 3)]), 
                 perm = c(2, 1, 3))
      dim(X)
      xdim3 <- dim(X)
      Xcporvar <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcporvar)<-individuos
      colnames(Xcporvar)<-variables
      Xcentrado<-Xcporvar
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[2], xdimA3D[1] * xdimA3D[3])
      
      
      #matrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
    }
    if (centering==3) {
      xdim <- dim(XN)
      #length(xdim) == 3L
      individuos<-rownames(XN)
      variables<-colnames(XN)
      X <- matrix(aperm(XN, perm = c(3, 1, 2)), xdim[3], 
                  xdim[1] * xdim[2])
      #dim(X) #2D por Col
      #X[1:10,1:10]
      X<-apply(X, 2, as.numeric)
      #X[1:10,1:10]
      xmns <- rowMeans(X, na.rm = T)
      xdim2<-dim(X)
      Y<- matrix(xmns, xdim2[1], xdim2[2])
      #Y[1:10,1:10]
      
      X <- X - Y
      #dim(X)
      X <- aperm(array(X, dim = xdim[c(3, 2, 1)]), 
                 perm = c(3, 2, 1))
      dim(X)
      xdim3 <- dim(X)
      Xcportime <- array(data = X, dim = c(xdim3[1],xdim3[2],xdim3[3]),dimnames = list(NULL, NULL, colnames(X[1,,])))
      rownames(Xcportime)<-individuos
      colnames(Xcportime)<-variables
      Xcentrado<-Xcportime
      A3D<-Xcentrado
      xdimA3D<-dim(A3D)
      #A<-matrix(A3D, xdimA3D[3], xdimA3D[1] * xdimA3D[2])
      
      
      #MatrixY
      ydim <- dim(YN)
      Ycentrado<-YN
      B3D<-Ycentrado
      ydimB3D<-dim(B3D)
      
    }
    
    xdims<-dim(A3D)
    ydims<-dim(B3D)
    B3D2 <- abind(B3D, B3D, along=2) 
    ydims2<-dim(B3D2)
    Tt <- U <- WsupraJ <- WsupraK <- QsupraJ <- QsupraK <- NULL
    B <- G <- matrix(0, ncol = factors, nrow = factors)
    Gu <- vector("list", factors)
    X <- matrix(A3D, xdims[1], xdims[2] * xdims[3])
    ordX <- length(dim(A3D))
    Y <- matrix(B3D, ydims[1], ydims[2] * ydims[3])
    Y2 <- matrix(B3D2, ydims[1], ydims[2] * ydims[3])
    ordY <- length(dim(B3D))
    Uf <- svd(Y)$u
    u <- Uf[, 1]
    f = 1
    #
    for (f in 1:factors) {
      it = 1
      while (it < max.iteration) {
        tX<-t(X)
        #tX<-tX[1:313]
        Zrow <- tX %*% u
        Z <- matrix(Zrow, nrow = dim(XN)[2], ncol = dim(XN)[3])
        svd.z <- svd(Z)
        wsupraj <- svd.z$u[, 1]
        wsuprak <- svd.z$v[, 1]
        tf <- X %*% kronecker(wsuprak, wsupraj)
        Vrow <- t(Y) %*% tf
        V <- matrix(Vrow, nrow = dim(YN)[2], ncol = dim(YN)[3])
        svd.v <- svd(V)
        qsupraj <- svd.v$u[, 1]
        qsuprak <- svd.v$v[, 1]
        uf <- Y %*% kronecker(qsuprak, qsupraj)
        if (sum((uf - u)^2) < conver) {
          print(paste("component number ", f))
          print(paste("number of iterations: ", it))
          it <- max.iteration
          Tt <- cbind(Tt, tf)
          WsupraJ <- cbind(WsupraJ, wsupraj)
          WsupraK <- cbind(WsupraK, wsuprak)
          QsupraJ <- cbind(QsupraJ, qsupraj)
          QsupraK <- cbind(QsupraK, qsuprak)
          bf <- ginv(t(Tt) %*% Tt) %*% t(Tt) %*% uf
          B[1:length(bf), f] <- bf
          U <- cbind(U, uf)
          TM <- ginv(t(Tt) %*% Tt) %*% t(Tt)
          WkM <- ginv(t(WsupraK) %*% WsupraK) %*% t(WsupraK)
          WjM <- ginv(t(WsupraJ) %*% WsupraJ) %*% t(WsupraJ)
          Gu[[f]] <- TM %*% X %*% kronecker(t(WkM), t(WjM))
          Y <- Y - Tt %*% bf %*% t(kronecker(qsuprak, qsupraj))
          f <- f + 1
          Uf <- svd(Y)$u
          u <- Uf[, 1]
        }
        else {
          u <- uf
          it <- it + 1
        }
      }
    }  
    rownames(Tt) <- dimnames(XN)[[1]]
    rownames(WsupraJ) <- dimnames(XN)[[2]]
    rownames(WsupraK) <- dimnames(XN)[[3]]
    rownames(U) <- dimnames(YN)[[1]]
    rownames(QsupraJ) <- dimnames(YN)[[2]]
    rownames(QsupraK) <- dimnames(YN)[[3]]
    G <- array(as.vector(Gu[[2]]), c(factors, factors,factors))
    FactorsX = list(Mode1 = Tt, Mode2 = WsupraJ, Mode3 = WsupraK)
    FactorsY = list(Mode1 = U, Mode2 = QsupraJ, Mode3 = QsupraK)
    
    #####
    # Evaluation of the model Vars
    PLSmodel<-plsreg(predictors = X, responses=Y2, comps = factors, crosval = TRUE)
    Ypred<-PLSmodel$y.pred
    explvar<-PLSmodel$expvar
    VIP2D<-PLSmodel$VIP 
    tucorr<-PLSmodel$cor.tu  
    residuals<-PLSmodel$resid
    Q22D<-PLSmodel$Q2
    Q2cum2D<-PLSmodel$Q2cum
    ##########################################################################################################
    
    explvar3d<-list()
    explvar3dtita<-NULL
    #
    VIP3Dsas<-NULL
    vip3d<-list()
    #
    Q23Dsas<-NULL
    Q23D<-list()
    #
    Q23Dcumsas<-NULL
    Q23Dcum<-list()
    
    for (i in 1:dim(A3D)[3]) {
      PLSmodel<-plsreg(predictors = A3D[,,i], responses=Y2, comps = factors, crosval = TRUE)
      explvar3dtita<-PLSmodel$expvar
      nameexplvar <- paste('item:',i,sep='')
      explvar3dsas<-list(explvar3dtita)
      explvar3d[[nameexplvar]]<-explvar3dsas
      
      VIP3Dtita<-PLSmodel$VIP 
      namevip <- paste('item:',i,sep='')
      vip3dsas<-list(VIP3Dtita)
      vip3d[[namevip]]<-vip3dsas
    #
      Q23Dtita<-PLSmodel$Q2
      nameQ23D <- paste('item:',i,sep='')
      Q23Dsas<-list(Q23Dtita)
      Q23D[[nameQ23D]]<-Q23Dsas
      #
      Q23Dcumtita<-PLSmodel$Q2cum
      nameQ23Dcum <- paste('item:',i,sep='')
      Q23Dcumsas<-list(Q23Dcumtita)
      Q23Dcum[[nameQ23Dcum]]<-Q23Dcumsas
    }
    R2Xmean<-NULL
    R2Ymean<-NULL
    for (i in 1:length(explvar3d)) {
      a=as.matrix(explvar3d[i][[1]][[1]][,1])
      R2Xmean<-cbind(R2Xmean,a)
      b=as.matrix(explvar3d[i][[1]][[1]][,3])
      R2Ymean<-cbind(R2Ymean,b)
    }
    R2Xmean3D<-as.matrix(rowMeans(R2Xmean)); colnames(R2Xmean3D)<-"R2Xmean3D"
    R2Ymean3D<-as.matrix(rowMeans(R2Ymean)); colnames(R2Ymean3D)<-"R2Ymean3D"
    Explvar3D<-cbind(R2Xmean3D,R2Ymean3D)
    ###############################################
    # Calculating Q2mean3D
    Q2mean<-NULL
    for (i in 1:length(Q23D)) {
      a=as.matrix(Q23D[i][[1]][[1]][,1])
      Q2mean<-cbind(Q2mean,a)
    }
    Q2mean3D<-as.matrix(rowMeans(Q2mean)); colnames(Q2mean3D)<-"Q2mean3D"
    ###############################################
    # Calculating Q2cummean3D
    Q2cummean<-NULL
    for (i in 1:length(Q23Dcum)) {
      a=as.matrix(Q23Dcum[i][[1]][[1]][,1])
      Q2cummean<-cbind(Q2cummean,a)
    }
    Q2cummean3D<-as.matrix(rowMeans(Q2cummean)); colnames(Q2cummean3D)<-"Q2cummean3D"
    
    #####################################################################################
    ####VIP3D model 1: Calculating the mean of all times per each t
    t1mean<-t2mean<-t3mean<-t4mean<-NULL
    if(factors==2) {
      for (i in 1:factors) {
        t1=as.matrix(vip3d[i][[1]][[1]][,1])
        t1mean<-cbind(t1mean,t1)
        t2=as.matrix(vip3d[i][[1]][[1]][,2])
        t2mean<-cbind(t2mean,t2)
        #t3=as.matrix(vip3d[i][[1]][[1]][,3])
        #t3mean<-cbind(t3mean,t3)
      }
      t1mean3D<-as.matrix(rowMeans(t1mean)); colnames(t1mean3D)<-"t1mean3D"
      t2mean3D<-as.matrix(rowMeans(t2mean)); colnames(t2mean3D)<-"t2mean3D"
      #t3mean3D<-as.matrix(rowMeans(t3mean)); colnames(t3mean3D)<-"t3mean3D"
      VIP3Dmodel1<-cbind(t1mean3D,t2mean3D)#,t3mean3D)
    }
    if(factors==3) {
      for (i in 1:factors) {
        t1=as.matrix(vip3d[i][[1]][[1]][,1])
        t1mean<-cbind(t1mean,t1)
        t2=as.matrix(vip3d[i][[1]][[1]][,2])
        t2mean<-cbind(t2mean,t2)
        t3=as.matrix(vip3d[i][[1]][[1]][,3])
        t3mean<-cbind(t3mean,t3)
      }
      t1mean3D<-as.matrix(rowMeans(t1mean)); colnames(t1mean3D)<-"t1mean3D"
      t2mean3D<-as.matrix(rowMeans(t2mean)); colnames(t2mean3D)<-"t2mean3D"
      t3mean3D<-as.matrix(rowMeans(t3mean)); colnames(t3mean3D)<-"t3mean3D"
      VIP3Dmodel1<-cbind(t1mean3D,t2mean3D,t3mean3D)
    } 
    if(factors==4) {
      for (i in 1:factors) {
        t1=as.matrix(vip3d[i][[1]][[1]][,1])
        t1mean<-cbind(t1mean,t1)
        t2=as.matrix(vip3d[i][[1]][[1]][,2])
        t2mean<-cbind(t2mean,t2)
        t3=as.matrix(vip3d[i][[1]][[1]][,3])
        t3mean<-cbind(t3mean,t3)
        t4=as.matrix(vip3d[i][[1]][[1]][,4])
        t4mean<-cbind(t4mean,t4)
      }
      t1mean3D<-as.matrix(rowMeans(t1mean)); colnames(t1mean3D)<-"t1mean3D"
      t2mean3D<-as.matrix(rowMeans(t2mean)); colnames(t2mean3D)<-"t2mean3D"
      t3mean3D<-as.matrix(rowMeans(t3mean)); colnames(t3mean3D)<-"t3mean3D"
      t4mean3D<-as.matrix(rowMeans(t4mean)); colnames(t4mean3D)<-"t4mean3D"
      VIP3Dmodel1<-cbind(t1mean3D,t2mean3D,t3mean3D,t4mean3D)
    } 
    ##########################################################################################################
    ####VIP3D model 2: Calculating VIP per time, sum all t's generating a VIPtable per timepoint
    
    time<-NULL
    timetable<-NULL
    #t1mean<-t2mean<-t3mean<-NULL
    for (i in 1:length(vip3d)) {
      time<- as.matrix(rowSums(vip3d[i][[1]][[1]]))
      timetable<-cbind(timetable,time)
    }
    #timetable
    VIP3Dmodel2<-timetable
    #########
    xdimA3D<-dim(A3D)
    A<-matrix(A3D, xdimA3D[1], xdimA3D[2] * xdimA3D[3])
    plsdaforplotting<-plsda (X=A, Y= YN[,1,1], ncomp=factors)
    NPLSDAvariates<-plsdaforplotting$variates       # *Para graficar plotIndiv
    NPLSDAloadings<-plsdaforplotting$loadings       # * Para graficar plotloadings
    
    ##########
    tabla<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,2))
    tabla4<-tabla[,c(4,1,2)]
    transformedNPLSDAloadings<-tabla4
    
    ##########
    NPvariates<-list()
    NPvariatestita<-NULL
    NPloadingstita<-NULL
    NPloadings<-list()
    
    for (i in 1:dim(A3D)[3]) {
      plsdaforplotting<-plsda (X=A3D[,,i], Y= YN[,1,1], ncomp=factors)
      
      NPvariatestita<-plsdaforplotting$variates       # *Para graficar plotIndiv
      namevariatestita<-paste("item:",i,sep = "")
      variatestitatemp<-list(NPvariatestita)
      NPvariates[[namevariatestita]]<-variatestitatemp
      
      NPloadingstita<-plsdaforplotting$loadings       # * Para graficar plotloadings
      nameloadingdtita<-paste("item:",i,sep = "")
      loadingstitatemp<-list(NPloadingstita)
      NPloadings[[nameloadingdtita]]<-loadingstitatemp
    }
    NPLSDAvariatesperMode3<-NPvariates
    NPLSDAloadingsperMode3<-NPloadings
   
    
     ###### Consensus configuration
    Xcomp1variatesmean<-Xcomp2variatesmean<-Xcomp3variatesmean<-Xcomp4variatesmean<-NULL
    Ycomp1variatesmean<-Ycomp2variatesmean<-Ycomp3variatesmean<-Ycomp4variatesmean<-NULL
    if(factors==2){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1variatesmean<-cbind(Xcomp1variatesmean,a)
        b=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2variatesmean<-cbind(Xcomp2variatesmean,b)
        #c=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,3])
        #Xcomp3variatesmean<-cbind(Xcomp3variatesmean,c)
        
        d=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1variatesmean<-cbind(Ycomp1variatesmean,d)
        e=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2variatesmean<-cbind(Ycomp2variatesmean,e)
        #f=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,3])
        #Ycomp3variatesmean<-cbind(Ycomp3variatesmean,f)
        
      }
      comp1mean3D<-as.matrix(rowMeans(Xcomp1variatesmean)); colnames(comp1mean3D)<-"Comp1mean3D"
      comp2mean3D<-as.matrix(rowMeans(Xcomp2variatesmean)); colnames(comp2mean3D)<-"Comp2mean3D"
      #comp3mean3D<-as.matrix(rowMeans(Xcomp3variatesmean)); colnames(comp3mean3D)<-"Comp3mean3D"
      Block.XConsensus<-cbind(comp1mean3D,comp2mean3D)#,comp3mean3D)
      
      Ycomp1mean3D<-as.matrix(rowMeans(Ycomp1variatesmean)); colnames(Ycomp1mean3D)<-"Comp1mean3D"
      Ycomp2mean3D<-as.matrix(rowMeans(Ycomp2variatesmean)); colnames(Ycomp2mean3D)<-"Comp2mean3D"
      #Ycomp3mean3D<-as.matrix(rowMeans(Ycomp3variatesmean)); colnames(Ycomp3mean3D)<-"Comp3mean3D"
      Block.YConsensus<-cbind(Ycomp1mean3D,Ycomp2mean3D)#,Ycomp3mean3D)
      
      NPLSDAConsensusvariates<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    if(factors==3){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1variatesmean<-cbind(Xcomp1variatesmean,a)
        b=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2variatesmean<-cbind(Xcomp2variatesmean,b)
        c=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3variatesmean<-cbind(Xcomp3variatesmean,c)
        
        d=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1variatesmean<-cbind(Ycomp1variatesmean,d)
        e=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2variatesmean<-cbind(Ycomp2variatesmean,e)
        f=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3variatesmean<-cbind(Ycomp3variatesmean,f)
        
      }
      comp1mean3D<-as.matrix(rowMeans(Xcomp1variatesmean)); colnames(comp1mean3D)<-"Comp1mean3D"
      comp2mean3D<-as.matrix(rowMeans(Xcomp2variatesmean)); colnames(comp2mean3D)<-"Comp2mean3D"
      comp3mean3D<-as.matrix(rowMeans(Xcomp3variatesmean)); colnames(comp3mean3D)<-"Comp3mean3D"
      Block.XConsensus<-cbind(comp1mean3D,comp2mean3D,comp3mean3D)
      
      Ycomp1mean3D<-as.matrix(rowMeans(Ycomp1variatesmean)); colnames(Ycomp1mean3D)<-"Comp1mean3D"
      Ycomp2mean3D<-as.matrix(rowMeans(Ycomp2variatesmean)); colnames(Ycomp2mean3D)<-"Comp2mean3D"
      Ycomp3mean3D<-as.matrix(rowMeans(Ycomp3variatesmean)); colnames(Ycomp3mean3D)<-"Comp3mean3D"
      Block.YConsensus<-cbind(Ycomp1mean3D,Ycomp2mean3D,Ycomp3mean3D)
      
      NPLSDAConsensusvariates<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    if(factors==4){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1variatesmean<-cbind(Xcomp1variatesmean,a)
        b=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2variatesmean<-cbind(Xcomp2variatesmean,b)
        c=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3variatesmean<-cbind(Xcomp3variatesmean,c)
        d=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[1]][,4])
        Xcomp4variatesmean<-cbind(Xcomp4variatesmean,d)
        
        e=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1variatesmean<-cbind(Ycomp1variatesmean,e)
        f=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2variatesmean<-cbind(Ycomp2variatesmean,f)
        g=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3variatesmean<-cbind(Ycomp3variatesmean,g)
        h=as.matrix(NPLSDAvariatesperMode3[i][[1]][[1]][[2]][,4])
        Ycomp4variatesmean<-cbind(Ycomp4variatesmean,h)
      }
      comp1mean3D<-as.matrix(rowMeans(Xcomp1variatesmean)); colnames(comp1mean3D)<-"Comp1mean3D"
      comp2mean3D<-as.matrix(rowMeans(Xcomp2variatesmean)); colnames(comp2mean3D)<-"Comp2mean3D"
      comp3mean3D<-as.matrix(rowMeans(Xcomp3variatesmean)); colnames(comp3mean3D)<-"Comp3mean3D"
      comp4mean3D<-as.matrix(rowMeans(Xcomp4variatesmean)); colnames(comp4mean3D)<-"Comp4mean3D"
      Block.XConsensus<-cbind(comp1mean3D,comp2mean3D,comp3mean3D,comp4mean3D)
      
      Ycomp1mean3D<-as.matrix(rowMeans(Ycomp1variatesmean)); colnames(Ycomp1mean3D)<-"Comp1mean3D"
      Ycomp2mean3D<-as.matrix(rowMeans(Ycomp2variatesmean)); colnames(Ycomp2mean3D)<-"Comp2mean3D"
      Ycomp3mean3D<-as.matrix(rowMeans(Ycomp3variatesmean)); colnames(Ycomp3mean3D)<-"Comp3mean3D"
      Ycomp4mean3D<-as.matrix(rowMeans(Ycomp4variatesmean)); colnames(Ycomp4mean3D)<-"Comp4mean3D"
      Block.YConsensus<-cbind(Ycomp1mean3D,Ycomp2mean3D,Ycomp3mean3D,Ycomp4mean3D)
      
      NPLSDAConsensusvariates<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    
    
    
    
    
    
    
    
    ### CONSENSUS LOADINGS
    Xcomp1loadingsmean<-Xcomp2loadingsmean<-Xcomp3loadingsmean<-Xcomp4loadingsmean<-NULL
    Ycomp1loadingsmean<-Ycomp2loadingsmean<-Ycomp3loadingsmean<-Ycomp4loadingsmean<-NULL
    
    if (factors==2){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1loadingsmean<-cbind(Xcomp1loadingsmean,a)
        b=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2loadingsmean<-cbind(Xcomp2loadingsmean,b)
        #c=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,3])
        #Xcomp3loadingsmean<-cbind(Xcomp3loadingsmean,c)
        
        d=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1loadingsmean<-cbind(Ycomp1loadingsmean,d)
        e=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2loadingsmean<-cbind(Ycomp2loadingsmean,e)
        #f=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,3])
        #Ycomp3loadingsmean<-cbind(Ycomp3loadingsmean,f)
        
      }
      comp1meanloadings3D<-as.matrix(rowMeans(Xcomp1loadingsmean)); colnames(comp1meanloadings3D)<-"Comp1loadingsmean3D"
      comp2meanloadings3D<-as.matrix(rowMeans(Xcomp2loadingsmean)); colnames(comp2meanloadings3D)<-"Comp2loadingsmean3D"
      #comp3meanloadings3D<-as.matrix(rowMeans(Xcomp3loadingsmean)); colnames(comp3meanloadings3D)<-"Comp3loadingsmean3D"
      Block.XConsensus<-cbind(comp1meanloadings3D,comp2meanloadings3D)#,comp3meanloadings3D)
      
      Ycomp1meanloadings3D<-as.matrix(rowMeans(Ycomp1loadingsmean)); colnames(Ycomp1meanloadings3D)<-"YComp1loadingsmean3D"
      Ycomp2meanloadings3D<-as.matrix(rowMeans(Ycomp2loadingsmean)); colnames(Ycomp2meanloadings3D)<-"YComp2loadingsmean3D"
      #Ycomp3meanloadings3D<-as.matrix(rowMeans(Ycomp3loadingsmean)); colnames(Ycomp3meanloadings3D)<-"YComp3loadingsmean3D"
      Block.YConsensus<-cbind(Ycomp1meanloadings3D,Ycomp2meanloadings3D)#,Ycomp3meanloadings3D)
      
      NPLSDAConsensusloadings<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    if (factors==3){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1loadingsmean<-cbind(Xcomp1loadingsmean,a)
        b=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2loadingsmean<-cbind(Xcomp2loadingsmean,b)
        c=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3loadingsmean<-cbind(Xcomp3loadingsmean,c)
        
        d=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1loadingsmean<-cbind(Ycomp1loadingsmean,d)
        e=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2loadingsmean<-cbind(Ycomp2loadingsmean,e)
        f=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3loadingsmean<-cbind(Ycomp3loadingsmean,f)
        
      }
      comp1meanloadings3D<-as.matrix(rowMeans(Xcomp1loadingsmean)); colnames(comp1meanloadings3D)<-"Comp1loadingsmean3D"
      comp2meanloadings3D<-as.matrix(rowMeans(Xcomp2loadingsmean)); colnames(comp2meanloadings3D)<-"Comp2loadingsmean3D"
      comp3meanloadings3D<-as.matrix(rowMeans(Xcomp3loadingsmean)); colnames(comp3meanloadings3D)<-"Comp3loadingsmean3D"
      Block.XConsensus<-cbind(comp1meanloadings3D,comp2meanloadings3D,comp3meanloadings3D)
      
      Ycomp1meanloadings3D<-as.matrix(rowMeans(Ycomp1loadingsmean)); colnames(Ycomp1meanloadings3D)<-"YComp1loadingsmean3D"
      Ycomp2meanloadings3D<-as.matrix(rowMeans(Ycomp2loadingsmean)); colnames(Ycomp2meanloadings3D)<-"YComp2loadingsmean3D"
      Ycomp3meanloadings3D<-as.matrix(rowMeans(Ycomp3loadingsmean)); colnames(Ycomp3meanloadings3D)<-"YComp3loadingsmean3D"
      Block.YConsensus<-cbind(Ycomp1meanloadings3D,Ycomp2meanloadings3D,Ycomp3meanloadings3D)
      
      NPLSDAConsensusloadings<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    if (factors==4){
      for (i in 1:factors) {
        a=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,1])
        Xcomp1loadingsmean<-cbind(Xcomp1loadingsmean,a)
        b=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,2])
        Xcomp2loadingsmean<-cbind(Xcomp2loadingsmean,b)
        c=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,3])
        Xcomp3loadingsmean<-cbind(Xcomp3loadingsmean,c)
        d=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[1]][,4])
        Xcomp4loadingsmean<-cbind(Xcomp4loadingsmean,d)
        
        e=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,1])
        Ycomp1loadingsmean<-cbind(Ycomp1loadingsmean,e)
        f=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,2])
        Ycomp2loadingsmean<-cbind(Ycomp2loadingsmean,f)
        g=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,3])
        Ycomp3loadingsmean<-cbind(Ycomp3loadingsmean,g)
        h=as.matrix(NPLSDAloadingsperMode3[i][[1]][[1]][[2]][,4])
        Ycomp4loadingsmean<-cbind(Ycomp4loadingsmean,h)
        
      }
      comp1meanloadings3D<-as.matrix(rowMeans(Xcomp1loadingsmean)); colnames(comp1meanloadings3D)<-"Comp1loadingsmean3D"
      comp2meanloadings3D<-as.matrix(rowMeans(Xcomp2loadingsmean)); colnames(comp2meanloadings3D)<-"Comp2loadingsmean3D"
      comp3meanloadings3D<-as.matrix(rowMeans(Xcomp3loadingsmean)); colnames(comp3meanloadings3D)<-"Comp3loadingsmean3D"
      comp4meanloadings3D<-as.matrix(rowMeans(Xcomp4loadingsmean)); colnames(comp4meanloadings3D)<-"Comp4loadingsmean3D"
      Block.XConsensus<-cbind(comp1meanloadings3D,comp2meanloadings3D,comp3meanloadings3D,comp4meanloadings3D)
      
      Ycomp1meanloadings3D<-as.matrix(rowMeans(Ycomp1loadingsmean)); colnames(Ycomp1meanloadings3D)<-"YComp1loadingsmean3D"
      Ycomp2meanloadings3D<-as.matrix(rowMeans(Ycomp2loadingsmean)); colnames(Ycomp2meanloadings3D)<-"YComp2loadingsmean3D"
      Ycomp3meanloadings3D<-as.matrix(rowMeans(Ycomp3loadingsmean)); colnames(Ycomp3meanloadings3D)<-"YComp3loadingsmean3D"
      Ycomp4meanloadings3D<-as.matrix(rowMeans(Ycomp4loadingsmean)); colnames(Ycomp4meanloadings3D)<-"YComp4loadingsmean3D"
      Block.YConsensus<-cbind(Ycomp1meanloadings3D,Ycomp2meanloadings3D,Ycomp3meanloadings3D,Ycomp4meanloadings3D)
      
      NPLSDAConsensusloadings<-list(Block.XConsensus=Block.XConsensus,Block.YConsensus=Block.YConsensus)
    }
    
    
    
   
    
    
    DAR2blockX<-plsdaforplotting$explained_variance$X
    blox<-sum(explvar[c(1:factors),1]) * DAR2blockX[1:factors]
    bloy<-explvar[c(1:factors),3] #BlockX tambien porque es que tanto X explica a Y
    NPLSDAexplVar<-t(rbind(R2.X=blox,R2.Y=bloy))
    
    #pch=c(0,1)
    #cex=c(1,2)
    #plotArrow(plsdaforplotting, ind.names=TRUE, abline = TRUE, pch=pch,cex=cex) # Esta si
    #plotDiablo(plsdaforplotting)
    # plotVar(plsdaforplotting) Mejorar esta
    # plotLoadings(plsdaforplotting)  Esta si
    #network(plsdaforplotting) # Muy pesada
    #selectVar ((plsdaforplotting))
    
    ######
    result <- list(FactorsX = FactorsX, FactorsY = FactorsY, 
                   T = Tt, WsupraJ = WsupraJ, WsupraK = WsupraK, U = U, QsupraJ = QsupraJ, 
                   QsupraK = QsupraK, B = B, Gu = Gu, G = G,Ypred =Ypred,
                   explvar=explvar, VIP2D=VIP2D, Explvar3D=Explvar3D, VIP3Dmodel1=VIP3Dmodel1,VIP3Dmodel2=VIP3Dmodel2,
                   tucorrelation=tucorr,residuals=residuals,
                   NPLSDAvariates=NPLSDAvariates, NPLSDAloadings=NPLSDAloadings, 
                   transformedNPLSDAloadings=transformedNPLSDAloadings,
                   NPLSDAexplVar=NPLSDAexplVar,
                   NPLSDAQ22D=Q22D,
                   NPLSDAQ22Dcum=Q2cum2D,
                   NPLSDAQ2mean3D=Q2mean3D,
                   NPLSDAQ2cummean3D=Q2cummean3D,
                   NPLSDAvariatesperMode3=NPLSDAvariatesperMode3, NPLSDAloadingsperMode3=NPLSDAloadingsperMode3,
                   NPLSDAConsensusvariates=NPLSDAConsensusvariates,NPLSDAConsensusloadings=NPLSDAConsensusloadings,
                   
                   Design=YN[,1,1]
                   
                   
    )
    return(result)
    
  }
}
  
  

  
#####################################################################################
# Plotting things
plotNPLSDAmod<-function (X, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                       cutoff = length(X), factors=2, penalty=1) 
  # Function to plot all necesary figures to interpret NPLSDA results
  # Input: NPLSDA result
  # Output:
  # Plot of the loadings per variable
  # Heatmap of the behavior of all variables in all individuals
  # Mode1 NPLS-DA relative to X explaining Y
  # Mode1 NPLS-DA relative to X explaining Y colored by AgeGroup
  # Mode1 NPLS-DA relative to X explaining Y colored by FirstAAb that appears
  # Mode3 NPLS-DA relative to X explaining Y 
  # Mode2 NPLS-DA relative to X explaining Y 
  # Analysis of the Core Array of the model to determine explained variances and the contributions of the components in each mode
# Here you can play with the script to obtain what you desire.

{
  library(sfsmisc)
  par(mar=c(5,4,4,2))
  
  #
  pc1 <- PCs[1]
  pc2 <- PCs[2]
  a = paste("Factors", "X", sep = "")
  label = NULL
  FactorsX <- X[[a]]
  x1 <- FactorsX[[1]][, pc1]
  y1 <- FactorsX[[1]][, pc2]
  x2 <- FactorsX[[2]][, pc1]
  y2 <- FactorsX[[2]][, pc2]
  x3 <- FactorsX[[3]][, pc1]
  y3 <- FactorsX[[3]][, pc2]
  
  b = paste("Factors", "Y", sep = "")
  label = NULL
  FactorsY <- X[[b]]
  x1Y <- FactorsY[[1]][, pc1]
  y1Y <- FactorsY[[1]][, pc2]
  x2Y <- FactorsY[[2]][, pc1]
  y2Y <- FactorsY[[2]][, pc2]
  x3Y <- FactorsY[[3]][, pc1]
  y3Y <- FactorsY[[3]][, pc2]
  par(mfrow=c(1,1))
  
  ######################################## Mode 1 Block X NPLSDA
  DVars<-NULL
  DVars$colourOutcome<- ""
  DVars$colourOutcome[X$Design == 0]<-"dodgerblue2"
  DVars$colourOutcome[X$Design == 1]<-"orange2"
  
  if (is.null(X$NPLSDAvariates$X)){   # BLOCK.NPLSDA
    plot(X$NPLSDAvariates$Block.X[,pc1],X$NPLSDAvariates$Block.X[,pc2], type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome, 
         xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
         ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
         main = paste("Block X Mode 1: NPLS-DA", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAvariates$Block.X[,pc1]), max(X$NPLSDAvariates$Block.X[,pc1])),
         ylim = c(min(X$NPLSDAvariates$Block.X[,pc2]), max(X$NPLSDAvariates$Block.X[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #text (X$NPLSDAvariates$Block.X[,pc1],X$NPLSDAvariates$Block.X[,pc2], labels=rownames(Factors$Mode1), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    ######################################## Mode 1 Block Y ########################################################
    plot(X$NPLSDAvariates$Block.Y[,pc1],X$NPLSDAvariates$Block.Y[,pc2], type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome,
         xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,2]*100,digits=2)," %", sep = ""),
         ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,2]*100,digits=2)," %", sep = ""), las=1,
         main = paste("Block Y Mode 1: NPLS-DA", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAvariates$Block.Y[,pc1]), max(X$NPLSDAvariates$Block.Y[,pc1])),
         ylim = c(min(X$NPLSDAvariates$Block.Y[,pc2]), max(X$NPLSDAvariates$Block.Y[,pc2])))
    #xlim = c((min(x1Y)+ (0.2*min(x1Y))), (max(x1Y)+ (0.2*max(x1Y)))),
    #ylim = c((min(y1Y)+ (0.2*min(y1Y))), (max(y1Y)+ (0.2*max(y1Y)))))
    abline(v=0,h=0)
   # text (X$NPLSDAvariates$Block.Y[,pc1],X$NPLSDAvariates$Block.Y[,pc2], labels=rownames(FactorsY$Mode1), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    
    ########################################################## Per time point
    for(i in 1:length(X$NPLSDAvariatesperMode3)){
      plot(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc1],X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc2], 
         type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome, 
         xlab = paste("Component ", pc1,sep = ""), 
         ylab = paste("Component ", pc2,sep = ""), las=1,
         main = paste("Block X Mode 1: NPLS-DA mode 3-",i,  sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc1]), max(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc1])),
         ylim = c(min(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc2]), max(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #text (X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc1],X$NPLSDAvariatesperMode3[[i]][[1]]$Block.X[,pc2], 
    #      labels=rownames(Factors$Mode1), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    }
    #####################################################################################
    # Consensus
    plot(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1],X$NPLSDAConsensusvariates$Block.XConsensus[,pc2], 
         type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome, 
         xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
         ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
         main = paste("Block X Mode 1: NPLS-DA Consensus", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1])),
         ylim = c(min(X$NPLSDAConsensusvariates$Block.XConsensus[,pc2]), max(X$NPLSDAConsensusvariates$Block.XConsensus[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #text (X$NPLSDAConsensusvariates$Block.XConsensus[,pc1],X$NPLSDAConsensusvariates$Block.XConsensus[,pc2], 
    #      labels=rownames(Factors$Mode1), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    ####################################################################################################################
    
    ########################################################## Block Y Per time point
    for(i in 1:length(X$NPLSDAvariatesperMode3)){
      plot(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc1],X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc2], 
           type = "p", col = "black", pch=21,
           bg=DVars$colourOutcome, 
           xlab = paste("Component ", pc1,sep = ""), 
           ylab = paste("Component ", pc2,sep = ""), las=1,
           main = paste("Block Y Mode 1: NPLS-DA mode 3-",i,  sep= " "), cex.main=2,
           cex.lab= 1.2,
           xlim = c(min(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc1]), max(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc1])),
           ylim = c(min(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc2]), max(X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc2])))
      #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
      #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
      abline(v=0,h=0)
     # text (X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc1],X$NPLSDAvariatesperMode3[[i]][[1]]$Block.Y[,pc2], 
    #        labels=rownames(Factors$Mode1), 
    #        col=DVars$colourOutcome,cex=1, pos=3)
    }
    #####################################################################################
    # Consensus Block Y mode 1
    #X$NPLSDAConsensusvariates$Block.YConsensus[,1]
    plot(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1],X$NPLSDAConsensusvariates$Block.YConsensus[,pc2], 
         type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome, 
         xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
         ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
         main = paste("Block Y Mode 1: NPLS-DA Consensus", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1])),
         ylim = c(min(X$NPLSDAConsensusvariates$Block.YConsensus[,pc2]), max(X$NPLSDAConsensusvariates$Block.YConsensus[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #text (X$NPLSDAConsensusvariates$Block.YConsensus[,pc1],X$NPLSDAConsensusvariates$Block.YConsensus[,pc2], 
    #      labels=rownames(Factors$Mode1), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    ####################################################################################################################
    ######################################## Mode 1 Joint Biplot 
    plot(X$NPLSDAvariates$Block.X[,pc1], X$NPLSDAvariates$Block.Y[,pc1], type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome,
         xlab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," % of Block X", sep = ""),
         ylab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,2]*100,digits=2)," % of Block Y", sep = ""), las=1,
         main = paste("Biplot Mode1 Blocks X & Y: NPLS-DA", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAvariates$Block.X[,pc1]), max(X$NPLSDAvariates$Block.X[,pc1])),
         ylim = c(min(X$NPLSDAvariates$Block.Y[,pc1]), max(X$NPLSDAvariates$Block.Y[,pc1])))
    abline(v=0,h=0)
    #text (X$NPLSDAvariates$Block.X[,pc1], X$NPLSDAvariates$Block.Y[,pc1], labels=rownames(Factors$Mode1), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    
    ####################################################################################################################
    ######################################## Consensus Mode 1 Joint Biplot 
    plot(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1], X$NPLSDAConsensusvariates$Block.YConsensus[,pc1], 
         type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome,
         xlab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," % of Block X", sep = ""),
         ylab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,2]*100,digits=2)," % of Block Y", sep = ""), las=1,
         main = paste("Biplot Mode1 Blocks X & Y: NPLS-DA Consensus", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1])),
         ylim = c(min(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1])))
    abline(v=0,h=0)
    #text (X$NPLSDAConsensusvariates$Block.XConsensus[,pc1], X$NPLSDAConsensusvariates$Block.YConsensus[,pc1], 
    #      labels=rownames(Factors$Mode1), 
     #     col=DVars$colourOutcome,cex=1, pos=3)
    #################################################################################################################
    ######################################## Just the Mode2 of X
    
    plot(X$NPLSDAloadings$Block.X[,pc1],X$NPLSDAloadings$Block.X[,pc2], type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Block X Mode 2: NPLS-DA", sep= " "), 
         xlim = c(min(X$NPLSDAloadings$Block.X[,pc1]), max(X$NPLSDAloadings$Block.X[,pc1])),
         ylim = c(min(X$NPLSDAloadings$Block.X[,pc2]), max(X$NPLSDAloadings$Block.X[,pc2])))
    #xlim = c(-1,1),
    #ylim = c(-1,1))
    
    abline(v=0,h=0)
    #points (X$NPLSDAloadings$BlockY[,1],X$NPLSDAloadings$BlockY[,2])
    ########################################################## Mode 2 of X Per time point
    for(i in 1:length(X$NPLSDAloadingsperMode3)){
      plot(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1],X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2], 
           type = "p", col = "black", pch=21,
           bg="dodgerblue2", 
           xlab = paste("Component ", pc1,sep = ""), 
           ylab = paste("Component ", pc2,sep = ""), las=1,
           main = paste("Block X Mode 1: NPLS-DA mode 3-",i,  sep= " "), cex.main=2,
           cex.lab= 1.2,
           xlim = c(min(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1]), max(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1])),
           ylim = c(min(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2]), max(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2])))
      #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
      #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
      abline(v=0,h=0)
      #text (X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,1],X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,2], 
       #     labels=rownames(Factors$Mode2), 
        #    col=DVars$colourOutcome,cex=1, pos=3)
    }
    #####################################################################################
    #######################   Consensus Block X Mode 2   ################################
    #X$NPLSDAConsensusloadings$Block.XConsensus
    plot(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1],X$NPLSDAConsensusloadings$Block.XConsensus[,pc2], 
         type = "p", col = "black", pch=21,
         bg="dodgerblue2", 
         xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
         ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
         main = paste("Block X Mode 2: NPLS-DA Consensus", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1]), max(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1])),
         ylim = c(min(X$NPLSDAConsensusloadings$Block.XConsensus[,pc2]), max(X$NPLSDAConsensusloadings$Block.XConsensus[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #labels=NULL
    #x<-NULL
    #y<-NULL
    #if (is.null(labels)) {
    #  labels <- names
    #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
    #             rownames(Factors$Mode3))
    #}
    #ABS <- abs(x2trans * y2trans)
    #r <- rank(-ABS)
    #value = match(c(1:length(ABS)), r)
    #text (x2,y2, 
    #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
    #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
    #      cex=0.8, pos=3)
    #for (n in 1:length(ABS)) {
    #  x[c(n)] <- x2trans[value[c(n)]]
    #  y[c(n)] <- y2trans[value[c(n)]]
    #  label[n] <- labels[value[c(n)]]
    #}
    #x <- x[-c((n + 1):(length(x) + 1))]
    #y <- y[-c((n + 1):(length(y) + 1))]
    #text(x, y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>0.05 | abs(y)>0.05),
    #                                                        "dodgerblue2", "transparent"), cex=0.7)
    
    
    
    #text (X$NPLSDAConsensusvariates$Block.XConsensus[,1],X$NPLSDAConsensusvariates$Block.XConsensus[,2], labels=rownames(FactorsY$Mode2), 
    #      col=DVars$colourOutcome,cex=1, pos=3)
    
    ######################################## Just the Mode2 of Y
    background<-NULL
    background<-c(background[rownames(X$NPLSDAloadings$Block.Y) == 0]<-"dodgerblue2",
                  background[rownames(X$NPLSDAloadings$Block.Y) == 1]<-"orange2")
    
    plot(X$NPLSDAloadings$Block.Y[,pc1],X$NPLSDAloadings$Block.Y[,pc2], type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), cex.main=2,
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Block Y Mode 2: NPLS-DA", sep= " "), 
         xlim = c(min(X$NPLSDAloadings$Block.Y[,pc1]), max(X$NPLSDAloadings$Block.Y[,pc1])),
         ylim = c(min(X$NPLSDAloadings$Block.Y[,pc2]), max(X$NPLSDAloadings$Block.Y[,pc2])))
    #xlim = c(-1,1),
    #ylim = c(-1,1))
    
    abline(v=0,h=0)
    #text(X$NPLSDAloadings$Block.Y[,1],X$NPLSDAloadings$Block.Y[,2], 
    #     labels= rownames(FactorsY$Mode2),
    #     cex= 1, pos=1, col="dodgerblue2")
    ##############################################################################################
    #######################   Consensus Block Y Mode 2   ################################
    #X$NPLSDAConsensusloadings$Block.XConsensus
    plot(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1],X$NPLSDAConsensusloadings$Block.YConsensus[,pc2], 
         type = "p", col = "black", pch=21,
         bg="dodgerblue2", 
         xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
         ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
         main = paste("Block Y Mode 2: NPLS-DA Consensus", sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1])),
         ylim = c(min(X$NPLSDAConsensusloadings$Block.YConsensus[,pc2]), max(X$NPLSDAConsensusloadings$Block.YConsensus[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #text (X$NPLSDAConsensusloadings$Block.YConsensus[,pc1],X$NPLSDAConsensusloadings$Block.YConsensus[,pc2], labels=rownames(FactorsY$Mode2), 
    #      col="dodgerblue2",cex=1, pos=3)
    ###########################################################################################
    ######################################## Biplot of Block X Mode 2
    plot(X$NPLSDAloadings$Block.X[,pc1]*penalty,X$NPLSDAloadings$Block.X[,pc2]*penalty, type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block X Modes 2 & 3: NPLS-DA", sep= " "), 
         xlim = c(-1,1),
         ylim = c(-1,1))
    abline(v=0,h=0)
    arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
    text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
         cex= 1, pos=1, col=background )
    ###########################################################################################
    ######################################## Consensus Biplot of Block X Mode 2
    plot(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1]*penalty,X$NPLSDAConsensusloadings$Block.XConsensus[,pc2]*penalty, 
         type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block X Modes 2 & 3: Consensus NPLS-DA", sep= " "), 
         xlim = c(-1,1),
         ylim = c(-1,1))
    abline(v=0,h=0)
    arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
    text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
         cex= 1, pos=1, col=background )
    ##########################################################################################
    ##########################################################################################
    # Biplot of Block X mode 2 per timepoint
    for(i in 1:length(X$NPLSDAloadingsperMode3)){
      plot(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1]*penalty,X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2]*penalty, 
           type = "p", col = "black", pch=21,
           bg="dodgerblue2", 
           xlab = paste("Component ", pc1,sep = ""), 
           ylab = paste("Component ", pc2,sep = ""), las=1,
           main = paste("Block X Mode 1: NPLS-DA mode 3-",i,  sep= " "), cex.main=2,
           cex.lab= 1.2,
           xlim = c(-1,1),
           ylim = c(-1,1))
           #xlim = c(min(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1]), max(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1])),
           #ylim = c(min(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2]), max(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2])))
      #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
      #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
      abline(v=0,h=0)
      arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
      text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
           cex= 1, pos=1, col=background )
      #text (X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc1],X$NPLSDAloadingsperMode3[[i]][[1]]$Block.X[,pc2], 
      #     labels=rownames(Factors$Mode2), 
      #    col=DVars$colourOutcome,cex=1, pos=3)
    }
   ##########################################################################################################
    ######################################## Biplot of Block Y Mode 2
    plot(X$NPLSDAloadings$Block.Y[,pc1]*penalty,X$NPLSDAloadings$Block.Y[,pc2]*penalty, type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block Y Modes 2 & 3: NPLS-DA", sep= " "), 
         xlim = c(-1,1),
         ylim = c(-1,1))
    abline(v=0,h=0)
    arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
    text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
         cex= 1, pos=1, col=background )
    #########################################################################################################
    ######################################## Consensus Biplot of Block Y Mode 2
    plot(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1]*penalty,X$NPLSDAConsensusloadings$Block.YConsensus[,pc2]*penalty, 
         type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block X Modes 2 & 3: Consensus NPLS-DA", sep= " "), 
         xlim = c(-1,1),
         ylim = c(-1,1))
    abline(v=0,h=0)
    arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
    text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
         cex= 1, pos=1, col=background )
    
    
    
    ##########################################################################################
    # Biplot of Block Y mode 2 per timepoint
    for(i in 1:length(X$NPLSDAvariatesperMode3)){
      plot(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc1]*penalty,X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc2]*penalty, 
           type = "p", col = "black", pch=21,
           bg="dodgerblue2", 
           xlab = paste("Component ", pc1,sep = ""), 
           ylab = paste("Component ", pc2,sep = ""), las=1,
           main = paste("Block Y Mode 1: NPLS-DA mode 3-",i,  sep= " "), cex.main=2,
           cex.lab= 1.2,
           xlim = c(-1,1),
           ylim = c(-1,1))
           #xlim = c(min(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc1]), max(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc1])),
           #ylim = c(min(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc2]), max(X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc2])))
           #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
           #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
      abline(v=0,h=0)
      arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
      text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
           cex= 1, pos=1, col=background )
      #text (X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc1],X$NPLSDAloadingsperMode3[[i]][[1]]$Block.Y[,pc2], 
      #      labels=rownames(FactorsY$Mode2), 
      #      col=DVars$colourOutcome,cex=1, pos=3)
    }
    
    
    ######################################## Triplot of Blocks X & Y Modes 2 & 3
    background2<-NULL
    background2<-c(background2[rownames(X$NPLSDAloadings$Block.Y) == 0]<-"transparent",
                  background2[rownames(X$NPLSDAloadings$Block.Y) == 1]<-"orange2")
    plot(X$NPLSDAloadings$Block.Y[,pc1],X$NPLSDAloadings$Block.Y[,pc2], type = "p", col = "black", pch=21,
         bg="dodgerblue2",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block Y Modes 2 & 3: NPLS-DA", sep= " "), 
         xlim = c(-1,1),
         ylim = c(-1,1))
    abline(v=0,h=0)
    points(X$NPLSDAloadings$Block.X[,pc1],X$NPLSDAloadings$Block.X[,pc2], type = "p", col = "black", pch=21,
         bg="orangered1",
         xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block X Modes 2 & 3: NPLS-DA", sep= " "), 
         xlim = c(-1,1),
         ylim = c(-1,1))
    arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background2, pch=21, bg = background)
    text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
         cex= 1, pos=1, col=background2 )
    legend("topright",  c("Block.X", "Block.Y", "Response"), col = c("orangered1", "dodgerblue2", "orange2"),
           pch = c(16, 16, 16))
    #
    
    #################################################################  Transformed Mode 2 of X
    #X$transformedNPLSDAloadings
    #x2trans<-X$transformedNPLSDAloadings[,2]
    #y2trans<-X$transformedNPLSDAloadings[,3]
    #names<-X$transformedNPLSDAloadings$names
    
    #plot(x2trans,y2trans, type = "p", col = "black", pch=21,
     #    bg="dodgerblue2", 
    #     cex=1,
         #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
         #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
    #     xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
    #     main = paste("Block X Mode 2: transformed NPLS-DA"), xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
    #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
    #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
    #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
    #abline(v=0,h=0)
    #labels=NULL
    #x<-NULL
    #y<-NULL
    #if (is.null(labels)) {
    #  labels <- names
      #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
      #             rownames(Factors$Mode3))
    #}
    #ABS <- abs(x2trans * y2trans)
    #r <- rank(-ABS)
    #value = match(c(1:length(ABS)), r)
    #text (x2,y2, 
    #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
    #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
    #      cex=0.8, pos=3)
    #for (n in 1:length(ABS)) {
    #  x[c(n)] <- x2trans[value[c(n)]]
    #  y[c(n)] <- y2trans[value[c(n)]]
    #  label[n] <- labels[value[c(n)]]
    #}
    #x <- x[-c((n + 1):(length(x) + 1))]
    #y <- y[-c((n + 1):(length(y) + 1))]
    #text(x, y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>0.05 | abs(y)>0.05),
    #                                                        "dodgerblue2", "transparent"), cex=0.7)
    
    ################################################################ Transformed Biplot
    
    #X$transformedNPLSDAloadings
    #x2trans<-X$transformedNPLSDAloadings[,2]
    #y2trans<-X$transformedNPLSDAloadings[,3]
    #names<-X$transformedNPLSDAloadings$names
    
    #plot(x2trans,y2trans, type = "p", col = "black", pch=21,
    #     bg="dodgerblue2", 
    #     cex=1,
         #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
         #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
    #    xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
    #     main = paste("Biplot Block X Modes 2 & 3: transformed NPLS-DA"), xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
    #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
    #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
    #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
    #abline(v=0,h=0)
    #labels=NULL
    #x<-NULL
    #y<-NULL
    #if (is.null(labels)) {
    #  labels <- names
      #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
      #             rownames(Factors$Mode3))
    #}
    #ABS <- abs(x2trans * y2trans)
    #r <- rank(-ABS)
    #value = match(c(1:length(ABS)), r)
    #text (x2,y2, 
    #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
    #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
    #      cex=0.8, pos=3)
    #for (n in 1:length(ABS)) {
    #  x[c(n)] <- x2trans[value[c(n)]]
    #  y[c(n)] <- y2trans[value[c(n)]]
    #  label[n] <- labels[value[c(n)]]
    #}
    #x <- x[-c((n + 1):(length(x) + 1))]
    #y <- y[-c((n + 1):(length(y) + 1))]
    #text(x, y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>0.05 | abs(y)>0.05),
    #                                                        "dodgerblue2", "transparent"), cex=0.7)
    
    #arrows (0,0,X$NPLSDAloadings$Y[,1],X$NPLSDAloadings$Y[,2], col = background, pch=21, bg = background)
    #text(X$NPLSDAloadings$Y[,1],X$NPLSDAloadings$Y[,2], labels= rownames(X$NPLSDAloadings$Y),
    #     cex= 1, pos=1, col=background )             
    
    
    
    
    ######################################## 
    ###########      Mode 3     ############
    ######################################## 
    # Without variable vectors
    plot(x3, y3, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Block X Mode 3: NPLS-DA"), 
         #xlim = c((min(x3)), (max(x3))),
         #ylim = c((min(y3)), (max(y3))))
         xlim = c((min(x3)), (max(x3)+ abs(0.2*max(x3)))),
         ylim = c((min(y3)), (max(y3)+ abs(0.2*max(y3)))))
    abline(v=0,h=0)
    text (x3,y3, labels=rownames(FactorsX$Mode3), col="dodgerblue2",cex=1, pos=3)
    
    ######################################## 
    plot(x3Y, y3Y, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Block Y Mode 3: NPLS-DA"), 
         xlim = c((min(x3Y)), (max(x3Y)+ abs(0.2*max(x3Y)))),
         ylim = c((min(y3Y)), (max(y3Y)+ abs(0.2*max(y3Y)))))
    abline(v=0,h=0)
    text (x3Y,y3Y, labels=rownames(FactorsY$Mode3), col="dodgerblue2",cex=1, pos=3)
    
    # With variable vectors
    plot(x3, y3, type = "p", col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1, sep = " "), las=1,
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block X Modes 2 and 3: NPLS-DA"), 
         xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
    abline(v=0,h=0)
    text (x3,y3, labels=rownames(FactorsX$Mode3), col="dodgerblue2",cex=1, pos=3)
    arrows(0,0,penalty*x2,penalty*y2, pch=16, lwd=2,
           col = ifelse ((abs(x2)>(penalty/2)*0.05 | abs(y2)>(penalty/2)*0.05),"orange2", "transparent") )
    labels=NULL
    x<-NULL
    y<-NULL
    if (is.null(labels)) {
      labels <- list(rownames(FactorsX$Mode1), rownames(FactorsX$Mode2), 
                     rownames(FactorsX$Mode3))
    }
    ABS <- abs(x2 * y2)
    r <- rank(-ABS)
    value = match(c(1:length(ABS)), r)
    for (n in 1:length(ABS)) {
      x[c(n)] <- FactorsX[[2]][value[c(n)], 1]
      y[c(n)] <- FactorsX[[2]][value[c(n)], 2]
      label[n] <- labels[[2]][value[c(n)]]
    }
    x <- x[-c((n + 1):(length(x) + 1))]
    y <- y[-c((n + 1):(length(y) + 1))]
    text(penalty*x, penalty*y, labels = label, adj = c(0, -0.5),col= ifelse ((abs(x)>(penalty/2)*0.05 | abs(y)>(penalty/2)*0.05),
                                                                             "orange2", "transparent"), cex=0.7)
    
    
    #text (x2,y2, labels= ifelse ((abs(x2)>0.05 | abs(y2)>0.05),names(x2), ""),
    #      col = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"orange2", "transparent"),    
    #      cex=1.5, pos=1)
    
    ######################################## 
    plot(x3Y, y3Y, type = "p",  col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1, sep = " "), las=1,
         ylab = paste("Component", pc2, sep = " "), 
         main = paste("Biplot Block Y Modes 2 and 3: NPLS-DA"), 
         xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
    abline(v=0,h=0)
    text (x3Y,y3Y, labels=rownames(FactorsY$Mode3), col="dodgerblue2",cex=1, pos=3)
    arrows(0,0,penalty*x2Y,penalty*y2Y, pch=16, lwd=2,
           col = ifelse ((abs(x2Y)>(penalty/2)*0.05 | abs(y2Y)>(penalty/2)*0.05),"orange2", "transparent"))
    labels=NULL
    label=NULL
    x<-NULL
    y<-NULL
    if (is.null(labels)) {
      labels <- list(rownames(FactorsY$Mode1), rownames(FactorsY$Mode2), 
                     rownames(FactorsY$Mode3))
    }
    ABS <- abs(x2Y * y2Y)
    r <- rank(-ABS)
    value = match(c(1:length(ABS)), r)
    for (n in 1:length(ABS)) {
      x[c(n)] <- FactorsY[[2]][value[c(n)], 1]
      y[c(n)] <- FactorsY[[2]][value[c(n)], 2]
      label[n] <- labels[[2]][value[c(n)]]
    }
    x <- x[-c((n + 1):(length(x) + 1))]
    y <- y[-c((n + 1):(length(y) + 1))]
    text(penalty*x, penalty*y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>(penalty/2)*0.05 | abs(y)>(penalty/2)*0.05),
                                                                            "orange2", "transparent"), cex=0.7)
    
    
    #text (x2Y,y2Y, labels= ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),names(x2Y), ""),
    #      col = ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),"orange2", "transparent"),    
    #      cex=1, pos=1)
    
    #    Joint Mode 3 ####################################### 
    plot(x3, x3Y, type = "p", col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1,"Block X", sep = " "), 
         ylab = paste("Component", pc1, "Block Y",sep = " "), las =1,
         main = paste("Joint plot: NPLS-DA Mode 3"), 
         xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
    #abline(v=0,h=0)
    abline(a=0,b=1,col= "gray47", lwd=1)
    text (x3,x3Y, labels=rownames(FactorsX$Mode3), col="dodgerblue2",cex=1, pos=3)
    
    ############################################################
    ex <- as.matrix(explcore(X$G, factors, n=5))
    plot.new()
    mtext(paste("Expl.var of X: ", round(sum(X$NPLSDAexplVar[c(1:factors),1])*100,digits=2),"%"))
    mtext(paste("Expl.var of Y by X: ", round(sum(X$NPLSDAexplVar[c(1:factors),2])*100,digits=2),"%"),line= -2)
    #mtext(paste("model expl.var: ", sum(as.numeric(ex[c(1:2),4]))))
    #mtext(paste("model expl.var: ", formatC(X$expl.var)))
    mtext(paste(colnames(ex), collapse = "                "), line = -4)
    par(cex = 0.8)
    for (i in 1:nrow(ex)) {
      mtext(paste(as.vector(ex[i, ]), collapse = "       "), line = -(i + 
                                                                        3) * 2, adj = 0)
    } 
  } else {
    ##################################################### Block X Mode 1
    plot(X$NPLSDAvariates$X[,pc1],X$NPLSDAvariates$X[,pc2], type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome, 
       xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
       ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
       main = paste("Block X Mode 1: NPLS-DA", sep= " "), cex.main=2,
       cex.lab= 1.2,
       xlim = c(min(X$NPLSDAvariates$X[,pc1]), max(X$NPLSDAvariates$X[,pc1])),
       ylim = c(min(X$NPLSDAvariates$X[,pc2]), max(X$NPLSDAvariates$X[,pc2])))
  #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
  #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
  abline(v=0,h=0)
  #text (X$NPLSDAvariates$X[,1],X$NPLSDAvariates$X[,2], labels=rownames(Factors$Mode1), 
  #      col=DVars$colourOutcome,cex=1, pos=3)
  ########################################################## Per time point
  for(i in 1:length(X$NPLSDAvariatesperMode3)){
    plot(X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc1],X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc2], 
         type = "p", col = "black", pch=21,
         bg=DVars$colourOutcome, 
         xlab = paste("Component ", pc1,sep = ""), 
         ylab = paste("Component ", pc2,sep = ""), las=1,
         main = paste("Block X Mode 1: NPLS-DA mode 3-",i,  sep= " "), cex.main=2,
         cex.lab= 1.2,
         xlim = c(min(X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc1]), max(X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc1])),
         ylim = c(min(X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc2]), max(X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc2])))
    #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
    #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
    abline(v=0,h=0)
    #text (X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc1],X$NPLSDAvariatesperMode3[[i]][[1]]$X[,pc2], 
    #      labels=rownames(Factors$Mode1), 
    #     col=DVars$colourOutcome,cex=1, pos=3)
  }
  #####################################################################################
  ##################################################### Consensus Block X Mode 1
  #plot(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1],X$NPLSDAConsensusvariates$Block.XConsensus[,pc2], 
  #     type = "p", col = "black", pch=21,
  #     bg=DVars$colourOutcome, 
  #     xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
  #     ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
  #     main = paste("Block X Mode 1: NPLS-DA Consensus", sep= " "), cex.main=2,
  #     cex.lab= 1.2,
  #     xlim = c(min(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1])),
  #     ylim = c(min(X$NPLSDAConsensusvariates$Block.XConsensus[,pc2]), max(X$NPLSDAConsensusvariates$Block.XConsensus[,pc2])))
  #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
  #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
  #abline(v=0,h=0)
 # text (X$NPLSDAConsensusvariates$Block.XConsensus[,1],X$NPLSDAConsensusvariates$Block.XConsensus[,2], 
#        labels=rownames(Factors$Mode1), 
 #       col=DVars$colourOutcome,cex=1, pos=3)
  
  
  ######################################## Mode 1 Block Y
  #plot(X$NPLSDAvariates$Y[,pc1],X$NPLSDAvariates$Y[,pc2], type = "p", col = "black", pch=21,
  #     bg=DVars$colourOutcome,
  #     xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,2]*100,digits=2)," %", sep = ""),
  #     ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,2]*100,digits=2)," %", sep = ""), las=1,
  #     main = paste("Block Y Mode 1: NPLS-DA", sep= " "), cex.main=2,
  #     cex.lab= 1.2,
  #     xlim = c(min(X$NPLSDAvariates$Y[,pc1]), max(X$NPLSDAvariates$Y[,pc1])),
  #     ylim = c(min(X$NPLSDAvariates$Y[,pc2]), max(X$NPLSDAvariates$Y[,pc2])))
  #xlim = c((min(x1Y)+ (0.2*min(x1Y))), (max(x1Y)+ (0.2*max(x1Y)))),
  #ylim = c((min(y1Y)+ (0.2*min(y1Y))), (max(y1Y)+ (0.2*max(y1Y)))))
  #abline(v=0,h=0)
  #text (X$NPLSDAvariates$Y[,pc1],X$NPLSDAvariates$Y[,pc2], labels=rownames(FactorsY$Mode1), 
   #     col=DVars$colourOutcome,cex=1, pos=3)
  ######################################## Consensus Mode 1 Block Y
  #plot(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1],X$NPLSDAConsensusvariates$Block.YConsensus[,pc2], 
  #     type = "p", col = "black", pch=21,
  #     bg=DVars$colourOutcome, 
  #     xlab = paste("Component ", pc1,": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," %",sep = ""), 
  #     ylab = paste("Component ", pc2, ": ",round(X$NPLSDAexplVar[pc2,1]*100,digits=2)," %", sep = ""), las=1,
  #     main = paste("Block X Mode 1: NPLS-DA Consensus", sep= " "), cex.main=2,
  #     cex.lab= 1.2,
  #     xlim = c(min(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1])),
  #     ylim = c(min(X$NPLSDAConsensusvariates$Block.YConsensus[,pc2]), max(X$NPLSDAConsensusvariates$Block.YConsensus[,pc2])))
  #xlim = c((min(x1)+ (0.2*min(x1))), (max(x1)+ (0.2*max(x1)))),
  #ylim = c((min(y1)+ (0.2*min(y1))), (max(y1)+ (0.2*max(y1)))))
  #abline(v=0,h=0)
  #text (X$NPLSDAConsensusvariates$Block.YConsensus[,pc1],X$NPLSDAConsensusvariates$Block.YConsensus[,pc2], 
  #      labels=rownames(Factors$Mode1), 
  #      col=DVars$colourOutcome,cex=1, pos=3)
  
  ######################################## Mode 1 Joint Biplot 
  #plot(X$NPLSDAvariates$X[,pc1], X$NPLSDAvariates$Y[,pc1], type = "p", col = "black", pch=21,
  #     bg=DVars$colourOutcome,
  #     xlab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," % of Block X", sep = ""),
  #     ylab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,2]*100,digits=2)," % of Block Y", sep = ""), las=1,
  #     main = paste("Biplot Mode1 Blocks X & Y: NPLS-DA", sep= " "), cex.main=2,
  #     cex.lab= 1.2,
  #     xlim = c(min(X$NPLSDAvariates$X[,pc1]), max(X$NPLSDAvariates$X[,pc1])),
  #     ylim = c(min(X$NPLSDAvariates$Y[,pc1]), max(X$NPLSDAvariates$Y[,pc1])))
  #abline(v=0,h=0)
 # text (X$NPLSDAvariates$X[,1], X$NPLSDAvariates$Y[,1], labels=rownames(Factors$Mode1), 
#        col=DVars$colourOutcome,cex=1, pos=3)
  ######################################## Mode 1 Consensus  Joint Biplot 
  #plot(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1], X$NPLSDAConsensusvariates$Block.YConsensus[,pc1], 
  #     type = "p", col = "black", pch=21,
  #     bg=DVars$colourOutcome,
  #     xlab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,1]*100,digits=2)," % of Block X", sep = ""),
  #     ylab = paste("Component ", pc1, ": ",round(X$NPLSDAexplVar[pc1,2]*100,digits=2)," % of Block Y", sep = ""), las=1,
  #     main = paste("Biplot Mode1 Blocks X & Y: Consensus NPLS-DA", sep= " "), cex.main=2,
  #     cex.lab= 1.2,
  #     xlim = c(min(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.XConsensus[,pc1])),
  #     ylim = c(min(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusvariates$Block.YConsensus[,pc1])))
  #abline(v=0,h=0)
  #text (X$NPLSDAConsensusvariates$Block.XConsensus[,1], X$NPLSDAConsensusvariates$Block.YConsensus[,1], 
  #      labels=rownames(Factors$Mode1), 
  #      col=DVars$colourOutcome,cex=1, pos=3)
  ######################################## Mode2 of X transformed
  #x2trans<-transformedNPLSDAloadings$x
  #y2trans<-transformedNPLSDAloadings$y
  #names<-transformedNPLSDAloadings$names
  #x2trans<-resultwithmatrizYnoNafactors2$transformedNPLSDAloadings$x
  #y2trans<-resultwithmatrizYnoNafactors2$transformedNPLSDAloadings$y
  #names<-resultwithmatrizYnoNafactors2$transformedNPLSDAloadings$names
  
  x2trans<-X$transformedNPLSDAloadings$x
  y2trans<-X$transformedNPLSDAloadings$y
  names<-X$transformedNPLSDAloadings$names
  
  plot(x2trans,y2trans, type = "p", col = "black", pch=21,
      bg="dodgerblue2", 
       cex=1,
  #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
  #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
       xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
       main = paste("Block X Mode 2: NPLS-DA"), 
  #xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
  xlim = c((min(x2trans)), (max(x2trans))), 
  ylim = c((min(y2trans)),(max(y2trans))))
  #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
  #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
  abline(v=0,h=0,col= "gray47", lwd=1)
  #labels=NULL
  #x<-NULL
  #y<-NULL
  #if (is.null(labels)) {
  #  labels <- names
  #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
  #             rownames(Factors$Mode3))
  #}
  #ABS <- abs(x2trans * y2trans)
  #r <- rank(-ABS)
  #value = match(c(1:length(ABS)), r)
  #text (x2,y2, 
  #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
  #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
  #      cex=0.8, pos=3)
  ##for (n in 1:length(ABS)) {
  #  x[c(n)] <- x2trans[value[c(n)]]
  #  y[c(n)] <- y2trans[value[c(n)]]
  #  label[n] <- labels[value[c(n)]]
  #}
  #x <- x[-c((n + 1):(length(x) + 1))]
  #y <- y[-c((n + 1):(length(y) + 1))]
  #text(x, y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>0.05 | abs(y)>0.05),
  #                                                      "dodgerblue2", "transparent"), cex=0.7)
  
  
  
  
  
  
  
  
  
  
  
  ######################################## Just the Mode2 of X
  #plot(X$NPLSDAloadings$X[,pc1],X$NPLSDAloadings$X[,pc2], type = "p", col = "black", pch=21,
  #     bg="dodgerblue2",
   #    xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Block X Mode 2: NPLS-DA", sep= " "), 
  #     xlim = c(min(X$NPLSDAloadings$X[,pc1]), max(X$NPLSDAloadings$X[,pc1])),
  #     ylim = c(min(X$NPLSDAloadings$X[,pc2]), max(X$NPLSDAloadings$X[,pc2])))
       #xlim = c(-1,1),
       #ylim = c(-1,1))
  
  #abline(v=0,h=0)
  #points (X$NPLSDAloadings$Y[,1],X$NPLSDAloadings$Y[,2])
  
  ######################################## Just the Mode2 of Y
  background<-NULL
  background<-c(background[rownames(X$NPLSDAloadings$Y) == 0]<-"dodgerblue2",
                background[rownames(X$NPLSDAloadings$Y) == 1]<-"orange2")
  
  plot(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], type = "p", col = "black", pch=21,
       bg=background,
       xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block Y Mode 2: NPLS-DA", sep= " "), 
       xlim = c(min(X$NPLSDAloadings$Y[,pc1]), max(X$NPLSDAloadings$Y[,pc1])),
       ylim = c(min(X$NPLSDAloadings$Y[,pc2]), max(X$NPLSDAloadings$Y[,pc2])))
  #xlim = c(-1,1),
  #ylim = c(-1,1))
  
  abline(v=0,h=0)
  text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], 
       labels= rownames(X$NPLSDAloadings$Y),
       cex= 1, pos=1, col=background )
  
  ######################################## Biplot Modes 2 and 3 of Block X
  #plot(X$NPLSDAloadings$X[,pc1]*penalty,X$NPLSDAloadings$X[,pc2]*penalty, type = "p", col = "black", pch=21,
  #     bg="dodgerblue2",
  #     xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Biplot Block X Modes 2 & 3: NPLS-DA", sep= " "), 
  #     xlim = c(-1,1),
  #     ylim = c(-1,1))
  #abline(v=0,h=0)
  #arrows (0,0,X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], col = background, pch=21, bg = background)
  #text(X$NPLSDAloadings$Y[,pc1],X$NPLSDAloadings$Y[,pc2], labels= rownames(X$NPLSDAloadings$Y),
  #     cex= 1, pos=1, col=background )
  #####################################################################################################################
  ######################################## Just the Mode2 of X Consensus
  #plot(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1],X$NPLSDAConsensusloadings$Block.XConsensus[,pc2], 
  #     type = "p", col = "black", pch=21,
  #     bg="dodgerblue2",
  #     xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Block X Mode 2: NPLS-DA", sep= " "), 
  #     xlim = c(min(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1]), max(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1])),
  #     ylim = c(min(X$NPLSDAConsensusloadings$Block.XConsensus[,pc2]), max(X$NPLSDAConsensusloadings$Block.XConsensus[,pc2])))
  #xlim = c(-1,1),
  #ylim = c(-1,1))
  
  #abline(v=0,h=0)
  #points (X$NPLSDAloadings$Y[,1],X$NPLSDAloadings$Y[,2])
  ######################################## Just the Mode2 of Y
  #background<-NULL
  #background<-c(background[rownames(X$NPLSDAloadings$Y) == 0]<-"dodgerblue2",
  #              background[rownames(X$NPLSDAloadings$Y) == 1]<-"orange2")
  
  #plot(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1],X$NPLSDAConsensusloadings$Block.YConsensus[,pc2], type = "p", col = "black", pch=21,
  #     bg=background,
  #     xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Block Y Mode 2: NPLS-DA", sep= " "), 
  #     xlim = c(min(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1])),
  #     ylim = c(min(X$NPLSDAConsensusloadings$Block.YConsensus[,pc2]), max(X$NPLSDAConsensusloadings$Block.YConsensus[,pc2])))
  #xlim = c(-1,1),
  #ylim = c(-1,1))
  
  #abline(v=0,h=0)
  #text(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1],X$NPLSDAConsensusloadings$Block.YConsensus[,pc2], 
  #     labels= rownames(X$NPLSDAloadings$Y),
  #     cex= 1, pos=1, col=background )
  
  ######################################## Consensus Biplot Modes 2 and 3 of Block X
  #lot(X$NPLSDAConsensusloadings$Block.XConsensus[,pc1]*penalty,X$NPLSDAConsensusloadings$Block.XConsensus[,pc2]*penalty,
  #     type = "p", col = "black", pch=21,
  #     bg="dodgerblue2",
  #     xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Biplot Block X Modes 2 & 3: Consensus NPLS-DA", sep= " "), 
       #xlim = c(-1,1),
       #ylim = c(-1,1))
  #     xlim = c(min(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1]), max(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1])),
  #     ylim = c(min(X$NPLSDAConsensusloadings$Block.YConsensus[,pc2]), max(X$NPLSDAConsensusloadings$Block.YConsensus[,pc2])))
  #abline(v=0,h=0)
  #arrows (0,0,X$NPLSDAConsensusloadings$Block.YConsensus[,pc1],X$NPLSDAConsensusloadings$Block.YConsensus[,pc2], 
  #        col = background, pch=21, bg = background)
  #text(X$NPLSDAConsensusloadings$Block.YConsensus[,pc1],X$NPLSDAConsensusloadings$Block.YConsensus[,pc2], 
  #     labels= rownames(X$NPLSDAloadings$Y),
  #     cex= 1, pos=1, col=background )
  
  
  
  
  #################################################################  Transformed Mode 2 of X
  #X$transformedNPLSDAloadings
  #x2trans<-X$transformedNPLSDAloadings[,2]
  #y2trans<-X$transformedNPLSDAloadings[,3]
  #names<-X$transformedNPLSDAloadings$names
  
  #plot(x2trans,y2trans, type = "p", col = "black", pch=21,
   #    bg="dodgerblue2", 
  #     cex=1,
       #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
       #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
  #     xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
  #     main = paste("Block X Mode 2: transformed NPLS-DA"), xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
  #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
  #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
  #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
  #abline(v=0,h=0)
  #labels=NULL
  #x<-NULL
  #y<-NULL
  #if (is.null(labels)) {
  #  labels <- names
    #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
    #             rownames(Factors$Mode3))
  #}
  #ABS <- abs(x2trans * y2trans)
  #r <- rank(-ABS)
  #value = match(c(1:length(ABS)), r)
  #text (x2,y2, 
  #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
  #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
  #      cex=0.8, pos=3)
  ##for (n in 1:length(ABS)) {
  #  x[c(n)] <- x2trans[value[c(n)]]
  #  y[c(n)] <- y2trans[value[c(n)]]
  #  label[n] <- labels[value[c(n)]]
  #}
  #x <- x[-c((n + 1):(length(x) + 1))]
  #y <- y[-c((n + 1):(length(y) + 1))]
  #text(x, y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>0.05 | abs(y)>0.05),
  #                                                      "dodgerblue2", "transparent"), cex=0.7)
  
################################################################ Transformed Biplot

  #X$transformedNPLSDAloadings
  #x2trans<-X$transformedNPLSDAloadings[,2]
  #y2trans<-X$transformedNPLSDAloadings[,3]
  #names<-X$transformedNPLSDAloadings$names
  
  #plot(x2trans,y2trans, type = "p", col = "black", pch=21,
  #     bg="dodgerblue2", 
  #     cex=1,
  #     #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
  #     #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
  #     xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
  #     main = paste("Biplot Block X Modes 2 & 3: transformed NPLS-DA"), xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
  #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
  #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
  #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
  #abline(v=0,h=0)
  #labels=NULL
  #x<-NULL
  #y<-NULL
  #if (is.null(labels)) {
  #  labels <- names
    #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
    #             rownames(Factors$Mode3))
  #}
  #ABS <- abs(x2trans * y2trans)
  #r <- rank(-ABS)
  #value = match(c(1:length(ABS)), r)
                #text (x2,y2, 
                #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
                #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
                #      cex=0.8, pos=3)
   #             for (n in 1:length(ABS)) {
  #                x[c(n)] <- x2trans[value[c(n)]]
   #               y[c(n)] <- y2trans[value[c(n)]]
  #                label[n] <- labels[value[c(n)]]
   #             }
    #            x <- x[-c((n + 1):(length(x) + 1))]
     #           y <- y[-c((n + 1):(length(y) + 1))]
    #            text(x, y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>0.05 | abs(y)>0.05),
    #                                                                    "dodgerblue2", "transparent"), cex=0.7)
                
     #           arrows (0,0,X$NPLSDAloadings$Y[,1],X$NPLSDAloadings$Y[,2], col = background, pch=21, bg = background)
    #            text(X$NPLSDAloadings$Y[,1],X$NPLSDAloadings$Y[,2], labels= rownames(X$NPLSDAloadings$Y),
     #                cex= 1, pos=1, col=background )             
  
  
  
  
  ######################################## 
  ###########      Mode 3     ############
  ######################################## 
  # Without variable vectors
  plot(x3, y3, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block X Mode 3: NPLS-DA"), 
       #xlim = c((min(x3)), (max(x3))),
       #ylim = c((min(y3)), (max(y3))))
       xlim = c((min(x3)), (max(x3)+ abs(0.2*max(x3)))),
       ylim = c((min(y3)), (max(y3)+ abs(0.2*max(y3)))))
  abline(v=0,h=0)
  text (x3,y3, labels=rownames(FactorsX$Mode3), col="dodgerblue2",cex=1, pos=3)
  
  ######################################## 
  #plot(x3Y, y3Y, type = "p", col = "black", pch=21,bg="dodgerblue2", xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Block Y Mode 3: NPLS-DA"), 
  #     xlim = c((min(x3Y)), (max(x3Y)+ abs(0.2*max(x3Y)))),
  #     ylim = c((min(y3Y)), (max(y3Y)+ abs(0.2*max(y3Y)))))
  #abline(v=0,h=0)
  #text (x3Y,y3Y, labels=rownames(FactorsY$Mode3), col="dodgerblue2",cex=1, pos=3)
  
  # With variable vectors
  #plot(x3, y3, type = "p", col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1, sep = " "), las=1,
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Biplot Block X Modes 2 and 3: NPLS-DA"), 
  #     xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
  #abline(v=0,h=0)
  #text (x3,y3, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1, pos=3)
  #arrows(0,0,penalty*x2,penalty*y2, pch=16, lwd=2,
  #       col = ifelse ((abs(x2)>(penalty/2)*0.05 | abs(y2)>(penalty/2)*0.05),"orange2", "transparent") )
  #labels=NULL
  #label<-NULL
  #x<-NULL
  #y<-NULL
  #if (is.null(labels)) {
  #  labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2), 
  #                 rownames(Factors$Mode3))
  #}
  #ABS <- abs(x2 * y2)
  #r <- rank(-ABS)
  #value = match(c(1:length(ABS)), r)
  #for (n in 1:length(ABS)) {
  #  x[c(n)] <- Factors[[2]][value[c(n)], 1]
  #  y[c(n)] <- Factors[[2]][value[c(n)], 2]
  #  label[n] <- labels[[2]][value[c(n)]]
  #}
  #x <- x[-c((n + 1):(length(x) + 1))]
  #y <- y[-c((n + 1):(length(y) + 1))]
  #text(penalty*x, penalty*y, labels = label, adj = c(0, -0.5),col= ifelse ((abs(x)>(penalty/2)*0.05 | abs(y)>(penalty/2)*0.05),
  #                                                       "orange2", "transparent"), cex=0.7)
  
  
  #text (x2,y2, labels= ifelse ((abs(x2)>0.05 | abs(y2)>0.05),names(x2), ""),
  #      col = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"orange2", "transparent"),    
  #      cex=1.5, pos=1)
  
  ######################################## 
  #plot(x3Y, y3Y, type = "p",  col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1, sep = " "), las=1,
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Biplot Block Y Modes 2 and 3: NPLS-DA"), 
  #     xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
  #abline(v=0,h=0)
  #text (x3Y,y3Y, labels=rownames(FactorsY$Mode3), col="dodgerblue2",cex=1, pos=3)
  #arrows(0,0,penalty*x2Y,penalty*y2Y, pch=16, lwd=2,
  #       col = ifelse ((abs(x2Y)>(penalty/2)*0.05 | abs(y2Y)>(penalty/2)*0.05),"orange2", "transparent"))
  #labels=NULL
  #x<-NULL
  #y<-NULL
  #if (is.null(labels)) {
  #  labels <- list(rownames(FactorsY$Mode1), rownames(FactorsY$Mode2), 
  #                 rownames(FactorsY$Mode3))
  #}
  #ABS <- abs(x2Y * y2Y)
  #r <- rank(-ABS)
  #value = match(c(1:length(ABS)), r)
  #for (n in 1:length(ABS)) {
  #  x[c(n)] <- FactorsY[[2]][value[c(n)], 1]
  #  y[c(n)] <- FactorsY[[2]][value[c(n)], 2]
  #  label[n] <- labels[[2]][value[c(n)]]
  #}
  #x <- x[-c((n + 1):(length(x) + 1))]
  #y <- y[-c((n + 1):(length(y) + 1))]
  #text(penalty*x, penalty*y, labels = label, adj = c(0, -0.5),col=ifelse ((abs(x)>(penalty/2)*0.05 | abs(y)>(penalty/2)*0.05),
  #                                                                        "orange2", "transparent"), cex=0.7)
  
  
  #text (x2Y,y2Y, labels= ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),names(x2Y), ""),
  #      col = ifelse ((abs(x2Y)>0.05 | abs(y2Y)>0.05),"orange2", "transparent"),    
  #      cex=1, pos=1)
  
  #    Joint Mode 3 ####################################### 
  #plot(x3, x3Y, type = "p", col = "black", pch=21,bg="dodgerblue2", cex=1, xlab = paste("Component", pc1,"Block X", sep = " "), 
  #     ylab = paste("Component", pc1, "Block Y",sep = " "), las =1,
  #     main = paste("Joint plot: NPLS-DA Mode 3"), 
  #     xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
  #abline(v=0,h=0)
  #abline(a=0,b=1,col= "gray47", lwd=1)
  #text (x3,x3Y, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1, pos=3)
  
  ############################################################
  ex <- as.matrix(explcore(X$G, factors, n=5))
  plot.new()
  mtext(paste("Expl.var of X: ", round(sum(X$NPLSDAexplVar[c(1:factors),1])*100,digits=2),"%"))
  mtext(paste("Expl.var of Y by X: ", round(sum(X$NPLSDAexplVar[c(1:factors),2])*100,digits=2),"%"),line= -2)
  #mtext(paste("model expl.var: ", sum(as.numeric(ex[c(1:2),4]))))
  #mtext(paste("model expl.var: ", formatC(X$expl.var)))
  mtext(paste(colnames(ex), collapse = "                "), line = -4)
  par(cex = 0.8)
  for (i in 1:nrow(ex)) {
    mtext(paste(as.vector(ex[i, ]), collapse = "       "), line = -(i + 
                                                                      3) * 2, adj = 0)
  } 
}
}

##################################################################################################################
# Plotting things
plotSpacemod3<-function (X, what, DVs, Bestfittedmodel, COMP, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                         cutoff = length(X)) 
  # Function to plot all necesary figures to interpret NPLS-DA results
  # Input: Array of the best model
  # DVs: Are the descriptive variables
  # Output:
  # Plot of the loadings per variable
  # Heatmap of the behavior of all variables in all individuals
  # Mode1 NPLS-DA relative to X explaining Y
  # Mode1 NPLS-DA relative to X explaining Y colored by AgeGroup
  # Mode1 NPLS-DA relative to X explaining Y colored by FirstAAb that appears
  # Mode3 NPLS-DA relative to X explaining Y 
  # Mode2 NPLS-DA relative to X explaining Y 
  # Analysis of the Core Array of the model to determine explained variances and the contributions of the components in each mode
# Here you can play with the script to obtain what you desire.

{
  library(sfsmisc)
  par(mar=c(9,9,9,9))
  A3D<-X$A3D
  xdims<-dim(A3D)
  XN2D <- matrix(A3D, xdims[1], xdims[2] * xdims[3])
  colnames(A3D)
  colp1<-paste(colnames(A3D),"-12", sep="_")
  colp2<-paste(colnames(A3D),"-9", sep="_")
  colp3<-paste(colnames(A3D),"-6", sep="_")
  colp4<-paste(colnames(A3D),"-3", sep="_")
  colp5<-paste(colnames(A3D),"0", sep="_")
  colnames(XN2D)<-c(colp1,colp2,colp3,colp4,colp5)
  plsdaforplotting<-plsda (X=XN2D, Y= as.vector(Outcomedummy),ncomp = 3)
  plotIndiv(plsdaforplotting)
  #plotVar(plsdaforplotting)
  #plotLoadings(plsdaforplotting)
  #cim(plsdaforplotting)
  #######
  IDs<- data.frame(Mask.Id=DVs$Mask.Id)
  DVars<-DVs[!duplicated(DVs$Mask.Id),]
  #dim(DVars)
  #length(unique(DVars$Mask.Id))
  #colnames(DVars)
  #plsdaforplotting$loadings
  par(mar=c(9,9,9,9))
  plotIndiv(plsdaforplotting,group = DVars$Outcome, centroid = FALSE, comp = c(1,2), 
            star=FALSE, ellipse = FALSE,col.per.group = c("dodgerblue2","orangered3"),
            title = "Block X Mode 1: OutcomeNPLSDA colored by Outcome",
            rep.space= 'X-variate', style = "graphics",
            legend = TRUE,ind.names=DVars$Outcome,cex=1
  )
  
  plotIndiv(plsdaforplotting,group = DVars$AgeGroup, centroid = FALSE, title = "Block X Mode 1: Individuals colored by AgeGroup",
            star=FALSE, ellipse = TRUE,col.per.group = c("mediumorchid4","seagreen","orangered3","dodgerblue3"),
            legend = TRUE, style="graphics",cex=0.7
  )
  
  plotIndiv(plsdaforplotting,group = DVars$FirstAAbCC, centroid = FALSE, title = "Block X Mode 1: Individuals colored by FirstAAb",
            star=FALSE, ellipse = TRUE,col.per.group = c("dodgerblue2", "darkolivegreen3", "darkorchid3","darkslategray3",
                                                         "mediumpurple1","lightsalmon2","orange1"),
            legend = TRUE, style="graphics",cex=0.7
  )
  
  #
  pc1 <- PCs[1]
  pc2 <- PCs[2]
  a = paste("Factors", what, sep = "")
  label = NULL
  Factors <- X[[a]]
  x1 <- Factors[[1]][, pc1]
  y1 <- Factors[[1]][, pc2]
  x2 <- Factors[[2]][, pc1]
  y2 <- Factors[[2]][, pc2]
  x3 <- Factors[[3]][, pc1]
  y3 <- Factors[[3]][, pc2]
  par(mfrow=c(1,1))
  
  DVars$colourOutcome<- ""
  DVars$colourOutcome[DVars$Outcome == 0]<-"dodgerblue2"
  DVars$colourOutcome[DVars$Outcome == 1]<-"orangered3"
  
  
  DVars$colourGender<- ""
  DVars$colourGender[DVars$Gender == "Male"]<-"dodgerblue2"
  DVars$colourGender[DVars$Gender == "Female"]<-"orangered3"
  
  
  
  DVars$colourAAb<- ""
  DVars$colourAAb[DVars$FirstAAbCC == 1]<-"dodgerblue2"
  DVars$colourAAb[DVars$FirstAAbCC == 2]<-"darkolivegreen3"#"transparent"    
  DVars$colourAAb[DVars$FirstAAbCC == 3]<-"darkorchid3"
  DVars$colourAAb[DVars$FirstAAbCC == 4]<-"darkslategray3"
  DVars$colourAAb[DVars$FirstAAbCC == 5]<-"mediumpurple1"
  DVars$colourAAb[DVars$FirstAAbCC == 6]<-"lightsalmon2"
  DVars$colourAAb[DVars$FirstAAbCC == 7]<-"orange1"
  
  DVars$colourAge<- ""
  DVars$colourAge[DVars$AgeGroup == 1]<-"mediumorchid4"
  DVars$colourAge[DVars$AgeGroup == 2]<-"seagreen"#"transparent"#
  DVars$colourAge[DVars$AgeGroup == 3]<-"orangered3"
  DVars$colourAge[DVars$AgeGroup == 4]<-"dodgerblue3"
  
  
  ######################################## 
  plot(x3, y3, type = "p", col = "dodgerblue2", pch=16, xlab = paste("Component", pc1, sep = " "), 
       ylab = paste("Component", pc2, sep = " "), 
       main = paste("Block X Mode 3: NPLSDA Time"), 
       # xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
       #xlim = c(-0.5, 0), ylim = c(-0.6, 0))
       #xlim = c((min(x3)), (max(x3))), ylim = c((min(y3)),(max(y3))))
       #ylim = c((min(y3)+ (0.2*min(y3))), 0.1))
       xlim = c((min(x3)+ (0.2*min(x3))), (max(x3)+ (0.2*max(x3)))),
       ylim = c((min(y3)+ (0.2*min(y3))), (max(y3)+ (0.2*max(y3)))))
  
  abline(v=0,h=0)
  text (x3,y3, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1.5, pos=3)
  ######################################## 
  tabla<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,2))
  head(tabla)
  x2<-tabla[,c(4,1)]
  y2<-tabla[,c(4,2)]
  
  #x2<-tabla[,1]
  #y2<-tabla[,2]
  #x2 <- Factors[[2]][, pc1]
  #y2 <- Factors[[2]][, pc2]
  names<-tabla$names
  
  plot(x2[,2],y2[,2], type = "p", pch=16, 
       col ="dodgerblue2",
       cex=1,
       #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
       #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
       xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
       main = paste("Mode 2: OutcomeNPLSDA"), xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
  #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
  #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
  #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
  abline(v=0,h=0)
  labels=NULL
  x<-NULL
  y<-NULL
  if (is.null(labels)) {
    labels <- names
    #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
     #             rownames(Factors$Mode3))
  }
  ABS <- abs(x2[,2] * y2[,2])
  r <- rank(-ABS)
  value = match(c(1:cutoff), r)
  #text (x2,y2, 
  #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
  #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
  #      cex=0.8, pos=3)
  for (n in 1:cutoff) {
    x[c(n)] <- x2[,2][value[c(n)]]
    y[c(n)] <- y2[,2][value[c(n)]]
    label[n] <- labels[value[c(n)]]
  }
  x <- x[-c((n + 1):(length(x) + 1))]
  y <- y[-c((n + 1):(length(y) + 1))]
  text(x, y, labels = label, adj = c(0, -1),col="dodgerblue2", cex=0.7)
  ############################################################
  
  P=COMP[1];Q=COMP[2];R=COMP[3]
  EV<-Bestfittedmodel$summarytable$Fitper[Bestfittedmodel$summarytable[c("Model")]==paste(c(P,Q,R),collapse="")]
  ex <- as.matrix(explcore(X$G, factors=3, n=5))
  plot.new()
  mtext(paste("Expl.var of X: ", round(EV,digits=2),"%"))
  mtext(paste("Expl.var of Y by X: ", round(sum(as.numeric(ex[c(pc1,pc2),4])),digits=2),"%"), line = -2)
  #mtext(paste("Expl.var of Y by X: ", (sum(X$explvar[,3])*100),"%"),line= -2)
  #mtext(paste("model expl.var: ", sum(as.numeric(ex[c(1:2),4]))))
  #mtext(paste("model expl.var: ", formatC(X$expl.var)))
  mtext(paste(colnames(ex), collapse = "                "), line = -4)
  par(cex = 0.8)
  for (i in 1:nrow(ex)) {
    mtext(paste(as.vector(ex[i, ]), collapse = "       "), line = -(i + 
                                                                      3) * 2, adj = 0)
  }
  ###############################################
  if(what == "Y"){
    pc1 <- PCs[1]
    pc2 <- PCs[2]
    a = paste("Factors", what, sep = "")
    label = NULL
    Factors <- X[[a]]
    x2 <- Factors[[2]][, pc1]
    y2 <- Factors[[2]][, pc2]
    par(mfrow=c(1,1))
    plot(x2, y2, type = "n", pch=16, 
         col ="dodgerblue2",
         cex = 1,
         #col =ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"dodgerblue2", "gray"),
         #cex = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),1, 0.5),
         xlab = paste("Component", pc1, sep = " "), ylab = paste("Component", pc2, sep = " "), 
         main = paste("Mode 2: OutcomeNPLS: Outcome"), 
         xlim = c(-1,1.2), ylim = c(-1.2, 1))
    #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
    #xlim = c((min(x2)+ (0.002*min(x2))), (max(x2)+ (0.002*max(x2)))), 
    #ylim = c((min(y2)+ (0.002*min(y2))), (max(y2)+ (0.8*max(y2)))))
    #xlim = c((min(x2)+ (0.002*min(x2))), (max(x2)+ 0.2)), 
    #ylim = c((min(y2)+ (0.2*min(y2))), 0.1))
    arrows(0,0,x2,y2 ,col = "dodgerblue2")
    
    abline(v=0,h=0)
    text (x2,y2, 
          #labels= ifelse ((abs(x2)>0.05 | abs(y2)>0.05),rownames(Factors$Mode2), ""),
          col = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"dodgerblue2", "transparent"),    
          labels= "IA Disease",
          #col = "dodgerblue2", 
          cex=1, pos=1)
    
    #
  }
}



##############################################################################################################
#To create figures by timepoint, meaning just 2D tables
plotSpacemod2D<-function (X, what="X", DVs, Bestfittedmodel, COMP, PCs = c(1, 2), labels = NULL, main = substitute(X), 
                         cutoff = length(X)) 
  # Function to plot all necesary figures to interpret NPLS-DA results
  # Input: Array of the best model
  # DVs: Are the descriptive variables
  # Output:
  # Plot of the loadings per variable
  # Heatmap of the behavior of all variables in all individuals
  # Mode1 NPLS-DA relative to X explaining Y
  # Mode1 NPLS-DA relative to X explaining Y colored by AgeGroup
  # Mode1 NPLS-DA relative to X explaining Y colored by FirstAAb that appears
  # Mode3 NPLS-DA relative to X explaining Y 
  # Mode2 NPLS-DA relative to X explaining Y 
# Analysis of the Core Array of the model to determine explained variances and the contributions of the components in each mode
# Here you can play with the script to obtain what you desire.

{
  library(sfsmisc)
  par(mar=c(9,9,9,9))
  A3D<-X
  A3Dm12<-A3D[,,1]; A3Dm9<-A3D[,,2];A3Dm6<-A3D[,,3]; A3Dm3<-A3D[,,4];A3Dm0<-A3D[,,5]
  A3Dto2D<-list(A3Dm12,A3Dm9,A3Dm6,A3Dm3,A3Dm0)
  for (i in 1:length(A3Dto2D)){
    
  plsdaforplotting<-plsda (X=A3Dto2D[[i]], Y= Outcomedummyarray[,,1],ncomp = 3)
  plotIndiv(plsdaforplotting)
  #plotVar(plsdaforplotting)
  #plotLoadings(plsdaforplotting)
  #cim(plsdaforplotting)
  #######
  IDs<- data.frame(Mask.Id=DVs$Mask.Id)
  DVars<-DVs[!duplicated(DVs$Mask.Id),]
  #dim(DVars)
  #length(unique(DVars$Mask.Id))
  #colnames(DVars)
  #plsdaforplotting$loadings
  par(mar=c(9,9,9,9))
  plotIndiv(plsdaforplotting,group = DVars$Outcome, centroid = FALSE, comp = c(1,2), 
            star=FALSE, ellipse = FALSE,col.per.group = c("dodgerblue2","orangered3"),
            title = "Block X Mode 1: OutcomeNPLSDA colored by Outcome",
            rep.space= 'X-variate', style = "graphics",
            legend = TRUE,ind.names=DVars$Outcome,cex=1
  )
  
  plotIndiv(plsdaforplotting,group = DVars$AgeGroup, centroid = FALSE, title = "Block X Mode 1: Individuals colored by AgeGroup",
            star=FALSE, ellipse = TRUE,col.per.group = c("mediumorchid4","seagreen","orangered3","dodgerblue3"),
            legend = TRUE, style="graphics",cex=0.7
  )
  
  plotIndiv(plsdaforplotting,group = DVars$FirstAAbCC, centroid = FALSE, title = "Block X Mode 1: Individuals colored by FirstAAb",
            star=FALSE, ellipse = TRUE,col.per.group = c("dodgerblue2", "darkolivegreen3", "darkorchid3","darkslategray3",
                                                         "mediumpurple1","lightsalmon2","orange1"),
            legend = TRUE, style="graphics",cex=0.7
  )
  
  #
  pc1 <- PCs[1]
  pc2 <- PCs[2]
  a = paste("Factors", what, sep = "")
  label = NULL
  Factors <- X[[a]]
  x1 <- Factors[[1]][, pc1]
  y1 <- Factors[[1]][, pc2]
  x2 <- Factors[[2]][, pc1]
  y2 <- Factors[[2]][, pc2]
  x3 <- Factors[[3]][, pc1]
  y3 <- Factors[[3]][, pc2]
  par(mfrow=c(1,1))
  
  DVars$colourOutcome<- ""
  DVars$colourOutcome[DVars$Outcome == 0]<-"dodgerblue2"
  DVars$colourOutcome[DVars$Outcome == 1]<-"orangered3"
  
  
  DVars$colourGender<- ""
  DVars$colourGender[DVars$Gender == "Male"]<-"dodgerblue2"
  DVars$colourGender[DVars$Gender == "Female"]<-"orangered3"
  
  
  
  DVars$colourAAb<- ""
  DVars$colourAAb[DVars$FirstAAbCC == 1]<-"dodgerblue2"
  DVars$colourAAb[DVars$FirstAAbCC == 2]<-"darkolivegreen3"#"transparent"    
  DVars$colourAAb[DVars$FirstAAbCC == 3]<-"darkorchid3"
  DVars$colourAAb[DVars$FirstAAbCC == 4]<-"darkslategray3"
  DVars$colourAAb[DVars$FirstAAbCC == 5]<-"mediumpurple1"
  DVars$colourAAb[DVars$FirstAAbCC == 6]<-"lightsalmon2"
  DVars$colourAAb[DVars$FirstAAbCC == 7]<-"orange1"
  
  DVars$colourAge<- ""
  DVars$colourAge[DVars$AgeGroup == 1]<-"mediumorchid4"
  DVars$colourAge[DVars$AgeGroup == 2]<-"seagreen"#"transparent"#
  DVars$colourAge[DVars$AgeGroup == 3]<-"orangered3"
  DVars$colourAge[DVars$AgeGroup == 4]<-"dodgerblue3"
  
  
  ######################################## 
  #plot(x3, y3, type = "p", col = "dodgerblue2", pch=16, xlab = paste("Component", pc1, sep = " "), 
  #     ylab = paste("Component", pc2, sep = " "), 
  #     main = paste("Block X Mode 3: NPLSDA Time"), 
       # xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1))
       #xlim = c(-0.5, 0), ylim = c(-0.6, 0))
       #xlim = c((min(x3)), (max(x3))), ylim = c((min(y3)),(max(y3))))
       #ylim = c((min(y3)+ (0.2*min(y3))), 0.1))
  #     xlim = c((min(x3)+ (0.2*min(x3))), (max(x3)+ (0.2*max(x3)))),
  #     ylim = c((min(y3)+ (0.2*min(y3))), (max(y3)+ (0.2*max(y3)))))
  
  #abline(v=0,h=0)
  #text (x3,y3, labels=rownames(Factors$Mode3), col="dodgerblue2",cex=1.5, pos=3)
  ######################################## 
  tabla<-plotVar(plsdaforplotting,plot=FALSE,comp=c(1,2))
  head(tabla)
  x2<-tabla[,c(4,1)]
  y2<-tabla[,c(4,2)]
  
  #x2<-tabla[,1]
  #y2<-tabla[,2]
  #x2 <- Factors[[2]][, pc1]
  #y2 <- Factors[[2]][, pc2]
  names<-tabla$names
  
  plot(x2[,2],y2[,2], type = "p", pch=16, 
       col ="dodgerblue2",
       cex=1,
       #col =ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "gray"),
       #cex = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),1, 0.5),
       xlab = paste("Component 1", sep = " "), ylab = paste("Component 2", sep = " "), 
       main = paste("Mode 2: OutcomeNPLSDA"), xlim = c(-1.2,1.2), ylim = c(-1.2, 1.2))
  #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
  #xlim = c((min(x2)+ (0.2*min(x2))), (max(x2)+ (0.2*max(x2)))), 
  #ylim = c((min(y2)+ (0.2*min(y2))), (max(y2)+ (0.2*max(y2)))))
  abline(v=0,h=0)
  labels=NULL
  x<-NULL
  y<-NULL
  if (is.null(labels)) {
    labels <- names
    #labels <- list(rownames(Factors$Mode1), rownames(Factors$Mode2),
    #             rownames(Factors$Mode3))
  }
  ABS <- abs(x2[,2] * y2[,2])
  r <- rank(-ABS)
  value = match(c(1:cutoff), r)
  #text (x2,y2, 
  #      labels= ifelse ((abs(x2)>0.9 | abs(y2)>0.9),rownames(tabla), ""),
  #      col = ifelse ((abs(x2)>0.9 | abs(y2)>0.9),"dodgerblue2", "transparent"),    
  #      cex=0.8, pos=3)
  for (n in 1:cutoff) {
    x[c(n)] <- x2[,2][value[c(n)]]
    y[c(n)] <- y2[,2][value[c(n)]]
    label[n] <- labels[value[c(n)]]
  }
  x <- x[-c((n + 1):(length(x) + 1))]
  y <- y[-c((n + 1):(length(y) + 1))]
  text(x, y, labels = label, adj = c(0, -1),col="dodgerblue2", cex=0.7)
  ############################################################
  
  #P=COMP[1];Q=COMP[2];R=COMP[3]
  #EV<-Bestfittedmodel$summarytable$Fitper[Bestfittedmodel$summarytable[c("Model")]==paste(c(P,Q,R),collapse="")]
  #ex <- as.matrix(explcore(X$G, factors=3, n=5))
  #plot.new()
  #mtext(paste("Expl.var of X: ", round(EV,digits=2),"%"))
  #mtext(paste("Expl.var of Y by X: ", round(sum(as.numeric(ex[c(pc1,pc2),4])),digits=2),"%"), line = -2)
  #mtext(paste("Expl.var of Y by X: ", (sum(X$explvar[,3])*100),"%"),line= -2)
  #mtext(paste("model expl.var: ", sum(as.numeric(ex[c(1:2),4]))))
  #mtext(paste("model expl.var: ", formatC(X$expl.var)))
  #mtext(paste(colnames(ex), collapse = "                "), line = -4)
  #par(cex = 0.8)
  #for (i in 1:nrow(ex)) {
   # mtext(paste(as.vector(ex[i, ]), collapse = "       "), line = -(i + 
  #                                                                    3) * 2, adj = 0)
  #}
  ###############################################
  if(what == "Y"){
    pc1 <- PCs[1]
    pc2 <- PCs[2]
    a = paste("Factors", what, sep = "")
    label = NULL
    Factors <- X[[a]]
    x2 <- Factors[[2]][, pc1]
    y2 <- Factors[[2]][, pc2]
    par(mfrow=c(1,1))
    plot(x2, y2, type = "n", pch=16, 
         col ="dodgerblue2",
         cex = 1,
         #col =ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"dodgerblue2", "gray"),
         #cex = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),1, 0.5),
         xlab = paste("Component", pc1, sep = " "), ylab = paste("Component", pc2, sep = " "), 
         main = paste("Mode 2: OutcomeNPLS: Outcome"), 
         xlim = c(-1,1.2), ylim = c(-1.2, 1))
    #xlim = c((min(x2)), (max(x2))), ylim = c((min(y2)),(max(y2))))
    #xlim = c((min(x2)+ (0.002*min(x2))), (max(x2)+ (0.002*max(x2)))), 
    #ylim = c((min(y2)+ (0.002*min(y2))), (max(y2)+ (0.8*max(y2)))))
    #xlim = c((min(x2)+ (0.002*min(x2))), (max(x2)+ 0.2)), 
    #ylim = c((min(y2)+ (0.2*min(y2))), 0.1))
    arrows(0,0,x2,y2 ,col = "dodgerblue2")
    
    abline(v=0,h=0)
    text (x2,y2, 
          #labels= ifelse ((abs(x2)>0.05 | abs(y2)>0.05),rownames(Factors$Mode2), ""),
          col = ifelse ((abs(x2)>0.05 | abs(y2)>0.05),"dodgerblue2", "transparent"),    
          labels= "IA Disease",
          #col = "dodgerblue2", 
          cex=1, pos=1)
    
    #
  }
}
}

