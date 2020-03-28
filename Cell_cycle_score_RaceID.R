
Scoring <- function(object, features, name=NULL, nbin = 24, ctrl= 100, n=sc@genes, st=NULL, nd=NULL, null=NULL ) {
  if ( is.null(features)) {
    stop("Input feature list is missing")
  }
  set.seed(1337)
  fdata <- object@ndata * min(object@counts)
  features.old <- features
  features <- lapply(seq_along(features), function(x, y, i) {
    z <- intersect(x[[i]], n)
    if ( length(z) == 0) {
      stop(paste("no feature of list element named ", y[[i]], " is element of n", sep = ""))
    } 
    return(z)
  }, y = names(features), x=features)
  names(features) <- names(features.old)
  feature.length <- length(features)
  universe <- n
  mean.universe <- apply(fdata[n,], 1, mean )
  mean.universe <- mean.universe[order(mean.universe)]
  mean.universe.binned <- ggplot2::cut_number(mean.universe + rnorm(length(mean.universe))/1e+30, n = nbin, labels=F, right=F)
  names(mean.universe.binned) <- names(mean.universe)
  
  ctrl.use <- vector(mode = "list", length = feature.length)
  for (i in 1:feature.length) {
    features.use <- features[[i]]   
    for (j in 1:length(x = features.use)) { 
      ctrl.use[[i]] <- c(ctrl.use[[i]], names(sample(mean.universe.binned[which( mean.universe.binned == 
                                                                                           as.numeric(mean.universe.binned[features.use[j]]))], size = ctrl, replace = FALSE))) ### samples for each feature gene 100 control genes from the same bin as feature gene  
    }
  }
  ctrl.use <- lapply(ctrl.use, unique)
  ctrl.scores <- matrix(numeric(length = 1L), length(ctrl.use), ncol(fdata))
  
  for (i in 1:length(ctrl.use)) {
    features.use <- ctrl.use[[i]]
    ctrl.scores[i, ] <- apply(fdata[features.use, ], 2, mean)    #### means across s.control or m.control genes for each cell
  }  
  features.scores <- matrix(numeric(length = 1L), feature.length, ncol(fdata))
  for (i in 1:feature.length) {
    features.use <- features[[i]]
    data.use <- fdata[features.use, , drop = FALSE]
    features.scores[i, ] <- apply(data.use, 2, mean) ### mean of feature genes across cells 
  }
  features.scores.use <- features.scores - ctrl.scores
  rownames(features.scores.use) <- paste0(names(features), "_factor")
  features.scores.use <- as.data.frame(t(features.scores.use))
  rownames(features.scores.use) <- colnames(fdata)
  
  if (length(features) == 1) { 
    assignments <- apply( features.scores.use, 1, function(scores) {
      if ( scores  < 0 ) { return(nd)}
      else { return(st)}
    })
  } 
  
  if ( length(features) > 1){
    if ( is.null(null)) {
      stop("set null")
    }
    assignments <- apply(features.scores.use, 1, FUN = function(scores) {
      if (all(scores < 0)) {
        return(null)
      }
      else {
        if (length(which(scores == max(scores))) > 1) {
          return("Undecided")
        }
        else {
          return(c(st, nd)[which(scores == max(scores))])
        }
      }
    })
  }
  features.scores.use <- cbind(features.scores.use, Phase=assignments)
  
  return(features.scores.use)
   
}
