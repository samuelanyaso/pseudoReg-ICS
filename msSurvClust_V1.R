## creates empty matrices to store the results of the dNs, Ys, and sum_dNs
## the dimension of each cluster specific matrix must conform
## to the dimension of the full matrix
FMat <- function(tree, tree0, Data, nt.states){
  # unique transition times from the full data set
  times <- sort(unique(Data$stop[!Data$stop == 0]))
  lng <- sapply(edges(tree0)[nodes(tree0) %in% nt.states], length) # gets the number of edges for each non-absorbing state
  # creates the possible transitions
  ds <- paste("dN", rep(nodes(tree0)[nodes(tree0) %in% nt.states],lng),  # replicates the non-absorbing states by the number of edges out of each state
              unlist(edges(tree0)[nodes(tree0) %in% nt.states])) # gets the edges proceeding out of each non-absorbing state
  nt.states2 <- nt.states[!nt.states == "LT"] # gets the non-absorbing states other than the LT state
  ys <- paste("y", nt.states2) # colnames for the risk set
  # matrix of possible transition, initialize to zeros
  dNs <- matrix(0, nrow = length(times), ncol = length(ds))
  # matrix of total possible transitions from each non-absorbing state, initialize to zeros
  sum_dNs <- matrix(0, nrow = length(times), ncol = length(nt.states))
  # matrix of at-risk sets for each stage at each time
  Ys <- matrix(NA, nrow = length(times), ncol = length(ys))
  rownames(dNs) <- rownames(sum_dNs) <- rownames(Ys) <- times
  colnames(dNs) <- ds
  colnames(Ys) <- ys
  colnames(sum_dNs) <- paste("dN", nt.states, ".")
  rownames(sum_dNs) <- rownames(dNs)
  colnames(sum_dNs) <- paste("dN", nt.states, ".")
  ans <- list(nt.states2=nt.states2, dNs = dNs, Ys = Ys, sum_dNs = sum_dNs)
  return(ans)
}

# adding start times - all subjects start from time 0
Add.start <- function (Data){
  Data$start <- 0  # here, we make all subject to start at time 0
  idx <- which(table(Data$id) > 1) 
  
  for (i in names(idx)) {
    ab <- Data[which(Data$id == i), ] 
    ab <- with(ab, ab[order(ab$stop), ]) 
    ab2 <- which(Data$id == i) 
    start2 <- vector(length = length(ab2)) 
    start2[1] <- 0 # initialize the start at 0
    start2[2:length(ab2)] <- ab$stop[1:length(ab2) - 1] 
    Data$start[ab2] <- start2 
  }
  new.data <- data.frame(cID=Data$cID, id = Data$id, csize=Data$csize, start = Data$start, 
                         stop = Data$stop, start.stage = Data$start.stage, end.stage = Data$end.stage)
  res <- new.data
}

## converting for censoring and/or LT node information to the tree
Add.States <-function (tree, LT){
  # adds censoring state (state 0) to nodes and edges
  Nodes <- c(nodes(tree), "0") # adds the censoring state (node)
  Edges <- edges(tree) # gets the edges in the tree
  Edges[["0"]] <- character(0) # adds 'edges' for the censoring node
  
  nt.states <- names(which(sapply(Edges, function(x) length(x) > 0))) # gets the non-absorbing states
  
  # adds the edge to the censoring state for the non-absorbing state. i.e. a subject can be censored from any non-absorbing state
  for (stage in nt.states) {
    Edges[[stage]] <- c("0", Edges[[stage]])
  }
  # constructs the tree for the censored data
  tree0 <- new("graphNEL", nodes = Nodes, edgeL = Edges, edgemode = "directed")
  
  # adds left truncated state
  if (LT) {
    Nodes <- c(nodes(tree0), "LT") # # adds the LT state (node)
    Edges[["LT"]] <- nt.states  # set the edges of the LT state to go into the non-absorbing states
    nt.states.LT <- names(which(sapply(Edges, function(x) length(x) > 0)))  # gets the non-absorbing states
    treeLT <- new("graphNEL", nodes = Nodes, edgeL = Edges, edgemode = "directed") # constructs the tree for the LT data
  }
  if (LT) {
    return(list(tree0 = tree0, nt.states = nt.states, nt.states.LT = nt.states.LT, 
                treeLT = treeLT))
  }
  else {
    return(list(tree0 = tree0, nt.states = nt.states))
  }
}

# adds an extra (LT) row information for the LT data
LT.Data <- function (Data){
  Data <- Data[order(Data$id), ]
  temp <- dplyr::distinct(Data,id,.keep_all = T)
  ids <- temp$id
  cIDs <- temp$cID
  csize <- temp$csize
  stop.time <- with(Data, tapply(start, id, min)) # creates stop times for dummy dataset
  enter.st <- by(Data, Data$id, function(x) x$start.stage[which.min(x$start)]) # state at which the subject entered the study
  dummy <- data.frame(cID=cIDs, id = ids, csize=csize, start = -1, stop = stop.time, 
                      start.stage = "LT", end.stage = as.character(enter.st)) 
  Data <- rbind(Data, dummy)
  Data <- with(Data, Data[order(id, stop), ])
  rownames(Data) <- NULL
  return(Data = Data)
}

# computes the counting process and the at-risk set
CP <- function (tree, tree0, Data, nt.states, nt.states2, dNs,Ys, sum_dNs){
  
  # gets all the nodes other than the initial node
  nodes.into <- nodes(tree0)[sapply(inEdges(tree0), function(x) length(x) > 0)]
  
  # counts the number of j to j' transitions for event times > 0
  for (i in nodes.into) {
    nodes.from <- inEdges(tree0)[[i]] 
    for (j in nodes.from) { 
      nam2 <- paste("dN", j, i) 
      idx <- which(Data$end.stage == i & Data$start.stage == j)
      tmp.tab <- table(Data$stop[idx][!Data$stop[idx] == 0]) 
      dNs[names(tmp.tab), nam2] <- tmp.tab 
    }
  }
  
  start.stages <- Data$start.stage[Data$start == 0]
  start.stages <- factor(start.stages, levels = nt.states2, 
                         labels = nt.states2)
  start.cnts <- table(start.stages)

  for (i in nt.states2) {
    n <- start.cnts[i] 
    nam <- paste("y", i) 
    if (length(inEdges(tree0)[[i]]) > 0) {
      into.node <- paste("dN", inEdges(tree0)[[i]], i)
    } else {
      into.node <- NULL
    }
    if (length(edges(tree0)[[i]]) > 0){
      from.node <- paste("dN", i, edges(tree0)[[i]])
    } else {
      from.node <- NULL
    }
    Ys[, nam] <- c(n, n + cumsum(rowSums(dNs[, into.node,drop = FALSE])) 
                   - cumsum(rowSums(dNs[, from.node,drop = FALSE])))[-(nrow(Ys) + 1)]
  }

  a <- strsplit(colnames(sum_dNs), " ")
  a2 <- strsplit(colnames(dNs), " ")
  uni <- unique(sapply(a, function(x) x[2]))
  
  for (i in uni) {
    b <- which(sapply(a, function(x) x[2] == i)) # gets the transition from
    b2 <- which(sapply(a2, function(x) x[2] == i & x[3] != 0)) # gets the transition to (except for the censoring state)
    sum_dNs[, b] <- rowSums(dNs[, b2, drop = FALSE])
  }
  list(dNs = dNs, Ys = Ys, sum_dNs = sum_dNs)
}

################################################################################
# Datta-Satten 2001 (IPCW) weighted quantities - for independent censoring
################################################################################
DS.ind <- function (nt.states, dNs, sum_dNs, Ys) {
  res <- strsplit(colnames(dNs), " ")
  res2 <- strsplit(colnames(Ys), " ")
  res3 <- strsplit(colnames(sum_dNs), " ")

  DS.col.idx <- which(sapply(res, function(x) x[3] == 0)) 
  DS2.col.idx <- which(sapply(res2, function(x) x[2] %in% nt.states)) 
  DS3.col.idx <- which(sapply(res3, function(x) x[2] %in% nt.states)) 
  K <- vector(length = nrow(dNs))

  dN0 <- rowSums(dNs[, DS.col.idx, drop = FALSE]) 
  Y0 <- rowSums(Ys[, DS2.col.idx, drop = FALSE]) 
  N.Y <- ifelse(dN0/Y0 == "NaN", 0, dN0/Y0) 
  colnames(N.Y) <- NULL
  H.t <- cumsum(N.Y) 
  k <- exp(-H.t) 
  K <- c(1, k[-length(k)])

  dNs.K <- dNs/K
  Ys.K <- Ys/K
  sum_dNs.K <- sum_dNs/K
  res <- list(dNs.K = dNs.K, Ys.K = Ys.K, sum_dNs.K = sum_dNs.K)
  return(res)
}

# reduces the CP & at-risk matrices for only timepoints where a transition occurs
Red <- function (tree, dNs, Ys, sum_dNs, dNs.K, Ys.K, sum_dNs.K){
  res <- strsplit(colnames(dNs), " ") 
  col.idx <- which(sapply(res, function(x) x[2] %in% nodes(tree) & x[3] %in% nodes(tree))) 
  row.idx <- which(apply(dNs[, col.idx, drop = FALSE], 1, function(x) any(x > 0))) 
  dNs.et <- dNs[row.idx, col.idx, drop = FALSE]  ## reduces dNs
  
  res2 <- strsplit(colnames(Ys), " ")
  nt.states.f <- names(which(sapply(edges(tree), function(x) length(x) > 0))) 
  col2.idx <- which(sapply(res2, function(x) x[2] %in% nt.states.f)) 
  Ys.et <- Ys[row.idx, col2.idx, drop = FALSE] 
  
  col3.idx <- which(sapply(strsplit(colnames(sum_dNs), " "), 
                           function(x) x[2] %in% nodes(tree))) 
  sum_dNs.et <- sum_dNs[row.idx, col3.idx, drop = FALSE]
  
  dNs.K.et <- dNs.K[row.idx, col.idx, drop = FALSE]
  Ys.K.et <- Ys.K[row.idx, col2.idx, drop = FALSE]
  sum_dNs.K.et <- sum_dNs.K[row.idx, col3.idx, drop = FALSE]
  ans <- list(dNs = dNs.et, Ys = Ys.et, sum_dNs = sum_dNs.et, 
              dNs.K = dNs.K.et, Ys.K = Ys.K.et, sum_dNs.K = sum_dNs.K.et)
  return(ans)
}

# computes the cumulative transition hazard and the stocc
AJ.estimator <- function (ns, tree, dNs.et, Ys.et, start.probs){
  cum.tm <- diag(ns) 
  colnames(cum.tm) <- rownames(cum.tm) <- nodes(tree) 
  ps <- matrix(NA, nrow = nrow(dNs.et), ncol = length(nodes(tree))) 
  rownames(ps) <- rownames(dNs.et)
  colnames(ps) <- paste("p", nodes(tree))

  all.dA <- all.I.dA <- all.AJs <- array(dim = c(ns, ns, nrow(dNs.et)), 
                                         dimnames = list(rows = nodes(tree), cols = nodes(tree), 
                                                         time = rownames(dNs.et)))
  
  for (i in 1:nrow(dNs.et)){
    I.dA <- diag(ns) 
    dA <- matrix(0, nrow = ns, ncol = ns) 
    colnames(I.dA) <- rownames(I.dA) <- colnames(dA) <- rownames(dA) <- nodes(tree)
    
    idx <- which(dNs.et[i, , drop = FALSE] > 0)
    t.nam <- colnames(dNs.et)[idx] 
    tmp <- strsplit(t.nam, " ") 
    start <- sapply(tmp, function(x) x[2])
    end <- sapply(tmp, function(x) x[3]) 
    idxs <- matrix(as.character(c(start, end)), ncol = 2) 
    idxs2 <- matrix(as.character(c(start, start)), ncol = 2) 
    dA[idxs] <- dNs.et[i, idx]/Ys.et[i, paste("y", start)]

    if (length(idx) == 1) {
      dA[start, start] <- -dNs.et[i, idx]/Ys.et[i, paste("y", start)]
    } else {
      dA[idxs2] <- -rowSums(dA[start, ])
    }
    I.dA <- I.dA + dA 
    all.dA[, , i] <- dA 
    all.I.dA[, , i] <- I.dA  
    cum.tm <- cum.tm %*% I.dA  
    all.AJs[, , i] <- cum.tm 
    ps[i, ] <- start.probs %*% all.AJs[, , i]  # the state occupation prob at time i
  }
  list(ps = ps, AJs = all.AJs, I.dA = all.I.dA)
}


##################################################################################
### The wrapper function to obtain the SOP, transition probabilities, 
### starting probabilities, cumulative hazard, and event times
### For clustered multistate data.
##################################################################################

msSurvClust <- function(Data, tree, LT = FALSE, weight=TRUE){
  if(any(!(c("id", "cID", "stop", "start.stage", "end.stage") %in% colnames(Data)))) 
    stop("Incorrect column names for 'Data'.  Column names should have 'id', 'cID', 'stop', 'start.stage', and 'end.stage'.")
  if (!("start" %in% colnames(Data)) & LT == TRUE) 
    stop("The 'start' times must be specified for left truncated data.")
  
  if ("start" %in% colnames(Data)) {
    init.start <- with(Data, tapply(start, id, min))
    if (max(init.start) > 0 & LT == FALSE) {
      warning("It appears left truncation is present.  If so please specify the 'LT=TRUE'\n  argument to 'msSurv'")
    }
  }
  if(min(Data$stop) == 0){
    tmp <- as.character(Data$id[which(Data$stop == 0)])
    stop(paste("For id =", paste(tmp, collapse = ", "), "Transition times between states should be greater than 0", sep = " "))
  }
  
  # creates the start times column
  if (!("start" %in% colnames(Data)) & LT == FALSE) 
    Data <- Add.start(Data)
  if (any(Data$start >= Data$stop)) 
    stop("'start' times must be < 'stop' times.")
  
  ## starting probabilities based on initial states for all individuals on study BEFORE
  ## 1st obs transition time
  idx <- which(Data$start < min(Data$stop[Data$end.stage != 0]))
  start.cnts <- table(factor(Data$start.stage[idx], levels = nodes(tree), labels = nodes(tree)))
  start.probs <- prop.table(start.cnts)
  
  if(sum(start.probs) != 1){
    warning("The sum of the starting probabilities is not equal to 1.")
  }
  
  n <- length(unique(Data$id))  ## total number of individuals in sample
  ns <- length(nodes(tree))  ## number of unique states
  
  # adds the censoring state, i.e. state 0 and LT state
  Cens <- Add.States(tree, LT) 
  
  ## creates empty matrices to store the results 
  ## the dimension of each cluster specific matrix must conform
  ## to the dimension of the full matrix
  if(LT){
    LT_Data <- LT.Data(Data)
    FM <- FMat(tree=tree, tree0=Cens$treeLT, Data=LT_Data, nt.states=Cens$nt.states.LT)
  } else {
    FM <- FMat(tree=tree, tree0=Cens$tree0, Data=Data, nt.states=Cens$nt.states)
  }
  
  M <- length(unique(Data$cID)) # length of unique clusters
  
  # objects for storing the CP and at-risk sets for each cluster
  all_dNs <- vector("list", M)
  names(all_dNs) <- as.character(sort(unique(Data$cID)))
  
  all_dNs.K <- all_Ys.K <- all_sum_dNs.K <- 
    all_Ys <- all_sum_dNs <- all_dNs
  
  for(m in sort(unique(Data$cID))){
    Cdat <- subset(Data, cID == m)
    
    # estimates the CP and at-risk quantities
    if (LT) { # investigate if the next line should be outside the m loop
      Cdat <- LT.Data(Cdat) # transforms data to accommodate LT
      cp <- CP(tree=tree, tree0=Cens$treeLT, Data=Cdat, nt.states=Cens$nt.states.LT, 
               nt.states2=FM$nt.states2, dNs=FM$dNs, Ys=FM$Ys, sum_dNs=FM$sum_dNs)
    } else {
      cp <- CP(tree=tree, tree0=Cens$tree0, Data=Cdat, nt.states=Cens$nt.states, 
               nt.states2=FM$nt.states2, dNs=FM$dNs, Ys=FM$Ys, sum_dNs=FM$sum_dNs)
    }
    
    # Update CP & at-risk quantities via (IPCW) Datta-Satten 2001 estimators 
    # actually, both conditions in the if statement below should give equivalent results
    if(LT){
      ds.est <- DS.ind(nt.states=Cens$nt.states.LT, dNs=cp$dNs, sum_dNs=cp$sum_dNs, Ys=cp$Ys)
    } else {
      ds.est <- DS.ind(nt.states=Cens$nt.states, dNs=cp$dNs, sum_dNs=cp$sum_dNs, Ys=cp$Ys)
    }
    
    # multiply each quantity by the inverse cluster size
    if(weight){
      icw <- 1/Cdat$csize[1]
    } else {
      icw <- 1
    }
    mm <- as.character(m)
    
    # unweighted (IPCW) quantities
    all_dNs[[mm]] <- cp$dNs*icw; 
    all_Ys[[mm]] <- cp$Ys*icw; 
    all_sum_dNs[[mm]] <- cp$sum_dNs*icw
    
    # weighted (IPCW) quantities
    all_dNs.K[[mm]] <- ds.est$dNs.K*icw; 
    all_Ys.K[[mm]] <- ds.est$Ys.K*icw; 
    all_sum_dNs.K[[mm]] <- ds.est$sum_dNs.K*icw
  }
  
  # take the sum for all clusters
  all_dNs <- Reduce("+", all_dNs); 
  all_Ys <- Reduce("+", all_Ys); 
  all_sum_dNs <- Reduce("+", all_sum_dNs)
  
  all_dNs.K <- Reduce("+", all_dNs.K); 
  all_Ys.K <- Reduce("+", all_Ys.K); 
  all_sum_dNs.K <- Reduce("+", all_sum_dNs.K)
  
  # reduce the estimators to event times where there were at least one transition
  cp.red <- Red(tree=tree, dNs=all_dNs, Ys=all_Ys, sum_dNs=all_sum_dNs,
                dNs.K=all_dNs.K, Ys.K=all_Ys.K, sum_dNs.K=all_sum_dNs.K)
  # event times for which a transition occurred
  et <- as.numeric(rownames(cp.red$dNs))
  
  # possible transitions other than transitions to the censoring state
  res.ci2 <- strsplit(colnames(cp.red$dNs), " ")
  a <- sapply(res.ci2, function(x) x[2])
  b <- sapply(res.ci2, function(x) x[3])
  pos.trans <- paste(a, b, sep = " ")
  
  stay <- paste(Cens$nt.states, Cens$nt.states, sep = " ")
  pos.trans <- sort(c(stay, pos.trans))
  # estimates the cumulative transition hazard and the stocc
  AJest <- AJ.estimator(ns=ns, tree=tree, dNs.et=cp.red$dNs.K,
                        Ys.et=cp.red$Ys.K, start.probs=start.probs)
  
  res <- list(tree = tree, ns = ns, et = et, pos.trans = pos.trans, 
              nt.states = Cens$nt.states, dNs = cp.red$dNs, Ys = cp.red$Ys, 
              sum_dNs = cp.red$sum_dNs, dNs.K = cp.red$dNs.K, Ys.K = cp.red$Ys.K, 
              sum_dNs.K = cp.red$sum_dNs.K, ps = AJest$ps, AJs = AJest$AJs, 
              I.dA = AJest$I.dA, start.probs=start.probs)
  return(res)
}

# estimates of occupation probabilities at specific time points
getProbs <- function(fit, timegrid=NULL, cutoffs){
  if(is.null(timegrid)){
    et <- as.numeric(rownames(fit$ps)) # time grid
  } else {
    et <- timegrid
  }
  
  if(length(et) != nrow(fit$ps)){
    stop("Length of the time grid should be equal to the number of rows of the estimated probabilities")
  }
  
  # set SOP at time 0 to be the initial (starting) SOP
  start.probs <- as.numeric(fit$start.probs)
  sops <- rbind(start.probs, fit$ps)
  et <- c(0.0000, et)
  rownames(sops) <- et
  
  # indices for the specified cutoffs
  indx <- sapply(cutoffs, function(x) max(which(et <= x)), simplify=TRUE)
  sops <- sops[indx, ]
  
  return(sops)
}
