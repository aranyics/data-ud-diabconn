balance = function(X){
  X = X-t(X)
  X[lower.tri(X)] = 0
  #return( (X-t(X))[upper.tri(X)] )
  return( X )
}

calc.stat = function (graph, types, mode, sumgraph = NULL, weighted = TRUE)
{
  graph <- as.matrix(graph)
  if (!weighted) {
    graph[graph > 0] <- 1
    graph[graph < 1] <- 0
  }
  
  if (is.null(sumgraph)){
    sumgraph = sum(graph)
  }
  
  if (mode == 'assortativity')
  {
    total_types <- unique(types)
    out <- matrix(0, nrow = length(total_types), ncol = length(total_types))
    for (i in 1:length(total_types)) {
      for (j in 1:length(total_types)) {
        out[i, j] <- sum(graph[which(types == total_types[i]), 
                               which(types == total_types[j])])/sumgraph#sum(graph)
      }
    }
    
    r <- (sum(diag(out)) - sum(rowSums(out) * colSums(out)))/(1 - sum(rowSums(out) * colSums(out)))
    
    nwout <- matrix(0, nrow = 1, ncol = length(total_types))
    for (i in 1:length(total_types)){
      nwout[i] <- (sum(c(graph[which(types == total_types[i]),], 
                      graph[,which(types == total_types[i])])) - 
                   sum(graph[which(types == total_types[i]), 
                        which(types == total_types[i])]) )   /  sumgraph#sum(graph)
    }
    
    
    out <- cbind(out, rowSums(out))
    out <- rbind(out, colSums(out))
    colnames(out) <- c(as.character(total_types), "ai")
    rownames(out) <- c(as.character(total_types), "bi")
    
    return(list(r = r, mixing_matrix = out, nwout = c(nwout)))
  }
  
  if (mode == 'strength')
  {
    total_types <- unique(types)
    out <- matrix(0, nrow = length(total_types), ncol = length(total_types))
    for (i in 1:length(total_types)) {
      for (j in 1:length(total_types)) {
        out[i, j] <- sum(graph[which(types == total_types[i]), 
                               which(types == total_types[j])])
      }
    }
    
    nwout <- matrix(0, nrow = 1, ncol = length(total_types))
    for (i in 1:length(total_types)){
      nwout[i] <- (sum(c(graph[which(types == total_types[i]),], 
                        graph[,which(types == total_types[i])])) - 
                   sum(graph[which(types == total_types[i]), 
                        which(types == total_types[i])]) )
    }
    
    out <- cbind(out, rowSums(out))
    out <- rbind(out, colSums(out))
    colnames(out) <- c(as.character(total_types), "ai")
    rownames(out) <- c(as.character(total_types), "bi")
    
    return(list( mixing_matrix = out, nwout = c(nwout)))
  }
  
  if (mode == 'balance')
  {
    graph = balance(graph)
    total_types <- unique(types)
    out <- matrix(0, nrow = length(total_types), ncol = length(total_types))
    outp <- matrix(0, nrow = length(total_types), ncol = length(total_types))
    for (i in 1:length(total_types)) {
      for (j in i:length(total_types)) {
        out[i, j] <- sum(graph[which(types == total_types[i]), 
                               which(types == total_types[j])])/sumgraph#sum(graph)
        outp[i, j] <- t.test(graph[which(types == total_types[i]), 
                                   which(types == total_types[j])])$p.value
      }
    }
    outp.full = outp
    outp[which(outp > 0.001)] = 0
    
    r <- mean( graph[upper.tri(graph)] )
    rp <- t.test( graph[upper.tri(graph)] )$p.value
    
    nwout <- matrix(0, nrow = 1, ncol = length(total_types))
    nwoutp <- matrix(0, nrow = 1, ncol = length(total_types))
    for (i in 1:length(total_types)){
      nwout[i] <- sum(c(graph[which(types == total_types[i]),], 
                        graph[,which(types == total_types[i])]))/sumgraph#sum(graph)
      nwoutp[i] <- t.test(c(graph[which(types == total_types[i]),], 
                        graph[,which(types == total_types[i])]))$p.value
    }
    #nwoutp[which(nwoutp > 0.001)] = 0
    
    regout <- matrix(0, nrow = 1, ncol = dim(graph)[1])
    regoutp <- matrix(0, nrow = 1, ncol = dim(graph)[1])
    for (i in 1:dim(graph)[1]){
      regout[i] <- sum(c(graph[i,], graph[,i])) / sumgraph#sum(graph)
      regoutp[i] <- t.test( c(graph[i,], graph[,i]) )$p.value
    }
    
    return(list(r = r, rp = rp, out = out, outp = outp, outp.full = outp.full, nwout = c(nwout), nwoutp = c(nwoutp), regout = regout, regoutp = regoutp))
  }
  
}
