#' Function to convert a set of branching times into a
#' phylogeny with random topology
#' This code is taken from the package TESS by Sebastian Hoehna, where the function
#' is called tess.create.phylo
#'
#' Converting a set of branching times to a phylogeny
#'
#'
#' @param times Set of branching times
#' @param root When root is FALSE, the largest branching time will be assumed to be
#' the crown age. When root is TRUE, it will be the stem age.
#' @param tip.label Tip labels. If set to NULL, the labels will be t1, t2, etc.
#' @return \item{ phy }{ A phylogeny of the phylo type }
#' @author Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @export brts2phylo
brts2phylo <- function(times,
                       root = FALSE,
                       tip.label = NULL)
{
  times = sort(times)
  n <- as.integer(length(times)) + 1
  if (root) {
    n <- n - 1
  }
  nbr <- 2 * n - 2
  
  # create the data types for edges and edge-lengths
  edge <- matrix(NA, nbr, 2)
  edge.length <- numeric(nbr)
  
  h <- numeric(2 * n - 1) # initialized with 0's
  pool <- 1:n
  # VERY VERY IMPORTANT: the root MUST have index n+1 !!!
  nextnode <- 2L * n - 1L
  if (n > 1) {
    for (i in 1:(n - 1)) {
      # sample two nodes that have no parent yet
      y <- sample(pool, size = 2)
      # compute the edge indices (we just order the edges from 1 to 2n-2)
      ind <- (i - 1) * 2 + 1:2
      # set the source node of the new edges (i.e. the new internal node)
      edge[ind, 1] <- nextnode
      # set the destination of the new edges (i.e. the two sampled nodes)
      edge[ind, 2] <- y
      # compute the edge length from the difference between the node heights (child <-> parent)
      edge.length[ind] <- times[i] - h[y]
      # store the node height of this new internal node
      # we cannot use times because then we would get into trouble with the indices and we would need to check for tip nodes ...
      h[nextnode] <- times[i]
      # reset the pool of available nodes to merge
      pool <- c(pool[!pool %in% y], nextnode)
      # increase the node index counter
      nextnode <- nextnode - 1L
    }
  }
  
  phy <- list(edge = edge, edge.length = edge.length)
  if (is.null(tip.label))
    tip.label <- paste("t", 1:n, sep = "")
  phy$tip.label <- sample(tip.label)
  phy$Nnode <- n - 1L
  
  if (root) {
    phy$root.edge <- times[n] - times[n - 1]
    phy$root <- times[n] - times[n - 1]
  }
  
  class(phy) <- "phylo"
  
  phy <- stats::reorder(phy)
  ## to avoid crossings when converting with as.hclust:
  phy$edge[phy$edge[, 2] <= n, 2] <- 1:n
  
  return(phy)
}

#' Function to do convolution of two vectors
#'
#' Convolution of two vectors
#'
#' @param x first vector
#' @param y second vector
#' @return vector that is the convolution of x and y
#' @author Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'
#' conv(1:10,1:10)
#'
#' @export conv
conv = function(x, y)
{
  lx = length(x)
  ly = length(y)
  lxy = length(x) + length(y)
  x = c(x, rep(0, lxy - lx))
  y = c(y, rep(0, lxy - ly))
  cvxy = rep(0, lxy)
  for (i in 2:lxy)
  {
    cvxy[i] = crossprod(x[(i - 1):1], y[1:(i - 1)])
  }
  return(cvxy[2:lxy])
}

flavec = function(ddep, la, mu, K, r, lx, kk, n0)
{
  nn = (0:(lx - 1)) + kk
  if (ddep == 1 | ddep == 5)
  {
    lavec = pmax(rep(0, lx), la - 1 / (r + 1) * (la - mu) / K * nn)
  }
  if (ddep == 1.3)
  {
    lavec = pmax(rep(0, lx), la * (1 - nn / K))
  }
  if (ddep == 1.4)
  {
    lavec = pmax(rep(0, lx), la * nn / (nn + K))
  }
  if (ddep == 1.5)
  {
    lavec = pmax(rep(0, lx), la * nn / K * (1 - nn / K))
  }
  if (ddep == 2 | ddep == 2.1 | ddep == 2.2)
  {
    x = -(log(la / mu) / log(K + n0)) ^ (ddep != 2.2)
    lavec = pmax(rep(0, lx), la * (nn + n0) ^ x)
  }
  if (ddep == 2.3)
  {
    x = K
    lavec = pmax(rep(0, lx), la * (nn + n0) ^ x)
  }
  if (ddep == 3 | ddep == 4 | ddep == 4.1 | ddep == 4.2)
  {
    lavec = la * rep(1, lx)
  }
  return(lavec)
}



#' Function to convert a table with speciation and extinction events to a
#' phylogeny
#'
#' Converting a table with speciation and extinction events to a phylogeny
#'
#'
#' @param L Matrix of events as produced by dd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param dropextinct Sets whether the phylogeny should drop species that are
#' extinct at the present
#' @return \item{ phy }{ A phylogeny of the phylo type }
#' @author Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'
#' sim = dd_sim(c(0.2,0.1,20),10)
#' phy = L2phylo(sim$L)
#' plot(phy)
#'
#' @export L2phylo
L2phylo = function(L, dropextinct = T)
  # makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
  L = L[order(abs(L[, 3])), 1:4]
  age = L[1, 1]
  L[, 1] = age - L[, 1]
  L[1, 1] = -1
  notmin1 = which(L[, 4] != -1)
  L[notmin1, 4] = age - L[notmin1, 4]
  if (dropextinct == T)
  {
    sall = which(L[, 4] == -1)
    tend = age
  } else {
    sall = which(L[, 4] >= -1)
    tend = (L[, 4] == -1) * age + (L[, 4] > -1) * L[, 4]
  }
  L = L[, -4]
  linlist = cbind(data.frame(L[sall, ]), paste("t", abs(L[sall, 3]), sep = ""), tend)
  linlist[, 4] = as.character(linlist[, 4])
  names(linlist) = 1:5
  done = 0
  while (done == 0)
  {
    j = which.max(linlist[, 1])
    daughter = linlist[j, 3]
    parent = linlist[j, 2]
    parentj = which(parent == linlist[, 3])
    parentinlist = length(parentj)
    if (parentinlist == 1)
    {
      spec1 = paste(linlist[parentj, 4], ":", linlist[parentj, 5] - linlist[j, 1], sep = "")
      spec2 = paste(linlist[j, 4], ":", linlist[j, 5] - linlist[j, 1], sep = "")
      linlist[parentj, 4] = paste("(", spec1, ",", spec2, ")", sep = "")
      linlist[parentj, 5] = linlist[j, 1]
      linlist = linlist[-j, ]
    } else {
      #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
      linlist[j, 1:3] = L[which(L[, 3] == parent), 1:3]
    }
    if (nrow(linlist) == 1) {
      done = 1
    }
  }
  linlist[4] = paste(linlist[4], ":", linlist[5], ";", sep = "")
  phy = ape::read.tree(text = linlist[1, 4])
  tree = ape::as.phylo(phy)
  return(tree)
}



#' Function to convert phylogeny to a table with speciation and extinction
#' events
#'
#' Converting a phylogeny to a table with speciation and extinction events
#'
#'
#' @param phy A phylogeny of the phylo type
#' @return \item{L}{Matrix of events as produced by dd_sim: \cr \cr - the first
#' column is the time at which a species is born in Mya\cr - the second column
#' is the label of the parent of the species; positive and negative values
#' indicate whether the species belongs to the left or right crown lineage \cr
#' - the third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' }
#' @author Liang Xu
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'
#' sim = dd_sim(c(0.2,0.1,20),10)
#' phy = sim$tas
#' L = phylo2L(phy)
#' phy2 = L2phylo(L, dropextinct = FALSE)
#' graphics::par(mfrow = c(1,3))
#' graphics::plot(phy)
#' graphics::plot(phy2)
#' graphics::plot(L2phylo(sim$L, dropextinct = FALSE))
#'
#' @export phylo2L
phylo2L = function(phy)
{
  emdata = phy
  # compute the relative branching times
  brt = ape::branching.times(emdata)
  if (min(brt) < 0)
  {
    brt = brt + abs(min(brt))
  }
  # number of species including extinct species.
  num.species = emdata$Nnode + 1
  brt_preL = c(brt[emdata$edge[, 1] - length(emdata$tip.label)])
  # check if the relative branching times are equal to the real branching times.
  # if not correct it to the real branching times.
  if (min(brt_preL) == 0)
  {
    correction = max(emdata$edge.length[which(brt_preL == 0)])
    brt_preL = brt_preL + correction
  }
  # preliminary L table
  pre.Ltable = cbind(brt_preL,
                     emdata$edge,
                     emdata$edge.length,
                     brt_preL - emdata$edge.length)
  # identify the extant species and the extinct species
  extantspecies.index = pre.Ltable[which(pre.Ltable[, 5] <= 1e-10), 3]
  tipsindex = c(1:num.species)
  extinct.index3 = subset(tipsindex, !(tipsindex %in% extantspecies.index))
  # assigen the extinct species with extinct times; the extant species with -1 and
  # the internal nodes with 0.
  eeindicator = matrix(0, length(emdata$edge.length), 1)
  eeindicator[match(extantspecies.index, pre.Ltable[, 3])] = -1
  ext.pos = match(extinct.index3, pre.Ltable[, 3])
  eeindicator[ext.pos] = pre.Ltable[ext.pos, 5]
  pre.Ltable = cbind(pre.Ltable, eeindicator)
  
  sort.L = pre.Ltable[order(pre.Ltable[, 1], decreasing = TRUE), ]
  nodesindex = unique(emdata$edge[, 1])
  L = sort.L
  realL = NULL
  do = 0
  while (do == 0) {
    j = which.min(L[, 3])
    daughter = L[j, 3]
    parent = L[j, 2]
    if (parent %in% nodesindex)
    {
      L[which(L[, 2] == parent), 2] = daughter
      if (length(which(L[, 3] == parent)) == 0) {
        realL = rbind(realL, L[j, ], row.names = NULL)
        L = L[-j, , drop = FALSE]
      } else {
        L[which(L[, 3] == parent), 6] = L[j, 6]
        L[which(L[, 3] == parent), 3] = daughter
        L = L[-j, , drop = FALSE]
      }
    } else {
      realL = rbind(realL, L[j, ], row.names = NULL)
      L = L[-j, , drop = FALSE]
    }
    
    if (nrow(L) == 0) {
      do = 1
    }
  }
  realL = realL[order(realL[, 1], decreasing = T), ]
  L = realL[, c(1, 2, 3, 6)]
  
  daughter.index = L[, 3]
  daughter.realindex = c(1:nrow(L))
  parent.index = L[, 2]
  parent.realindex = match(parent.index, daughter.index)
  
  L[, 2] = parent.realindex
  L[, 3] = daughter.realindex
  L[1, 2] = 0
  L[1, 3] = -1
  L[2, 2] = -1
  for (i in c(2:nrow(L)))
  {
    if (L[i - 1, 3] < 0) {
      mrows = which(L[, 2] == abs(L[i - 1, 3]))
      L[mrows, 2] = L[i - 1, 3]
      L[mrows, 3] = -1 * L[mrows, 3]
    }
  }
  dimnames(L) = NULL
  return(L)
}



#' Function to convert a table with speciation and extinction events to a set
#' of branching times
#'
#' Converting a table with speciation and extinction events to a set of
#' branching times
#'
#'
#' @param L Matrix of events as produced by dd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param dropextinct Sets whether the phylogeny should drop species that are
#' extinct at the present
#' @return \item{ brts }{ A set of branching times }
#' @author Rampal S. Etienne
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @examples
#'
#' sim = dd_sim(c(0.2,0.1,20),10)
#' phy = L2brts(sim$L)
#' plot(phy)
#'
#' @export L2brts
L2brts = function(L, dropextinct = T)
  # makes a phylogeny out of a matrix with branching times, parent and daughter species, and extinction times
{
  brts = NULL
  L = L[order(abs(L[, 3])), 1:4]
  age = L[1, 1]
  L[, 1] = age - L[, 1]
  L[1, 1] = -1
  notmin1 = which(L[, 4] != -1)
  L[notmin1, 4] = age - L[notmin1, 4]
  if (dropextinct == T)
  {
    sall = which(L[, 4] == -1)
    tend = age
  } else {
    sall = which(L[, 4] >= -1)
    tend = (L[, 4] == -1) * age + (L[, 4] > -1) * L[, 4]
  }
  L = L[, -4]
  linlist = cbind(data.frame(L[sall, ]), paste("t", abs(L[sall, 3]), sep = ""), tend)
  linlist[, 4] = as.character(linlist[, 4])
  names(linlist) = 1:5
  done = 0
  while (done == 0)
  {
    j = which.max(linlist[, 1])
    daughter = linlist[j, 3]
    parent = linlist[j, 2]
    parentj = which(parent == linlist[, 3])
    parentinlist = length(parentj)
    if (parentinlist == 1)
    {
      spec1 = paste(linlist[parentj, 4], ":", linlist[parentj, 5] - linlist[j, 1], sep = "")
      spec2 = paste(linlist[j, 4], ":", linlist[j, 5] - linlist[j, 1], sep = "")
      linlist[parentj, 4] = paste("(", spec1, ",", spec2, ")", sep = "")
      linlist[parentj, 5] = linlist[j, 1]
      brts = c(brts, linlist[j, 1])
      linlist = linlist[-j, ]
    } else {
      #linlist[j,1:3] = L[abs(as.numeric(parent)),1:3]
      linlist[j, 1:3] = L[which(L[, 3] == parent), 1:3]
    }
    if (nrow(linlist) == 1) {
      done = 1
    }
  }
  #linlist[4] = paste(linlist[4],":",linlist[5],";",sep = "")
  brts = rev(sort(age - brts))
  return(brts)
}


#' Rounds up in the usual manner
#'
#' The standard round function in R rounds x.5 to the nearest even integer.
#' This is odd behavior that is corrected in roundn
#'
#'
#' @param x Number to be rounded
#' @param digits Sets the number of decimals in rounding.
#' @return \item{n}{ A number }
#' @author Rampal S. Etienne
#' @keywords models
#' @examples
#'
#' round(2.5)
#' roundn(2.5)
#' round(3.5)
#' roundn(3.5)
#' round(2.65,digits = 1)
#' roundn(2.65,digits = 1)
#' round(2.75,digits = 1)
#' roundn(2.75,digits = 1)
#'
#' @export roundn
roundn = function(x, digits = 0)
{
  fac = 10 ^ digits
  n = trunc(fac * x + 0.5) / fac
  return(n)
}



#' Takes samples in the usual manner
#'
#' The standard sample function in R samples from n numbers when x = n. This is
#' unwanted behavior when the size of the vector to sample from changes
#' dynamically. This is corrected in sample2
#'
#'
#' @param x A vector of one or more elements
#' @param size A non-negative integer giving the number of items to choose.
#' @param replace Should sampling be with replacement?
#' @param prob A vector of probability weights for obtaining the elements of
#' the vector being sampled.
#' @return \item{sam}{A vector of length \code{size} that is sampled from
#' \code{x}. }
#' @author Rampal S. Etienne
#' @keywords models
#' @examples
#'
#' sample(x = 10,size = 5,replace = TRUE)
#' sample2(x = 10,size = 5,replace = TRUE)
#'
#' @export sample2
sample2 = function(x,
                   size,
                   replace = FALSE,
                   prob = NULL)
{
  if (length(x) == 1)
  {
    x = c(x, x)
    prob = c(prob, prob)
    if (is.null(size))
    {
      size = 1
    }
    if (replace == FALSE & size > 1)
    {
      stop('It is not possible to sample without replacement multiple times from a single item.')
    }
  }
  sam = sample(x, size, replace, prob)
  return(sam)
}

#' Carries out optimization using a simplex algorithm (finding a minimum)
#'
#' Function to optimize target function using a
#' simplex method adopted from Matlab
#'
#' @param fun Function to be optimized
#' @param trparsopt Initial guess of the parameters to be optimized
#' @param ... Any other arguments of the function to be optimimzed, or settings
#' of the optimization routine
#' @param optimpars Parameters of the optimization: relative tolerance in
#' function arguments, relative tolerance in function value, absolute tolerance
#' in function arguments, and maximum number of iterations
#' @return \item{out}{ A list containing optimal function arguments
#' (\code{par}, the optimal function value (\code{fvalues}) and whether the
#' optimization converged (\code{conv})}.
#' @author Rampal S. Etienne
#' @keywords models
#' @examples
#'
#' cat("No examples")
#'
#' @export simplex
simplex = function(fun, trparsopt, optimpars, ...)
{
  numpar = length(trparsopt)
  reltolx = optimpars[1]
  reltolf = optimpars[2]
  abstolx = optimpars[3]
  maxiter = optimpars[4]
  
  ## Setting up initial simplex
  v = t(matrix(rep(trparsopt, each = numpar + 1), nrow = numpar + 1))
  for (i in 1:numpar)
  {
    parsoptff = 1.05 * untransform_pars(trparsopt[i])
    trparsoptff = transform_pars(parsoptff)
    fac = trparsoptff / trparsopt[i]
    if (v[i, i + 1] == 0)
    {
      v[i, i + 1] = 0.00025
    } else {
      v[i, i + 1] = v[i, i + 1] * min(1.05, fac)
    }
  }
  
  fv = rep(0, numpar + 1)
  for (i in 1:(numpar + 1))
  {
    fv[i] = -fun(trparsopt = v[, i], ...)
  }
  
  how = "initial"
  itercount = 1
  string = itercount
  for (i in 1:numpar)
  {
    string = paste(string, untransform_pars(v[i, 1]), sep = " ")
  }
  string = paste(string, -fv[1], how, "\n", sep = " ")
  cat(string)
  utils::flush.console()
  
  tmp = order(fv)
  if (numpar == 1)
  {
    v = matrix(v[tmp], nrow = 1, ncol = 2)
  } else {
    v = v[, tmp]
  }
  fv = fv[tmp]
  
  ## Iterate until stopping criterion is reached
  rh = 1
  ch = 2
  ps = 0.5
  si = 0.5
  
  v2 = t(matrix(rep(v[, 1], each = numpar + 1), nrow = numpar + 1))
  
  while (itercount <= maxiter &
         ((is.nan(max(abs(
           fv - fv[1]
         ))) |
         (max(abs(
           fv - fv[1]
         )) - reltolf * abs(fv[1]) > 0)) + ((max(
           abs(v - v2) - reltolx * abs(v2)
         ) > 0) | (max(abs(
           v - v2
         )) - abstolx > 0))))
  {
    ## Calculate reflection point
    
    if (numpar == 1)
    {
      xbar = v[1]
    } else {
      xbar = rowSums(v[, 1:numpar]) / numpar
    }
    xr = (1 + rh) * xbar - rh * v[, numpar + 1]
    fxr = -fun(trparsopt = xr, ...)
    
    if (fxr < fv[1])
    {
      ## Calculate expansion point
      xe = (1 + rh * ch) * xbar - rh * ch * v[, numpar + 1]
      fxe = -fun(trparsopt = xe, ...)
      if (fxe < fxr)
      {
        v[, numpar + 1] = xe
        fv[numpar + 1] = fxe
        how = "expand"
      } else {
        v[, numpar + 1] = xr
        fv[numpar + 1] = fxr
        how = "reflect"
      }
    } else {
      if (fxr < fv[numpar])
      {
        v[, numpar + 1] = xr
        fv[numpar + 1] = fxr
        how = "reflect"
      } else {
        if (fxr < fv[numpar + 1])
        {
          ## Calculate outside contraction point
          xco = (1 + ps * rh) * xbar - ps * rh * v[, numpar + 1]
          fxco = -fun(trparsopt = xco, ...)
          if (fxco <= fxr)
          {
            v[, numpar + 1] = xco
            fv[numpar + 1] = fxco
            how = "contract outside"
          } else {
            how = "shrink"
          }
        } else {
          ## Calculate inside contraction point
          xci = (1 - ps) * xbar + ps * v[, numpar + 1]
          fxci = -fun(trparsopt = xci, ...)
          if (fxci < fv[numpar + 1])
          {
            v[, numpar + 1] = xci
            fv[numpar + 1] = fxci
            how = "contract inside"
          } else {
            how = "shrink"
          }
        }
        if (how == "shrink")
        {
          for (j in 2:(numpar + 1))
          {
            v[, j] = v[, 1] + si * (v[, j] - v[, 1])
            fv[j] = -fun(trparsopt = v[, j], ...)
          }
        }
      }
    }
    tmp = order(fv)
    if (numpar == 1)
    {
      v = matrix(v[tmp], nrow = 1, ncol = 2)
    } else {
      v = v[, tmp]
    }
    fv = fv[tmp]
    itercount = itercount + 1
    string = itercount
    
    for (i in 1:numpar)
    {
      string = paste(string, untransform_pars(v[i, 1]), sep = " ")
    }
    string = paste(string, -fv[1], how, "\n", sep = " ")
    cat(string)
    utils::flush.console()
    v2 = t(matrix(rep(v[, 1], each = numpar + 1), nrow = numpar + 1))
  }
  if (itercount < maxiter)
  {
    cat("Optimization has terminated successfully.", "\n")
  } else {
    cat("Maximum number of iterations has been exceeded.", "\n")
  }
  out = list(
    par = v[, 1],
    fvalues = -fv[1],
    conv = as.numeric(itercount > maxiter)
  )
  invisible(out)
}

#' Carries out optimization (finding a minimum)
#'
#' A wrapper to use several optimization routines, currently only 'simplex' (a
#' method adopted from Matlab, or 'subplex', from the R package subplex). The
#' function is called from several packages by the same author.
#'
#'
#' @param optimmethod The method to use for optimization, either 'simplex' or
#' 'subplex'
#' @param optimpars Parameters of the optimization: relative tolerance in
#' function arguments, relative tolerance in function value, absolute tolerance
#' in function arguments, and maximum number of iterations
#' @param num_cycles Number of cycles of the optimization. When set to Inf, the
#' optimization will be repeated until the result is, within the tolerance,
#' equal to the starting values, with a maximum of 10 cycles.
#' @param fun Function to be optimized
#' @param trparsopt Initial guess of the parameters to be optimized
#' @param ... Any other arguments of the function to be optimimzed, or settings
#' of the optimization routine
#' @return \item{out}{ A list containing optimal function arguments
#' (\code{par}, the optimal function value (\code{fvalues}) and whether the
#' optimization converged (\code{conv})}.
#' @author Rampal S. Etienne
#' @keywords models
#' @examples
#'
#' cat("No examples")
#'
#' @export optimizer
optimizer = function(optimmethod = 'simplex',
                     optimpars = c(1E-4, 1E-4, 1E-6, 1000),
                     num_cycles = 1,
                     fun,
                     trparsopt,
                     ...)
{
  if (num_cycles == Inf)
  {
    max_cycles <- 10
  } else
  {
    max_cycles <- num_cycles
  }
  cy <- 1
  fvalue <- rep(-Inf, max_cycles)
  out <- NULL
  while (cy <= max_cycles)
  {
    if (num_cycles > 1)
      cat(paste('Cycle ', cy, '\n', sep = ''))
    if (optimmethod == 'simplex')
    {
      outnew = suppressWarnings(simplex(
        fun = fun,
        trparsopt = trparsopt,
        optimpars = optimpars,
        ...
      ))
    }
    if (optimmethod == 'subplex')
    {
      minfun = function(fun, trparsopt, ...)
      {
        return(-fun(trparsopt = trparsopt, ...))
      }
      outnew = suppressWarnings(subplex::subplex(
        par = trparsopt,
        fn = minfun,
        control = list(
          abstol = optimpars[3],
          reltol = optimpars[1],
          maxit = optimpars[4]
        ),
        fun = fun,
        ...
      ))
      outnew = list(
        par = outnew$par,
        fvalues = -outnew$value,
        conv = outnew$convergence
      )
    }
    if (cy > 1 &
        (
          any(is.na(outnew$par)) |
          any(is.nan(outnew$par)) |
          is.na(outnew$fvalues) | is.nan(outnew$fvalues) | outnew$conv != 0
        ))
    {
      cat('The last cycle failed; second last cycle result is returned.\n')
      return(out)
    } else
    {
      out <- outnew
      trparsopt <- out$par
      fvalue[cy] <- out$fvalues
    }
    if (cy > 1)
    {
      if (abs(fvalue[cy] - fvalue[cy - 1]) < optimpars[3])
      {
        if (cy < max_cycles)
          cat('No more cycles needed.\n')
        cy <- max_cycles
      } else if (cy == max_cycles)
      {
        cat('More cycles in optimization recommended.\n')
      }
    }
    cy <- cy + 1
  }
  return(out)
}

#' @name transform_pars
#' @title Transforming parameters from -Inf to Inf into parameters
#' from -1 to 1
#' @description Function to transform pars in a way that is more
#' useful for optimization: trpars <- sign(pars) * pars/(sign(pars) + pars);
#' @param pars Parameters to be transformed
#' @return Transformed parameters
#' @author Rampal S. Etienne
#' @export transform_pars
transform_pars <- function(pars)
{
  trpars1 <- sign(pars) * pars / (sign(pars) + pars)
  
  trpars1[which(pars == 0)] <- 0
  
  trpars1[which(pars == -Inf)] <- -1
  
  trpars1[which(pars == Inf)] <- 1
  
  return(trpars1)
  
}

#' @name untransform_pars
#' @title Untransforming parameters from -1 to 1 into parameters
#' from -Inf to Inf.
#' @description Function to untransform pars after optimization:
#' pars <- sign(trpars) * trpars/(sign(trpars) - trpars);
#' @param trpars Parameters to be untransformed
#' @return Untransformed parameters
#' @author Rampal S. Etienne
#' @export untransform_pars
untransform_pars <- function(trpars)
{
  pars <- sign(trpars) * trpars / (sign(trpars) - trpars)
  
  pars[which(trpars == 0)] <- 0
  
  pars[which(trpars == 1)] <- Inf
  
  pars[which(trpars == -1)] <- -Inf
  
  return(pars)
}

check_probs <- function(loglik, probs, verbose)
{
  probs <- probs * (probs > 0)
  if (is.na(sum(probs)) | is.nan(sum(probs)))
  {
    if (verbose)
      cat('NA or NaN issues encountered.\n')
    loglik <- -Inf
    probs <- rep(-Inf, length(probs))
  } else if (sum(probs) <= 0)
  {
    if (verbose)
      cat('Probabilities smaller than 0 encountered\n')
    loglik <- -Inf
    probs <- rep(-Inf, length(probs))
  } else {
    loglik <- loglik + log(sum(probs))
    probs <- probs / sum(probs)
  }
  return(list(loglik, probs))
}

#' Sampling in which zero probabilities are removed
#' @inheritParams base::sample
#' @return a vector of length size with elements drawn
#'   from either x or from the integers 1:x.
#' @examples
#'   # Number of draws
#'   n <- 1000
#'
#'   # Do normal sampling
#'   set.seed(42)
#'   draws_1 <- DDD:::rng_respecting_sample(
#'     1:3, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0)
#'   )
#'
#'   # Do a sampling with one element of probability zero
#'   set.seed(42)
#'   draws_2 <- DDD:::rng_respecting_sample(
#'     1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
#'   )
#'   testit::assert(sum(draws_2 == 4) == 0)
#'   testit::assert(draws_1 == draws_2)
#'
#'   # Use base sampling will give different results,
#'   # as it results in different RNG values
#'   set.seed(42)
#'   draws_3 <- sample(
#'     1:4, size = n, replace = TRUE, prob = c(1.0, 1.0, 1.0, 0.0)
#'   )
#'   testit::assert(sum(draws_3 == 4) == 0)
#'   testit::assert(!all(draws_1 == draws_3))
#'
#' @author Richel J.C. Bilderbeek
#' @seealso See \code{\link[base]{sample}} for more details
#' @note thanks to Pedro Neves for finding this feature in base::sample
#' @export
rng_respecting_sample <- function(x, size, replace, prob) {
  which_non_zero <- prob > 0.0
  non_zero_prob <- prob[which_non_zero]
  non_zero_x <- x[which_non_zero]
  return(sample(
    x = non_zero_x,
    size = size,
    replace = replace,
    prob = non_zero_prob
  ))
}

#' @name L2Phi
#' @title Converting a table with speciation and extinction events to a phylogenetic
#' diversity metric
#' @description Function to convert a table with speciation and extinction events to a
#' phylogenetic diversity metric
#' @param L Matrix of events as produced by pdd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @param t Sets whether the phylogeny should drop species that are
#' extinct at the present
#' @param metric specifies which phylogenetic diversity metric should be used
#' @return a value of one of the phylogenetic diversity metrices
#' @author Tianjian Qin
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @export L2Phi
L2Phi <- function(L, t, metric) {
  # reverse time scale
  L[, 1] <- t - c(L[, 1])
  notmin1 <- which(L[, 4] != -1)
  L[notmin1, 4] <- t - c(L[notmin1, 4])
  L[which(L[, 4] == t + 1), 4] <- -1
  
  # metrics
  if (metric == "pd") {
    return(sum(DDD::L2phylo(L, dropextinct = T)$edge.length))
  } else if (metric == "mpd") {
    phy <- DDD::L2phylo(L, dropextinct = T)
    n <- length(phy$tip.label)
    dist <- dist.nodes(phy)[1:n, 1:n]
    return(mean(dist[lower.tri(dist)]))
  }
}

#' @name L2ED
#' @title Converting a table with speciation and extinction events to evolutionary 
#' distances
#' @description Function to convert a table with speciation and extinction events to 
#' mean evolutionary distances between each species and the rest of the community
#' @param L Matrix of events as produced by pdd_sim: \cr \cr - the first column
#' is the time at which a species is born in Mya\cr - the second column is the
#' label of the parent of the species; positive and negative values indicate
#' whether the species belongs to the left or right crown lineage \cr - the
#' third column is the label of the daughter species itself; positive and
#' negative values indicate whether the species belongs to the left or right
#' crown lineage \cr - the fourth column is the time of extinction of the
#' species; if the fourth element equals -1, then the species is still extant.
#' @return a named vector of evolutionary distances
#' @author Tianjian Qin
#' @references - Etienne, R.S. et al. 2012, Proc. Roy. Soc. B 279: 1300-1309,
#' doi: 10.1098/rspb.2011.1439 \cr - Etienne, R.S. & B. Haegeman 2012. Am. Nat.
#' 180: E75-E89, doi: 10.1086/667574
#' @keywords models
#' @export L2ED
L2ED <- function(L, t) {
  require(gtools)
  require(ape)

  # reverse time scale
  L[, 1] <- t - c(L[, 1])
  notmin1 <- which(L[, 4] != -1)
  L[notmin1, 4] <- t - c(L[notmin1, 4])
  L[which(L[, 4] == t + 1), 4] <- -1
  
  dist.tips <- cophenetic.phylo(L2phylo(L, dropextinct = T))
  diag(dist.tips) <- NA
  dist.means <- rowMeans(dist.tips, na.rm = T)
  dist.means.sorted <- dist.means[mixedorder(names(dist.means))]
  return(range01(dist.means.sorted))
}

#' @name range01
#' @title standardizing data to the range from 0 to 1
#' @description Function to standardize data to the range from 0 to 1 while retaining 
#' the relative size
#' @param x Data to be standardized
#' @author Tianjian Qin
#' @keywords standardize
#' @export range01
range01 <- function(x) {
  if (max(x) < 0) {
    stop("ED < 0")
  } else if (max(x) == 0) {
    sapply(x, function(x) {
      1
    })
  } else{
    maxval <- max(x)
    sapply(x, function(x) {
      x / maxval
    })
  }
}

#' @name phylo2mpd
#' @title Converting a phylogenetic tree of phylo object to mean pairwise distance
#' among all tips
#' @description Function to convert a phylogenetic tree with extant species to mean pairwise
#' distance
#' @param phy a phylo object that only contains extant species
#' extinct at the present
#' @return a numeric value of mean pairwise distance
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export phylo2mpd
phylo2mpd <- function(phy) {
  n <- length(phy$tip.label)
  dist <- dist.nodes(phy)[1:n, 1:n]
  return(mean(dist[lower.tri(dist)]))
}

#' @name phylo2pd
#' @title Converting a phylogenetic tree of phylo object to total lengths of all
#' branches in the tree
#' @description Function to convert a phylogenetic tree with extant species to
#' phylogenetic diversity
#' @param phy
#' @param a phylo object that only contains extant species
#' extinct at the present
#' @return a numeric value of phylogenetic diversity
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export phylo2pd
phylo2pd <- function(phy) {
  return(sum(phy$edge.length))
}

#' @name extract_pdd_result
#' @title Extracting all results generated in simulation
#' @description Function to extract all results generated in pdd simulation
#' @param result a list generated by pdd simulation function
#' @param nrep the number of replication used in the simulation
#' @return a list of extracted results
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export extract_pdd_result
extract_pdd_result <- function(result, nrep, which, nlist = 5) {
  out <-
    return(result[seq(which, nlist * nrep, by = nlist)])
}

#' @name bind_result
#' @title Binding extracted results in simulation
#' @description Function to bind all results extracted from pdd simulation
#' @param result a list generated by pdd simulation function
#' @return a binded list of extracted results
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export bind_result
bind_result <- function(result) {
  out <- lapply(result, as.data.frame)
  n <- length(out)
  
  for (i in 1:n) {
    out[[i]]["rep"] <- rep(i, times = nrow(out[[i]]))
  }
  
  return(dplyr::bind_rows(out))
}

#' @name pdd_simulation_replicated
#' @title Running a replicated pdd simulation
#' @description Function to automatically run a replicated pdd simulation
#' @param rep the number of replication
#' @param pars parameters used in a pdd simulation
#' @param age the total running time of a pdd simulation
#' @param model the type of model used in a pdd simulation
#' @param metric the metric used in a pdd simulation
#' @param offset the offset method for Phi used in a pdd simulation
#' @return a list of simulation results
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export pdd_simulation_replicated
pdd_simulation_replicated <-
  function(rep, pars, age, model, metric, offset) {
    require(tidyverse)
    results <- replicate(rep, {
      pdd_sim(pars,
              age,
              model,
              metric,
              offset)
    })
    result_pars <- extract_pdd_result(results, 10, 5)
    result_pars_binded <- bind_result(result_pars)
    
    return(list(result_pars = result_pars, result_pars_binded = result_pars_binded))
  }

#' @name pdd_simulation_plot
#' @title Plotting results of a replicated pdd simulation
#' @description Function to automatically produce several plots from the results
#' of a replicated pdd simulation
#' @param result a list of extracted results
#' @param pars the parameters used in simulation being plotted
#' @param age the total simulation time
#' @param model the model used in the simulation
#' @return a plot grid generated by cow_plot
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export pdd_simulation_plot
pdd_simulation_plot <- function(result, pars, age, model, offset) {
  require(tidyverse)
  require(cowplot)
  
  plotN <-
    ggplot(result[[2]], aes(time, N, group = rep, color = rep)) + geom_line() +
    theme(legend.position = "none") + xlab("age")
  
  if (model == "dsce1" | model == "dsde1") {
    plotPhi <-
      ggplot(result[[2]], aes(time, Phi, group = rep, color = rep)) + geom_line() +
      theme(legend.position = "none") + xlab("age") + geom_hline(yintercept = pars[3])
    title <-
      ggdraw() + draw_label(
        paste(
          model,
          ", ",
          "lambda = ",
          pars[1],
          ", ",
          "mu = ",
          pars[2],
          ", ",
          "K = ",
          pars[3],
          ", ",
          "age = ",
          age,
          ", ",
          "offset = ",
          offset
        ),
        fontface = 'bold'
      )
  } else if (model == "dsde2") {
    plotPhi <-
      ggplot(result[[2]], aes(time, Phi, group = rep, color = rep)) + geom_line() +
      theme(legend.position = "none") + xlab("age")
    title <-
      ggdraw() + draw_label(
        paste(
          model,
          ", ",
          "lambda = ",
          pars[1],
          ", ",
          "mu = ",
          pars[2],
          ", ",
          "betaN = ",
          pars[3],
          ", ",
          "betaPhi = ",
          pars[4],
          ", ",
          "gammaN = ",
          pars[5],
          ", ",
          "gammaPhi = ",
          pars[6],
          ", ",
          "age = ",
          age,
          ", ",
          "offset = ",
          offset
        ),
        fontface = 'bold'
      )
  } else if (model == "dsce2"){
    plotPhi <-
      ggplot(result[[2]], aes(time, Phi, group = rep, color = rep)) + geom_line() +
      theme(legend.position = "none") + xlab("age")
    title <-
      ggdraw() + draw_label(
        paste(
          model,
          ", ",
          "lambda = ",
          pars[1],
          ", ",
          "mu = ",
          pars[2],
          ", ",
          "betaN = ",
          pars[3],
          ", ",
          "betaPhi = ",
          pars[4],
          ", ",
          "age = ",
          age,
          ", ",
          "offset = ",
          offset
        ),
        fontface = 'bold'
      )
  }
  
  plotla <-
    ggplot(result[[2]], aes(time, lambda, group = rep, color = rep)) + geom_line() +
    theme(legend.position = "none") + xlab("age")
  plotmu <-
    ggplot(result[[2]], aes(time, mu, group = rep, color = rep)) + geom_line() +
    theme(legend.position = "none") + xlab("age")
  plotall <- cowplot::plot_grid(plotN, plotPhi, plotla, plotmu)
  
  plotwithtitle <-
    plot_grid(title,
              plotall,
              ncol = 1,
              rel_heights = c(0.1, 1))
  
  return(plotwithtitle)
}

#' @name pdd_simulation_wrapper
#' @title Conducting a simulation and plotting the results
#' @description Function to automatically run a replicated pdd simulation and then
#' produce several plots from the results
#' @param rep the number of replication
#' @param pars parameters used in a pdd simulation
#' @param age the total running time of a pdd simulation
#' @param model the type of model used in a pdd simulation
#' @param metric the metric used in a pdd simulation
#' @param offset the offset method for Phi used in a pdd simulation
#' @author Tianjian Qin
#' @keywords phylogenetics
#' @export pdd_simulation_wrapper
pdd_simulation_wrapper <- function(rep, pars, age, model, metric, offset) {
  result <-
    pdd_simulation_replicated(
      rep = rep,
      pars = pars,
      age = age,
      model = model,
      metric = metric,
      offset = offset
    )
  
  if (model == "dsce1" | model == "dsde1") {
    if (is.null(pars[3])) {
      stop("missing K")
    }
  }
  pdd_simulation_plot(result, pars, age, model, offset)
}

edd_simulation_replicated <-
  function(rep, pars, age, model, metric, offset) {
    require(tidyverse)
    results <- replicate(rep, {
      edd_sim(pars,
              age,
              model,
              metric,
              offset)
    })
    result_pars <- extract_pdd_result(results, 10, 5)
    result_pars_binded <- bind_result(result_pars)
    
    return(list(result_pars = result_pars, result_pars_binded = result_pars_binded))
  }

edd_simulation_wrapper <- function(rep, pars, age, model, metric, offset) {
  result <-
    edd_simulation_replicated(
      rep = rep,
      pars = pars,
      age = age,
      model = model,
      metric = metric,
      offset = offset
    )
  
  if (model == "dsce1" | model == "dsde1") {
    if (is.null(pars[3])) {
      stop("missing K")
    }
  }
  pdd_simulation_plot(result, pars, age, model, offset)
}
