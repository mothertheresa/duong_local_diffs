verify_inputs <- function(dat, dim, m, alpha, plot){
  # throw an error if dat is null
  
  # throw an error if dim is not an integer
  
  # throw an error if m is not an integer
  
  # throw an error if alpha is not numeric
  
  # throw an error if plot is not a boolean
  
  return("done verifying inputs")
}

get_m_from_d <- function(d){
  # Wand 1994 recommends (d, m) = (1, 401), (2, 151^2), (3, 51^3)
  # and a fixed uniform grid of points
  if (d == 1){
    return(401)
  } else if (d == 2){
    return(151)
  } else if (d == 3) {
    return(51)
  } else {
    stop("invalid input for d: not 1, 2, or 3")
  }
}

density.est <- function(dat, m, from, to, give.Rkern){
  # bandwidth from Sheather & Jones (1991)
  # using pilot estimation of derivatives
  
  # SJ is better than nrd0 according to Venables and Ripley (2002)
  return(density(dat, kernel="gaussian", n = m, bw="SJ",
                 from = from, to = to, give.Rkern = give.Rkern))
}

f.hat <- function(dat, m, from, to){
  d.est = density.est(dat, m, from, to, give.Rkern = FALSE)
  
  function(x){
    return(d.est$y[d.est$x == x])}
}

U.stat <- function(x, f1.hat, f2.hat){
  return((f1.hat(x) - f2.hat(x))^2)
}

sigma.sq.hat <- function(x, group1, group2, f1.hat, f2.hat, m, from, to){
  d.est1 = density.est(group1, m, from, to, give.Rkern = FALSE)
  d.est2 = density.est(group2, m, from, to, give.Rkern = FALSE)
  
  h1 = d.est1$bw
  h2 = d.est2$bw
  
  n.val1 = d.est1$n
  n.val2 = d.est2$n
  
  # R(K) same for both groups, so we can just from the first one
  Rkern = density.est(group1, m, from, to, give.Rkern = TRUE)
  
  part1 = f1.hat(x) / (n.val1 * sqrt(h1))
  part2 = f2.hat(x) / (n.val2 * sqrt(h2))
  return(Rkern * (part1 + part2))
}

t.stat <- function(x, group1, group2, f1.hat, f2.hat, m, from, to){
  sigma.sq = sigma.sq.hat(x, group1, group2, f1.hat, f2.hat, m, from, to)
  U = U.stat(x, f1.hat, f2.hat)
  return(U^2/sigma.sq) #should U be squared?
}

p.val <- function(t.stat){
  return(pchisq(q = t.stat, df = 1, lower.tail=F)) # P(T.stat >= chi_sq_1)
}

get.j.star <- function(p.vals, alpha, m){
  idx = rep(NA, m)
  for (j in 1:m){
    criterion = alpha/(m - j + 1)
    idx[j] = p.vals[j] <= criterion
  }
  max.idx = max(which(idx))
  return(p.vals[max.idx])
}

