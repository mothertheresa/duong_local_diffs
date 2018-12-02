# import util.R

#' Calculate regions of locally significant differences for two samples.
#' 
#' @param dat data frame with the two groups to compare
#' @param dim an integer for the dimension of the comparison
#' @param m an integer for the number of estimation points
#' @param alpha alpha level for the hypothesis testing
#' @param plot boolean whether the results should be plotted
#' @return list of locally significant differences
#' @examples
#' test_local_hypotheses(data.frame(x1 = rnorm(100), x2 = rnorm(100)))
test_local_hypotheses <- function(dat,
                                  dim = 1,
                                  m = NA,
                                  alpha = 0.05,
                                  plot=F){
  # check inputs
  print("verifying inputs")
  verify_inputs(dat, dim, m, alpha, plot)
  
  # set m
  if (is.na(m)){
    m <- get_m_from_d(dim)
  }
  
  # extract data
  group_names <- names(dat)
  group1 <- dat[,1]
  group2 <- dat[,2]
  
  # set parameters based on all data
  d = density(c(group1, group2),
              kernel="gaussian", n = m, bw="SJ")
  from = min(d$x)
  to = max(d$x)
  xs = d$x
  
  # compute the test statistic
  print("computing f1.hat and f2.hat")
  f1.hat = f.hat(group1, m, from, to)
  f2.hat = f.hat(group2, m, from, to)
  
  # compute p-values
  print("computing p-values")
  p.vals = sapply(xs, function(x) {
    p.val(t.stat(x, group1, group2, f1.hat, f2.hat, m, from, to))
  }
  )
  
  # sort the p-values into ascending order
  sorted.p.vals = sort(p.vals[!is.nan(p.vals)])
  
  # find j.star = argmax {P(j) <= alpha/(m - j + 1)}
  print("computing j.star")
  j.star = get.j.star(sorted.p.vals, alpha, m)
  
  # compute rejection region {x: P(j) <= j.star}
  print("computing rejection region")
  rej.region = xs[p.vals <= j.star]
  
  if (plot){
    main = "TODO"
    plot_results(from, to, m, main, group_names,
                 group1, group2, f1.hat, f2.hat,
                 density.est, xs, sorted.p.vals, j.star)
  }
  
  return(rej.region)
}

plot_results <- function(from, to, m, main, names, group1, group2,
                         f1.hat, f2.hat, density.est, xs, p.vals, j.star){
  print("beginning plotting")
  
  # calculate the density diff for plotting purposes
  print("calculating density difference")
  h1 = hist(group1,
            breaks=seq(from=from, to=to, length.out=400),
            plot=F)
  h2 = hist(group2,
            breaks=h1$breaks, plot=F)
  
  density.diff = h1$density - h2$density
  f1_above_f2 = sapply(xs, function(x) {f1.hat(x) > f2.hat(x)})
  
  # standardize xlims and ylims for all plots
  xlims = c(from+1, to-1)
  ylims = c(0, 1) # todo: make this dynamic
  
  # put plots tighter together
  par(mfrow=c(3,1), omi=c(0.5,0,0,0), plt=c(0.1,0.9,0,0.8))
  print("plot generation")
  
  # plot group 1
  hist(group1, col="lightblue", xlab="",
       ylab=expression("Density"*" (f"[1]*")"), main=main, cex.main=1.5,
       breaks=seq(from=from, to=to, length.out=200),
       freq=F, xlim=xlims, ylim=ylims, lty="blank")
  lines(density.est(group1, m, from, to, give.Rkern = FALSE), lwd=2)
  legend("topright",
         legend=c(paste("kernel density estimation for \n", names[1])),
         lty=c(1), lwd=c(2), box.lty=0)
  
  # plot group 2
  hist(group2, col="pink", xlab="",
       ylab=expression("Density"*" (f"[2]*")"), main="",
       breaks=seq(from=from, to=to, length.out=200),
       freq=F, xlim=xlims, ylim=ylims, lty="blank")
  lines(density.est(group2, m, from, to, give.Rkern = FALSE), lwd=2)
  legend("topright",
         legend=c(paste("kernel density estimation for \n", names[2])),
         lty=c(1), lwd=c(2), box.lty=0)
  
  # plot the density difference
  plot(y = density.diff, x = h1$mids, main="", col="white", xlab="", xlim=xlims, ylim=c(-(ylims[2]+0.5), ylims[2]), bty="n", ylab=expression("f"[1]*" - "*"f"[2]))
  lines(y = density.diff, x = h1$mids, lwd=2)
  
  # and overlay the significantly different region
  pval.colors = rep("white", length(p.vals))
  pval.colors[(p.vals <= j.star) & f1_above_f2] <- "blue"
  pval.colors[(p.vals <= j.star) & !f1_above_f2] <- "red"
  points(y = (p.vals <= j.star)-(ylims[2]+0.5), x = xs, col=pval.colors, xlab="", pch=15, cex=0.7)
  
  print("done")
}

dat <- data.frame(Mutant_1 = rnorm(10000, 0, 1),
                  Mutant_2 = rnorm(10000, 1, 2))
rej.region <- test_local_hypotheses(dat, plot=T)
head(rej.region, 20)