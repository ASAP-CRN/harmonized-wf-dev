
NormalMixCutoff <- function(mix, lower=0, upper=1) {

  f <- function(x, m1, sd1, m2, sd2, p1, p2) {
    dnorm(x, m1, sd1) * p1 - dnorm(x, m2, sd2) * p2 
  }

  cutoff <- rootSolve::uniroot.all(f, lower=lower, upper=upper, 
                                    m1=mix$mu[1], sd1=mix$sigma[1], 
                                    m2=mix$mu[2], sd2=mix$sigma[2], 
                                    p1=mix$lambda[1], p2=mix$lambda[2])

  return(cutoff)

}

