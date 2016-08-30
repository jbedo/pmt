# Copyright © 2016 by Justin Bedő
#
# Permission to use, copy, modify, and/or distribute this software for any
# purpose with or without fee is hereby granted, provided that the above
# copyright notice and this permission notice appear in all copies.
#
# THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
# WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
# ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
# WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
# ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
# OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

pmt <- function(a, b) {
  if(any(a == 0) | any(b == 0))
    stop("cannot handle 0 counts")

  ce <- function(mu){
    pt1 <- mapply(function(xa, mu, rho) logsumexp(cumsum(log(((1 + rho) * (xa - (0:xa-1))/(2 * mu))))), a, mu, rho)
    pt2 <- mapply(function(xb, mu, rho) logsumexp(cumsum(log((2*mu*rho) / ((1 + rho)*(xb + (1:100)))))), b, mu, rho)
    rho + rho * (exp(pt1) - exp(pt2))
  }

  logsumexp <- function(x) {
    x <- sort(x, decreasing=T)
    log1p(sum(exp(x[-1] - x[1]))) + x[1]
  }

  Ra <- sum(a)
  Rb <- sum(b)
  chi <- a/b
  rho <- rep(Ra/Rb, length(a))

  msk <- a/Ra > b/Rb
  rho[msk] <- 1/rho[msk]
  tmp <- a[msk]
  a[msk] <- b[msk]
  b[msk] <- tmp
  chi[msk] <- 1/chi[msk]

  lb <- (rho - chi + sqrt(((rho - chi)^2 + 4 * chi * rho * (1 + rho)))) / (4 * rho) * b * 0.75
  ub <- (1 + chi) * (1 + rho) / (2 * rho + 1) * b * 1.25

  while(max(ub - lb) > 1e-6){
    mu <- (ub + lb) / 2
    scmu <- ce(mu) >= 0
    ub <- ub * scmu + mu * (1 - scmu)
    lb <- lb * (1 - scmu) + mu * scmu
  }

  t1 <- mapply(function(xa, mu, rho) { i <- 0:xa; logsumexp(i * log(2 * mu) - lfactorial(i) - i * log(1 + rho))}, a, mu, rho)
  t2 <- mapply(function(xb, mu, rho) { j <- xb:(xb*100); logsumexp(j * log(2 * mu * rho) - lfactorial(j) - j * log(1 + rho))}, b, mu, rho)

  (2*mu-t1-t2) / log(10)
}
