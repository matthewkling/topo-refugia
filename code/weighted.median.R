weighted.median <- function(x, w){
      ord <- order(x)
      x <- x[ord]
      w <- w[ord]
      w <- w/sum(w)
      cs <- cumsum(w)
      x[cs > .5][1]
}
