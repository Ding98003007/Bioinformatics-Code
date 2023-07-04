fpkmToTpm <- function(fpkm) { # nolint
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
