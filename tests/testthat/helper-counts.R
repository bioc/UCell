# Shared fixtures for the UCell test suite.

# Random counts, with ties, as a plain matrix
make_counts <- function(nfeat = 60, ncell = 40, seed = 42) {
  set.seed(seed)
  m <- matrix(rpois(nfeat * ncell, lambda = 3), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("gene", seq_len(nfeat)),
                      paste0("cell", seq_len(ncell)))
  m
}

# Counts with a known, tie-free ranking: gene1 is the top-expressed gene in
# every cell, gene2 the second, and so on. This makes UCell scores exactly
# predictable, and makes ranks independent of 'ties.method'.
make_graded <- function(nfeat = 50, ncell = 5) {
  m <- matrix(rep(rev(seq_len(nfeat)), times = ncell),
              nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("gene", seq_len(nfeat)),
                      paste0("cell", seq_len(ncell)))
  m
}

sig <- function(m) list(sigA = rownames(m)[1:10], sigB = rownames(m)[11:20])
