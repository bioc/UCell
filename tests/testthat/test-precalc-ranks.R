# StoreRankings_UCell and the rankings2Uscore() scoring path

test_that("StoreRankings_UCell returns a gene x cell sparse rank matrix", {
  m <- make_counts()
  ranks <- StoreRankings_UCell(m, ncores = 1)

  expect_s4_class(ranks, "dgCMatrix")
  expect_equal(dim(ranks), dim(m))
  expect_equal(rownames(ranks), rownames(m))
  expect_equal(colnames(ranks), colnames(m))
  # maxRank is capped at nrow(m), and ranks at or above it are stored as 0
  expect_true(all(ranks < nrow(m)))
})

test_that("scoring from pre-calculated ranks matches scoring from the matrix", {
  m <- make_graded(nfeat = 50)  # tie-free, so ranks are integers
  feats <- list(top = rownames(m)[1:5], mid = rownames(m)[20:25])

  ranks <- StoreRankings_UCell(m, ncores = 1)
  expect_equal(
    ScoreSignatures_UCell(features = feats, precalc.ranks = ranks, ncores = 1),
    ScoreSignatures_UCell(m, features = feats, ncores = 1)
  )
})

test_that("chunk.size does not change the scores", {
  m <- make_counts(ncell = 40)
  ref <- ScoreSignatures_UCell(m, features = sig(m), ncores = 1)

  for (cs in c(1, 7, 40, 1000)) {
    expect_equal(
      ScoreSignatures_UCell(m, features = sig(m), ncores = 1, chunk.size = cs),
      ref,
      info = sprintf("chunk.size = %d", cs)
    )
  }
})
