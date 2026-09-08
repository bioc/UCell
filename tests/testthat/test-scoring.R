# Scoring semantics of ScoreSignatures_UCell / u_stat. The graded matrix has
# no ties, so the expected scores below are exact rather than approximate.

test_that("UCell scores follow the normalized U statistic", {
  m <- make_graded(nfeat = 50)
  # maxRank is capped at nrow(m), so the top 5 genes hold ranks 1:5 and reach
  # the maximum score, while the bottom 5 genes sit just above the minimum
  feats <- list(top = rownames(m)[1:5], bottom = rownames(m)[46:50])
  u <- ScoreSignatures_UCell(m, features = feats, ncores = 1)

  expect_equal(unname(u[, "top_UCell"]), rep(1, ncol(m)))
  expect_true(all(u >= 0 & u <= 1))
  expect_true(all(u[, "top_UCell"] > u[, "bottom_UCell"]))
})

test_that("negative gene sets are subtracted and weighted by w_neg", {
  m <- make_graded(nfeat = 50)
  pos <- list(s = rownames(m)[1:5])
  mixed <- list(s = c(paste0(rownames(m)[1:5], "+"),
                      paste0(rownames(m)[6:10], "-")))

  ref <- ScoreSignatures_UCell(m, features = pos, ncores = 1)
  expect_equal(ScoreSignatures_UCell(m, features = mixed, w_neg = 0,
                                     ncores = 1), ref)
  expect_true(all(ScoreSignatures_UCell(m, features = mixed, w_neg = 1,
                                        ncores = 1) < ref))
  # a heavily weighted negative set drives the score below zero; it is clipped
  expect_equal(unname(ScoreSignatures_UCell(m, features = mixed, w_neg = 5,
                                            ncores = 1)[, 1]),
               rep(0, ncol(m)))
  expect_error(ScoreSignatures_UCell(m, features = pos, w_neg = -1), "w_neg")
})

test_that("missing_genes controls how absent genes are handled", {
  m <- make_graded(nfeat = 50)
  absent <- paste0("absent", 1:5)
  feats <- list(s = c(rownames(m)[1:5], absent))

  # 'impute' treats absent genes as ranked at maxRank, penalizing the score
  imputed <- ScoreSignatures_UCell(m, features = feats, ncores = 1,
                                   missing_genes = "impute")
  # 'skip' drops them, leaving the score of the 5 genes that are present
  skipped <- ScoreSignatures_UCell(m, features = feats, ncores = 1,
                                   missing_genes = "skip")

  expect_true(all(skipped > imputed))
  expect_equal(skipped,
               ScoreSignatures_UCell(m, features = list(s = rownames(m)[1:5]),
                                     ncores = 1))
  # a signature with nothing but absent genes scores 0 either way
  only.absent <- list(s = absent)
  expect_equal(unname(ScoreSignatures_UCell(m, features = only.absent,
                                            ncores = 1)[, 1]),
               rep(0, ncol(m)))
  expect_equal(unname(ScoreSignatures_UCell(m, features = only.absent,
                                            ncores = 1,
                                            missing_genes = "skip")[, 1]),
               rep(0, ncol(m)))
})

test_that("signature names and the name suffix are reflected in colnames", {
  m <- make_counts()
  expect_equal(colnames(ScoreSignatures_UCell(m, features = sig(m),
                                              ncores = 1)),
               c("sigA_UCell", "sigB_UCell"))
  expect_equal(colnames(ScoreSignatures_UCell(m, features = sig(m),
                                              ncores = 1, name = "")),
               c("sigA", "sigB"))
  # unnamed signatures get default names from check_signature_names()
  expect_equal(colnames(ScoreSignatures_UCell(m, features = unname(sig(m)),
                                              ncores = 1)),
               c("signature_1_UCell", "signature_2_UCell"))
})

test_that("signatures longer than maxRank are rejected", {
  m <- make_counts()
  expect_error(
    ScoreSignatures_UCell(m, features = list(s = rownames(m)), maxRank = 5),
    "maxRank")
})
