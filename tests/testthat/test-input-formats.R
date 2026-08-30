set.seed(42)

make_counts <- function(nfeat = 60, ncell = 40) {
  m <- matrix(rpois(nfeat * ncell, lambda = 3), nrow = nfeat, ncol = ncell)
  dimnames(m) <- list(paste0("gene", seq_len(nfeat)), paste0("cell", seq_len(ncell)))
  m
}

sig <- function(m) list(sigA = rownames(m)[1:10], sigB = rownames(m)[11:20])

test_that("ScoreSignatures_UCell accepts the classic input formats", {
  m <- make_counts()
  ref <- ScoreSignatures_UCell(m, features = sig(m), ncores = 1)
  expect_true(is.matrix(ref))
  expect_identical(rownames(ref), colnames(m))
  expect_equal(
    ScoreSignatures_UCell(Matrix::Matrix(m, sparse = TRUE), features = sig(m), ncores = 1),
    ref
  )
  expect_equal(
    ScoreSignatures_UCell(as.data.frame(m), features = sig(m), ncores = 1),
    ref
  )
})

test_that("ScoreSignatures_UCell accepts other matrix-like classes", {
  skip_if_not_installed("DelayedArray")
  m <- make_counts()
  ref <- ScoreSignatures_UCell(m, features = sig(m), ncores = 1)
  # A DelayedMatrix used to be rejected with "Unrecognized input format",
  # even though calculate_Uscore() coerces whatever it is given
  d <- DelayedArray::DelayedArray(Matrix::Matrix(m, sparse = TRUE))
  expect_equal(ScoreSignatures_UCell(d, features = sig(m), ncores = 1), ref)
  # dgeMatrix, a dense Matrix class, was rejected for the same reason
  expect_equal(
    ScoreSignatures_UCell(Matrix::Matrix(m, sparse = FALSE), features = sig(m), ncores = 1),
    ref
  )
})

test_that("StoreRankings_UCell accepts other matrix-like classes", {
  skip_if_not_installed("DelayedArray")
  m <- make_counts()
  ref <- StoreRankings_UCell(m, ncores = 1)
  d <- DelayedArray::DelayedArray(Matrix::Matrix(m, sparse = TRUE))
  expect_equal(StoreRankings_UCell(d, ncores = 1), ref)
})

test_that("input that is not matrix-like is still rejected", {
  m <- make_counts()
  expect_error(ScoreSignatures_UCell(letters, features = sig(m)), "Unrecognized input format")
  expect_error(ScoreSignatures_UCell(NULL, features = sig(m)), "Unrecognized input format")
  expect_error(StoreRankings_UCell(letters), "Unrecognized input format")
})

test_that("a Seurat object gets a pointed error rather than a generic one", {
  skip_if_not_installed("SeuratObject")
  m <- make_counts()
  obj <- suppressWarnings(SeuratObject::CreateSeuratObject(counts = Matrix::Matrix(m, sparse = TRUE)))
  expect_error(ScoreSignatures_UCell(obj, features = sig(m)), "AddModuleScore_UCell")
  expect_error(StoreRankings_UCell(obj), "LayerData")
})

test_that("scores are unchanged for a sparse input after the coercion change", {
  m <- make_counts()
  sp <- Matrix::Matrix(m, sparse = TRUE)
  expect_equal(
    ScoreSignatures_UCell(sp, features = sig(m), ncores = 1),
    ScoreSignatures_UCell(as.matrix(sp), features = sig(m), ncores = 1)
  )
})
