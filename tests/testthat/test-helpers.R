# Internal helpers, tested directly

test_that("check_signature_names fills in missing and duplicate names", {
  f <- list(1, 2, 3)
  expect_equal(names(UCell:::check_signature_names(f)),
               c("signature_1", "signature_2", "signature_3"))

  names(f) <- c("a", "", "a")
  # the empty name and the duplicate are replaced, the first "a" is kept
  expect_equal(names(UCell:::check_signature_names(f)),
               c("a", "signature_2", "signature_3"))
})

test_that("get_gene_idx maps names to indices, in signature order", {
  all.genes <- paste0("gene", 1:5)

  expect_equal(UCell:::get_gene_idx(all.genes, c("gene3", "gene1")), c(3, 1))
  # 'impute' marks absent genes with -1, 'skip' drops them
  expect_equal(UCell:::get_gene_idx(all.genes, c("gene2", "nope"),
                                    missing_genes = "impute"), c(2, -1))
  expect_equal(UCell:::get_gene_idx(all.genes, c("gene2", "nope"),
                                    missing_genes = "skip"), 2)
  expect_length(UCell:::get_gene_idx(all.genes, "nope",
                                     missing_genes = "skip"), 0)
})

test_that("split_data.matrix partitions cells without loss or reordering", {
  m <- make_counts(nfeat = 10, ncell = 25)

  for (cs in c(1, 4, 10, 25, 100)) {
    chunks <- UCell:::split_data.matrix(m, chunk.size = cs)
    expect_equal(sum(vapply(chunks, ncol, integer(1))), ncol(m),
                 info = sprintf("chunk.size = %d", cs))
    expect_equal(Reduce(cbind, chunks), m,
                 info = sprintf("chunk.size = %d", cs))
  }
})

test_that("knn_smooth_scores averages over neighbors with an exponential decay", {
  scores <- matrix(c(0, 1, 1, 0), ncol = 1,
                   dimnames = list(paste0("cell", 1:4), "sigA"))
  # each cell's single neighbor is the next one, cyclically
  nn <- list(index = matrix(c(2, 3, 4, 1), ncol = 1))

  # decay=1 gives all the weight to the cell itself: scores are unchanged
  expect_equal(UCell:::knn_smooth_scores(scores, nn, decay = 1)$sigA,
               as.numeric(scores))
  # decay=0 weighs cell and neighbor equally
  expect_equal(UCell:::knn_smooth_scores(scores, nn, decay = 0)$sigA,
               c(0.5, 1, 0.5, 0))
  # up.only keeps the original score wherever smoothing would lower it
  expect_equal(UCell:::knn_smooth_scores(scores, nn, decay = 0,
                                         up.only = TRUE)$sigA,
               c(0.5, 1, 1, 0))

  smoothed <- UCell:::knn_smooth_scores(scores, nn, decay = 0.1)
  expect_s3_class(smoothed, "data.frame")
  expect_equal(rownames(smoothed), rownames(scores))
})
