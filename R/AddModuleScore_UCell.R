#' Calculate module enrichment scores from single-cell data (Seurat interface)
#'
#' Given a Seurat object, calculates module/signature enrichment scores at
#' single-cell level using the Mann-Whitney U statistic.
#' UCell scores are normalized U statistics (between 0 and 1), and they are
#' mathematically related to the Area under the ROC curve
#' (see [Mason and Graham](https://doi.org/10.1256/003590002320603584))
#' 
#' In contrast to Seurat's AddModuleScore, which is normalized by binning genes
#' of similar expression at the population level, UCell scores depend 
#' only on the gene expression ranks of individual cell, and therefore they are
#' robust across datasets regardless of dataset composition.
#'
#' @param obj Seurat object
#' @param features A list of signatures, for example:
#'     \code{list( Tcell_signature = c("CD2","CD3E","CD3D"),
#'     Myeloid_signature = c("SPI1","FCER1G","CSF1R"))}
#'     You can also specify positive and negative gene sets by adding a + or -
#'     sign to genes in the signature; see an example below
#' @param maxRank Maximum number of genes to rank per cell; above this rank,
#'     a given gene is considered as not expressed.
#' @param storeRanks Store ranks matrix in Seurat object ('UCellRanks' assay)
#'     for fast subsequent computations. This option may demand large
#'     amounts of RAM.
#' @param w_neg Weight on negative genes in signature. e.g. `w_neg=1` weighs
#'     equally up- and down-regulated genes, `w_neg=0.5` gives 50% less
#'     importance to negative genes
#' @param assay Pull out data from this assay of the Seurat object
#'     (if NULL, use \code{DefaultAssay(obj)})
#' @param slot Pull out data from this slot (layer in v5) of the Seurat object
#' @param chunk.size Number of cells to be processed simultaneously (lower
#'     size requires slightly more computation but reduces memory demands)
#' @param BPPARAM A [BiocParallel::bpparam()] object that tells UCell
#'     how to parallelize. If provided, it overrides the `ncores` parameter.     
#' @param ncores Number of processors to parallelize computation. If
#'     \code{BPPARAM = NULL}, the function uses
#'     \code{BiocParallel::MulticoreParam(workers=ncores)}
#' @param ties.method How ranking ties should be resolved -
#'      passed on to [data.table::frank]
#' @param force.gc Explicitly call garbage collector to reduce memory footprint
#' @param name Name tag that will be appended at the end of each signature
#'    name, "_UCell" by default (e.g. signature score in meta data will be
#'    named: Myeloid_signature_UCell)
#' @return Returns a Seurat object with module/signature enrichment scores
#'     added to object meta data; each score is stored as the corresponding
#'     signature name provided in \code{features} followed by the tag given
#'     in \code{name} (or "_UCell" by default )
#' @examples
#' library(UCell)
#' gene.sets <- list(Tcell = c("CD2","CD3E","CD3D"),
#'                 Myeloid = c("SPI1","FCER1G","CSF1R"))
#' data(sample.matrix)
#' obj <- Seurat::CreateSeuratObject(sample.matrix)                
#' 
#' obj <- AddModuleScore_UCell(obj,features = gene.sets)
#' head(obj[[]])
#'
#' ## Using positive and negative gene sets
#' gene.sets <- list()
#' gene.sets$Tcell_gd <- c("TRDC+","TRGC1+","TRGC2+","TRDV1+",
#'     "TRAC-","TRBC1-","TRBC2-")
#' gene.sets$NKcell <- c("FGFBP2+", "SPON2+", "KLRF1+",
#'     "FCGR3A+", "CD3E-","CD3G-")
#' obj <- AddModuleScore_UCell(obj, features = gene.sets, name=NULL)
#' head(obj$NKcell)
#'
#' @export
AddModuleScore_UCell <- function(obj, features, maxRank=1500,
        chunk.size=100, BPPARAM=NULL, ncores=1, storeRanks=FALSE,
        w_neg=1, assay=NULL, slot="counts", ties.method="average",
        force.gc=FALSE, name="_UCell") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Function 'AddModuleScore_UCell' requires the Seurat package.
            Please install it.", call. = FALSE)
    }  
    features <- check_signature_names(features)
    if (is.null(assay)) {
        assay <- Seurat::DefaultAssay(obj)
    }
    # If rank matrix was pre-computed, evaluate the new signatures
    # from these ranks.
    if ("UCellRanks" %in% Seurat::Assays(obj)) {
        meta.list <- rankings2Uscore(
            Seurat::GetAssayData(obj, layer="counts", assay="UCellRanks"),
            features=features, chunk.size=chunk.size, w_neg=w_neg,
            ncores=ncores, BPPARAM=BPPARAM, force.gc=force.gc, name=name)
    } else {
        layers <- SeuratObject::Layers(obj, assay=assay, search = slot)
        if (is.null(layers)) {
          stop(sprintf("Cannot find layer %s in assay %s", slot, assay))
        }

        meta.list <- lapply(layers, function(x) {
          calculate_Uscore(
            Seurat::GetAssayData(obj, layer=x, assay=assay),
            features=features, maxRank=maxRank,
            chunk.size=chunk.size, w_neg=w_neg,
            ncores=ncores, BPPARAM=BPPARAM, ties.method=ties.method,
            force.gc=force.gc, storeRanks=storeRanks, name=name)
        })
        
        meta.list <- unlist(meta.list, recursive = FALSE)
        #store ranks matrix?
        if (storeRanks==TRUE){
            cells_rankings.merge <- lapply(meta.list,
                function(x) rbind(x[["cells_rankings"]]))
            cells_rankings.merge <- Reduce(cbind, cells_rankings.merge)
            obj[["UCellRanks"]] <- Seurat::CreateAssayObject(
                cells_rankings.merge)
        }
    }
    meta.merge <- lapply(meta.list,function(x) rbind(x[["cells_U"]]))
    meta.merge <- Reduce(rbind, meta.merge)
    obj <- Seurat::AddMetaData(obj, as.data.frame(meta.merge))
    return(obj)
}
