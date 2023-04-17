

#' Title
#'
#' @slot sdf SDFset.
#' @slot smiles SMIset.
#'
#' @return
#' @export
#'
#' @examples
chemmol <- setClass("chemmol",
                    slots=list(
                      sdf="SDFset",
                      smiles="SMIset"
                      ),
                    prototype=list(
                      sdf=NULL,
                      smiles=NULL
                    )
)

setMethod(initialize, "chemmol", function(.Object, sdf=NULL, smiles=NULL, ...) {
  if(is.null(sdf) && is.null(smiles))
      stop("at least one of sdf or smiles must be non-NULL")

  if(is.null(smiles)){
    smiles <- list()
    for(i in seq_along(sdf)) smiles[[i]] <- sdf2smiles(sdf[i])
    .Object@smiles <- new("SMIset", smilist=smiles)
  }

  if(is.null(sdf)){
    print("version check")
    sdflist <- list()
    for(i in seq_along(smiles)) {
      sdf <- tryCatch(
        smiles2sdf(smiles[i])[[1]],
        error = function(e) NA
        )
      sdflist[[i]] <- sdf
    }
    .Object@sdf <- new("SDFset", SDF=sdflist)
  }

  #callNextMethod(.Object, sdf=.Object@sdf, smiles=.Object@smiles, ...)
  callNextMethod(...)
})

setMethod("dist",
          signature(x = "chemmol"),
          function (x, method = "FP2")
          {
            M <- fingerprintOB(chemm@sdf, "FP2")@fpma
            dist(M, "manhattan")
          }
)



# smiles <- read.csv(
#   paste0("/Users/francesco/Dropbox/git/Cell exhaustion/data/in/GSE70138_Broad_LINCS_pert_info.txt"),
#   sep="\t")$canonical_smiles[1:10]
# names(smiles) <- as.character(1:length(smiles))
# smiles[[3]] <- "reserved"
# smi <- new("SMIset", smilist=as.list(smiles))
# sdf <- smiles2sdf(smi)
# chemmol(sdf, smi)
# chemmol(sdf=sdf)
# chemm <- chemmol(smiles = smi)
