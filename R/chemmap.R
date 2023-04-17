

generate_chem_glyphs <- function(smiles, outdir=tempdir(),
                                 bondcol="black", showatoms, ...){
  s <- list()
  fnames <- list()
  for(i in 1:length(smiles)) {
    s <- tryCatch(
      {smiles2sdf(smiles[i])},
      error=function(e){
        warning(paste(
          "Could not generate image for SMILES:",
          smiles[i]
        ))
        return(NULL)
      }
    )

    fnames[[i]] <- NA
    if(!is.null(s)) {
      if(length(bondcol) > 1){
        col <- bondcol[i]} else {col <- bondcol}
      fnames[[i]] <- paste0(tempfile("chemmap_glyph_"), ".png")
      png(width = 200, height = 200, fnames[[i]])
      par(bg=NA, mar = rep(0, 4))
      nbonds <- strsplit(s[[1]]@header["Counts_Line"], "\\s+")[[1]][3]
      atomnames <- unique(gsub("_.*", "", rownames(atomblock(s[[1]]))))
      if(!showatoms) {
        ChemmineR::plot(s, colbonds=1:nbonds,
           bondcol=col, asp=1, no_print_atoms=atomnames, ...)
      } else {
        ChemmineR::plot(s, colbonds=1:nbonds,
                        bondcol=col, asp=1, ...)
      }
      dev.off()
    }
  }
  return(unlist(fnames))
}

GeomChemMap <- ggproto(
  "ChemMap", Geom,
  default_aes = aes(colour = "black", fill="black", alpha=1,
                    chemsize=.05, showatoms=T),
  required_aes = c("x", "y", "smiles"),
  draw_key = draw_key_image,
  draw_group = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)
    image_files <- generate_chem_glyphs(data$smiles,
                                        bondcol = coords$colour,
                                        showatoms = data$showatoms)
    myImage <- png::readPNG(image_files)
    rasterGrob(myImage, x=coords$x, y=coords$y,
               width=data$chemsize, height=data$chemsize)
  }
)


#' Title
#'
#' @param mapping
#' @param data
#' @param stat
#' @param position
#' @param na.rm
#' @param show.legend
#' @param inherit.aes
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
geom_chemmap <- function(mapping = NULL, data = NULL, stat = "identity",
                          position = "identity", na.rm = FALSE, show.legend = NA,
                          inherit.aes = TRUE, ...) {

  file.remove(list.files(tempdir(), pattern = "chemmap_glyph.*.png", full.names = T))
  layer(
    geom = GeomChemMap, mapping = mapping,  data = data, stat = stat,
    position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}



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
      sdflist <- list()
      for(i in seq_along(smiles)) {
        sdf <- tryCatch(
          smiles2sdf(smiles[i])[[1]],
          error = function(e) stop(paste0("Open Babel returned an error while converting",
          "the following SMILES to SDF: ", smiles[i]))
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
            M <- fingerprintOB(x@sdf, "FP2")@fpma
            dist(M, "manhattan")
          }
)


