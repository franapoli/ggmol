#' Generate Chemical Glyphs
#'
#' This function generates PNG images of chemical structures from SMILES strings using the `ChemmineR` package.
#'
#' @param smiles A character vector of SMILES strings.
#' @param outdir Directory where the output PNG files will be saved. Default is a temporary directory.
#' @param bond_col Color of the bonds in the chemical structure. Default is "black".
#' @param show_atoms Logical value indicating whether to show atom labels. Default is TRUE.
#' @param resolution Resolution of the output images. Default is 200.
#' @param ... Additional arguments passed to `ChemmineR::plot`.
#' @return A character vector of file paths to the generated PNG images.
#' @export
#' @examples
#' smiles <- c("CCO", "CCC")
#' generate_chem_glyphs(smiles)
generate_chem_glyphs <- function(smiles, outdir = tempdir(), bond_col = "black", show_atoms = TRUE, resolution = 200, ...) {
  s <- list()
  fnames <- list()
  for (i in 1:length(smiles)) {
    s <- tryCatch({
      smiles2sdf(smiles[i])
    }, error = function(e) {
      warning(paste("Could not generate image for SMILES:", smiles[i]))
      return(NULL)
    })

    fnames[[i]] <- NA
    if (!is.null(s)) {
      if (length(bond_col) > 1) {
        col <- bond_col[i]
      } else {
        col <- bond_col
      }
      fnames[[i]] <- paste0(tempfile("chemmap_glyph_"), ".png")
      png(width = resolution, height = resolution, fnames[[i]])
      par(bg = NA, mar = rep(0, 4))
      nbonds <- strsplit(s[[1]]@header["Counts_Line"], "\\s+")[[1]][3]
      atom_names <- unique(gsub("_.*", "", rownames(atomblock(s[[1]]))))
      if (!show_atoms) {
        ChemmineR::plot(s, colbonds = 1:nbonds, bondcol = col, asp = 1, no_print_atoms = atom_names, regenCoords=T, ...)
      } else {
        ChemmineR::plot(s, colbonds = 1:nbonds, bondcol = col, asp = 1, regenCoords=T, ...)
      }
      dev.off()
    }
  }
  return(unlist(fnames))
}

#' GeomMol: A ggplot2 Geom for Plotting Chemical Structures
#'
#' This `ggplot2` geom adds chemical structure images to a plot based on SMILES strings.
#'
#' @inheritParams ggplot2::layer
#' @param mapping Set of aesthetic mappings created by `aes()`.
#' @param data The data to be displayed in this layer.
#' @param stat The statistical transformation to use on the data for this layer. Default is "identity".
#' @param position Position adjustment. Default is "identity".
#' @param na.rm If `FALSE`, removes missing values with a warning. If `TRUE`, silently removes missing values. Default is FALSE.
#' @param show.legend Logical. Should this layer be included in the legends? Default is `NA`.
#' @param inherit.aes If `FALSE`, overrides the default aesthetics. Default is `TRUE`.
#' @param ... Other arguments passed on to `layer()`.
#' @return A `ggplot2` layer.
#' @export
#' @examples
#' library(ggplot2)
#' smiles <- c("CCO", "CCC")
#' df <- data.frame(x = 1:2, y = 1:2, smiles = smiles)
#' ggplot(df, aes(x, y, smiles = smiles)) +
#'   geom_mol()
geom_mol <- function(mapping = NULL, data = NULL, stat = "identity", position = "identity", na.rm = FALSE, show.legend = NA, inherit.aes = TRUE, ...) {
  #file.remove(list.files(tempdir(), pattern = "chemmap_glyph.*.png", full.names = TRUE))
  layer(
    geom = GeomMol, mapping = mapping, data = data, stat = stat, position = position, show.legend = show.legend, inherit.aes = inherit.aes,
    params = list(na.rm = na.rm, ...)
  )
}

#' GeomMol ggproto Object
#'
#' The ggproto object defining the GeomMol for ggplot2.
GeomMol <- ggproto(
  "GeomMol", Geom,
  default_aes = aes(colour = "black", fill = "black", image_files=NULL, alpha = 1, chemsize = .05, show_atoms = TRUE, resolution = 200),
  required_aes = c("x", "y"),
  draw_key = draw_key_image,
  draw_group = function(data, panel_params, coord) {
    coords <- coord$transform(data, panel_params)
    if(is.null(data$image_files)){
    image_files <- generate_chem_glyphs(data$smiles, bond_col = coords$colour, show_atoms = data$show_atoms, resolution = data$resolution)
    } else image_files <- data$image_files
    myImage <- png::readPNG(image_files)
    rasterGrob(myImage, x = coords$x, y = coords$y, width = data$chemsize, height = data$chemsize)
  }
)

#' Chemmol Class
#'
#' An S4 class to represent chemical molecules using SDF and SMILES formats.
#'
#' @slot sdf An object of class `SDFset` containing the SDF representation.
#' @slot smiles An object of class `SMIset` containing the SMILES representation.
#' @export
#' @examples
#' sdf <- read.SDFset("path_to_sdf_file")
#' smiles <- c("CCO", "CCC")
#' chem <- Chemmol(sdf = sdf, smiles = smiles)
Chemmol <- setClass("Chemmol",
                    slots = list(
                      sdf = "SDFset",
                      smiles = "SMIset"
                    ),
                    prototype = list(
                      sdf = NULL,
                      smiles = NULL
                    ))

#' Initialize Chemmol Object
#'
#' Initialize a new object of class `Chemmol`.
#'
#' @param .Object The `Chemmol` object to initialize.
#' @param sdf An optional `SDFset` object.
#' @param smiles An optional character vector of SMILES strings.
#' @param ... Additional arguments passed to `callNextMethod()`.
#' @export
setMethod(initialize, "Chemmol", function(.Object, sdf = NULL, smiles = NULL, ...) {
  if (is.null(sdf) && is.null(smiles))
    stop("at least one of sdf or smiles must be non-NULL")

  if (is.null(smiles)) {
    smiles <- list()
    for (i in seq_along(sdf)) smiles[[i]] <- sdf2smiles(sdf[i])
    .Object@smiles <- new("SMIset", smilist = smiles)
  }

  if (is.null(sdf)) {
    sdflist <- list()
    for (i in seq_along(smiles)) {
      sdf <- tryCatch(
        smiles2sdf(smiles[i])[[1]],
        error = function(e) stop(paste0("Open Babel returned an error while converting",
                                        "the following SMILES to SDF: ", smiles[i]))
      )
      sdflist[[i]] <- sdf
    }
    .Object@sdf <- new("SDFset", SDF = sdflist)
  }

  callNextMethod(...)
})

#' Calculate Distance Between Chemmol Objects
#'
#' Compute a distance matrix for objects of class `Chemmol` using specified fingerprint method.
#'
#' @param x An object of class `Chemmol`.
#' @param method The fingerprint method to use. Default is "FP2".
#' @return A distance matrix.
#' @export
#' @examples
#' chem <- Chemmol(sdf = sdf, smiles = smiles)
#' dist_matrix <- dist(chem, method = "FP2")
setMethod("dist",
          signature(x = "Chemmol"),
          function(x, method = "FP2") {
            print("running dist")
            M <- fingerprintOB(x@sdf, "FP2")@fpma
            stats::dist(M, "manhattan")
          })
