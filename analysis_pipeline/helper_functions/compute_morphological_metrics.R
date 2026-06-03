# Required packages are Seurat, sf, sp

# extract cell area from Seurat object
extract_cell_area <- function(obj) {
  samples <- names(obj@images)
  cell_area <- unlist(
    lapply(samples, function(fov) {
      polygons <- obj@images[[fov]]$segmentation@polygons
      vapply(polygons, function(poly) poly@area, numeric(1))
    })
  )
  return(cell_area)
}

compute_aspect_ratio <- function(obj) {
  samples <- names(obj@images)
  ratio_list <- lapply(samples, function(fov) {
    polys <- obj[[fov]]$segmentation@polygons
    if (length(polys) == 0) return(NULL)
    
    # extract the bounding box from Polygons
    bb_matrix <- sapply(polys, function(p) {
      bb <- sp::bbox(p)
      c(xmin = bb["x", "min"],
        xmax = bb["x", "max"],
        ymin = bb["y", "min"], 
        ymax = bb["y", "max"])
    })
    
    # calculate the aspect ratio
    apply(bb_matrix, 2, function(bb) {
      bbx <- bb["xmax"] - bb["xmin"]
      bby <- bb["ymax"] - bb["ymin"]
      # standardise aspect ratio so that it's between 0 and 1 and higher value indicates more squared shape
      min(bbx, bby) / max(bbx, bby)
    })
  })
  
  # combine all FOVs into a single names vector
  return(unlist(ratio_list, use.names = TRUE))
}

computeSpatialOutlier <- function(
    object,
    computeBy = NULL,
    method = c("mc", "scuttle", "both"),
    mcDoScale = FALSE,
    scuttleType = c("both", "lower", "higher")
) {
  
  method <- match.arg(method)
  scuttleType <- match.arg(scuttleType)
  options(mc_doScale_quiet = TRUE)
  
  ## ---- Determine object type & extract metadata ----
  if (inherits(object, "SpatialExperiment")) {
    md <- as.data.frame(colData(object))
    cell_names <- colnames(object)
    updateMetaData <- function(object, md) {
      colData(object) <- md
      return(object)
    }
    
  } else if (inherits(object, "Seurat")) {
    md <- object@meta.data
    cell_names <- rownames(md)
    updateMetaData <- function(object, md) {
      object@meta.data <- md
      return(object)
    }
    
  } else {
    stop("object must be a SpatialExperiment or Seurat object")
  }
  
  stopifnot(!is.null(computeBy))
  stopifnot(computeBy %in% colnames(md))
  
  cdcol <- md[[computeBy]]
  
  mcfl <- scuttlefl <- FALSE
  switch(method,
         both    = { mcfl <- scuttlefl <- TRUE },
         mc      = { mcfl <- TRUE },
         scuttle = { scuttlefl <- TRUE }
  )
  
  ## ---- Medcouple method ----
  if (mcfl) {
    skw <- e1071::skewness(cdcol, na.rm = TRUE)
    if (skw > -1 & skw < 1)
      warning("Distribution is symmetric: mc is intended for asymmetric distributions.")
    
    mcval <- robustbase::mc(cdcol, doScale = mcDoScale, na.rm = TRUE)
    if (mcval <= -0.6 | mcval >= 0.6)
      stop("mc is: ", round(mcval, 4), " — adjusted boxplot assumptions violated")
    
    names(cdcol) <- cell_names
    outl <- robustbase::adjbox(cdcol, plot = FALSE)
    
    outsmc <- rep("NO", nrow(md))
    outsmc[rownames(md) %in% names(outl$out)] <-
      ifelse(outl$out <= outl$fence[1], "LOW", "HIGH")
    
    outlier_mc <- scuttle::outlier.filter(outsmc)
    
    thrs <- as.numeric(outl$fence)
    names(thrs) <- c("lower", "higher")
    attr(outlier_mc, "thresholds") <- thrs
    
    md[[paste0(computeBy, "_outlier_mc")]] <- outlier_mc
  }
  
  ## ---- Scuttle MAD method ----
  if (scuttlefl) {
    outssc <- scuttle::isOutlier(cdcol, type = scuttleType)
    
    sctri <- rep("NO", nrow(md))
    sctri[outssc & cdcol <= attr(outssc, "thresholds")[1]] <- "LOW"
    sctri[outssc & cdcol >= attr(outssc, "thresholds")[2]] <- "HIGH"
    
    outlier_sc <- scuttle::outlier.filter(sctri)
    attr(outlier_sc, "thresholds") <- attr(outssc, "thresholds")
    
    md[[paste0(computeBy, "_outlier_sc")]] <- outlier_sc
  }
  
  ## ---- Write metadata back ----
  object <- updateMetaData(object, md)
  return(object)
}
