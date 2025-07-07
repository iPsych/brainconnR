#' Plots brains
#'
#' \code{brainconn} plots and returns a ggraph object of plotted brain connectivity matrix..
#' @author Sidhant Chopra
#'
#' @param atlas Either a string of one of the included atlases \code{brainconn::list_atlases()} or a \code{data.frame()} that meets specifications, see \code{vignette("brainconn")}
#' @param background 'ICBM152', currently the only background option
#' @param background.alpha Number between 0-1 to  set the transparency of the background.
#' @param view A sting to choose the view. Can be any of these: c("ortho", "top", "bottom", "left", "right")
#' @param conmat A adjacency matrix. Can be binary, weights, directed or undirected. see example_* data.
#' @param node.size A integer that determines the diameter of the nodes. Can also be a vector of integers with a length equal to the number of ROIs in the atlas
#' @param node.color A string that sets the node color. e.g. "blue". If set to "network", then nodes will be colored according to the network column of the atlas
#' @param edge.color.weighted A boolean that applies when the conmat is weighted. if \code{TRUE}, then edges will be colored according to the weight \code{FALSE}, the edges will be sized according to weight.
#' @param all.nodes if \code{TRUE}, then all nodes will be shown be hemisphere without ticks. If \code{FALSE}, then only nodes with connecting edges will be shown.
#' @param edge.color A string that sets the edge color. e.g. "blue".
#' @param edge.alpha Number between 0-1 to  set the transparency of the edges.
#' @param edge.width Number to set the width of the edges.
#' @param labels if \code{TRUE}, ROI labels for all visible nodes will be shown. If \code{FALSE}, then no labes will be shown.
#' @param show.legend if \code{TRUE}, legend will be shown. If \code{FALSE}, then no legend will be shown.
#' @param thr a optional value to set a threshold on the conmat (e.g. edges with a weighted value lower than the one set here will not be shown)
#' @param uthr a optional value to set a upper threshold on the conmat (e.g. edges with a weighted value higher than the one set here will not be shown)
#' @param scale.edge.width If \code{edge.color.weighted=FALSE}, you can use this rescale the edge.width according to weight. e.g. \code{scale.edge.width = c(1,3)}
#' @param label.size If labels=TRUE then, \code{label.size} can can be set as in integer to control the size of the labels.
#' @param label.edge.weight if \code{TRUE}, then the edge weight will be labels along the edge.
#' @param use_mesh_background if \code{TRUE}, generates background from 3D mesh instead of using pre-generated .rda files. Note: Pre-generated backgrounds provide higher quality; mesh generation is primarily a fallback option.
#'
#' @return a ggraph object
#'
#' @import ggraph
#' @import cowplot
#' @import grid
#' @importFrom grDevices rgb
#' @importFrom scales colour_ramp brewer_pal rescale
#' @examples
#' library(brainconn)
#' x <- example_unweighted_undirected
#' # Use pre-generated backgrounds (default - highest quality)
#' brainconn(atlas ="schaefer300_n7", conmat=x, node.size = 3, view="ortho")
#' # Use 3D mesh backgrounds (fallback option, lower quality but more flexible)
#' brainconn(atlas ="schaefer300_n7", conmat=x, node.size = 3, view="ortho", use_mesh_background = TRUE)
#' @export
brainconn <- function(atlas,
                      background='ICBM152',
                      view ="top",
                      conmat=NULL,
                      # interactive = F,
                      node.size=4,
                      node.color="network",
                      all.nodes=FALSE,
                      edge.color="black",
                      edge.alpha=0.8,
                      edge.width=1,
                      edge.color.weighted=FALSE,
                      labels=FALSE,
                      show.legend=TRUE,
                      thr=NULL,
                      uthr=NULL,
                      scale.edge.width = NULL,
                      label.size=1.5,
                      label.edge.weight = FALSE,
                      background.alpha = 1,
                      bg_xmax=0,
                      bg_xmin=0,
                      bg_ymax=0,
                      bg_ymin=0,
                      use_mesh_background = FALSE) {

  # Helper function to generate 2D background from 3D mesh
  generate_background <- function(view, background.alpha, bg_xmin, bg_ymin, bg_xmax, bg_ymax, use_mesh = FALSE) {
    
    # First try to load existing background if available and not forcing mesh
    if (!use_mesh) {
      bg_name <- paste0("ICBM152_", view)
      if (exists(bg_name)) {
        m <- get(bg_name)
        w <- matrix(rgb(m[,,1], m[,,2], m[,,3], m[,,4] * background.alpha), nrow=dim(m)[1])
        return(rasterGrob(w))
      }
    }
    # Load 3D mesh data once
    vb <- get("ICBM152_mesh_vb")
    it <- get("ICBM152_mesh_it")
    
    # Extract 3D coordinates
    x_3d <- vb[1,]
    y_3d <- vb[2,]
    z_3d <- vb[3,]
    
    # Apply view-specific coordinate transformations
    if (view == "top") {
      x_2d <- x_3d
      y_2d <- y_3d
      depth <- z_3d
    } else if (view == "bottom") {
      x_2d <- x_3d * -1
      y_2d <- y_3d
      depth <- z_3d * -1
    } else if (view == "front") {
      x_2d <- x_3d
      y_2d <- z_3d
      depth <- y_3d
    } else if (view == "back") {
      x_2d <- x_3d * -1
      y_2d <- z_3d
      depth <- y_3d * -1
    } else if (view == "left") {
      x_2d <- y_3d * -1
      y_2d <- z_3d
      depth <- x_3d
    } else if (view == "right") {
      x_2d <- y_3d
      y_2d <- z_3d
      depth <- x_3d * -1
    }
    
    # Set view-specific limits (2% expansion for top/bottom, 12-15% for others)
    if (view %in% c("top", "bottom")) {
      xmax = 73 + bg_xmax; xmin = -78 + bg_xmin  # 2% expansion for top/bottom
      ymax = 77 + bg_ymax; ymin = -111 + bg_ymin
    } else if (view %in% c("front", "back")) {
      xmax = 78 + bg_xmax; xmin = -78 + bg_xmin  # 12% expansion
      ymax = 88 + bg_ymax; ymin = -56 + bg_ymin
    } else if (view == "left") {
      xmax = 116 + bg_xmax; xmin = -85 + bg_xmin  # 15% expansion
      ymax = 87 + bg_ymax; ymin = -60 + bg_ymin
    } else if (view == "right") {
      xmax = 121 + bg_xmax; xmin = -158 + bg_xmin  # 15% expansion
      ymax = 87 + bg_ymax; ymin = -60 + bg_ymin
    }
    
    # Create brain outline by rasterizing mesh triangles
    resolution <- 1024  # Higher resolution for better quality
    
    # Create coordinate matrices for the 2D grid
    x_grid <- seq(xmin, xmax, length.out = resolution)
    y_grid <- seq(ymin, ymax, length.out = resolution)
    
    # Initialize the image matrix with white background
    brain_img <- array(1, dim = c(resolution, resolution, 4))
    brain_mask <- array(0, dim = c(resolution, resolution))  # Binary mask for brain regions
    
    # Calculate triangle face depths for z-buffering and coloring
    triangle_depths <- apply(it, 2, function(face_indices) {
      mean(depth[face_indices])
    })
    
    # Create brain silhouette first by marking all brain-covered pixels
    for (i in 1:ncol(it)) {
      # Get triangle vertices in 2D
      v1_idx <- it[1, i]
      v2_idx <- it[2, i]
      v3_idx <- it[3, i]
      
      tri_x <- c(x_2d[v1_idx], x_2d[v2_idx], x_2d[v3_idx])
      tri_y <- c(y_2d[v1_idx], y_2d[v2_idx], y_2d[v3_idx])
      
      # Skip triangles outside our view bounds
      if (max(tri_x) < xmin || min(tri_x) > xmax || 
          max(tri_y) < ymin || min(tri_y) > ymax) {
        next
      }
      
      # Convert triangle vertices to pixel coordinates
      tri_px <- round((tri_x - xmin) / (xmax - xmin) * (resolution - 1)) + 1
      tri_py <- round((tri_y - ymin) / (ymax - ymin) * (resolution - 1)) + 1
      
      # Ensure coordinates are within bounds
      tri_px <- pmax(1, pmin(resolution, tri_px))
      tri_py <- pmax(1, pmin(resolution, tri_py))
      
      # Find bounding box
      x_min_px <- max(1, min(tri_px) - 2)
      x_max_px <- min(resolution, max(tri_px) + 2)
      y_min_px <- max(1, min(tri_py) - 2)
      y_max_px <- min(resolution, max(tri_py) + 2)
      
      # Fill triangle using barycentric coordinates
      for (px in x_min_px:x_max_px) {
        for (py in y_min_px:y_max_px) {
          # Convert back to world coordinates for triangle test
          world_x <- xmin + (px - 1) / (resolution - 1) * (xmax - xmin)
          world_y <- ymin + (py - 1) / (resolution - 1) * (ymax - ymin)
          
          # Barycentric coordinate test
          denom <- (tri_y[2] - tri_y[3]) * (tri_x[1] - tri_x[3]) + 
                   (tri_x[3] - tri_x[2]) * (tri_y[1] - tri_y[3])
          
          if (abs(denom) > 1e-10) {
            a <- ((tri_y[2] - tri_y[3]) * (world_x - tri_x[3]) + 
                  (tri_x[3] - tri_x[2]) * (world_y - tri_y[3])) / denom
            b <- ((tri_y[3] - tri_y[1]) * (world_x - tri_x[3]) + 
                  (tri_x[1] - tri_x[3]) * (world_y - tri_y[3])) / denom
            c <- 1 - a - b
            
            # Point is inside triangle
            if (a >= -0.01 && b >= -0.01 && c >= -0.01) {  # Small tolerance for edge pixels
              brain_mask[py, px] <- 1
            }
          }
        }
      }
    }
    
    # Apply morphological operations to smooth the brain outline
    # Simple dilation followed by erosion to smooth edges
    kernel_size <- 3
    smoothed_mask <- brain_mask
    
    # Simple smoothing by averaging neighborhood
    for (px in (kernel_size+1):(resolution-kernel_size)) {
      for (py in (kernel_size+1):(resolution-kernel_size)) {
        if (brain_mask[py, px] == 1) {
          # Check if we're near an edge and smooth
          neighborhood <- brain_mask[(py-kernel_size):(py+kernel_size), 
                                   (px-kernel_size):(px+kernel_size)]
          if (mean(neighborhood) < 0.8) {  # Near edge
            smoothed_mask[py, px] <- 0.7  # Slight anti-aliasing
          }
        }
      }
    }
    
    # Create realistic brain coloring based on distance from center and depth variation
    center_x <- resolution / 2
    center_y <- resolution / 2
    
    for (px in 1:resolution) {
      for (py in 1:resolution) {
        if (smoothed_mask[py, px] > 0) {
          # Distance from center for radial shading
          dist_from_center <- sqrt((px - center_x)^2 + (py - center_y)^2) / (resolution/2)
          
          # Add some noise for texture
          noise <- sin(px * 0.1) * cos(py * 0.1) * 0.1
          
          # Create realistic gray matter coloring
          base_gray <- 0.4 + 0.3 * (1 - dist_from_center) + noise
          base_gray <- pmax(0.2, pmin(0.8, base_gray))
          
          if (smoothed_mask[py, px] < 1) {
            # Edge pixels - create anti-aliasing
            alpha_val <- smoothed_mask[py, px] * background.alpha
            brain_img[py, px, 1] <- base_gray  # R
            brain_img[py, px, 2] <- base_gray  # G  
            brain_img[py, px, 3] <- base_gray  # B
            brain_img[py, px, 4] <- alpha_val  # Alpha for anti-aliasing
          } else {
            brain_img[py, px, 1] <- base_gray  # R
            brain_img[py, px, 2] <- base_gray  # G  
            brain_img[py, px, 3] <- base_gray  # B
            brain_img[py, px, 4] <- background.alpha  # Alpha
          }
        }
      }
    }
    
    # Convert to rasterGrob
    w <- matrix(rgb(brain_img[,,1], brain_img[,,2], brain_img[,,3], 
                    brain_img[,,4]), nrow = resolution)
    return(rasterGrob(w))
  }

  ifelse(is.character(atlas), data <- get(atlas), data <- atlas)

  #set background (add ability to add custom background image)
  if(background != "ICBM152"  && view == "ortho") {
    stop("Custom background image detected, view cannot be 'ortho', please select top,
        bottom, left, right, front or back.")}

  if (!is.null(thr)) {conmat[conmat < thr] <- 0} #lower threshold graph
  if (!is.null(uthr)) {conmat[conmat > thr] <- 0} #upper threshold graph
  #loop three times for the three vies that make ortho view

  if (view == "ortho") {
    ortho_list <- list()
    ortho_views  <- c("top", "left", "front")
    for (v in 1:3) {
      current_view <- ortho_views[v]
      
      # Generate background from 3D mesh
      background <- generate_background(current_view, background.alpha, 
                                      bg_xmin, bg_ymin, bg_xmax, bg_ymax, use_mesh_background)

      #if no conmat is provided, build nparc x  nparc empty one
      nparc <- dim(data)[1]
      if (!exists("conmat")){conmat <- matrix(0L, nrow=nparc, ncol=nparc)
      }

      #convert conmat to matrix
      conmat <- as.matrix(conmat)

      #Remove nodes with no edges
      rownames(conmat) <- colnames(conmat) #this needs to be same same if is.sym to work
      ifelse(isSymmetric.matrix(conmat)==TRUE,
             directed <- FALSE,
             directed <- TRUE)

      if(all.nodes == FALSE && directed == FALSE) {
        include.vec <- vector(length=dim(data)[1])
        for (i in 1:dim(conmat)[1]){
          ifelse(any(conmat[i, ] != 0), include.vec[i] <- 1, include.vec[i] <- 0)
        }
        data <- data[as.logical(include.vec), ,drop=F]
        conmat <- conmat[which(rowSums(conmat, na.rm = T) != 0), which(colSums(conmat, na.rm = T) != 0), drop = F]
      }

      if(all.nodes==FALSE && directed == TRUE) {
        include.vec <- vector(length=dim(data)[1])
        for (i in 1:dim(conmat)[1]){
          ifelse(any(conmat[i, ] != 0) | any(conmat[, i] != 0), include.vec[i] <- 1, include.vec[i] <- 0)
        }
      }

      if(all.nodes==TRUE) {
        include.vec <- vector(length=dim(data)[1])
        include.vec <- rep(1, length=dim(data)[1])
      }

      #in ortho view, only show legend for top view to avoid redundancy
      ifelse(v == 1, show.legend <- T, show.legend <- F)

      ortho_list[[v]] <- build_plot(conmat=conmat,
                                    data=data,
                                    background=background,
                                    node.size=node.size,
                                    view=current_view,
                                    node.color=node.color,
                                    thr=thr,
                                    uthr=uthr,
                                    edge.color=edge.color,
                                    edge.alpha=edge.alpha,
                                    edge.width=edge.width,
                                    scale.edge.width=scale.edge.width,
                                    show.legend=show.legend,
                                    labels=labels, label.size=label.size,
                                    include.vec=include.vec,
                                    edge.color.weighted=edge.color.weighted,
                                    label.edge.weight=label.edge.weight)
      if(is.environment(edge.color) == T) {
        ortho_list[[v]] <- ortho_list[[v]] + edge.color
      }
    }

    right_col <- plot_grid(ortho_list[[2]],
                           ortho_list[[3]],
                           nrow=2,
                           rel_heights = c(1.3, 1.6))  # Give much more space to both views
    p <- plot_grid(ortho_list[[1]], right_col, rel_widths = c(1.5, 1.5))  # Equal space allocation
    return(p)
  }

  # If not ortho, then do the below:
  if(background=='ICBM152') {
    background <- generate_background(view, background.alpha, 
                                    bg_xmin, bg_ymin, bg_xmax, bg_ymax, use_mesh_background)
  }

  if(background!='ICBM152') {
    m <- OpenImageR::readImage(background)
    w <- matrix(rgb(m[,,1],m[,,2],m[,,3], m[,,4] * background.alpha), nrow=dim(m)[1])
    background <- rasterGrob(w)
  }

  #if no conmat is provided, build nparc x  nparc empty one
  nparc <- dim(data)[1]
  if (!exists("conmat")){conmat <- matrix(0L, nrow=nparc, ncol=nparc)}

  #convert conmat to matrix
  conmat <- as.matrix(conmat)

  #Remove nodes with no edges
  rownames(conmat) <- colnames(conmat) #this needs to be same same if is.sym to work
  ifelse(isSymmetric.matrix(conmat)==TRUE,
         directed <- FALSE,
         directed <- TRUE)

  if(all.nodes == FALSE && directed == FALSE) {
    include.vec <- vector(length=dim(data)[1])
    for (i in 1:dim(conmat)[1]){
      ifelse(any(conmat[i, ] != 0), include.vec[i] <- 1, include.vec[i] <- 0)
    }
    data <- data[as.logical(include.vec), ,drop=F]
    conmat <- conmat[which(rowSums(conmat, na.rm = T) != 0), which(colSums(conmat, na.rm = T) != 0), drop = F]
  }

  if(all.nodes==FALSE && directed == TRUE) {
    include.vec <- vector(length=dim(data)[1])
    for (i in 1:dim(conmat)[1]){
      ifelse(any(conmat[i, ] != 0) | any(conmat[, i] != 0), include.vec[i] <- 1, include.vec[i] <- 0)
    }
  }

  if(all.nodes==TRUE) {
    include.vec <- vector(length=dim(data)[1])
    include.vec <- rep(1, length=dim(data)[1])
  }

  p <- build_plot(conmat=conmat,
                  data=data,
                  background=background,
                  node.size=node.size,
                  view=view,
                  node.color=node.color,
                  thr=thr,
                  uthr=uthr,
                  edge.color=edge.color,
                  edge.alpha=edge.alpha,
                  edge.width=edge.width,
                  scale.edge.width=scale.edge.width,
                  show.legend=show.legend,
                  labels=labels,
                  label.size=label.size,
                  include.vec=include.vec,
                  edge.color.weighted=edge.color.weighted,
                  label.edge.weight=label.edge.weight,
                  bg_xmax=bg_xmax,
                  bg_xmin=bg_xmin,
                  bg_ymax=bg_ymax,
                  bg_ymin=bg_ymin)

  return(p)
}
