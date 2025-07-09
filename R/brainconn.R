brainconn <- function(atlas,
                      background='ICBM152',
                      view ="top",
                      conmat=NULL,
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
                      background.alpha = 0.5,
                      bg_xmax=0,
                      bg_xmin=0,
                      bg_ymax=0,
                      bg_ymin=0,
                      d.factor=1.05,  # Add distance factor like in 3D
                      use_mesh_background = FALSE) {

  # Helper function to create rotation matrix from theta and phi
  create_rotation_matrix <- function(theta_deg, phi_deg) {
    # Convert degrees to radians
    theta <- theta_deg * pi / 180
    phi <- phi_deg * pi / 180
    
    # Rotation around X-axis (phi)
    Rx <- matrix(c(1,        0,         0,
                   0, cos(phi), -sin(phi),
                   0, sin(phi),  cos(phi)), nrow = 3, byrow = TRUE)
    
    # Rotation around Z-axis (theta)  
    Rz <- matrix(c(cos(theta), -sin(theta), 0,
                   sin(theta),  cos(theta), 0,
                   0,           0,          1), nrow = 3, byrow = TRUE)
    
    # Combined rotation: apply phi first, then theta
    return(Rz %*% Rx)
  }

  # Helper function to project 3D coordinates to 2D using proper rotation matrices
  project_3d_to_2d <- function(coords_3d, view) {
    # Get view parameters (same as used in background generation)
    view_params <- switch(view,
      "top"    = list(theta = 0,   phi = 0),
      "bottom" = list(theta = 180, phi = -180),
      "left"   = list(theta = 90,  phi = -90),
      "right"  = list(theta = -90, phi = -90),
      "front"  = list(theta = -180, phi = -90),
      "back"   = list(theta = 0,   phi = -90)
    )
    
    # Create rotation matrix
    R <- create_rotation_matrix(view_params$theta, view_params$phi)
    
    # Apply rotation to each point
    coords_matrix <- as.matrix(coords_3d)  # Ensure it's a matrix
    rotated_coords <- t(R %*% t(coords_matrix))  # Apply rotation
    
    # For 2D projection, we typically use the first two dimensions after rotation
    # The third dimension becomes depth
    proj_x <- rotated_coords[, 1]
    proj_y <- rotated_coords[, 2] 
    depth <- rotated_coords[, 3]
    
    return(data.frame(x = proj_x, y = proj_y, depth = depth))
  }

  # Helper function to generate accurate 2D background from 3D mesh
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
    
    # Load 3D mesh data
    vb <- get("ICBM152_mesh_vb")
    it <- get("ICBM152_mesh_it")
    
    # Calculate accurate view bounds from mesh
    mesh_coords_3d <- t(vb[1:3, ])  # Only x, y, z coordinates
    mesh_2d <- project_3d_to_2d(mesh_coords_3d, view)
    
    # Calculate bounds with small margin for background
    x_range <- range(mesh_2d$x, na.rm = TRUE)
    y_range <- range(mesh_2d$y, na.rm = TRUE)
    margin_factor <- 0.05
    
    x_margin <- diff(x_range) * margin_factor
    y_margin <- diff(y_range) * margin_factor
    
    xmax <- x_range[2] + x_margin + bg_xmax
    xmin <- x_range[1] - x_margin + bg_xmin
    ymax <- y_range[2] + y_margin + bg_ymax
    ymin <- y_range[1] - y_margin + bg_ymin
    
    # Generate brain background using rgl approach
    if (!requireNamespace("rgl", quietly = TRUE)) {
      stop("rgl package is required for mesh background generation")
    }
    library(rgl)
    
    # Create mesh3d object for rgl
    vertices <- t(vb)  # transpose to get n x 3 matrix for rgl
    indices <- t(it)   # transpose triangle indices
    
    # Create mesh3d object
    brain_mesh <- tmesh3d(
      vertices = t(vertices),  # rgl expects 3 x n format
      indices = t(indices),    # rgl expects 3 x m format
      homogeneous = FALSE
    )
    
    # Helper function to render and capture views
    capture_view <- function(mesh, theta=0, phi=0, zoom=1, width=400, height=400, 
                             bg_color="white", mesh_color="gray85", alpha=0.3) {
      
      # Open 3D device
      open3d(windowRect = c(0, 0, width, height))
      
      # Set background
      bg3d(color = bg_color)
      
      # Render the mesh
      shade3d(mesh, color = mesh_color, alpha = alpha, smooth = TRUE)
      
      # Set camera position and orientation
      user_matrix <- rotationMatrix(phi * pi/180, 1, 0, 0) %*%
        rotationMatrix(theta * pi/180, 0, 0, 1)
      
      # Apply view transformation
      view3d(userMatrix = user_matrix, zoom = zoom)
      
      # Capture the image
      tmpf <- tempfile(fileext = ".png")
      snapshot3d(tmpf, fmt = "png", width = width, height = height)
      
      # Read the PNG back as array [h x w x 4] (RGBA)
      if (!requireNamespace("png", quietly = TRUE)) {
        stop("png package is required for mesh background generation")
      }
      library(png)
      img <- readPNG(tmpf)
      
      # Clean up
      rgl.close()
      unlink(tmpf)  # Remove temporary file
      
      return(img)
    }
    
    # Get view parameters (same as used in coordinate projection)
    view_params <- switch(view,
      "top"    = list(theta = 0,   phi = 0,   zoom = 0.48),
      "bottom" = list(theta = 180, phi = -180, zoom = 0.49),
      "left"   = list(theta = 90,  phi = -90,  zoom = 0.56),
      "right"  = list(theta = -90, phi = -90,  zoom = 0.56),
      "front"  = list(theta = -180, phi = -90, zoom = 0.50),
      "back"   = list(theta = 0,   phi = -90,  zoom = 0.50)
    )
    
    # Calculate appropriate dimensions based on view limits
    width_range <- xmax - xmin
    height_range <- ymax - ymin
    aspect_ratio <- width_range / height_range
    
    # Base resolution
    base_height <- 400
    render_width <- round(base_height * aspect_ratio)
    render_height <- base_height
    
    # Capture the view using rgl
    brain_img <- capture_view(
      mesh = brain_mesh,
      theta = view_params$theta,
      phi = view_params$phi,
      zoom = view_params$zoom,
      width = render_width,
      height = render_height,
      mesh_color = "gray85",
      alpha = background.alpha
    )
    
    # Ensure RGBA format
    if (length(dim(brain_img)) == 3 && dim(brain_img)[3] == 3) {
      # only RGB present, add alpha channel
      library(abind)
      alpha_chan <- array(background.alpha, dim = c(dim(brain_img)[1], dim(brain_img)[2], 1))
      brain_img <- abind(brain_img, alpha_chan, along = 3)
    }
    
    # Convert to rasterGrob
    w <- matrix(rgb(brain_img[,,1], brain_img[,,2], brain_img[,,3], 
                    brain_img[,,4]), nrow = dim(brain_img)[1])
    return(rasterGrob(w))
  }

  ifelse(is.character(atlas), data <- get(atlas), data <- atlas)

  #set background (add ability to add custom background image)
  if(background != "ICBM152"  && view == "ortho") {
    stop("Custom background image detected, view cannot be 'ortho', please select top,
        bottom, left, right, front or back.")}

  if (!is.null(thr)) {conmat[conmat < thr] <- 0} #lower threshold graph
  if (!is.null(uthr)) {conmat[conmat > thr] <- 0} #upper threshold graph
  
  #loop three times for the three views that make ortho view
  if (view == "ortho") {
    ortho_list <- list()
    ortho_views  <- c("top", "left", "front")
    for (v in 1:3) {
      current_view <- ortho_views[v]
      
      # Generate background from 3D mesh with accurate coordinates
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
                                    label.edge.weight=label.edge.weight,
                                    d.factor=d.factor)  # Pass d.factor
      if(is.environment(edge.color) == T) {
        ortho_list[[v]] <- ortho_list[[v]] + edge.color
      }
    }

    right_col <- plot_grid(ortho_list[[2]],
                           ortho_list[[3]],
                           nrow=2,
                           rel_heights = c(1.6, 1.5))  # Give much more space to both views
    p <- plot_grid(ortho_list[[1]], right_col, rel_widths = c(1.8, 1.8))  # Equal space allocation
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
                  bg_ymin=bg_ymin,
                  d.factor=d.factor)  # Pass d.factor

  return(p)
}