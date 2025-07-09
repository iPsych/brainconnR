build_plot <- function(conmat,
                       data,
                       data.row=NULL,
                       data.col=NULL,
                       background,
                       node.size,
                       node.color="network",
                       thr=NULL,
                       uthr=NULL,
                       view,
                       edge.color,
                       edge.alpha,
                       edge.width,
                       show.legend,
                       label.size,
                       labels,
                       include.vec=NULL,
                       scale.edge.width,
                       edge.color.weighted,
                       label.edge.weight,
                       bg_xmin=0,
                       bg_ymin=0,
                       bg_xmax=0,
                       bg_ymax=0,
                       d.factor=1.05,
                       ...) {

  # Helper function to project 3D coordinates to 2D based on view
  # Using the same transformations as the original working code
  project_3d_to_2d <- function(coords_3d, view) {
    x <- coords_3d[,1]  # x.mni
    y <- coords_3d[,2]  # y.mni  
    z <- coords_3d[,3]  # z.mni
    
    if (view == "top") {
      proj_x <- x
      proj_y <- y
      depth <- z
    } else if (view == "bottom") {
      proj_x <- x * -1
      proj_y <- y
      depth <- z * -1
    } else if (view == "front") {
      proj_x <- x
      proj_y <- z
      depth <- y
    } else if (view == "back") {
      proj_x <- x * -1
      proj_y <- z
      depth <- y * -1
    } else if (view == "left") {
      # Original: x.mni <- data$y.mni*-1, y.mni <- data$z.mni, depth <- data$x.mni
      proj_x <- y * -1
      proj_y <- z
      depth <- x
    } else if (view == "right") {
      # Original: x.mni <- data$y.mni, y.mni <- data$z.mni, depth <- data$x.mni*-1
      proj_x <- y
      proj_y <- z
      depth <- x * -1
    }
    
    return(data.frame(x = proj_x, y = proj_y, depth = depth))
  }

  # Get actual 3D mesh coordinates to calculate accurate bounds
  calculate_view_bounds <- function(view, margin_factor = 0.05) {
    # Load mesh vertices
    vb <- get("ICBM152_mesh_vb")
    
    # Create 3D coordinate matrix (transpose to get n x 3)
    mesh_coords_3d <- t(vb[1:3, ])  # Only x, y, z coordinates
    
    # Project mesh coordinates to 2D
    mesh_2d <- project_3d_to_2d(mesh_coords_3d, view)
    
    # Calculate bounds with margin
    x_range <- range(mesh_2d$x, na.rm = TRUE)
    y_range <- range(mesh_2d$y, na.rm = TRUE)
    
    x_margin <- diff(x_range) * margin_factor
    y_margin <- diff(y_range) * margin_factor
    
    bounds <- list(
      x_min = x_range[1] - x_margin,
      x_max = x_range[2] + x_margin,
      y_min = y_range[1] - y_margin,
      y_max = y_range[2] + y_margin
    )
    
    return(bounds)
  }

  # Calculate accurate bounds for this view
  bounds <- calculate_view_bounds(view)
  
  # Apply user-specified background adjustments
  bg_bounds <- list(
    x_min = bounds$x_min + bg_xmin,
    x_max = bounds$x_max + bg_xmax,
    y_min = bounds$y_min + bg_ymin,
    y_max = bounds$y_max + bg_ymax
  )

  # Project node coordinates to 2D
  node_coords_3d <- data.frame(
    x.mni = data$x.mni * d.factor,  # Apply distance factor like in 3D
    y.mni = data$y.mni * d.factor,
    z.mni = data$z.mni * d.factor
  )
  
  node_2d <- project_3d_to_2d(node_coords_3d, view)
  x.mni <- node_2d$x
  y.mni <- node_2d$y
  depth <- node_2d$depth

  if (!exists("conmat")) stop(print("Please enter a valid connectivity matrix"))

  # Determine graph properties
  directed <- !isSymmetric.matrix(conmat)
  weighted <- !all(conmat %in% c(0,1))

  # Create layout
  if (!directed) {
    conmat[upper.tri(conmat)] <- 0 # only take bottom tri to avoid duplicate edges
  }
  
  layout <- create_layout(graph = conmat, layout = "stress", circular = TRUE)
  layout$x <- x.mni
  layout$y <- y.mni
  
  if (directed) {
    layout$facet <- include.vec
  }

  # Build the base plot with background using accurate bounds
  p <- ggraph(layout, circular = FALSE) +
    annotation_custom(background, 
                     xmax = bg_bounds$x_max, xmin = bg_bounds$x_min, 
                     ymax = bg_bounds$y_max, ymin = bg_bounds$y_min) +
    coord_fixed(xlim = c(bounds$x_min, bounds$x_max), 
                ylim = c(bounds$y_min, bounds$y_max))

  # Add edges based on graph type
  if (directed && !weighted) {
    p <- p + geom_edge_parallel(color = edge.color,
                               edge_width = edge.width,
                               edge_alpha = edge.alpha,
                               arrow = arrow(length = unit(3, 'mm')),
                               end_cap = circle((node.size/2) + 0.6, 'mm'))
  } else if (directed && weighted) {
    if (!edge.color.weighted && !label.edge.weight) {
      p <- p + geom_edge_parallel(aes(width = weight),
                                 color = edge.color,
                                 edge_alpha = edge.alpha,
                                 arrow = arrow(length = unit(3, 'mm')),
                                 end_cap = circle(node.size/2, 'mm')) +
        geom_edge_loop0(aes(strength = node.size * 3, width = weight),
                       color = edge.color,
                       edge_alpha = edge.alpha,
                       arrow = arrow(length = unit(3, 'mm')))
    } else if (!edge.color.weighted && label.edge.weight) {
      p <- p + geom_edge_parallel(aes(width = weight, label = round(weight, 3)),
                                 color = edge.color,
                                 edge_alpha = edge.alpha,
                                 arrow = arrow(length = unit(3, 'mm')),
                                 end_cap = circle(node.size/2, 'mm'),
                                 angle_calc = 'along',
                                 alpha = 0,
                                 label_dodge = unit(2.5, 'mm'),
                                 label_size = 2,
                                 fontface = "bold") +
        geom_edge_loop0(aes(strength = node.size * 3, width = weight, label = round(weight, 3)),
                       color = edge.color,
                       edge_alpha = edge.alpha,
                       arrow = arrow(length = unit(3, 'mm')),
                       angle_calc = 'none',
                       alpha = 0,
                       label_dodge = unit(6, 'mm'),
                       label_size = 2,
                       vjust = -1,
                       fontface = "bold")
    } else if (edge.color.weighted && !label.edge.weight) {
      p <- p + geom_edge_parallel(aes(color = weight),
                                 edge_alpha = edge.alpha,
                                 edge_width = edge.width,
                                 arrow = arrow(length = unit(3, 'mm')),
                                 end_cap = circle(node.size/2, 'mm')) +
        geom_edge_loop(aes(strength = node.size * 3, color = weight),
                      edge_width = edge.width,
                      edge_alpha = edge.alpha,
                      arrow = arrow(length = unit(3, 'mm')))
    } else if (edge.color.weighted && label.edge.weight) {
      p <- p + geom_edge_parallel(aes(color = weight, label = round(weight, 3)),
                                 edge_alpha = edge.alpha,
                                 edge_width = edge.width,
                                 arrow = arrow(length = unit(3, 'mm')),
                                 end_cap = circle(node.size/2, 'mm'),
                                 angle_calc = 'along',
                                 alpha = 0,
                                 label_dodge = unit(2.5, 'mm'),
                                 label_size = 2,
                                 fontface = "bold") +
        geom_edge_loop(aes(strength = node.size * 3, color = weight, label = round(weight, 3)),
                      edge_width = edge.width,
                      edge_alpha = edge.alpha,
                      arrow = arrow(length = unit(3, 'mm')),
                      angle_calc = 'none',
                      alpha = 0,
                      label_dodge = unit(6, 'mm'),
                      label_size = 2,
                      vjust = -1,
                      fontface = "bold")
    }
  } else if (!directed && !weighted) {
    p <- p + geom_edge_link(color = edge.color,
                           edge_width = edge.width,
                           edge_alpha = edge.alpha)
  } else if (!directed && weighted) {
    if (!edge.color.weighted && !label.edge.weight) {
      p <- p + geom_edge_link(aes(width = weight),
                             color = edge.color,
                             edge_alpha = edge.alpha)
    } else if (!edge.color.weighted && label.edge.weight) {
      p <- p + geom_edge_link(aes(width = weight, label = round(weight, 3)),
                             color = edge.color,
                             edge_alpha = edge.alpha,
                             angle_calc = 'along',
                             alpha = 0,
                             label_dodge = unit(2.5, 'mm'),
                             label_size = 2,
                             fontface = "bold")
    } else if (edge.color.weighted && !label.edge.weight) {
      p <- p + geom_edge_link(aes(colour = weight),
                             edge_width = edge.width,
                             edge_alpha = edge.alpha)
    } else if (edge.color.weighted && label.edge.weight) {
      p <- p + geom_edge_link(aes(colour = weight, label = round(weight, 3)),
                             edge_width = edge.width,
                             edge_alpha = edge.alpha,
                             angle_calc = 'along',
                             alpha = 0,
                             label_dodge = unit(2.5, 'mm'),
                             label_size = 2,
                             fontface = "bold")
    }
  }

  # Scale edge width for weighted networks
  if (weighted && !is.null(scale.edge.width)) {
    p <- p + scale_edge_width(range = scale.edge.width)
  }

  # Add nodes
  if (directed) {
    if (node.color == "network") {
      p <- p + geom_node_point(size = node.size, 
                              aes(colour = as.factor(data$network), 
                                  filter = as.logical(facet)))
    } else {
      p <- p + geom_node_point(size = node.size, colour = node.color,
                              aes(filter = as.logical(facet)))
    }
  } else {
    if (node.color == "network") {
      p <- p + geom_node_point(size = node.size, 
                              aes(colour = as.factor(data$network)))
    } else {
      p <- p + geom_node_point(size = node.size, colour = node.color)
    }
  }

  # Add labels
  if (labels) {
    if (directed) {
      p <- p + geom_node_text(aes(label = data$ROI.Name, 
                                 filter = as.logical(facet)),
                             size = label.size, repel = TRUE,
                             nudge_x = node.size + 2, 
                             nudge_y = node.size)
    } else {
      p <- p + geom_node_text(aes(label = data$ROI.Name),
                             size = label.size, repel = TRUE,
                             nudge_x = node.size + 2, 
                             nudge_y = node.size)
    }
  }

  # Apply theme
  p <- p + theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank())

  # Handle legend
  if (!show.legend) {
    p <- p + theme(legend.position = "none")
  } else {
    p <- p + scale_color_discrete(name = "Network")
  }

  return(p)
}