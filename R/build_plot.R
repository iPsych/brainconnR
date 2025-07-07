#' workhorse function for \code{brainconn()}
#'
#' returns a ggraph object of plotted brain connectivity matrix
#' @author Sidhant Chopra
#' @import ggraph
#' @import ggplot2
#' @import grid

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
                       ...) {

  # Helper function to get view-specific coordinates and limits
  get_view_coords <- function(view, data, bg_xmin, bg_ymin, bg_xmax, bg_ymax) {
    coords <- list()
    
    if (view == "top") {
      coords$x.mni <- data$x.mni
      coords$y.mni <- data$y.mni
      coords$depth <- data$z.mni
      coords$xlim <- c(-75 + bg_xmin, 70 + bg_xmax)  # Node coordinate limits (unchanged)
      coords$ylim <- c(-107 + bg_ymin, 73 + bg_ymax)
      coords$bg_xlim <- c(-78 + bg_xmin, 73 + bg_xmax)  # Background limits (2% expanded)
      coords$bg_ylim <- c(-111 + bg_ymin, 77 + bg_ymax)
    } else if (view == "bottom") {
      coords$x.mni <- data$x.mni * -1
      coords$y.mni <- data$y.mni
      coords$depth <- data$z.mni * -1
      coords$xlim <- c(-70 + bg_xmin, 70 + bg_xmax)
      coords$ylim <- c(-107 + bg_ymin, 73 + bg_ymax)
      coords$bg_xlim <- c(-78 + bg_xmin, 73 + bg_xmax)
      coords$bg_ylim <- c(-111 + bg_ymin, 77 + bg_ymax)
    } else if (view == "front") {
      coords$x.mni <- data$x.mni
      coords$y.mni <- data$z.mni
      coords$depth <- data$y.mni
      coords$xlim <- c(-70 + bg_xmin, 70 + bg_xmax)
      coords$ylim <- c(-54 + bg_ymin, 80 + bg_ymax)  # Expanded bottom limit
      coords$bg_xlim <- c(-78 + bg_xmin, 78 + bg_xmax)  # 12% expansion
      coords$bg_ylim <- c(-56 + bg_ymin, 88 + bg_ymax)
    } else if (view == "back") {
      coords$x.mni <- data$x.mni * -1
      coords$y.mni <- data$z.mni
      coords$depth <- data$y.mni * -1
      coords$xlim <- c(-70 + bg_xmin, 70 + bg_xmax)
      coords$ylim <- c(-54 + bg_ymin, 80 + bg_ymax)  # Expanded bottom limit
      coords$bg_xlim <- c(-78 + bg_xmin, 78 + bg_xmax)  # 12% expansion
      coords$bg_ylim <- c(-56 + bg_ymin, 88 + bg_ymax)
    } else if (view == "left") {
      coords$x.mni <- data$y.mni * -1
      coords$y.mni <- data$z.mni
      coords$depth <- data$x.mni
      coords$xlim <- c(-64, 98)  # Use the adjusted limits directly
      coords$ylim <- c(-50, 80)  # Expanded bottom limit to prevent cropping
      coords$bg_xlim <- c(-85 + bg_xmin, 116 + bg_xmax)  # 15% expansion
      coords$bg_ylim <- c(-60 + bg_ymin, 87 + bg_ymax)
    } else if (view == "right") {
      coords$x.mni <- data$y.mni
      coords$y.mni <- data$z.mni
      coords$depth <- data$x.mni * -1
      coords$xlim <- c(-98, 64)  # Use the adjusted limits directly
      coords$ylim <- c(-50, 80)  # Expanded bottom limit to prevent cropping
      coords$bg_xlim <- c(-158 + bg_xmin, 121 + bg_xmax)  # 15% expansion
      coords$bg_ylim <- c(-60 + bg_ymin, 87 + bg_ymax)
    }
    
    return(coords)
  }

  # Get coordinates for the specified view
  view_coords <- get_view_coords(view, data, bg_xmin, bg_ymin, bg_xmax, bg_ymax)
  x.mni <- view_coords$x.mni
  y.mni <- view_coords$y.mni
  depth <- view_coords$depth
  xlim <- view_coords$xlim  # Node coordinate limits (unchanged)
  ylim <- view_coords$ylim
  bg_xlim <- view_coords$bg_xlim  # Background limits (expanded)
  bg_ylim <- view_coords$bg_ylim

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

  # Build the base plot with background
  p <- ggraph(layout, circular = FALSE) +
    annotation_custom(background, 
                     xmax = bg_xlim[2], xmin = bg_xlim[1], 
                     ymax = bg_ylim[2], ymin = bg_ylim[1]) +
    coord_fixed(xlim = xlim, ylim = ylim)

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
