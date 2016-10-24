#' Function to add color to facet strip
#'
#' @param d the data
#' @param p the base plot
#' @param statepalette character list where each element is a color and
#'                     the name of each element is the state number
plot_colorstrip <- function(d, p, statepalette){
  # Create new strips of color (to be updated)
  dummy <- p
  dummy$layers <- NULL
  dummy <- dummy +
    geom_rect(data=d, xmin=-Inf, ymin=-Inf, xmax=Inf, ymax=Inf,
              aes(fill = as.character(state))) +
    theme_minimal() +
    scale_fill_manual(values = statepalette, breaks=NULL)

  library(gtable)

  g1 <- ggplotGrob(p)
  g2 <- ggplotGrob(dummy)

  gtable_select <- function (x, ...)
  {
    matches <- c(...)
    x$layout <- x$layout[matches, , drop = FALSE]
    x$grobs <- x$grobs[matches]
    x
  }

  # Update plot grobs/new strips
  panels <- grepl(pattern="panel", g2$layout$name)
  strips <- grepl(pattern="strip-top|strip_t", g2$layout$name)
  g2$layout$t[panels] <- g2$layout$t[panels] - 1
  g2$layout$b[panels] <- g2$layout$b[panels] - 1
  new_strips <- gtable_select(g2, panels | strips)
  if (!grepl("panel", new_strips$layout$name[1])) {
    new_strips$layout$z <- rev(new_strips$layout$z)
  }
  #grid::grid.newpage()
  #grid::grid.draw(new_strips)

  gtable_stack <- function(g1, g2){
    g1$grobs <- c(g1$grobs, g2$grobs)
    g1$layout <- transform(g1$layout, z= z-max(z), name="g2")
    g1$layout <- rbind(g1$layout, g2$layout)
    g1
  }
  ## ideally you'd remove the old strips, for now they're just covered
  new_plot <- gtable_stack(g1, new_strips)
  grid::grid.newpage()
  grid::grid.draw(new_plot)
}
