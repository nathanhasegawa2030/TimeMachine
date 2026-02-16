###plotting.R

#Written by Nathan Hasegawa

#Functions to generate plots

angle_plot <- function(circData, ccol = "#c9cfd4", pcol = "#007fff", 
                       title=NULL, loczero=0, mode="none",
                       edgelabs = c("0","","","")) {
  #Generates an angle plot of every angle in circData on the unit circle.
  #circData: an object of type circular that includes the angles to be plotted.
  #Must be in RADIANS.
  #ccol: optional argument for the color of the circle.
  #pcol: optional argument for the color of the points.
  #title: optional argument for the plot title.
  #loczero: optional argument for the location of the zero angle. By default, it
  #is at an angle of 0 (i.e. the right). If loczero isn't 0, the location of the
  #zero value, and the line accompanying it, is rotated counterclockwise by
  #loczero degrees.
  #mode: optional argument for the mode of the display. If mode = "clock", then
  #loczero is pi/2.
  #edgelabs: optional argument for what should be written at the top, right,
  #bottom, and left of the circle. Do not use with "mode."
  if (mode == "clock") {
    loczero = (pi/2)
    edgelabs = c("12 AM", "6 AM", "12 PM", "6 PM")
  }
  plotx <- cos(circData + loczero)
  ploty <- sin(circData + loczero)
  ggplot() +
    geom_segment(aes(x = 0, y = 0, xend = cos(loczero), yend = sin(loczero)), 
                 color = "black", size = 2.5, alpha=1) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = NA, color = ccol,
                linewidth=4, alpha=0.3) +
    geom_point(aes(x = plotx, y = ploty), color = pcol, size=3, alpha=0.2) +
    xlim(-1.1,1.1) +
    ylim(-1.1,1.1) +
    labs(title = title) +
    theme_void() +
    theme(plot.title=element_text(hjust = 0.5,family="Trebuchet MS", size = 22,
                                  margin = margin(b = -15, unit = "pt")),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          aspect.ratio=1,
          plot.margin = unit(c(.2,.2,.2,.2), "inches"))
}

donut_plot <- function(estData, obsData, ccol="#c9cfd4", 
                       pcol="#007fff", ncol="#f40000", title=NULL) {
  #Generates a "donut plot" (see Jha & Biswas, 2017) that plots both the
  #residuals and estimated data. In this way, we can view predictive power and
  #accuracy by region.
  #estData: the data estimated by a model.
  #obsData: the actual data.
  #pcol: optional argument for the color used when the residuals are POSITIVE.
  #That is, the true angle is greater than the estimated angle (resid > 0).
  #ncol: optional argument for the color used when the residuals are NEGATIVE.
  #title: optional argument for the plot title.
  
  #Reorder the rows at random. This is done for viewability.
  order <- sample(x=1:length(estData), replace=FALSE)
  estData <- estData[order]
  obsData <- obsData[order]
  
  resid <- sapply(obsData - estData, rescale_angle)
  mag <- 2 + cos(resid)
  plotx <- mag * cos(estData)
  ploty <- mag * sin(estData)
  
  ggplot() +
    geom_circle(aes(x0 = 0, y0 = 0, r = 3), fill = NA, color = ccol,
                linewidth=4, alpha=0.3) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 1), fill = NA, color = ccol,
                linewidth=4, alpha=0.3) +
    geom_circle(aes(x0 = 0, y0 = 0, r = 2), fill = NA, color = ccol,
                linewidth=4, alpha=0.01) +
    geom_segment(aes(x = 0, y = 0, xend = 3, yend = 0), 
                 color = "black", size = 2.5, alpha=1) +
    geom_point(aes(x = plotx, y = ploty, 
                   color = ifelse(resid > 0, "Positive", "Negative")), 
                   size=3, alpha=0.1) +
    scale_color_manual(values = c("Positive"=pcol,"Negative"=ncol),
                       guide = "none") +
    xlim(-3.1,3.1) +
    ylim(-3.1,3.1) +
    labs(title = title) +
    theme_void() +
    theme(plot.title=element_text(hjust = 0.5,family="Trebuchet MS", size = 22,
                                  margin = margin(b = -15, unit = "pt")),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          aspect.ratio=1,
          plot.margin = unit(c(.2,.2,.2,.2), "inches"))
}

#This helper function was written by Google Gemini.
rescale_angle <- function(angle) {
  angle <- angle %% (2 * pi) # Wrap to [0, 2*pi)
  angle[angle >= pi] <- angle[angle >= pi] - (2 * pi) # Shift angles > pi to the negative range
  angle
}