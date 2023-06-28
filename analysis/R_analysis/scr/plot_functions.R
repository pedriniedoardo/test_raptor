## Useful Plot functions

#### Prettify functions ####
# Make a color transparent
transp <- function(col, a=.5){
  colr <- col2rgb(col)
  return(rgb(colr[1,], colr[2,], colr[3,], a*255, maxColorValue = 255))
}

# Generates 2 (edge) ticks on given axes of current plot
twoTicks <- function(side = c(1,2), col = NA, col.ticks = 1, las = c(1,2), log=NULL, ...){
  # col = NA & col.ticks = 1 makes the axis line dissapear, but keeps ticks
  for(i in seq_along(side)){
    axt <- axTicks(side[i], log=log)
    axis(side[i], at = axt[c(1, length(axt))], col = col, col.ticks = col.ticks, las = las[i], ...)
  }
}


#### Custom plot functions ####
# Plot pca components
plot_pca <- function(form, data, pca, nc = ncol(pca$x), legend = 1, 
                     lwd = 2, pch = 1, col = NULL, main = NULL, 
                     ylab = "PC", xlab = NULL, nc.only = F,
                     col.pal = col.palette, 
                     expr_per_comp = "", expr_legend = "", ...){
  if(length(nc) == 1 & !nc.only)
    nc <- seq_len(nc)
  
  form <- update(form, NULL~.)
  evars <- all.vars(form)
  if(!all(evars%in%colnames(data)))
    stop("All terms of formula must be in data.")
  
  x <- data[,evars[1]]
  if(!is.numeric(x))
    stop("First term of formula must be numeric variable (x-axis).")
  
  
  expr_per_comp <- as.character(expr_per_comp)
  
  xs <- lapply(evars[-1], function(xn) factor(data[, xn]))
  
  if(missing(col)){
    col <- rep(col.palette[1], length(x))
    if(length(evars) > 1){
      col <- col.pal[as.numeric(xs[[1]])]
    }
  }
  if(missing(pch)){
    if(length(evars) > 2){
      pch <- as.numeric(xs[[2]])
    }
  }
  
  if(missing(xlab))
    xlab <- evars[1]
  if(missing(main)) 
    main <- paste0(ylab, nc) 
  else main <- paste0(main, " (", ylab, nc,")")
  
  invisible(lapply(seq_along(nc), function(i){
    plot(x, pca$x[,nc[i]], 
         ylab = ylab, xlab = xlab, 
         col = col, pch = pch, lwd = lwd, 
         main = main[i], ...)
    eval(expr = parse(text = expr_per_comp))
    if(i %in% legend)
      eval(expr = parse(text = expr_legend))
    
  }))
}

# plot dots with mean bar
plot_dotmean <- function(x, y, ats=1:length(levels(x)), cols=1:length(levels(x)), 
                         ml=.5, mlwd=2, 
                         method= "center", ...){
  beeswarm::beeswarm(y~x, at=ats, col = cols, method = method, ...)
  l <- lm(y~x)
  segments(x0 = ats-(ml/2), x1= ats+(ml/2), y0=c(l$coefficients[1], l$coefficients[1] + l$coefficients[-1]), col= cols,
           lwd=mlwd)
}

# 2D density plot
gg_hex2d <- function(x, y=NULL, rg = range(xy, na.rm = T)*1.02,
                     l.breaks = log1p(c(0, 1, 5, 10, 50, 100, Inf)),
                     l.labels = c('1', '5', '10', '50', '100', '100+'),
                     xlab = "log2(FC) in x",
                     ylab = "log2(FC) in y",
                     main = "", get.r =F, add.vd = F,
                     DEgsel = NULL,
                     ...){
  require(ggplot2)
  if(is.null(y))
    xy <- x
  else
    xy <- cbind(x, y)
  xy <- as.data.frame(xy)
  
  colnames(xy) <- c("x", "y") 
  
  g <- ggplot(data = xy, mapping = aes(x=x, y=y)) + 
    geom_hex(aes(fill = stat(cut(log(count), breaks = l.breaks, 
                                 labels = F, right = T, include.lowest = T))), bins=100) +
    scale_fill_gradientn(colors = viridisLite::inferno(length(l.breaks)), name = 'count', labels = l.labels) + 
    theme_classic() + xlim(rg) + ylim(rg) +
    xlab(xlab) + ylab(ylab) + ggtitle(main) + coord_fixed()
  if(get.r){
    if(any(is.na(xy))){
      message("Warning: removed NA values to compute correlation coefs")
      rmsel <- apply(xy, 1, function(r) any(is.na(r)))
      xy <- xy[!rmsel,]
      if(!is.null(DEgsel))
        DEgsel <- DEgsel[!rmsel]
    }
    
    cc <- cor(xy)[1,2]
    cctxt <- paste0("r = ", round(cc, 3))
    if(add.vd)
      cctxt <- paste0(cctxt, "\n", round(100*(cc^2), 1), "% VarDev")
    if(!is.null(DEgsel)){
      ccDE <- cor(xy[DEgsel,])[1,2]
      cctxt <- paste0(cctxt, "\nr(DE) = ", round(ccDE, 3))
      if(add.vd)
        cctxt <- paste0(cctxt, "\n", round(100*(ccDE^2), 1), "% VarDev(DE)")
    }
    
    g <- g + annotate(geom="text", hjust=1, vjust=0,
                      x = rg[2],
                      y = rg[1],
                      label = cctxt)
  }
  return(g)
}



#### Figure saving functions ####

# Wrapper function for png figure generation
png_custom <- function(figname, path = "figs/", 
                       fig.width = 7, fig.height = 5, res = 150, ...){
  png(filename = paste0(path, figname, ".png"), 
      width = fig.width, height = fig.height, res = res, units = "in")
}

# Wrapper function for pdf figure generation
pdf_custom <- function(figname, path = "figs/", 
                       fig.width = 7, fig.height = 5, ...){
  pdf(file = paste0(path, figname, ".pdf"), 
      width = fig.width, height = fig.height, ...)
}

# Wrapper of wrappers
fig_custom <- function(figname, output = c("png", "pdf"), path = "figs/", fig.width = 7, fig.height = 5, ...){
  output <- match.arg(output)
  if("png" == output)
    png_custom(figname = figname, path = path, fig.width = fig.width, fig.height = fig.height, ...)
  if("pdf" == output)
    pdf_custom(figname = figname, path = path, fig.width = fig.width, fig.height = fig.height, ...)
}


