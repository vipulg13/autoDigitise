df1 <- data.frame(id=letters[1:10], x=1:10, y=rep(5,10),mou=paste("point",letters[1:10]),
                  link=file.path(tempdir(),paste0(LETTERS[1:10],".html")), stringsAsFactors=FALSE)
## Typically one wants to get pixel-coordinates for plots written to file.
## Here we'll use R's tempdir, later you may want to choose other locations
pngFile <- file.path(tempdir(),"test01.png")
png(pngFile, width=800, height=600, res=72)
## here we'll just plot a set of horiontal points at default parameters ...
plot(df1[,2:3], las=1)
dev.off()
## Note: Special characters should be converted for proper display in html during mouse-over
library(wrMisc)
df1$mou <- htmlSpecCharConv(df1$mou)
## Let's add the x- and y-coordiates of the points in pixels to the data.frame
df1 <- cbind(df1,convertPlotCoordPix(x=df1[,2], y=df1[,3], plotD=c(800,600),plotRes=72))
head(df1)

#userMar = (b,l,t,r)

img <- EBImage::readImage(pngFile)
EBImage::display(img)

df1 <- data.frame(id=letters[1:10], x=1:10, y=rep(5,10))
## Typically one wants to get pixel-coordinates for plots written to file.
## Here we'll use R's tempdir, later you may want to choose other locations
pngFile <- file.path(tempdir(),"test01.png")
png(pngFile, width=800, height=600, res=72)
## here we'll just plot a set of horiontal points at default parameters ...
plot(df1[,2:3], las=1, main="test01")
dev.off()
library(wrMisc)
df1$mou <- htmlSpecCharConv(df1$mou)
## Let's add the x- and y-coordiates of the points in pixels to the data.frame
df1 <- cbind(df1,convertPlotCoordPix(x=df1[,2], y=df1[,3], useMar = c(6,4.3,5,2), plotD=c(800,600),plotRes=72))
head(df1)

imgr <- imager::load.image(pngFile)
draw_circle(imgr,85,293,1,"darkgreen") %>% plot

## using mouseOverHtmlFile() one could now make an html document with interactive
## display of names and clockable links to the coordinates we determined here ...



library(data.table)
library(wrGraph)
library(ggplot2)
library(EBImage)
library(egg)

pngFile <- file.path(tempdir(),"test01.png")
png(pngFile, width=640, height=640, res=100)
rand_size <- sample(3:4, 1, replace = TRUE)
rand_theme <- sample(c("theme_minimal", "theme_gray", "theme_bw", 
                       "theme_light", "theme_dark", "theme_classic",
                       "theme_test", "theme_void", "theme_linedraw",
                       "theme_grey"), 1)
rand_sample_size <- sample(1:15, 1)
rand_shapes_total <- sample(1:6, 1)
rand_shapes <- sample(1:25, rand_shapes_total)
sample_shapes <- function(x, ...) x[sample(length(x), ...)]

rand_x_scale <- sample(c("identity", "log10"), 1)
rand_y_scale <- sample(c("identity", "log10"), 1)
getSampleData <- function(range_type) {
  if (rand_scale_range_type == "tens") {
    rand_min_scale <- sample(seq(10,50,5), 1)
    rand_max_scale <- sample(seq(rand_min_scale+25,90,5), 1)
    rand_range <- seq(rand_min_scale, rand_max_scale, by = 1)
  }
  if (rand_scale_range_type == "hundreds") {
    rand_min_scale <- sample(seq(0,500,50), 1)
    rand_max_scale <- sample(seq(rand_min_scale+250,900,50), 1)
    rand_range <- seq(rand_min_scale, rand_max_scale, by = 1)
  }
  if (rand_scale_range_type == "miniscule") {
    rand_min_scale <- sample(0.1 * 10^(seq(-12,0,1)), 1)
    rand_max_scale <- sample(10 * 10^(seq(log10(rand_min_scale),log10(rand_min_scale)+5,1)), 1)
    rand_range <- 10^(seq(log10(rand_min_scale),log10(rand_max_scale),0.01))
  }
  return(rand_range)
}

rand_scale_range_type <- sample(c("tens", "hundreds", "miniscule"), 1)
rand_x_data <- getSampleData(rand_scale_range_type)
rand_scale_range_type <- sample(c("tens", "hundreds", "miniscule"), 1)
rand_y_data <-  getSampleData(rand_scale_range_type)
rand_label_size <- sample(5:15,1)

dt <- data.table(x=sample(rand_x_data, rand_sample_size, replace=TRUE), 
                 y=sample(rand_y_data, rand_sample_size, replace=TRUE), 
                 shapes=factor(sample_shapes(rand_shapes, rand_sample_size, replace=TRUE)))

legendPosVec <- c("none", "top", "bottom", "right", "c(0,1)", "c(0,0)", "c(1,0)", "c(1,1)")
if (length(unique(dt[[3]])) == 1) {
  rand_legend_pos <- "none"
} else {
  rand_legend_pos <- sample(legendPosVec[head(2):length(legendPosVec)], 1)
}
legendPosIndex <- which(legendPosVec %in% rand_legend_pos)

#plot(df2[,1:2], las=1, main="test01")

p <- ggplot(dt, aes(x, y)) +
  geom_point(aes(shape = shapes), size = rand_size) +
  get(rand_theme)()

if(legendPosIndex > 4) {
  p = p + theme(legend.position = eval(parse(text=rand_legend_pos)),  
                legend.justification = eval(parse(text=rand_legend_pos)),
                #plot.margin = margin(1, 2, 1, 1, "cm"),
                axis.title=element_text(size=rand_label_size)) # (t,r,b,l)
} else {
  p = p + theme(legend.position = rand_legend_pos,
                #plot.margin = margin(1, 2, 1, 1, "cm"),
                axis.title=element_text(size=rand_label_size)) 
}
p <- p + scale_x_continuous(trans = rand_x_scale, expand = c(0, 0)) + 
  scale_y_continuous(trans = rand_y_scale, expand = c(0, 0))

p <- p + scale_x_continuous(trans = rand_x_scale, expand = c(0, 0)) + 
  scale_y_continuous(trans = rand_y_scale, expand = c(0, 0)) +
  coord_fixed()
p
#p
grid::grid.draw(egg::set_panel_size(p, height = grid::unit(4, "in"), width = grid::unit(5, "in")))
dev.off()



#dt1 <- cbind(dt,convertPlotCoordPix(x=dt[["x"]], 
#                                    y=dt[["y"]], 
#                                    useMar=getMargin(pxVec), 
#                                    plotDim=c(640,640), 
#                                    plotRes=100))

#ggsave(pngFile, p, height = 640, width = 640, units = "px", dpi = 100)
img <- EBImage::readImage(pngFile)
EBImage::display(img)
#useMar = c(6,4.3,5,2)
# useMar=c(4.97,4.12,4.03,2.12) works for plot() function with (640,640) and res = 100

#(b,l,t,r)
getMargin<- function(pxVec) {
  if (!length(pxVec) == 4) stop("four sides pixel margins not provided")  
  unit_px_to_margin <- 0.0575
  return(pxVec*unit_px_to_margin)
}

theme_to_margin_px <- function(theme_name){
  if (theme_name == "theme_classic" && rand_legend_pos == "none") pxVec <- c(45,48,13,13)
  if (theme_name == "theme_gray" && rand_legend_pos == "none") pxVec <- c(45,48,13,13)
}

#p_built <- ggplot_build(p)
#p_gtable <- ggplot_gtable(p_built)