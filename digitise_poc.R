library(EBImage)

#binary_threshold <- 0.7
img <- EBImage::readImage("C:/Users/guptav/Pictures/test2/gr2.png")

# convert to binary, remove where axis unlikely, extract
img_gray <- EBImage::channel(img, mode = "gray") < binary_threshold
hb_img_gray <- img_gray[, round(dim(img_gray)[2]/2):dim(img_gray)[2]]
vl_img_gray <- img_gray[1:round(dim(img_gray)[1]/2), ]
lineBrush <- EBImage::makeBrush((dim(hb_img_gray)[2]), shape = "line", angle = 0)
lines_on_hb <- EBImage::opening(EBImage::distmap(hb_img_gray), lineBrush)
x_axis <- EBImage::watershed(EBImage::distmap(lines_on_hb), tolerance = 1, ext = 1)

# eliminate all but the lowermost
allLines <- EBImage::computeFeatures.shape(x_axis)
exclusionList <- which(allLines[, "s.area"] != max(allLines[, "s.area"]))
allDetectedX <- EBImage::rmObjects(allDetectedX, exclusionList)
theCoordinates <- EBImage::computeFeatures.moment(allDetectedX)
exclusionList <- which(theCoordinates[, "m.cy"] != max(theCoordinates[, "m.cy"]))
allDetectedX <- EBImage::rmObjects(allDetectedX, exclusionList)

# erase everything outside detected axis
Xcontr <- EBImage::ocontour(x_axis)
Xmax <- max(Xcontr[[1]][, 1]); Xmin <- min(Xcontr[[1]][, 1])
img_gray[c(1:(Xmin + 3), Xmax:dim(img_gray)[1]), ] <- 0


theAxis <- EBImage::ocontour(allDetectedX)
coordX1 <- min(theAxis[[1]][, 1]); coordY1 <- min(theAxis[[1]][, 2]);
coordX2 <- max(theAxis[[1]][, 1]); coordY2 <- max(theAxis[[1]][, 2]);



# works with EBIMage package
# input: an EBImage
# output: a list containing coordinate values (width, height) of x1, x2, y1, y2  
getXYCoordinates <- function(img) {
  # find the coordinates of axes intersection point (x1 and y1)
  img_gray <- channel(img, "gray")
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  sum_width <- apply(img_gray, 1, sum)
  sum_width[round(width/2):width] <- NA
  width_axes_meet <- which.min(sum_width)
  sum_height <- apply(img_gray, 2, sum)
  sum_height[1:round(height/2)] <- NA
  height_axes_meet <- which.min(sum_height)
  
  # find the coordinates of x2 and y2
  white_y_indices <- which(img_gray[width_axes_meet,] > 0.8)
  height_y_axis_end <- max(white_y_indices[white_y_indices < height_axes_meet]) + 1
  white_x_indices <- which(img_gray[, height_axes_meet] > 0.8)
  width_x_axis_end <- min(white_x_indices[white_x_indices > width_axes_meet]) - 1
  
  lst <- list()
  lst$X1_coord <- c(width_axes_meet, height_axes_meet)
  lst$Y1_coord <- c(width_axes_meet, height_axes_meet)
  lst$X2_coord <- c(width_x_axis_end, height_axes_meet)
  lst$Y2_coord <- c(width_axes_meet, height_y_axis_end)
  return(lst)
}

displayXYCoordinates <- function(img, lst) {
  X1Y1 <- EBImage::drawCircle(img, lst$X1_coord[1], lst$X1_coord[2], 3, col = 1, fill = T)
  X2 <- EBImage::drawCircle(X1Y1, lst$X2_coord[1], lst$X2_coord[2], 3, col = 1, fill = T)
  Y2 <- EBImage::drawCircle(X2, lst$Y2_coord[1], lst$Y2_coord[2], 3, col = 1, fill = T)
  display(Y2)
}



# works with png package
getXYCoordinates_PNG <- function(img_obj) {
  height <- dim(img_obj)[1]
  width <- dim(img_obj)[2]
  sum_row <- apply(img_obj[,,1:3], 1, sum)
  sum_row[1:round(height/2)] <- NA
  x_axis_start <- which.min(sum_row)
  sum_column <- apply(img_obj[,,1:3], 2, sum)
  sum_column[round(width/2):width] <- NA
  y_axis_start <- which.min(sum_column)
  white_column_index <- which(img_obj[x_axis_start, , 1] > 0.8)
  x_axis_end <- min(white_column_index[white_column_index > y_axis_start]) - 1
  white_row_index <- which(img_obj[, y_axis_start, 1] > 0.8)
  y_axis_end <- max(white_row_index[white_row_index < x_axis_start]) + 1
  lst <- list()
  lst$X1_coord <- c(x_axis_start, y_axis_start)
  lst$Y1_coord <- c(x_axis_start, y_axis_start)
  lst$X2_coord <- c(x_axis_start, x_axis_end)
  lst$Y2_coord <- c(y_axis_end, y_axis_start)
  return(lst)
}


######## post processing ########
dtt <- data.table(dt)
setnames(dtt, c("x", "y"), c("x_variable", "y_variable"))
dtt <- dtt[, -c("x_variable", "y_variable")]
View(dtt)
dtt <- data.table(dt)
dtt <- data.table(dt)
setnames(dtt, c("x", "y"), c(unique(dtt$x_variable), unique(dtt$y_variable)))
View(dtt)
dtt <- dtt[, -c("group", "col", "pch", "x_variable", "y_variable")]