install.packages("pracma")
library(pracma)
library(EBImage)


# works with EBIMage package
# input: an EBImage
# output: a list containing coordinate values (width, height) of x1, x2, y1, y2  
getXYCoordinates <- function(img) {
  # find the coordinates of axes intersection point (x1 and y1)
  img_gray <- channel(img, "gray")
  img_gray <- img_gray > 0.6
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  #px_fact <- ceiling(max(width, height)/512)
  px_fact <- round(max(width, height)/512)
  px_fact <- ifelse(px_fact == 0, 1, px_fact)
  sum_row <- apply(img_gray, 1, sum)
  sum_row[round(width/2):width] <- NA
  x_axis_start <- which.min(sum_row)
  # logic to make decision for axis by checking 
  # the presence of the sequence of black pixels
  for (w in 1:(round(width/2)-1)) {
    p_x_start <- which(img_gray[w,] == FALSE)
    if (length(p_x_start)) {
      seq_black_px <- rle(diff(p_x_start))$lengths
      if (length(seq_black_px) && max(seq_black_px) > height/5) {
        if (abs(x_axis_start-w) > 3*px_fact) {
          x_axis_start <- w
        } else {
          x_axis_start <- ceiling((x_axis_start+w)/2)
        }
        break()
      }
    }
  }
  #rle(test)$lengths[rle(test)$values==1]
  #diff(unique(cumsum(test == 1)[test != 1]))
  sum_column <- apply(img_gray, 2, sum)
  sum_column[1:round(height/2)] <- NA
  y_axis_start <- which.min(sum_column)
  #y_axis_start <- last(which(sum_column > 0.75))
  # this supports the cases where the logarithmic ticks are span 
  # over the plot, which hallucinates the algorithm to select the x axis
  for (h in height:(round(height/2))) {
    p_y_start <- which(img_gray[,h] == FALSE)
    if (length(p_y_start)) {
      seq_black_px <- rle(diff(p_y_start))$lengths
      if (length(seq_black_px) && max(seq_black_px) > width/5) {
        if (abs(y_axis_start-h) > 3*px_fact) {
          y_axis_start <- h
        } else {
          y_axis_start <- floor((y_axis_start+h)/2)
        }
        break()
      }
    }
  }
  
  # find where x and y axes end
  #white_y_indices <- which(img_gray[x_axis_start,] > 0.8)
  #max.ids <- cumsum(p_y_start)[rle(diff(p_y_start))$lengths == max(seq_black_px)] - c(max(seq_black_px) - 1, 0)
  #max.ids
  white_y_indices <- which(img_gray[x_axis_start,] == TRUE)
  p_y_axis_end <- white_y_indices[white_y_indices < y_axis_start & white_y_indices < height/3]
  y_axis_end <- 1
  if (length(p_y_axis_end))
    y_axis_end <- max(p_y_axis_end) + 1
  #white_x_indices <- which(img_gray[, y_axis_start] > 0.8)
  white_x_indices <- which(img_gray[, y_axis_start] == TRUE)
  p_x_axis_end <- white_x_indices[white_x_indices > x_axis_start & white_x_indices > 2*(width/3)]
  x_axis_end <- width
  if (length(p_x_axis_end))
    x_axis_end <- min(p_x_axis_end) - 1
  initAxesVec <- c(x_axis_start, x_axis_end, y_axis_start, y_axis_end)
  lst <- list()
  lst$x1_coord <- locateXCoordinates(img_gray, initAxesVec, axis_pos = "start", px_fact = px_fact)
  lst$x2_coord <- locateXCoordinates(img_gray, initAxesVec, axis_pos = "end", px_fact = px_fact)
  lst$y1_coord <- locateYCoordinates(img_gray, initAxesVec, axis_pos = "start", px_fact = px_fact)
  lst$y2_coord <- locateYCoordinates(img_gray, initAxesVec, axis_pos = "end", px_fact = px_fact)
  return(lst)
}

# find the y1 and y2
locateYCoordinates <- function(img_gray, initAxesVec, axis_pos = "start", ax_val_th = 0.99, img_gry_th = 0.6, px_fact = 1) {
  loc <- 2
  if (!axis_pos == "start")
    loc <- 3
  axesVec <- c()
  
  # logic to check ticks outside (cartesian coordinate 2)
  #p_coord_points <- which(img_gray[(initAxesVec[1]-2),] < gray_threshold)
  p_coord_points <- which(img_gray[(initAxesVec[1]-(3*px_fact)),] == FALSE)
  p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1+3]+(2*px_fact)) & (p_coord_points<initAxesVec[1+2]-(2*px_fact))]
  if (!length(p_coord_points)) {
    # logic to check ticks inside (cartesian coordinate 1)
    #p_coord_points <- which(img_gray[(initAxesVec[1]+2),] < gray_threshold)
    p_coord_points <- which(img_gray[(initAxesVec[1]+(2*px_fact)),] == FALSE)
    p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1+3]+(2*px_fact)) & (p_coord_points<initAxesVec[1+2]-(2*px_fact))]
  } 
  if (!length(p_coord_points))
    return(c(initAxesVec[1],initAxesVec[1+loc]))
  # logic for linear or exponential
  #p_coord_points <- p_coord_points[which(seq_diff > 1)] # has to be evaluated
  tick_chunks <- split(p_coord_points, cumsum(c(1, diff(p_coord_points) != 1)))
  p_coord_points <- unname(sapply(tick_chunks, function(x) ceiling(mean(x))))
  
  # locate start point
  crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[1+loc], axis = "y", px_fact = px_fact)
  total_area <- dim(crop_area)[1]*dim(crop_area)[2]
  if (sum(crop_area)/total_area > ax_val_th) {
    seq_diff <- diff(p_coord_points)
    # logarithmic with increamental ticks
    if (length(unique(seq_diff)) > 7) {
      p_peak_pos <- findpeaks(seq_diff, threshold = 10)[,2]
      p_peak <- median(seq_diff[p_peak_pos])
      peak_pos_index <- which(seq_diff[p_peak_pos] >= peak-(3*px_fact) & seq_diff[p_peak_pos] <= peak+(3*px_fact))
      peak_pos <- p_peak_pos[peak_pos_index]
      offset_val <- 0
      range <- (length(peak_pos)-1):1
      if (!axis_pos == "start") {
        range <- 1:length(peak_pos)
        offset_val <- 1
      }
      for(len in range) {
        p_coord_point <- p_coord_points[peak_pos[len]+offset_val]
        crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", px_fact = px_fact)
        total_area <- dim(crop_area)[1]*dim(crop_area)[2]
        if (sum(crop_area)/total_area < ax_val_th) {
          axesVec <- c(initAxesVec[1],p_coord_point)
          break()
        }
      }
    } else {
      # linear with equally or nearly equally distributed ticks
      range <- length(p_coord_points):1
      if (!axis_pos == "start")
        range <- 1:length(p_coord_points)
      for(len in range) {
        p_coord_point <- p_coord_points[len]
        crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", px_fact = px_fact)
        total_area <- dim(crop_area)[1]*dim(crop_area)[2]
        if (sum(crop_area)/total_area < ax_val_th) {
          axesVec <- c(initAxesVec[1],p_coord_point)
          break()
        }
      }
    }
  } else {
    if (axis_pos == "start")
      prox_tick <- min(p_coord_points)
    else
      prox_tick <- max(p_coord_points)
    if (abs(initAxesVec[1+loc]-prox_tick) <= 15*px_fact)
      axesVec <- c(initAxesVec[1], prox_tick)
    else
      axesVec <- c(initAxesVec[1],initAxesVec[1+loc])
  }
  return(axesVec)
}


# find the x1 and x2
locateXCoordinates <- function(img_gray, initAxesVec, axis_pos = "start", ax_val_th = 0.99, img_gry_th = 0.6, px_fact = 1) {
  loc <- 1
  if (!axis_pos == "start")
    loc <- 2
  axesVec <- c()
  
  # locate initial area
  crop_area <- getAxisValueCropArea(img_gray, initAxesVec[loc], initAxesVec[3], axis = "x", px_fact = px_fact)
  total_area <- dim(crop_area)[1]*dim(crop_area)[2]
  
  # logic to check ticks outside (cartesian coordinate 4)
  #p_coord_points <- which(img_gray[,(initAxesVec[3]+2)] < gray_threshold)
  p_coord_points <- which(img_gray[,(initAxesVec[3]+(3*px_fact))] == FALSE)
  p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1]+(2*px_fact)) & (p_coord_points<initAxesVec[2]-(2*px_fact))]
  if (!length(p_coord_points)) {
    # logic to check ticks inside (cartesian coordinate 1)
    #p_coord_points <- which(img_gray[,(initAxesVec[3]-2)] < gray_threshold)
    p_coord_points <- which(img_gray[,(initAxesVec[3]-(2*px_fact))] == FALSE)
    p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1]+(2*px_fact)) & (p_coord_points<initAxesVec[2]-(2*px_fact))]
  }
  if (!length(p_coord_points))
    return(c(initAxesVec[loc],initAxesVec[3]))
  # logic for finding normally distributed or logarithmic increamental peaks
  tick_chunks <- split(p_coord_points, cumsum(c(1, diff(p_coord_points) != 1)))
  p_coord_points <- unname(sapply(tick_chunks, function(x) ceiling(mean(x)))) # USE.NAMES = FALSE not working
  
  if (sum(crop_area)/total_area > ax_val_th) {
    seq_diff <- diff(p_coord_points)
    # logarithmic with increamental ticks
    if (length(unique(seq_diff)) > 7) {
      p_peak_pos <- findpeaks(seq_diff, threshold = 10)[,2]
      p_peak <- median(seq_diff[p_peak_pos])
      peak_pos_index <- which(seq_diff[p_peak_pos] >= peak-(3*px_fact) & seq_diff[p_peak_pos] <= peak+(3*px_fact))
      peak_pos <- p_peak_pos[peak_pos_index]
      offset_val <- 1
      range <- 1:(length(peak_pos))
      if (!axis_pos == "start") {
        range <- length(peak_pos):1
        offset_val <- 0
      }
      for(len in range) {
        p_coord_point <- p_coord_points[peak_pos[len]+offset_val]
        crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", px_fact = px_fact)
        total_area <- dim(crop_area)[1]*dim(crop_area)[2]
        if (sum(crop_area)/total_area < ax_val_th) {
          axesVec <- c(p_coord_point,initAxesVec[3])
          break()
        }
      }
    } else {
      # linear with equally or nearly equally distributed ticks
      range <- 1:length(p_coord_points)
      if (!axis_pos == "start")
        range <- length(p_coord_points):1
      for(len in range) {
        p_coord_point <- p_coord_points[len]
        crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", px_fact = px_fact)
        total_area <- dim(crop_area)[1]*dim(crop_area)[2]
        if (sum(crop_area)/total_area < ax_val_th) {
          axesVec <- c(p_coord_point,initAxesVec[3])
          break()
        }
      }
    }
  } else {
    if (axis_pos == "start")
      prox_tick <- min(p_coord_points)
    else
      prox_tick <- max(p_coord_points)
    if (abs(initAxesVec[loc]-prox_tick) <= 15*px_fact)
      axesVec <- c(prox_tick,initAxesVec[3])
    else
      axesVec <- c(initAxesVec[loc],initAxesVec[3])
  }
  return(axesVec)
}

# get axis crop area to estimate presence of value
getAxisValueCropArea <- function(img_gray, w, h, axis = "x", section = "left", px_fact = 1) {
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  wVec <- c(1,width)
  hVec <- c(1,height)
  
  if (section == "left") {
    if (axis == "y") {
      width_min <- w-(17*px_fact)
      width_max <- w-(6*px_fact)
      height_min <- h-px_fact
      height_max <- h+(2*px_fact)
    } else if (axis == "x") {
      width_min <- w-(3*px_fact)
      width_max <- w+px_fact
      height_min <- h+(6*px_fact)
      height_max <- h+(17*px_fact)
    }
  } else {
    if (axis == "y") {
      width_min <- w-(17*px_fact)
      width_max <- w-(6*px_fact)
      height_min <- h-px_fact
      height_max <- h+(2*px_fact)
    } else if (axis == "x") {
      width_min <- w-(3*px_fact)
      width_max <- w+px_fact
      height_min <- h+(6*px_fact)
      height_max <- h+(17*px_fact)
    }
  }
  #condition to prevent values outside of image
  if (width_min < 1 || width_min > width)
    width_min <- wVec[which.min(c(width_min-wVec[1],wVec[2]-width_min))]
  if (width_max < 1 || width_max > width)
    width_max <- wVec[which.min(c(width_max-wVec[1],wVec[2]-width_max))]
  if (height_min < 1 || height_min > height)
    height_min <- hVec[which.min(c(height_min-hVec[1],hVec[2]-height_min))]
  if (height_max < 1 || height_max > height)
    height_max <- hVec[which.min(c(height_max-hVec[1],hVec[2]-height_max))]
  
  crop_area <- img_gray[width_min:width_max,height_min:height_max]
  return(crop_area)
}

displayXYCoordinates <- function(img, lst) {
  px_coords <- list()
  img_draw <- EBImage::drawCircle(img, lst$x1_coord[1], lst$x1_coord[2], 3, col = 1, fill = T)
  img_draw <- EBImage::drawCircle(img_draw, lst$y1_coord[1], lst$y1_coord[2], 3, col = 1, fill = T)
  img_draw <- EBImage::drawCircle(img_draw, lst$x2_coord[1], lst$x2_coord[2], 3, col = 1, fill = T)
  img_draw <- EBImage::drawCircle(img_draw, lst$y2_coord[1], lst$y2_coord[2], 3, col = 1, fill = T)
  display(img_draw)
}


img <- EBImage::readImage("C:/Users/guptav/Pictures/test2/gr1.png")
display(img)
coords <- getXYCoordinates(img)
displayXYCoordinates(img, coords)


a_width_min <- 60/2
58 <- 59
a_height_min <- 193 - 2
a_height_max <- 193 + 2
plot(img_gray[a_width_min:a_width_max,a_height_min:a_height_max])
sum_row <- apply(img_gray[a_width_min:a_width_max,a_height_min:a_height_max], 1, sum)
plot(sum_row, type = "l")
dim <- dim(img_gray[a_width_min:a_width_max,a_height_min:a_height_max])
area <- dim[1]*dim[2]
sum(img_gray[a_width_min:a_width_max,a_height_min:a_height_max])


# logic for linear or exponential
#breaks <- c(0, which(diff(p_coord_points) != 1), length(p_coord_points))
#cons_list <- sapply(seq(length(breaks) - 1), function(i) p_coord_points[(breaks[i] + 1):breaks[i+1]])



fit1 <- fitdistr(a, "exponential") 
fit2 <- fitdistr(control, "exponential")

# goodness of fit test
ks.test(ex, "pexp", fit1$estimate) # p-value > 0.05 -> distribution not refused
ks.test(control, "pexp", fit2$estimate) #  significant p-value -> distribution refused

# plot a graph
hist(ex, freq = FALSE, breaks = 100, xlim = c(0, quantile(ex, 0.99)))
curve(dexp(x, rate = fit1$estimate), from = 0, col = "red", add = TRUE)