install.packages("pracma")
library(pracma)
library(EBImage)


# works with EBIMage package, <to be changed to magick>
# input: an EBImage object
# output: a list containing coordinate values (width, height) of x1, x2, y1, y2  
getXYCoordinates <- function(img) {
  
  # find the coordinates of axes intersection point (x1 and y1)
  img_gray <- channel(img, "gray")
  img_gray <- img_gray > 0.8
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  px_fact <- round(max(width, height)/512)
  px_fact <- ifelse(px_fact == 0, 1, px_fact)
  sum_row <- apply(img_gray, 1, sum)
  sum_row[round(width/2):width] <- NA
  axes_intersect_width <- which.min(sum_row)
  
  # logic to make decision for axis by checking 
  # the presence of the sequence of black pixels
  for (w in 1:(round(width/2)-1)) {
    p_axes_intersect_width <- which(img_gray[w,] == FALSE)
    if (length(p_axes_intersect_width)) {
      seq_black_px <- rle(diff(p_axes_intersect_width))$lengths
      if (length(seq_black_px) && max(seq_black_px) > height/4) {
        if (abs(axes_intersect_width-w) > 4*px_fact) {
          axes_intersect_width <- w
        }
        if (px_fact > 1) {
          axes_intersect_width <- w + 2*(px_fact-(1/2))
        }
        
        # find the y2 height
        y_black_px <- which(img_gray[axes_intersect_width,] == FALSE)
        rle_y <- rle(diff(y_black_px))
        y_seq_indices <- which(rle_y$values == 1)
        max_value_y_index <- which.max(rle_y$lengths[y_seq_indices])
        max_seq_y_index <- y_seq_indices[max_value_y_index]
        cum_sum_y_index <- cumsum(rle_y$lengths)
        y_axis_end_height_index <- ifelse(max_seq_y_index == 1, 1, cum_sum_y_index[max_seq_y_index-1]+1)
        y_axis_end_height <- y_black_px[y_axis_end_height_index] - 1 + 2*(px_fact-(1/2))
        break
      }
    }
  }
  sum_column <- apply(img_gray, 2, sum)
  sum_column[1:round(height/2)] <- NA
  axes_intersect_height <- which.min(sum_column)
  
  # this supports the cases where the logarithmic ticks are span 
  # over the plot, which hallucinates the algorithm to select the x axis
  for (h in height:(round(height/2))) {
    p_axes_intersection_height <- which(img_gray[,h] == FALSE)
    if (length(p_axes_intersection_height)) {
      seq_black_px <- rle(diff(p_axes_intersection_height))$lengths
      
      # check if a line-like object is detected
      if (length(seq_black_px) && max(seq_black_px) > width/4) {
        
        # check the proximity of the detected line object with the initial line object
        if (abs(axes_intersect_height-h) > 4*px_fact) {
          axes_intersect_height <- h #- 2*(px_fact-(1/2)) 
        }
        
        # find the center point of the line object in cases when pixel factor is more than 1
        if (px_fact > 1) {
          axes_intersect_height <- h - 2*(px_fact-(1/2)) 
        }
        
        # find the x2 width
        x_black_px <- which(img_gray[,axes_intersect_height] == FALSE)
        rle_x <- rle(diff(x_black_px))
        x_seq_indices <- which(rle_x$values == 1)
        max_value_x_index <- which.max(rle_x$lengths[x_seq_indices])
        max_seq_x_index <- x_seq_indices[max_value_x_index]
        cum_sum_x_index <- cumsum(rle_x$lengths)
        x_axis_end_width_index <- ifelse(max_seq_x_index == length(rle_x$lengths), length(x_black_px), cum_sum_x_index[max_seq_x_index-1])
        x_axis_end_width <-  x_black_px[x_axis_end_width_index] + 1 - 2*(px_fact-(1/2))
        break
      }
    }
  }
  
  # special case: when axes are not intersecting
  # if (img_gray[axes_intersect_width,axes_intersect_height]) {
  #   x_axis_start_width_index <- ifelse(max_seq_x_index == 1, 1, cum_sum_x_index[max_seq_x_index-1]) + 1
  #   axes_intersect_width <- x_black_px[x_axis_start_width_index] + 1 - 2*(px_fact-(1/2))
  #   y_axis_start_height_index <- ifelse(max_seq_y_index == 1, rle_y$lengths+1, cum_sum_y_index[max_seq_y_index]+1)
  #   axes_intersect_height <- y_black_px[y_axis_start_height_index] - 1 + 2*(px_fact-(1/2))
  # }
  
  
  initAxesVec <- c(axes_intersect_width, x_axis_end_width, axes_intersect_height, y_axis_end_height)
  if (any(is.na(initAxesVec)) || !length(initAxesVec) == 4)
    stop("all initial axes points could not be located")
  lst <- list()
  lst$x1_coord <- locateXCoordinates(img_gray, initAxesVec, axis_pos = "start", px_fact = px_fact)
  lst$x2_coord <- locateXCoordinates(img_gray, initAxesVec, axis_pos = "end", px_fact = px_fact)
  lst$y1_coord <- locateYCoordinates(img_gray, initAxesVec, axis_pos = "start", px_fact = px_fact)
  lst$y2_coord <- locateYCoordinates(img_gray, initAxesVec, axis_pos = "end", px_fact = px_fact)
  
  # special case: when axes are not intersecting
  if (all(img_gray[lst$x1_coord[1]:(lst$x1_coord[1]+1), lst$x1_coord[2]:(lst$x1_coord[2]+1)]))
    lst$x1_coord[1] <- lst$x1_coord[1] + which.min(img_gray[lst$x1_coord[1]:lst$x2_coord[1], lst$x1_coord[2]]) - 1 + 2*(px_fact-(1/2))
  if (all(img_gray[lst$y1_coord[1]:(lst$y1_coord[1]+1), lst$y1_coord[2]:(lst$y1_coord[2]+1)]))
    lst$y1_coord[2] <- lst$y1_coord[2] - which.min(img_gray[lst$y1_coord[1], lst$y1_coord[2]:lst$y2_coord[2]]) + 2*(px_fact-(1/2))
  return(lst)
}

# find the y1 and y2
locateYCoordinates <- function(img_gray, initAxesVec, axis_pos = "start", ax_val_th = 0.99, px_fact = 1) {
  loc <- 2
  if (!axis_pos == "start")
    loc <- 3
  axesVec <- c()
  
  # check the presence of the start point
  top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[1+loc], axis = "y", section = "top", px_fact = px_fact)
  bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[1+loc], axis = "y", section = "bottom", px_fact = px_fact)
  total_top_area <- dim(top_crop_area)[1]*dim(top_crop_area)[2]
  total_bottom_area <- dim(bottom_crop_area)[1]*dim(bottom_crop_area)[2]
  
  # logic to check whether axis value is found or not
  if (sum(top_crop_area)/total_top_area > ax_val_th || sum(bottom_crop_area)/total_bottom_area > ax_val_th) {
    
    # logic to check ticks outside (cartesian coordinate 2)
    p_coord_points <- which(img_gray[(initAxesVec[1]-(3*px_fact)),] == FALSE)
    p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1+3]+(2*px_fact)) & (p_coord_points<initAxesVec[1+2]-(2*px_fact))]
    if (!length(p_coord_points)) {
      
      # logic to check ticks inside (cartesian coordinate 1)
      p_coord_points <- which(img_gray[(initAxesVec[1]+(3*px_fact)),] == FALSE)
      p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1+3]+(2*px_fact)) & (p_coord_points<initAxesVec[1+2]-(2*px_fact))]
    } 
    if (!length(p_coord_points))
      return(c(initAxesVec[1],initAxesVec[1+loc]))
    
    # logic for linear or exponential
    tick_chunks <- split(p_coord_points, cumsum(c(1, diff(p_coord_points) != 1)))
    p_coord_points <- unname(sapply(tick_chunks, function(x) ceiling(mean(x))))
    seq_diff <- diff(p_coord_points)
    
    # logarithmic with increamental ticks
    if (length(unique(seq_diff)) > 7) {
      p_peak_pos <- findpeaks(seq_diff, threshold = 5)[,2]
      p_peak <- median(seq_diff[p_peak_pos])
      peak_pos_index <- which(seq_diff[p_peak_pos] >= p_peak-(3*px_fact) & seq_diff[p_peak_pos] <= p_peak+(3*px_fact))
      peak_pos <- p_peak_pos[peak_pos_index]
      range <- (length(peak_pos)-1):1
      if (!axis_pos == "start")
        range <- 1:length(peak_pos)
      for(len in range) {
        p_coord_point <- p_coord_points[peak_pos[len]+1]
        top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "top", px_fact = px_fact)
        bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "bottom", px_fact = px_fact)
        total_top_area <- dim(top_crop_area)[1]*dim(top_crop_area)[2]
        total_bottom_area <- dim(bottom_crop_area)[1]*dim(bottom_crop_area)[2]
        if (sum(top_crop_area)/total_top_area < ax_val_th && sum(bottom_crop_area)/total_bottom_area < ax_val_th) {
          axesVec <- c(initAxesVec[1],p_coord_point)
          break
        }
      }
    } else {
      
      # linear with equally or nearly equally distributed ticks
      range <- length(p_coord_points):1
      if (!axis_pos == "start")
        range <- 1:length(p_coord_points)
      for(len in range) {
        p_coord_point <- p_coord_points[len]
        top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "top", px_fact = px_fact)
        bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "bottom", px_fact = px_fact)
        total_top_area <- dim(top_crop_area)[1]*dim(top_crop_area)[2]
        total_bottom_area <- dim(bottom_crop_area)[1]*dim(bottom_crop_area)[2]
        if (sum(top_crop_area)/total_top_area < ax_val_th && sum(bottom_crop_area)/total_bottom_area < ax_val_th) {
          axesVec <- c(initAxesVec[1],p_coord_point)
          break
        }
      }
    }
  } else {
    # if axis value found
    axesVec <- c(initAxesVec[1],initAxesVec[1+loc])
  }
  # if axis value couldn't be located then return initial vector
  if (any(is.na(axesVec)) || !length(axesVec) == 2)
    axesVec <- c(initAxesVec[1],initAxesVec[1+loc])

  return(axesVec)
}


# find the x1 and x2
locateXCoordinates <- function(img_gray, initAxesVec, axis_pos = "start", ax_val_th = 0.99, px_fact = 1) {
  loc <- 1
  if (!axis_pos == "start")
    loc <- 2
  axesVec <- c()
  
  # locate initial area
  left_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[loc]-1, initAxesVec[3], axis = "x", section = "left", px_fact = px_fact)
  right_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[loc]-1, initAxesVec[3], axis = "x", section = "right", px_fact = px_fact)
  total_left_area <- dim(left_crop_area)[1]*dim(left_crop_area)[2]
  total_right_area <- dim(right_crop_area)[1]*dim(right_crop_area)[2]  
  
  # logic to check whether axis value is found or not
  if (sum(left_crop_area)/total_left_area > ax_val_th || sum(right_crop_area)/total_right_area > ax_val_th) {
    
    # if axis value not found
    # logic to check ticks outside (cartesian coordinate 4)
    p_coord_points <- which(img_gray[,(initAxesVec[3]+(3*px_fact))] == FALSE)
    p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1]+(2*px_fact)) & (p_coord_points<initAxesVec[2]-(2*px_fact))]
    if (!length(p_coord_points)) {
      
      # logic to check ticks inside (cartesian coordinate 1)
      p_coord_points <- which(img_gray[,(initAxesVec[3]-(3*px_fact))] == FALSE)
      p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1]+(2*px_fact)) & (p_coord_points<initAxesVec[2]-(2*px_fact))]
    }
    if (!length(p_coord_points))
      return(c(initAxesVec[loc],initAxesVec[3]))
    
    # logic for finding normally distributed or logarithmic increamental peaks
    tick_chunks <- split(p_coord_points, cumsum(c(1, diff(p_coord_points) != 1)))
    p_coord_points <- unname(sapply(tick_chunks, function(x) ceiling(mean(x)))) # USE.NAMES = FALSE not working
    seq_diff <- diff(p_coord_points)
    
    # logarithmic with increamental ticks
    if (length(unique(seq_diff)) > 7) {
      p_peak_pos <- findpeaks(seq_diff, threshold = 5)[,2]
      p_peak <- median(seq_diff[p_peak_pos])
      peak_pos_index <- which(seq_diff[p_peak_pos] >= p_peak-(3*px_fact) & seq_diff[p_peak_pos] <= p_peak+(3*px_fact))
      peak_pos <- p_peak_pos[peak_pos_index]
      range <- 1:(length(peak_pos))
      if (!axis_pos == "start")
        range <- length(peak_pos):1
      for(len in range) {
        p_coord_point <- p_coord_points[peak_pos[len]]
        left_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "left", px_fact = px_fact)
        right_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "right", px_fact = px_fact)
        total_left_area <- dim(left_crop_area)[1]*dim(left_crop_area)[2]
        total_right_area <- dim(right_crop_area)[1]*dim(right_crop_area)[2]
        if (sum(left_crop_area)/total_left_area < ax_val_th && sum(right_crop_area)/total_right_area < ax_val_th) {
          axesVec <- c(p_coord_point,initAxesVec[3])
          break
        }
      }
    } else {
      
      # linear with equally or nearly equally distributed ticks
      range <- 1:length(p_coord_points)
      if (!axis_pos == "start")
        range <- length(p_coord_points):1
      for(len in range) {
        p_coord_point <- p_coord_points[len]
        left_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "left", px_fact = px_fact)
        right_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "right", px_fact = px_fact)
        total_left_area <- dim(left_crop_area)[1]*dim(left_crop_area)[2]
        total_right_area <- dim(right_crop_area)[1]*dim(right_crop_area)[2]
        if (sum(left_crop_area)/total_left_area < ax_val_th && sum(right_crop_area)/total_right_area < ax_val_th) {
          axesVec <- c(p_coord_point,initAxesVec[3])
          break
        }
      }
    }
  } else {
    
    # if axis value found
    axesVec <- c(initAxesVec[loc],initAxesVec[3])
  }
  
  # if axis value couldn't be located then return initial vector
  if (any(is.na(axesVec)) || !length(axesVec) == 2)
    axesVec <- c(initAxesVec[loc],initAxesVec[3])
  return(axesVec)
}

# get axis crop area to estimate presence of tick value
getAxisValueCropArea <- function(img_gray, w, h, axis = "x", section = "top", px_fact = 1) {
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  wVec <- c(1,width)
  hVec <- c(1,height)
  
  if (!section %in% c("bottom", "top", "left", "right"))
    stop("section value is incorrect")
  
  if (is.null(w) || is.null(h))
    return(0)
  
  if (section == "top" || section == "left") {
    if (axis == "y") {
      prox_dist_fact <- ceiling((w/px_fact)/64)
      #prox_dist_fact <- ifelse(prox_dist_fact == 0, 1, prox_dist_fact)
      width_min <- w-(20*px_fact*prox_dist_fact)
      width_max <- w-(6*px_fact*prox_dist_fact)
      height_min <- h-(3*px_fact)
      height_max <- h-px_fact
    } else if (axis == "x") {
      prox_dist_fact <- ceiling(((height-h)/px_fact)/64)
      #prox_dist_fact <- ifelse(prox_dist_fact == 0, 1, prox_dist_fact)
      width_min <- w-(5*px_fact)
      width_max <- w
      height_min <- h+(6*px_fact*prox_dist_fact)
      height_max <- h+(18*px_fact*prox_dist_fact)
    }
  } else if (section == "bottom" || section == "right") {
    if (axis == "y") {
      prox_dist_fact <- ceiling((w/px_fact)/64)
      #prox_dist_fact <- ifelse(prox_dist_fact == 0, 1, prox_dist_fact)
      width_min <- w-(20*px_fact*prox_dist_fact)
      width_max <- w-(6*px_fact*prox_dist_fact)
      height_min <- h+px_fact
      height_max <- h+(3*px_fact)
    } else if (axis == "x") {
      prox_dist_fact <- ceiling(((height-h)/px_fact)/64)
      #prox_dist_fact <- ifelse(prox_dist_fact == 0, 1, prox_dist_fact)
      width_min <- w+px_fact
      width_max <- w+(5*px_fact)
      height_min <- h+(6*px_fact*prox_dist_fact)
      height_max <- h+(18*px_fact*prox_dist_fact)
    }
  }
  
  # condition to prevent values outside of image
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

