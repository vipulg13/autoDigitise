install.packages("pracma")
library(pracma)
library(EBImage)


# works with EBIMage package at the moment, <to be changed to magick>
# input: an EBImage object
# output: a list containing coordinate values (width, height) of x1, x2, y1, y2  
getXYCoordinates <- function(img) {
  
  # find the coordinates of axes intersection point (x1 and y1)
  img_gray <- channel(img, "gray")
  img_gray <- img_gray > 0.85
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  px_fact_width <- round(width/512)
  px_fact_width <- ifelse(px_fact_width == 0, 1, px_fact_width)
  px_fact_height <- round(height/512)
  px_fact_height <- ifelse(px_fact_height == 0, 1, px_fact_height)
  
  # find the axes intersection width
  sum_row <- apply(img_gray, 1, sum)
  sum_row[round(width/2):width] <- NA
  p_axes_intersect_width <- which.min(sum_row)
  y_black_px <- which(!img_gray[p_axes_intersect_width,])
  init_y_seq_black_px <- rle(diff(y_black_px))$lengths
  if (length(init_y_seq_black_px) && max(init_y_seq_black_px) > height/5) {
    if (!px_fact_width == 1) {
      width_min_range <- p_axes_intersect_width - 6*px_fact_width
      if (width_min_range < 1) width_min_range <- 1
      width_max_range <- p_axes_intersect_width + 6*px_fact_width
      width_vec <- c()
      for (w in width_min_range:width_max_range) {
        y_black_px <- which(!img_gray[w,])
        y_seq_black_px <- rle(diff(y_black_px))$lengths
        if (!length(y_seq_black_px)) y_seq_black_px <- 0
        if (abs(max(init_y_seq_black_px) - max(y_seq_black_px))/max(init_y_seq_black_px) < .1) {
          width_vec <- c(width_vec, w)
        }
      }
      axes_intersect_width <- round(mean(width_vec))
    } else {
      axes_intersect_width <- p_axes_intersect_width
    }
  }

  # find the y2 height
  y_black_px <- which(!img_gray[axes_intersect_width,])
  rle_y <- rle(diff(y_black_px))
  y_seq_indices <- which(rle_y$values == 1)
  max_value_y_index <- which.max(rle_y$lengths[y_seq_indices])
  max_seq_y_index <- y_seq_indices[max_value_y_index]
  cum_sum_y_index <- cumsum(rle_y$lengths)
  y_axis_end_height_index <- ifelse(max_seq_y_index == 1, 1, cum_sum_y_index[max_seq_y_index-1]+1)
  y_axis_end_height <- y_black_px[y_axis_end_height_index] + px_fact_height - 1
  
  # find the y1 height
  sum_column <- apply(img_gray, 2, sum)
  sum_column[1:round(height/2)] <- NA
  init_p_axes_intersect_height <- which.min(sum_column)

  # find the axes intersection height
  # this supports the cases where the logarithmic ticks are span 
  # over the plot, which hallucinates the algorithm to select the x axis
  for (h in height:(init_p_axes_intersect_height - 6*px_fact_height)) {
    init_x_black_px <- which(!img_gray[,h])
    init_x_seq_black_px <- rle(diff(init_x_black_px))$lengths
    if (length(init_x_seq_black_px) && max(init_x_seq_black_px) > width/5) {
      if (px_fact_height == 1 && abs(h - init_p_axes_intersect_height) < 3) {
        axes_intersect_height <- init_p_axes_intersect_height
      } else if (abs(h - init_p_axes_intersect_height) >= 3) {
        height_max_range <- h - 1
        height_min_range <- h - 7*px_fact_height
        height_vec <- c(h)
        # loop to find the center value
        for (h1 in height_max_range:height_min_range) {
          x_black_px <- which(!img_gray[,h1])
          x_seq_black_px <- rle(diff(x_black_px))$lengths
          if (!length(x_seq_black_px)) x_seq_black_px <- 0
          if (abs(max(init_x_seq_black_px) - max(x_seq_black_px))/max(init_x_seq_black_px) < .1) {
            height_vec <- c(height_vec, h1)
          } else {
            break
          }
        }
        axes_intersect_height <- round(mean(height_vec))
      } else {
        axes_intersect_height <- h
      }
      break
    }
  }
  
  # find the x2 width
  x_black_px <- which(!img_gray[,axes_intersect_height])
  rle_x <- rle(diff(x_black_px))
  x_seq_indices <- which(rle_x$values == 1)
  max_value_x_index <- which.max(rle_x$lengths[x_seq_indices])
  max_seq_x_index <- x_seq_indices[max_value_x_index]
  cum_sum_x_index <- cumsum(rle_x$lengths)
  x_axis_end_width_index <- ifelse(max_seq_x_index == length(rle_x$lengths), length(x_black_px), cum_sum_x_index[max_seq_x_index-1])
  x_axis_end_width <-  x_black_px[x_axis_end_width_index] - px_fact_width + 1
  
  initAxesVec <- c(axes_intersect_width, x_axis_end_width, axes_intersect_height, y_axis_end_height)
  if (any(is.na(initAxesVec)) || !length(initAxesVec) == 4)
    stop("all initial axes points could not be located")
  lst <- list()
  
  # approximate the pixel distance of values from axes
  empty_flag <- FALSE
  prox_dist <- NA
  for (i in initAxesVec[1]:1) {
    black_px <-  which(!img_gray[i,(initAxesVec[4]+2):(initAxesVec[3]-2)])
    if (!length(black_px)) empty_flag <- TRUE
    if (length(black_px) && empty_flag) {
      prox_dist <- initAxesVec[1] - i -1
      break 
    }
  }
  
  # case when x1 and y1 have a common value point
  top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[3], axis = "y", section = "top", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[3], axis = "y", section = "bottom", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  left_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[3], axis = "x", section = "left", common_point = TRUE, prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  right_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[3], axis = "x", section = "right", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  total_top_area <- dim(top_crop_area)[1]*dim(top_crop_area)[2]
  total_bottom_area <- dim(bottom_crop_area)[1]*dim(bottom_crop_area)[2]
  total_left_area <- dim(left_crop_area)[1]*dim(left_crop_area)[2]
  total_right_area <- dim(right_crop_area)[1]*dim(right_crop_area)[2]
  top_crit <- sum(top_crop_area)/total_top_area > 0.99
  bottom_crit <- sum(bottom_crop_area)/total_bottom_area > 0.99
  left_crit <- sum(left_crop_area)/total_left_area > 0.99
  right_crit <- sum(right_crop_area)/total_right_area > 0.99
  
  tryCatch({
    # if value is common
    if ((right_crit && top_crit) && !(left_crit && bottom_crit)) {
      lst$x1_coord <- c(initAxesVec[1],initAxesVec[3])
      lst$y1_coord <- c(initAxesVec[1],initAxesVec[3])
    } else {
      # if value is not common
      lst$x1_coord <- locateXCoordinates(img_gray, initAxesVec, axis_pos = "start", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
      lst$y1_coord <- locateYCoordinates(img_gray, initAxesVec, axis_pos = "start", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
    }
    # get x2 and y2
    lst$x2_coord <- locateXCoordinates(img_gray, initAxesVec, axis_pos = "end", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
    lst$y2_coord <- locateYCoordinates(img_gray, initAxesVec, axis_pos = "end", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  }, 
  error = function(e) {
    if(length(lst$x1_coord) < 2) lst$x1_coord <- c(initAxesVec[1],initAxesVec[3])
    if(length(lst$y1_coord) < 2) lst$y1_coord <- c(initAxesVec[1],initAxesVec[3])
    if(length(lst$x2_coord) < 2) lst$x2_coord <- c(initAxesVec[2],initAxesVec[3])
    if(length(lst$y2_coord) < 2) lst$y2_coord <- c(initAxesVec[1],initAxesVec[4])
  })
  
  # special case: when axes are not intersecting
  if (all(img_gray[lst$x1_coord[1]:(lst$x1_coord[1]+1), (lst$x1_coord[2]-1):(lst$x1_coord[2])]))
    lst$x1_coord[1] <- lst$x1_coord[1] + which.min(img_gray[lst$x1_coord[1]:lst$x2_coord[1], lst$x1_coord[2]]) - 1 + 2*(px_fact_width-(1/2))
  if (all(img_gray[lst$y1_coord[1]:(lst$y1_coord[1]+1), lst$y1_coord[2]:(lst$y1_coord[2]-1)]))
    lst$y1_coord[2] <- lst$y1_coord[2] - which.min(img_gray[lst$y1_coord[1], lst$y1_coord[2]:lst$y2_coord[2]]) + 2*(px_fact_height-(1/2))
  return(lst)
  
}

# find the y1 and y2
locateYCoordinates <- function(img_gray, initAxesVec, axis_pos = "start", ax_val_th = 0.99, prox_dist = prox_dist, px_fact_width = 1, px_fact_height = 1) {
  loc <- 2
  if (!axis_pos == "start")
    loc <- 3
  axesVec <- c()
  
  # check the presence of the start point
  top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[1+loc], axis = "y", section = "top", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], initAxesVec[1+loc], axis = "y", section = "bottom", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  total_top_area <- dim(top_crop_area)[1]*dim(top_crop_area)[2]
  total_bottom_area <- dim(bottom_crop_area)[1]*dim(bottom_crop_area)[2]
  
  # logic to check whether axis value is found or not
  if (sum(top_crop_area)/total_top_area > ax_val_th || sum(bottom_crop_area)/total_bottom_area > ax_val_th) {
    
    #logic to find ticks
    for (i in 3:1) {
      # logic to check ticks outside (cartesian coordinate 2)
      if (initAxesVec[1]-(i*px_fact_width) >= 1) {
        p_coord_points <- which(!img_gray[(initAxesVec[1]-(i*px_fact_width)),])
        p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1+3]+(2*px_fact_height)) & (p_coord_points<initAxesVec[1+2]-(2*px_fact_height))] 
      }
      # logic to check ticks inside (cartesian coordinate 1)
      if (!length(p_coord_points)) {
        p_coord_points <- which(!img_gray[(initAxesVec[1]+(i*px_fact_width)),])
        p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1+3]+(2*px_fact_height)) & (p_coord_points<initAxesVec[1+2]-(2*px_fact_height))]
      }
      if (length(p_coord_points) && length(p_coord_points[!diff(p_coord_points) == 1]))
        break
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
      peak_pos_index <- which(seq_diff[p_peak_pos] >= p_peak-(3*px_fact_height) & seq_diff[p_peak_pos] <= p_peak+(3*px_fact_height))
      peak_pos <- p_peak_pos[peak_pos_index]
      range <- (length(peak_pos)-1):1
      if (!axis_pos == "start")
        range <- 1:length(peak_pos)
      for(len in range) {
        p_coord_point <- p_coord_points[peak_pos[len]+1]
        top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "top", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
        bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "bottom", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
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
        top_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "top", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
        bottom_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[1], p_coord_point, axis = "y", section = "bottom", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
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
locateXCoordinates <- function(img_gray, initAxesVec, axis_pos = "start", ax_val_th = 0.99, prox_dist = prox_dist, px_fact_width = 1, px_fact_height = 1) {
  loc <- 1
  if (!axis_pos == "start")
    loc <- 2
  axesVec <- c()
  
  # locate initial area
  left_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[loc]-1, initAxesVec[3], axis = "x", section = "left", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  right_crop_area <- getAxisValueCropArea(img_gray, initAxesVec[loc]-1, initAxesVec[3], axis = "x", section = "right", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
  total_left_area <- dim(left_crop_area)[1]*dim(left_crop_area)[2]
  total_right_area <- dim(right_crop_area)[1]*dim(right_crop_area)[2]  
  
  # logic to check whether axis value is found or not
  if (sum(left_crop_area)/total_left_area > ax_val_th || sum(right_crop_area)/total_right_area > ax_val_th) {
    
    #logic to find ticks
    for (i in 3:1) {
      # logic to check ticks outside (cartesian coordinate 4)
      if (initAxesVec[3]+i*px_fact_height <= dim(img_gray)[2]) {
        p_coord_points <- which(!img_gray[,(initAxesVec[3]+(i*px_fact_height))])
        p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1]+(2*px_fact_width)) & (p_coord_points<initAxesVec[2]-(2*px_fact_width))]
      }
      if (!length(p_coord_points)) {
        # logic to check ticks inside (cartesian coordinate 1)
        p_coord_points <- which(!img_gray[,(initAxesVec[3]-(i*px_fact_height))])
        p_coord_points <- p_coord_points[(p_coord_points>initAxesVec[1]+(2*px_fact_width)) & (p_coord_points<initAxesVec[2]-(2*px_fact_width))]
      }
      if (length(p_coord_points) && length(p_coord_points[!diff(p_coord_points) == 1]))
        break
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
      peak_pos_index <- which(seq_diff[p_peak_pos] >= p_peak-(3*px_fact_width) & seq_diff[p_peak_pos] <= p_peak+(3*px_fact_width))
      peak_pos <- p_peak_pos[peak_pos_index]
      range <- 1:(length(peak_pos))
      if (!axis_pos == "start")
        range <- length(peak_pos):1
      for(len in range) {
        p_coord_point <- p_coord_points[peak_pos[len]]
        left_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "left", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
        right_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "right", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
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
        left_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "left", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
        right_crop_area <- getAxisValueCropArea(img_gray, p_coord_point, initAxesVec[3], axis = "x", section = "right", prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height)
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
getAxisValueCropArea <- function(img_gray, w, h, axis = "x", section = "top", common_point = FALSE, prox_dist = prox_dist, px_fact_width = px_fact_width, px_fact_height = px_fact_height) {
  width <- dim(img_gray)[1]
  height <- dim(img_gray)[2]
  wVec <- c(1,width)
  hVec <- c(1,height)
  
  if (!section %in% c("bottom", "top", "left", "right"))
    stop("section value is incorrect")
  
  if (is.null(w) || is.null(h))
    return(0)
  
  if (is.na(prox_dist))
    prox_dist <- 5
  
  if (section == "top" || section == "left") {
    if (axis == "y") {
      width_min <- w - prox_dist - 15*px_fact_width
      width_max <- w-(prox_dist)
      height_min <- h-(3*px_fact_height)
      height_max <- h-px_fact_height
    } else if (axis == "x") {
      lb <- 5
      if (common_point)
        lb <- 10
      width_min <- w-(lb*px_fact_width)
      width_max <- w
      height_min <- h+(prox_dist)
      height_max <- h+(prox_dist + 15)
    }
  } else if (section == "bottom" || section == "right") {
    if (axis == "y") {
      width_min <- w - prox_dist - 15*px_fact_width
      width_max <- w-(prox_dist)
      height_min <- h+px_fact_height
      height_max <- h+(3*px_fact_height)
    } else if (axis == "x") {
      width_min <- w+px_fact_width
      width_max <- w+(5*px_fact_width)
      height_min <- h+(prox_dist)
      height_max <- h+(prox_dist + 15*px_fact_height)
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
  img_draw <- EBImage::drawCircle(img, lst$x1_coord[1], lst$x1_coord[2], 3, col = 2, fill = T)
  img_draw <- EBImage::drawCircle(img_draw, lst$y1_coord[1], lst$y1_coord[2], 3, col = 2, fill = T)
  img_draw <- EBImage::drawCircle(img_draw, lst$x2_coord[1], lst$x2_coord[2], 3, col = 2, fill = T)
  img_draw <- EBImage::drawCircle(img_draw, lst$y2_coord[1], lst$y2_coord[2], 3, col = 2, fill = T)
  display(img_draw)
}

