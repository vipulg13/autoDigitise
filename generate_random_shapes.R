library(ggplot2)
library(ggrepel)

set.seed(2)
n_pixels <- 512
# create a random  x and y value in the range (0:100)
df <- data.frame(
  x <- as.integer(runif(2, 0, n_pixels)),
  y <- as.integer(runif(2, 0, n_pixels))
)

# ppi = 72
# image_width_in_inches = n_pixels/72
# image_width_in_mm = 25.4*n_pixels/72
# pixels_per_mm = n_pixels/image_width_in_mm
# pixel_size = 1/pixels_per_mm  

ggplot(df, aes(x, y)) +
  geom_point(size = 1, color = "black") +
  # ggplot has (0,0) in the bottom left corner, so reverse the y value
  geom_text_repel(label = paste0("(", x, ",", n_pixels - y, ")"), size = 4) +
  scale_x_continuous(limits = c(0, n_pixels), expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, n_pixels), expand = c(0, 0)) +
  coord_fixed() +
  theme_void() +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# save a 100 x 100 pixel image with ggsave
ggsave(
  "point.png",
  plot = last_plot(),
  units = "px",
  width = n_pixels,
  height = n_pixels,
  dpi = 72
)
img <- EBImage::readImage("C:/Users/guptav/autodigitise/autodigitise/point.png")
EBImage::display(img)
