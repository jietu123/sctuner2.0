
suppressPackageStartupMessages({
  library(Seurat)
  library(jsonlite)
})
args <- commandArgs(trailingOnly = TRUE)
rdata <- args[[1]]
out_dir <- args[[2]]
object_name <- args[[3]]
image_name <- args[[4]]
env <- new.env()
load(rdata, envir = env)
obj <- get(object_name, envir = env)
img <- obj@images[[image_name]]
arr <- img@image
h <- dim(arr)[1]
w <- dim(arr)[2]
sf <- as.list(img@scale.factors)
lowres <- as.numeric(sf$lowres)
hires <- as.numeric(sf$hires)
coords <- img@coordinates
x <- coords$imagecol * lowres
y <- coords$imagerow * lowres
bounds <- list(
  embedded_image_source = paste0("ST_data.RData::", object_name, "@images[['", image_name, "']]@image"),
  image_slot_class = class(img),
  embedded_image_dimensions = dim(arr),
  lowres_scale_used = lowres,
  hires_scale_available = hires,
  coordinate_columns_used = c("imagecol", "imagerow"),
  x_min = min(x, na.rm = TRUE),
  x_max = max(x, na.rm = TRUE),
  y_min = min(y, na.rm = TRUE),
  y_max = max(y, na.rm = TRUE),
  image_width = w,
  image_height = h,
  bounds_check_passed = all(x >= 0, x <= w, y >= 0, y <= h, na.rm = TRUE)
)
png(file.path(out_dir, "fig_bioapp_v3_embedded_tissue_background.png"), width = w, height = h)
par(mar = c(0, 0, 0, 0))
plot.new()
plot.window(xlim = c(0, w), ylim = c(h, 0), asp = 1)
rasterImage(as.raster(arr), 0, h, w, 0)
dev.off()
writeLines(jsonlite::toJSON(bounds, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "fig_bioapp_v3_embedded_image_metadata.json"), useBytes = TRUE)
