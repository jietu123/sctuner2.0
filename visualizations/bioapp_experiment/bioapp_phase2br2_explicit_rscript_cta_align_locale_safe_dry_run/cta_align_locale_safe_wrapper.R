
# Locale-safe CTA_align wrapper generated for Phase 2B-R2.
# Original CTA_align.R is not modified.
# Only the non-ASCII Centroid.X/Centroid.Y column access is generalized.
# The geometric alignment remains equivalent:
#   CTA x = Centroid X / pixel_size
#   CTA y_reverse = raw_image_height - Centroid Y / pixel_size
#   spot square = full-resolution imagecol/imagerow +/- spot radius
# For performance, Phase 2B-R2 computes this geometry with a vectorized
# Python KD-tree implementation after exporting Seurat spot coordinates.
