# BioApp Coordinate Registration Recovery Decision

## Final Decision

`REGISTRATION_RECOVERED_EMBEDDED_IMAGE`

## Embedded Image

BCSA2TumB1 embedded image found: `True`

Embedded image dimensions: `[600, 594, 3]`

The embedded image is stored inside:

`ST_data.RData::bcsa@images[['BCSA2TumB1']]@image`

## Referenced Local Tissue/H&E Image

Referenced image found from path-like metadata: `False`

Referenced image path: `None`

## Scale Factors

Scale factors found: `True`

Scale factor names/values are recorded in `bioapp_image_slot_scale_factor_audit.csv`.

The object contains Seurat scale factors such as `spot`, `fiducial`, `hires`, and `lowres`. The CTA-specific factor `0.3` is documented in README/CTA_align and used to convert Seurat image coordinates to raw microscope coordinates; it is not itself one of the Seurat image-slot scale factors.

## Transform / Inverse Transform

Transform metadata found: `False`

Inverse transform found: `False`

No external inverse transform file was found. However, the embedded Seurat image can be used with the processed BCSA2TumB1 coordinate system using the image slot's lowres scale factor.

## H&E Overlay Recoverable

`True`

## Reason

The Seurat `VisiumV1` image slot contains an embedded raster image and scale factors. The BioApp frozen coordinates are explicitly linked to this image slot. The audit-only bounds check confirms that `imagecol/imagerow` scaled by the image slot lowres factor falls within the embedded image dimensions.

This recovers a Seurat image-slot tissue-background route. It does not recover the original raw microscope TIFF or a full Space Ranger spatial bundle.

## Biological Results

No biological result, endpoint, metric, threshold, label, or conclusion was changed.

## Recommended Next Step

`BioApp Main Figure V3 - Seurat image-slot tissue-background overlay`
