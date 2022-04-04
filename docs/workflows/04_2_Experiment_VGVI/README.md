Experiment - Viewshed
================
Sebastian Brinkmann
31 3 2022

## Libraries

For this analysis we used our AGILE R package, which has been built for
this paper specifically. The functions used in this analysis are the
same as in the GVI R package, that is being used for calculating VGVI.

``` r
library(protoVS)
library(terra)
library(sf)
library(tidyverse)
```

## Elevation Data and Observer Locations

For the viewshed experiment we used LiDAR derived a Digital Surface
Model (DSM) and Digital Terrain Model (DTM) with 1m spatial resolution
of the Vancouver metropolitan area (Natural Resources Canada 2019).

``` r
DTM <- rast("data/DTM_Vancouver_1m.tif")
DSM <- rast("data/DSM_Vancouver_1m.tif")
```

Publicly available Land Use and Land Cove (LULC) data has been acquired
by Metro Vancouver (31.11.2019) at 2m resolution and reclassified to a
binary greenness raster (0=no-green; 1=green).

``` r
LULC <- rast("data/LULC_2014.tif")

rcl_mat <- matrix(c(1, 6, 0,    # no vegetation
                    6, 11, 1,   # vegetation
                    11, 14, 0), # no vegetation
                  ncol = 3, byrow = TRUE)

GS <- terra::classify(LULC, rcl = rcl_mat, include.lowest = TRUE, right = FALSE)
```

To demonstrate runtime for a large-scale greenness visibility study we
estimated city-wide VGVI at 5 m intervals, resulting in 17,329,345
observer locations. We will remove observer locations that are located
inside buildings or on water based on the LULC raster.

``` r
# 5 m intervals
DSM_5 <- aggregate(DSM, 5)
observer <- xyFromCell(DSM_5, which(!is.na(DSM_5[]))) %>% 
  sfheaders::sf_point() %>% 
  st_sf(crs = crs(DSM_5)) %>%
  dplyr::rename(geom = geometry)

# Buildings and waters
lulc_build_water <- (LULC %in% c(1, 12)) %>% 
  terra::mask(LULC)

# intersect with LULC
observer_lulc_vals <- unlist(terra::extract(lulc_build_water, y = st_coordinates(observer)),
                             use.names = FALSE)
valid_obs <- which(observer_lulc_vals == 0)
observer <- observer[valid_obs,]
```

## VGVI calculation

As suggested by Labib et al. (2021) we assumed eye-level observer height
of 1.7 m and used an exponential distance decay weight function with a
viewing distance threshold of 800 m. Total runtime time using 20 CPU
threads was 20.3 hours, computation time per observer was 84
milliseconds.

``` r
vgvi_1 <- vgvi_from_sf(observer = observer,
                       dsm_rast = DSM, dtm_rast = DTM, greenspace_rast = GS,
                       max_distance = 800, observer_height = 1.7,
                       m = 1, b = 3, mode = "exponential", cores = 20, progress = TRUE)
```

To compare our algorithm against other results by Labib et al. (2021)
and Cimburova & Blumentrath (2022), respectively, we also computed VGVI
with 5 m resolution elevation models and 100 m viewing distance, and 1 m
resolution elevation models and 100 m viewing distance.

### Labib et al. (2021)

Total runtime: 73 minutes  
Computation time per point: 5 milliseconds

``` r
DSM_5 <- aggregate(DSM, 5)
DTM_5 <- aggregate(DTM, 5)

vgvi_2 <- vgvi_from_sf(observer = observer,
                       dsm_rast = DSM_5, dtm_rast = DTM_5, greenspace_rast = GS,
                       max_distance = 800, observer_height = 1.7,
                       m = 1, b = 3, mode = "exponential", cores = 20, progress = TRUE)
```

### Cimburova & Blumentrath (2022)

Total runtime: 33 minutes  
Computation time per point: 2 milliseconds

``` r
vgvi_3 <- vgvi_from_sf(observer = observer,
                       dsm_rast = DSM, dtm_rast = DTM, greenspace_rast = GS,
                       max_distance = 100, observer_height = 1.7,
                       m = 1, b = 3, mode = "exponential", cores = 20, progress = TRUE)
```

## Bibliography

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-MetroVancouver.31.11.2019" class="csl-entry">

Metro Vancouver. 31.11.2019. “Land Cover Classification 2014 - 2m LiDAR
(Raster).” <http://www.metrovancouver.org/data>.

</div>

<div id="ref-NaturalResourcesCanada.2019" class="csl-entry">

Natural Resources Canada. 2019. “High Resolution Digital Elevation Model
(HRDEM) - CanElevation Serie.”
<https://open.canada.ca/data/en/dataset/957782bf-847c-4644-a757-e383c0057995>.

</div>

</div>