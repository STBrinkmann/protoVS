#' @title Viewshed Benchmark
#' @description Compares computation times of the traditional against our novel viewshed algorithm.
#' At each distance level a binary viewshed for each observer point will be calculated based on a Digital Surface Model (DSM).
#'
#' @param observer object of class \code{sf}; Observer locations from which viewsheds will be calculated
#' @param dsm_rast object of class \code{\link[terra]{rast}}; \code{\link[terra]{rast}} of the DSM
#' @param dtm_rast object of class \code{\link[terra]{rast}}; \code{\link[terra]{rast}} of the DTM
#' @param max_distance_vec numeric vector; Buffer distances to calculate the viewsheds
#' @param observer_height numeric > 0; Height of the observer (e.g. 1.7 meters)
#' @param sample_size numeric > 0; Number of times to evaluate the viewshed calculation at each distance level. 
#' @param ncores integer; he number of cores to use, i.e. at most how many child processes will be run simultaneously
#' @param progress logical; Show progress bar and computation time?
#'
#' @return object of class \code{\link[tibble]{tibble}};
#' @export
#' 
#' @importFrom magrittr %>%
#' @importFrom sf st_buffer
#' @importFrom sf st_coordinates
#' @importFrom terra extract
#' @importFrom terra crop
#' @importFrom terra mask
#' @importFrom terra vect
#' @importFrom terra rowFromY
#' @importFrom terra colFromX
#' @importFrom terra values
#' @importFrom terra rast
#' @importFrom raster raster
#' @importFrom dplyr tibble
#' @importFrom dplyr add_row
#' @useDynLib protoVS, .registration = TRUE

viewshed_comparison <- function(observer, dsm_rast, dtm_rast,
                                max_distance_vec = 1:5, observer_height = 1.7,
                                sample_size = 10, ncores = 1, progress = FALSE) {
  #### 1. Check input ####
  # observer
  valid_sf_types <- c("POINT", "MULTIPOINT")
  if (!is(observer, "sf")) {
    stop("observer must be a sf object")
  } else if (is.null(sf::st_crs(observer)$units) || sf::st_crs(observer)$units_gdal == "degree") {
    stop("observer CRS unit must not be degree")
  } else if (!as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) %in% valid_sf_types) {
    stop("observer has no valid geometry")
  } else if (as.character(sf::st_geometry_type(observer, by_geometry = FALSE)) == "MULTIPOINT") {
    observer <- sf::st_cast(observer, "POINT")
  }
  rm(valid_sf_types)
  
  # dsm_rast
  if (!is(dsm_rast, "SpatRaster")) {
    stop("dsm_rast needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(dsm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dsm_rast needs to have the same CRS as observer")
  } else if(dsm_rast@ptr$res[1] != dsm_rast@ptr$res[2]) {
    stop("dsm_rast: x and y resolution must be equal.\nSee https://github.com/STBrinkmann/GVI/issues/1")
  }
  
  # dtm_rast
  if (!is(dtm_rast, "SpatRaster")) {
    stop("dtm_rast needs to be a SpatRaster object!")
  } else if (sf::st_crs(terra::crs(dtm_rast))$epsg != sf::st_crs(observer)$epsg) {
    stop("dtm_rast needs to have the same CRS as observer")
  }
  
  # max_distance_vec
  max_distance_vec <- round(max_distance_vec, digits = 0)
  if (min(max_distance_vec) < min(res(dsm_rast))) {
    stop("all distances from max_distance_vec must be greater than the DSM resolution")
  }
  
  
  if(progress) {
    message("Preprocessing:")
    pb = txtProgressBar(min = 0, max = 4, initial = 0, style = 3)
  }
  if (progress) setTxtProgressBar(pb, 1)
  
  output <- dplyr::tibble(
    method = as.character(),
    distance = as.numeric(),
    time_mean = as.numeric(),
    time_sd = as.numeric()
  )
  
  # Coordinates of start point
  x0 <- sf::st_coordinates(observer)[,1]
  y0 <- sf::st_coordinates(observer)[,2]
  
  terraOptions(progress = 0)
  height0 <- unlist(terra::extract(dtm_rast, cbind(x0, y0)), use.names = FALSE) + observer_height
  
  if (progress) setTxtProgressBar(pb, 2)
  
  # Mask to AOI
  aoi <- observer %>% 
    sf::st_bbox() %>% 
    sf::st_as_sfc() %>% 
    sf::st_buffer(max(max_distance_vec)) %>% 
    terra::vect()
  
  dsm_masked <- dsm_rast %>% 
    terra::crop(aoi) %>% 
    terra::mask(aoi)
  
  if (progress) setTxtProgressBar(pb, 3)
  
  # Start row/col
  r0 <- terra::rowFromY(dsm_masked, y0)
  c0 <- terra::colFromX(dsm_masked, x0)
  
  # Convert dsm to vector and raster
  dsm_vec <- terra::values(dsm_masked, mat = FALSE)
  dsm_cpp_rast <- terra::rast(dsm_masked) %>% raster::raster()
  
  if (progress) setTxtProgressBar(pb, 4)
  if (progress) close(pb)
  
  terraOptions(progress = 3)
  invisible(gc())
  
  if(progress) {
    message("Calculating Viewsheds:")
    pb = txtProgressBar(min = 0, max = (sample_size*length(max_distance_vec))^2, initial = 0, style = 3) 
  }
  for(i in seq_along(max_distance_vec)){
    max_distance <- max_distance_vec[i]
    
    
    # Run viewshed analysis
    old <- list()
    new <- list()
    new_par <- list()
    for(t in 1:sample_size){
      old[[t]] <- viewshed_test_old(dsm = dsm_cpp_rast, dsm_values = dsm_vec,
                                    x0 = c0, y0 = r0, h0 = height0,
                                    radius = max_distance)
      
      new[[t]] <- viewshed_test_new(dsm = dsm_cpp_rast, dsm_values = dsm_vec,
                                    x0 = c0, y0 = r0, h0 = height0,
                                    radius = max_distance, ncores = 1)
      
      if(ncores > 1){
        new_par[[t]] <- viewshed_test_new(dsm = dsm_cpp_rast, dsm_values = dsm_vec,
                                          x0 = c0, y0 = r0, h0 = height0,
                                          radius = max_distance, ncores = ncores)
      }
      
      if (progress) setTxtProgressBar(pb, (((i-1)*sample_size) + t)^2)
    }
    
    # Old Method
    output <- output %>% 
      dplyr::add_row(
        dplyr::tibble(
          method = "Old",
          distance = max_distance,
          time_mean = mean(unlist(old, use.names = FALSE)),
          time_sd = sd(unlist(old, use.names = FALSE))
        )
      )
    
    # New Method
    output <- output %>% 
      dplyr::add_row(
        dplyr::tibble(
          method = "New",
          distance = max_distance,
          time_mean = mean(unlist(new, use.names = FALSE)),
          time_sd = sd(unlist(new, use.names = FALSE))
        )
      )
    
    if (ncores > 1){
      # New Method - Parallel
      output <- output %>% 
        dplyr::add_row(
          dplyr::tibble(
            method = paste0("New (", ncores, " cores)"),
            distance = max_distance,
            time_mean = mean(unlist(new_par, use.names = FALSE)),
            time_sd = sd(unlist(new_par, use.names = FALSE))
          )
        )
    }
  }
  rm(dsm_masked, dsm_vec, dsm_cpp_rast); invisible(gc())
  if (progress) close(pb)
  
  return(output)
}