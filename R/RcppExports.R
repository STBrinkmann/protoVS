# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

IDW_cpp <- function(rast, x, sf_x, sf_y, sf_z, n, b, radius, mode = 1L, na_only = FALSE, ncores = 1L, display_progress = FALSE) {
    .Call(`_protoVS_IDW_cpp`, rast, x, sf_x, sf_y, sf_z, n, b, radius, mode, na_only, ncores, display_progress)
}

LoS_reference <- function(x0_ref, y0_ref, r, nc_ref) {
    .Call(`_protoVS_LoS_reference`, x0_ref, y0_ref, r, nc_ref)
}

VGVI_cpp <- function(dsm, dsm_values, greenspace, greenspace_values, x0, y0, h0, radius, fun, m, b, ncores = 1L, display_progress = FALSE) {
    .Call(`_protoVS_VGVI_cpp`, dsm, dsm_values, greenspace, greenspace_values, x0, y0, h0, radius, fun, m, b, ncores, display_progress)
}

viewshed_test_old <- function(dsm, dsm_values, x0, y0, h0, radius) {
    .Call(`_protoVS_viewshed_test_old`, dsm, dsm_values, x0, y0, h0, radius)
}

viewshed_test_new <- function(dsm, dsm_values, x0, y0, h0, radius, ncores = 1L) {
    .Call(`_protoVS_viewshed_test_new`, dsm, dsm_values, x0, y0, h0, radius, ncores)
}

viewshed_cpp <- function(dsm, dsm_values, x0, y0, h0, radius, ncores = 1L, display_progress = FALSE) {
    .Call(`_protoVS_viewshed_cpp`, dsm, dsm_values, x0, y0, h0, radius, ncores, display_progress)
}

viewshed_cpp2 <- function(dsm, dsm_values, x0, y0, h0, radius, ncores = 1L, display_progress = FALSE) {
    .Call(`_protoVS_viewshed_cpp2`, dsm, dsm_values, x0, y0, h0, radius, ncores, display_progress)
}

