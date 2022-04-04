#include <Rcpp.h>
#include <chrono>
#include "rsinfo.h"
#include "rasterutils.h"
#include "bresenham.h"

using namespace Rcpp;


#ifdef _OPENMP
#include <omp.h>
#endif
// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>


// [[Rcpp::export]]
double viewshed_test_old(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                         const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                         const Rcpp::NumericVector &h0, const int radius)
{
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Cells from x0,y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int l = r+1;
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  
  // Main loop
  for(int k = 0; k < input_cells.size(); k++){
    if(h0[k] > dsm_values[input_cells[k]]){
      // Viewshed: x0/y0 is always visible
      Rcpp::IntegerVector viewshed(nc_ref*nr_ref, NA_INTEGER);
      viewshed[c0_ref] = input_cells[k]+1;
      
      // 1. Reference Bresenham Lines for the first eigth of the circle
      const Rcpp::IntegerMatrix bh_mat = bresenham_map(x0_ref, y0_ref, r, nc_ref);
      
      // 2. Reference matrix (los_ref_mat): Project reference Bresenham Lines to all perimeter cells 
      Rcpp::IntegerMatrix los_ref_mat = Rcpp::IntegerMatrix(l*8 - 8, r);
      for(int i=0; i<l; i++){
        const Rcpp::IntegerMatrix bs_xy = colRowFromCell(bh_mat( i,_ ), nc_ref);
        
        // Fill los_ref_mat with projected Bresenham Lines, using Bresenham's Midpoint algorithm
        for(int j=0; j<r; j++){
          
          if ( Rcpp::IntegerVector::is_na(bs_xy( j,0 )) || Rcpp::IntegerVector::is_na(bs_xy( j,1 )) ) {
            los_ref_mat( i,j )             = NA_INTEGER; 
            los_ref_mat( 2*(l-1) + i,j )   = NA_INTEGER;
            los_ref_mat( 4*(l-1) + i,j )   = NA_INTEGER;
            los_ref_mat( 6*(l-1) + i,j )   = NA_INTEGER;
            
            if( i != 0 && i != (l-1) ) {
              los_ref_mat( 2*(l-1) - i,j ) = NA_INTEGER;
              los_ref_mat( 4*(l-1) - i,j ) = NA_INTEGER;
              los_ref_mat( 6*(l-1) - i,j ) = NA_INTEGER;
              los_ref_mat( 8*(l-1) - i,j ) = NA_INTEGER;
            }
          } else {
            const int x = bs_xy( j,0 )-x0_ref;
            const int y = bs_xy( j,1 )-x0_ref;
            
            los_ref_mat( i,j )             = (y+y0_ref)*nc_ref + (x+x0_ref); 
            los_ref_mat( 2*(l-1) + i,j )   = (-x+y0_ref)*nc_ref + (y+x0_ref);
            los_ref_mat( 4*(l-1) + i,j )   = (-y+y0_ref)*nc_ref + (-x+x0_ref);
            los_ref_mat( 6*(l-1) + i,j )   = (x+y0_ref)*nc_ref + (-y+x0_ref);
            
            if( i != 0 && i != (l-1) ) {
              los_ref_mat( 2*(l-1) - i,j ) = (x+y0_ref)*nc_ref + (y+x0_ref);
              los_ref_mat( 4*(l-1) - i,j ) = (-y+y0_ref)*nc_ref + (x+x0_ref);
              los_ref_mat( 6*(l-1) - i,j ) = (-x+y0_ref)*nc_ref + (-y+x0_ref);
              los_ref_mat( 8*(l-1) - i,j ) = (y+y0_ref)*nc_ref + (-x+x0_ref);
            }
          }
        }
      }
      
      
      // 3: Tangents Map
      const Rcpp::NumericVector tan_vec = tangentsMap(x0_o[k], y0_o[k], h0[k], ras.nrow, ras.ncol, r, dsm_values);
      
      // 4: Visibility
      Rcpp::LogicalVector visibility_vec(nc_ref*nr_ref); 
      
      for(int j = 0; j < (l*8 - 8); j++){ 
        // Compute tangents of LoS_vec and update visibility_vec
        LoS(los_ref_mat, tan_vec, visibility_vec, j, r);
      }
      
      // 5: Project reference cells to actual cell value
      int x = input_cells[k] - c0_ref - r*(ras.ncol-nc_ref);
      for(int j = 0; j < visibility_vec.size(); j++){
        if(visibility_vec[j]){
          const int cell = x + j + trunc(j/nc_ref)*(ras.ncol-nc_ref);
          const int row = cell - (trunc(cell/ras.ncol) * ras.ncol);
          const int drow = abs(row-x0_o[0]);
          
          viewshed[j] = (cell<0 || cell > ras.ncell ||
            Rcpp::NumericVector::is_na(dsm_values[cell]) ||
            drow>r) ? NA_INTEGER : (cell+1);
        }
      }
      
      Rcpp::IntegerVector viewshed_na = Rcpp::na_omit(viewshed);
    }
  }
    
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto time = end_time - start_time;
  
  return ((time/std::chrono::microseconds(1)) / 1000.0) / input_cells.size();
}


// [[Rcpp::export]]
double viewshed_test_new(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                         const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                         const Rcpp::NumericVector &h0, const int radius,
                         const int ncores=1)
{
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Cells from x0, y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  //Rcpp::IntegerVector output(input_cells.size(), NA_INTEGER);
  
  // Protoptype of Line of Sight (LoS) paths:
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  const Rcpp::IntegerVector los_start = shared_LoS(r, los_ref_vec);
  
  // Main loop
#if defined(_OPENMP)
  omp_set_num_threads(ncores);
#pragma omp parallel for //shared(output)
#endif
  for(int k = 0; k < input_cells.size(); ++k){
    // Eye-level height at x0/y0 > DSM height?
    const int this_input_cell = input_cells[k];
    if(h0[k] > dsm_values[this_input_cell]){
      
      // Viewshed
      std::vector<int> viewshed;
      viewshed.reserve(nc_ref*nr_ref);
      
      // Viewshed at x0/y0 is always visible
      viewshed.push_back(this_input_cell+1);
      
      // Counter for visible points
      size_t n = 1;
      
      // Parameter for projecting reference cell to true cell values
      const int x = this_input_cell - c0_ref - r*(ras.ncol-nc_ref);
      
      // Vector for storing max tangent from previous LoS
      std::vector<double> max_tan_vec(r, -9999.0);
      
      
      // Iterate over all LoS paths
      for(int i = 0; i < (r*8); ++i){
        // Re-use tangents calculated from prior LoS or assign -9999.0
        const int k_i = los_start[i];
        double max_tan = (k_i > 1) ? max_tan_vec[k_i-1] : -9999.0;
        
        // Iterate over all cells of this LoS path starting at k_i
        for(int j = k_i; j < r; ++j){
          
          // This LoS path cell
          const int los_ref_cell = los_ref_vec[i*r + j];
          
          if(!Rcpp::IntegerVector::is_na(los_ref_cell)){ //!Rcpp::IntegerVector::is_na(los_ref_cell)
            
            // Project reference cell to true cell value
            const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
            
            // DSM height at this LoS path cell
            const double h_cell = dsm_values[cell];
            
            // Test if this LoS path cell is within DSM raster
            const int row = trunc(cell/ras.ncol);
            const int col = cell - (row * ras.ncol);
            const int dcol = abs(col-x0_o[k]);
            
            if(!(cell<0 || cell > ras.ncell || h_cell != h_cell || dcol>r)){
              // Compute tangent of x0/y0 (observer location) and this LoS path cell
              const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
              const double this_tan = (h_cell - h0[k]) / (distance_traveled);
              
              // Update viewshed and max tangent
              if(this_tan > max_tan){
                max_tan = this_tan;
                
                viewshed.push_back(cell+1);
                ++n;
              }
            }
            
            max_tan_vec[j] = max_tan;
          } else {
            break;
          }
        }
      }
      viewshed.resize(n);
      //output[k] = viewshed.size();
    }
  }
    
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto time = end_time - start_time;
  
  return ((time/std::chrono::microseconds(1)) / 1000.0) / input_cells.size();
}


std::vector<double> LoS_tangent(const std::vector<int> &shared_LoS,
                                const Rcpp::NumericVector &dsm_values,
                                const int p, const int nc_ref, const int nc_ras, const int ncell_ras,
                                const int x0, const int y0, const int h0, const int r){
  const int n = shared_LoS.size();
  std::vector<double> out(n, -9999.0);
  
  for(int i = 0; i < n; ++i){
    // 1. Get LoS cell
    const int los_cell = shared_LoS[i];
    
    // 2. Project LoS cell
    const int cell = p + los_cell + trunc(los_cell/nc_ref)*(nc_ras-nc_ref);
    
    // 3. DSM height at this LoS path cell
    const double h_cell = dsm_values[cell];
    
    // 4. Test if this LoS path cell is within DSM raster
    const int row = trunc(cell/nc_ras);
    const int col = cell - (row * nc_ras);
    const int dcol = abs(col-x0);
    
    if(!(cell<0 || cell > ncell_ras || h_cell != h_cell || dcol>r)){
      // Compute tangent of x0/y0 (observer location) and this LoS path cell
      const double distance_traveled = sqrt((x0 - col)*(x0 - col) + (y0 - row)*(y0 - row));
      const double this_tan = (h_cell - h0) / (distance_traveled);
      
      // Update viewshed and max tangent
      out[i] = this_tan;
    }
  }
  
  return out;
}



// [[Rcpp::export]]
std::vector<int> viewshed_cpp(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                              const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                              const Rcpp::NumericVector &h0, const int radius,
                              const int ncores=1, const bool display_progress=false)
{
  // Cells from x0, y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)round(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  std::vector<int> output;
  
  // Protoptype of Line of Sight (LoS) paths:
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  const Rcpp::IntegerVector los_start = shared_LoS(r, los_ref_vec);
  
  // Progress bar
  Progress pb(input_cells.size(), display_progress);
  
  // Main loop
  for(int k = 0; k < input_cells.size(); ++k){
    
    const int this_input_cell = input_cells[k];
    
    // Viewshed
    std::vector<int> viewshed(nc_ref*nr_ref, NA_INTEGER);
    
    // Viewshed at x0/y0 is always visible
    viewshed[c0_ref] = input_cells[k]+1;
    
    if (!pb.is_aborted()) {
      Progress::check_abort();
      
      // Eye-level height at x0/y0 > DSM height?
      if(h0[k] > dsm_values[this_input_cell]){
        
        // Parameter for projecting reference cell to true cell values
        const int x = this_input_cell - c0_ref - r*(ras.ncol-nc_ref);
        
        // Vector for storing max tangent from previous LoS
        std::vector<double> max_tan_vec(r, -9999.0);
        
        // Iterate over all LoS paths
        for(int i = 0; i < (r*8); ++i){
          // Re-use tangents calculated from prior LoS or assign -9999.0
          const int k_i = los_start[i];
          double max_tan = (k_i > 1) ? max_tan_vec[k_i-1] : -9999.0;
          
          // Iterate over all cells of this LoS path starting at k_i
          for(int j = k_i; j < r; ++j){
            
            // This LoS path cell
            const int los_ref_cell = los_ref_vec[i*r + j];
            
            if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
              
              // Project reference cell to true cell value
              const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
              
              // DSM height at this LoS path cell
              const double h_cell = dsm_values[cell];
              
              // Test if this LoS path cell is within DSM raster
              const int row = trunc(cell/ras.ncol);
              const int col = cell - (row * ras.ncol);
              const int dcol = abs(col-x0_o[k]);
              
              if(!(cell<0 || cell > ras.ncell || h_cell != h_cell || dcol>r)){
                // Compute tangent of x0/y0 (observer location) and this LoS path cell
                const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
                const double this_tan = (h_cell - h0[k]) / (distance_traveled);
                
                // Update viewshed and max tangent
                if(this_tan > max_tan){
                  max_tan = this_tan;
                  viewshed[los_ref_cell] = cell+1;
                }
              }
              
              max_tan_vec[j] = max_tan;
            } else {
              break;
            }
          }
        }
      }
    }
    
    viewshed.erase(std::remove_if(std::begin(viewshed),
                                  std::end(viewshed),
                                  [](const int& value) { return Rcpp::IntegerVector::is_na(value); }),
                                  std::end(viewshed));
    output = viewshed;
    pb.increment();
  }
  
  return output;
}
// [[Rcpp::export]]
std::vector<double> viewshed_cpp2(Rcpp::S4 &dsm, const Rcpp::NumericVector &dsm_values,
                               const Rcpp::IntegerVector &x0, const Rcpp::IntegerVector &y0,
                               const Rcpp::NumericVector &h0, const int radius,
                               const int ncores=1, const bool display_progress=false)
{
  auto start_time = std::chrono::high_resolution_clock::now();
  
  // Cells from x0, y0
  Rcpp::IntegerVector x0_o = x0-1;
  Rcpp::IntegerVector y0_o = y0-1;
  const Rcpp::IntegerVector input_cells = Rcpp::na_omit(cellFromColRowSensitive(dsm, x0_o, y0_o));
  
  // Basic raster information
  const RasterInfo ras(dsm);
  
  // Parameters
  const int r = (int)round(radius/ras.res);
  const int nc_ref = (2*r)+1, nr_ref = (2*r)+1;
  const int x0_ref = r, y0_ref = x0_ref;
  const int c0_ref = y0_ref*nc_ref + x0_ref;
  
  // Output vector
  std::vector<double> output;
  
  // Protoptype of Line of Sight (LoS) paths:
  // Will be used as a reference for all input points
  const Rcpp::IntegerVector los_ref_vec = LoS_reference(x0_ref, y0_ref, r, nc_ref);
  const Rcpp::IntegerVector los_start = shared_LoS(r, los_ref_vec);
  
  auto middle_time = std::chrono::high_resolution_clock::now();
  
  // Progress bar
  Progress pb(input_cells.size(), display_progress);
  
  // Main loop
  for(int k = 0; k < input_cells.size(); ++k){
    
    const int this_input_cell = input_cells[k];
    
    // Viewshed
    std::vector<int> viewshed(nc_ref*nr_ref, NA_INTEGER);
    
    // Viewshed at x0/y0 is always visible
    viewshed[c0_ref] = input_cells[k]+1;
    
    if (!pb.is_aborted()) {
      Progress::check_abort();
      
      // Eye-level height at x0/y0 > DSM height?
      if(h0[k] > dsm_values[this_input_cell]){
        
        // Parameter for projecting reference cell to true cell values
        const int x = this_input_cell - c0_ref - r*(ras.ncol-nc_ref);
        
        // Vector for storing max tangent from previous LoS
        std::vector<double> max_tan_vec(r, -9999.0);
        
        // Iterate over all LoS paths
        for(int i = 0; i < (r*8); ++i){
          // Re-use tangents calculated from prior LoS or assign -9999.0
          const int k_i = los_start[i];
          double max_tan = (k_i > 1) ? max_tan_vec[k_i-1] : -9999.0;
          
          // Iterate over all cells of this LoS path starting at k_i
          for(int j = k_i; j < r; ++j){
            
            // This LoS path cell
            const int los_ref_cell = los_ref_vec[i*r + j];
            
            if(!Rcpp::IntegerVector::is_na(los_ref_cell)){
              
              // Project reference cell to true cell value
              const int cell = x + los_ref_cell + trunc(los_ref_cell/nc_ref)*(ras.ncol-nc_ref);
              
              // DSM height at this LoS path cell
              const double h_cell = dsm_values[cell];
              
              // Test if this LoS path cell is within DSM raster
              const int row = trunc(cell/ras.ncol);
              const int col = cell - (row * ras.ncol);
              const int dcol = abs(col-x0_o[k]);
              
              if(!(cell<0 || cell > ras.ncell || h_cell != h_cell || dcol>r)){
                // Compute tangent of x0/y0 (observer location) and this LoS path cell
                const double distance_traveled = sqrt((x0_o[k] - col)*(x0_o[k] - col) + (y0_o[k] - row)*(y0_o[k] - row));
                const double this_tan = (h_cell - h0[k]) / (distance_traveled);
                
                // Update viewshed and max tangent
                if(this_tan > max_tan){
                  max_tan = this_tan;
                  viewshed[los_ref_cell] = cell+1;
                }
              }
              
              max_tan_vec[j] = max_tan;
            } else {
              break;
            }
          }
        }
      }
    }
    
    viewshed.erase(std::remove_if(std::begin(viewshed),
                                  std::end(viewshed),
                                  [](const int& value) { return Rcpp::IntegerVector::is_na(value); }),
                                  std::end(viewshed));
    pb.increment();
  }
  
  auto end_time = std::chrono::high_resolution_clock::now();
  
  
  auto time_prototyping = middle_time - start_time;
  auto time_viewshed = end_time - middle_time;
  auto time_total = end_time - start_time;
  
  output.push_back((time_prototyping/std::chrono::microseconds(1)) / 1000.0);
  output.push_back((time_viewshed/std::chrono::microseconds(1)) / 1000.0);
  output.push_back((time_total/std::chrono::microseconds(1)) / 1000.0);

  return output;
}