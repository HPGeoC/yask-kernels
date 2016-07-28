/*
 * 4-D loop code.
 * Generated automatically from the following pseudo-code:
 *
 * omp loop(nv,xv,yv,zv) { calc(halo(t)); }
 *
 */

 // Number of iterations to get from begin_nv to (but not including) end_nv, stepping by step_nv.
 const idx_t num_nv = ((end_nv - begin_nv) + (step_nv - 1)) / step_nv;

 // Number of iterations to get from begin_xv to (but not including) end_xv, stepping by step_xv.
 const idx_t num_xv = ((end_xv - begin_xv) + (step_xv - 1)) / step_xv;

 // Number of iterations to get from begin_yv to (but not including) end_yv, stepping by step_yv.
 const idx_t num_yv = ((end_yv - begin_yv) + (step_yv - 1)) / step_yv;

 // Number of iterations to get from begin_zv to (but not including) end_zv, stepping by step_zv.
 const idx_t num_zv = ((end_zv - begin_zv) + (step_zv - 1)) / step_zv;

 // Number of iterations in loop collapsed across nv, xv, yv, zv dimensions.
 const idx_t num_nv_xv_yv_zv = (idx_t)num_nv * (idx_t)num_xv * (idx_t)num_yv * (idx_t)num_zv;

 // Computation loop.

 // Distribute iterations among OpenMP threads.
_Pragma("omp parallel for schedule(dynamic,1)")
 for (idx_t loop_index_nv_xv_yv_zv = 0; loop_index_nv_xv_yv_zv < num_nv_xv_yv_zv; loop_index_nv_xv_yv_zv++) {

 // Zero-based, unit-stride index var for nv.
 idx_t index_nv = loop_index_nv_xv_yv_zv / (num_xv*num_yv*num_zv);

 // Zero-based, unit-stride index var for xv.
 idx_t index_xv = (loop_index_nv_xv_yv_zv / (num_yv*num_zv)) % num_xv;

 // Zero-based, unit-stride index var for yv.
 idx_t index_yv = (loop_index_nv_xv_yv_zv / (num_zv)) % num_yv;

 // Zero-based, unit-stride index var for zv.
 idx_t index_zv = (loop_index_nv_xv_yv_zv) % num_zv;

 // This value of index_nv covers nv from start_nv to stop_nv-1.
 const idx_t start_nv = begin_nv + (index_nv * step_nv);
 const idx_t stop_nv = std::min(start_nv + step_nv, end_nv);

 // This value of index_xv covers xv from start_xv to stop_xv-1.
 const idx_t start_xv = begin_xv + (index_xv * step_xv);
 const idx_t stop_xv = std::min(start_xv + step_xv, end_xv);

 // This value of index_yv covers yv from start_yv to stop_yv-1.
 const idx_t start_yv = begin_yv + (index_yv * step_yv);
 const idx_t stop_yv = std::min(start_yv + step_yv, end_yv);

 // This value of index_zv covers zv from start_zv to stop_zv-1.
 const idx_t start_zv = begin_zv + (index_zv * step_zv);
 const idx_t stop_zv = std::min(start_zv + step_zv, end_zv);

 // Computation.
  calc_halo(context, t, start_nv, start_xv, start_yv, start_zv, stop_nv, stop_xv, stop_yv, stop_zv);
 }
// End of generated code.
