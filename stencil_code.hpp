// Automatically generated code; do not edit.

////// Implementation of the 'awp' stencil //////

namespace yask {

 ////// Overall stencil-context class //////
struct StencilContext_awp : public StencilContext {

 // Grids.
 const idx_t vel_x_halo_x = 2;
 const idx_t vel_x_halo_y = 2;
 const idx_t vel_x_halo_z = 2;
 Grid_TXYZ* vel_x; // updated by stencil.
 const idx_t vel_y_halo_x = 2;
 const idx_t vel_y_halo_y = 2;
 const idx_t vel_y_halo_z = 2;
 Grid_TXYZ* vel_y; // updated by stencil.
 const idx_t vel_z_halo_x = 2;
 const idx_t vel_z_halo_y = 2;
 const idx_t vel_z_halo_z = 2;
 Grid_TXYZ* vel_z; // updated by stencil.
 const idx_t stress_xx_halo_x = 2;
 const idx_t stress_xx_halo_y = 0;
 const idx_t stress_xx_halo_z = 0;
 Grid_TXYZ* stress_xx; // updated by stencil.
 const idx_t stress_yy_halo_x = 0;
 const idx_t stress_yy_halo_y = 2;
 const idx_t stress_yy_halo_z = 0;
 Grid_TXYZ* stress_yy; // updated by stencil.
 const idx_t stress_zz_halo_x = 0;
 const idx_t stress_zz_halo_y = 0;
 const idx_t stress_zz_halo_z = 2;
 Grid_TXYZ* stress_zz; // updated by stencil.
 const idx_t stress_xy_halo_x = 2;
 const idx_t stress_xy_halo_y = 2;
 const idx_t stress_xy_halo_z = 0;
 Grid_TXYZ* stress_xy; // updated by stencil.
 const idx_t stress_xz_halo_x = 2;
 const idx_t stress_xz_halo_y = 0;
 const idx_t stress_xz_halo_z = 2;
 Grid_TXYZ* stress_xz; // updated by stencil.
 const idx_t stress_yz_halo_x = 0;
 const idx_t stress_yz_halo_y = 2;
 const idx_t stress_yz_halo_z = 2;
 Grid_TXYZ* stress_yz; // updated by stencil.
 const idx_t stress_mem_xx_halo_x = 0;
 const idx_t stress_mem_xx_halo_y = 0;
 const idx_t stress_mem_xx_halo_z = 0;
 Grid_TXYZ* stress_mem_xx; // updated by stencil.
 const idx_t stress_mem_yy_halo_x = 0;
 const idx_t stress_mem_yy_halo_y = 0;
 const idx_t stress_mem_yy_halo_z = 0;
 Grid_TXYZ* stress_mem_yy; // updated by stencil.
 const idx_t stress_mem_zz_halo_x = 0;
 const idx_t stress_mem_zz_halo_y = 0;
 const idx_t stress_mem_zz_halo_z = 0;
 Grid_TXYZ* stress_mem_zz; // updated by stencil.
 const idx_t stress_mem_xy_halo_x = 0;
 const idx_t stress_mem_xy_halo_y = 0;
 const idx_t stress_mem_xy_halo_z = 0;
 Grid_TXYZ* stress_mem_xy; // updated by stencil.
 const idx_t stress_mem_xz_halo_x = 0;
 const idx_t stress_mem_xz_halo_y = 0;
 const idx_t stress_mem_xz_halo_z = 0;
 Grid_TXYZ* stress_mem_xz; // updated by stencil.
 const idx_t stress_mem_yz_halo_x = 0;
 const idx_t stress_mem_yz_halo_y = 0;
 const idx_t stress_mem_yz_halo_z = 0;
 Grid_TXYZ* stress_mem_yz; // updated by stencil.
 const idx_t weight_halo_x = 0;
 const idx_t weight_halo_y = 0;
 const idx_t weight_halo_z = 0;
 Grid_XYZ* weight; // not updated by stencil.
 const idx_t tau2_halo_x = 0;
 const idx_t tau2_halo_y = 0;
 const idx_t tau2_halo_z = 0;
 Grid_XYZ* tau2; // not updated by stencil.
 const idx_t anelastic_ap_halo_x = 0;
 const idx_t anelastic_ap_halo_y = 0;
 const idx_t anelastic_ap_halo_z = 0;
 Grid_XYZ* anelastic_ap; // not updated by stencil.
 const idx_t anelastic_as_diag_halo_x = 0;
 const idx_t anelastic_as_diag_halo_y = 0;
 const idx_t anelastic_as_diag_halo_z = 0;
 Grid_XYZ* anelastic_as_diag; // not updated by stencil.
 const idx_t anelastic_xy_halo_x = 0;
 const idx_t anelastic_xy_halo_y = 0;
 const idx_t anelastic_xy_halo_z = 0;
 Grid_XYZ* anelastic_xy; // not updated by stencil.
 const idx_t anelastic_xz_halo_x = 0;
 const idx_t anelastic_xz_halo_y = 0;
 const idx_t anelastic_xz_halo_z = 0;
 Grid_XYZ* anelastic_xz; // not updated by stencil.
 const idx_t anelastic_yz_halo_x = 0;
 const idx_t anelastic_yz_halo_y = 0;
 const idx_t anelastic_yz_halo_z = 0;
 Grid_XYZ* anelastic_yz; // not updated by stencil.
 const idx_t lambda_halo_x = 1;
 const idx_t lambda_halo_y = 1;
 const idx_t lambda_halo_z = 1;
 Grid_XYZ* lambda; // not updated by stencil.
 const idx_t rho_halo_x = 1;
 const idx_t rho_halo_y = 1;
 const idx_t rho_halo_z = 1;
 Grid_XYZ* rho; // not updated by stencil.
 const idx_t mu_halo_x = 1;
 const idx_t mu_halo_y = 1;
 const idx_t mu_halo_z = 1;
 Grid_XYZ* mu; // not updated by stencil.
 const idx_t sponge_halo_x = 0;
 const idx_t sponge_halo_y = 0;
 const idx_t sponge_halo_z = 0;
 Grid_XYZ* sponge; // not updated by stencil.

 // Max halos across all grids.
 const idx_t max_halo_x = 2;
 const idx_t max_halo_y = 2;
 const idx_t max_halo_z = 2;

 // Parameters.
 GenericGrid0d<real_t>* delta_t;
 GenericGrid0d<real_t>* h;

 StencilContext_awp() {
  name = "awp";
  vel_x = 0;
  vel_y = 0;
  vel_z = 0;
  stress_xx = 0;
  stress_yy = 0;
  stress_zz = 0;
  stress_xy = 0;
  stress_xz = 0;
  stress_yz = 0;
  stress_mem_xx = 0;
  stress_mem_yy = 0;
  stress_mem_zz = 0;
  stress_mem_xy = 0;
  stress_mem_xz = 0;
  stress_mem_yz = 0;
  weight = 0;
  tau2 = 0;
  anelastic_ap = 0;
  anelastic_as_diag = 0;
  anelastic_xy = 0;
  anelastic_xz = 0;
  anelastic_yz = 0;
  lambda = 0;
  rho = 0;
  mu = 0;
  sponge = 0;
  delta_t = 0;
  h = 0;
 }

 virtual void allocGrids() {
  gridPtrs.clear();
  eqGridPtrs.clear();
  vel_x = new Grid_TXYZ(dx, dy, dz, vel_x_halo_x + px, vel_x_halo_y + py, vel_x_halo_z + pz, "vel_x");
  gridPtrs.push_back(vel_x);
  eqGridPtrs.push_back(vel_x);
  vel_y = new Grid_TXYZ(dx, dy, dz, vel_y_halo_x + px, vel_y_halo_y + py, vel_y_halo_z + pz, "vel_y");
  gridPtrs.push_back(vel_y);
  eqGridPtrs.push_back(vel_y);
  vel_z = new Grid_TXYZ(dx, dy, dz, vel_z_halo_x + px, vel_z_halo_y + py, vel_z_halo_z + pz, "vel_z");
  gridPtrs.push_back(vel_z);
  eqGridPtrs.push_back(vel_z);
  stress_xx = new Grid_TXYZ(dx, dy, dz, stress_xx_halo_x + px, stress_xx_halo_y + py, stress_xx_halo_z + pz, "stress_xx");
  gridPtrs.push_back(stress_xx);
  eqGridPtrs.push_back(stress_xx);
  stress_yy = new Grid_TXYZ(dx, dy, dz, stress_yy_halo_x + px, stress_yy_halo_y + py, stress_yy_halo_z + pz, "stress_yy");
  gridPtrs.push_back(stress_yy);
  eqGridPtrs.push_back(stress_yy);
  stress_zz = new Grid_TXYZ(dx, dy, dz, stress_zz_halo_x + px, stress_zz_halo_y + py, stress_zz_halo_z + pz, "stress_zz");
  gridPtrs.push_back(stress_zz);
  eqGridPtrs.push_back(stress_zz);
  stress_xy = new Grid_TXYZ(dx, dy, dz, stress_xy_halo_x + px, stress_xy_halo_y + py, stress_xy_halo_z + pz, "stress_xy");
  gridPtrs.push_back(stress_xy);
  eqGridPtrs.push_back(stress_xy);
  stress_xz = new Grid_TXYZ(dx, dy, dz, stress_xz_halo_x + px, stress_xz_halo_y + py, stress_xz_halo_z + pz, "stress_xz");
  gridPtrs.push_back(stress_xz);
  eqGridPtrs.push_back(stress_xz);
  stress_yz = new Grid_TXYZ(dx, dy, dz, stress_yz_halo_x + px, stress_yz_halo_y + py, stress_yz_halo_z + pz, "stress_yz");
  gridPtrs.push_back(stress_yz);
  eqGridPtrs.push_back(stress_yz);
  stress_mem_xx = new Grid_TXYZ(dx, dy, dz, stress_mem_xx_halo_x + px, stress_mem_xx_halo_y + py, stress_mem_xx_halo_z + pz, "stress_mem_xx");
  gridPtrs.push_back(stress_mem_xx);
  eqGridPtrs.push_back(stress_mem_xx);
  stress_mem_yy = new Grid_TXYZ(dx, dy, dz, stress_mem_yy_halo_x + px, stress_mem_yy_halo_y + py, stress_mem_yy_halo_z + pz, "stress_mem_yy");
  gridPtrs.push_back(stress_mem_yy);
  eqGridPtrs.push_back(stress_mem_yy);
  stress_mem_zz = new Grid_TXYZ(dx, dy, dz, stress_mem_zz_halo_x + px, stress_mem_zz_halo_y + py, stress_mem_zz_halo_z + pz, "stress_mem_zz");
  gridPtrs.push_back(stress_mem_zz);
  eqGridPtrs.push_back(stress_mem_zz);
  stress_mem_xy = new Grid_TXYZ(dx, dy, dz, stress_mem_xy_halo_x + px, stress_mem_xy_halo_y + py, stress_mem_xy_halo_z + pz, "stress_mem_xy");
  gridPtrs.push_back(stress_mem_xy);
  eqGridPtrs.push_back(stress_mem_xy);
  stress_mem_xz = new Grid_TXYZ(dx, dy, dz, stress_mem_xz_halo_x + px, stress_mem_xz_halo_y + py, stress_mem_xz_halo_z + pz, "stress_mem_xz");
  gridPtrs.push_back(stress_mem_xz);
  eqGridPtrs.push_back(stress_mem_xz);
  stress_mem_yz = new Grid_TXYZ(dx, dy, dz, stress_mem_yz_halo_x + px, stress_mem_yz_halo_y + py, stress_mem_yz_halo_z + pz, "stress_mem_yz");
  gridPtrs.push_back(stress_mem_yz);
  eqGridPtrs.push_back(stress_mem_yz);
  weight = new Grid_XYZ(dx, dy, dz, weight_halo_x + px, weight_halo_y + py, weight_halo_z + pz, "weight");
  gridPtrs.push_back(weight);
  tau2 = new Grid_XYZ(dx, dy, dz, tau2_halo_x + px, tau2_halo_y + py, tau2_halo_z + pz, "tau2");
  gridPtrs.push_back(tau2);
  anelastic_ap = new Grid_XYZ(dx, dy, dz, anelastic_ap_halo_x + px, anelastic_ap_halo_y + py, anelastic_ap_halo_z + pz, "anelastic_ap");
  gridPtrs.push_back(anelastic_ap);
  anelastic_as_diag = new Grid_XYZ(dx, dy, dz, anelastic_as_diag_halo_x + px, anelastic_as_diag_halo_y + py, anelastic_as_diag_halo_z + pz, "anelastic_as_diag");
  gridPtrs.push_back(anelastic_as_diag);
  anelastic_xy = new Grid_XYZ(dx, dy, dz, anelastic_xy_halo_x + px, anelastic_xy_halo_y + py, anelastic_xy_halo_z + pz, "anelastic_xy");
  gridPtrs.push_back(anelastic_xy);
  anelastic_xz = new Grid_XYZ(dx, dy, dz, anelastic_xz_halo_x + px, anelastic_xz_halo_y + py, anelastic_xz_halo_z + pz, "anelastic_xz");
  gridPtrs.push_back(anelastic_xz);
  anelastic_yz = new Grid_XYZ(dx, dy, dz, anelastic_yz_halo_x + px, anelastic_yz_halo_y + py, anelastic_yz_halo_z + pz, "anelastic_yz");
  gridPtrs.push_back(anelastic_yz);
  lambda = new Grid_XYZ(dx, dy, dz, lambda_halo_x + px, lambda_halo_y + py, lambda_halo_z + pz, "lambda");
  gridPtrs.push_back(lambda);
  rho = new Grid_XYZ(dx, dy, dz, rho_halo_x + px, rho_halo_y + py, rho_halo_z + pz, "rho");
  gridPtrs.push_back(rho);
  mu = new Grid_XYZ(dx, dy, dz, mu_halo_x + px, mu_halo_y + py, mu_halo_z + pz, "mu");
  gridPtrs.push_back(mu);
  sponge = new Grid_XYZ(dx, dy, dz, sponge_halo_x + px, sponge_halo_y + py, sponge_halo_z + pz, "sponge");
  gridPtrs.push_back(sponge);
 }

 virtual void allocParams() {
  paramPtrs.clear();
  delta_t = new GenericGrid0d<real_t>();
  paramPtrs.push_back(delta_t);
  h = new GenericGrid0d<real_t>();
  paramPtrs.push_back(h);
 }
};

 ////// Stencil equation 'velocity' //////

struct Stencil_velocity {
 std::string name = "velocity";

 // 78 FP operation(s) per point:
 // vel_x(t+1, x, y, z) = ((vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))))) * sponge(x, y, z)).
 // vel_y(t+1, x, y, z) = ((vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))))) * sponge(x, y, z)).
 // vel_z(t+1, x, y, z) = ((vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))))) * sponge(x, y, z)).
 const int scalar_fp_ops = 78;

 // All grids updated by this equation.
 std::vector<RealVecGridBase*> eqGridPtrs;
 void init(StencilContext_awp& context) {
  eqGridPtrs.clear();
  eqGridPtrs.push_back(context.vel_x);
  eqGridPtrs.push_back(context.vel_y);
  eqGridPtrs.push_back(context.vel_z);
 }

 // Calculate one scalar result relative to indices t, x, y, z.
 void calc_scalar(StencilContext_awp& context, idx_t t, idx_t x, idx_t y, idx_t z) {

 // temp1 = delta_t().
 real_t temp1 = (*context.delta_t)();

 // temp2 = h().
 real_t temp2 = (*context.h)();

 // temp3 = rho(x, y, z).
 real_t temp3 = context.rho->readElem(x, y, z, __LINE__);

 // temp4 = rho(x, y-1, z).
 real_t temp4 = context.rho->readElem(x, y-1, z, __LINE__);

 // temp5 = rho(x, y, z) + rho(x, y-1, z).
 real_t temp5 = temp3 + temp4;

 // temp6 = rho(x, y, z-1).
 real_t temp6 = context.rho->readElem(x, y, z-1, __LINE__);

 // temp7 = rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1).
 real_t temp7 = temp5 + temp6;

 // temp8 = rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1).
 real_t temp8 = temp7 + context.rho->readElem(x, y-1, z-1, __LINE__);

 // temp9 = h() * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)).
 real_t temp9 = temp2 * temp8;

 // temp10 = 0.25.
 real_t temp10 = 2.50000000000000000e-01;

 // temp11 = h() * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25.
 real_t temp11 = temp9 * temp10;

 // temp12 = (delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)).
 real_t temp12 = temp1 / temp11;

 // temp13 = 1.125.
 real_t temp13 = 1.12500000000000000e+00;

 // temp14 = 1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z)).
 real_t temp14 = temp13 * (context.stress_xx->readElem(t, x, y, z, __LINE__) - context.stress_xx->readElem(t, x-1, y, z, __LINE__));

 // temp15 = -0.0416667.
 real_t temp15 = -4.16666666666666644e-02;

 // temp16 = -0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z)).
 real_t temp16 = temp15 * (context.stress_xx->readElem(t, x+1, y, z, __LINE__) - context.stress_xx->readElem(t, x-2, y, z, __LINE__));

 // temp17 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))).
 real_t temp17 = temp14 + temp16;

 // temp18 = stress_xy(t, x, y, z).
 real_t temp18 = context.stress_xy->readElem(t, x, y, z, __LINE__);

 // temp19 = (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z)).
 real_t temp19 = temp18 - context.stress_xy->readElem(t, x, y-1, z, __LINE__);

 // temp20 = 1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z)).
 real_t temp20 = temp13 * temp19;

 // temp21 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))).
 real_t temp21 = temp17 + temp20;

 // temp22 = -0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z)).
 real_t temp22 = temp15 * (context.stress_xy->readElem(t, x, y+1, z, __LINE__) - context.stress_xy->readElem(t, x, y-2, z, __LINE__));

 // temp23 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))).
 real_t temp23 = temp21 + temp22;

 // temp24 = stress_xz(t, x, y, z).
 real_t temp24 = context.stress_xz->readElem(t, x, y, z, __LINE__);

 // temp25 = (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1)).
 real_t temp25 = temp24 - context.stress_xz->readElem(t, x, y, z-1, __LINE__);

 // temp26 = 1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1)).
 real_t temp26 = temp13 * temp25;

 // temp27 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))).
 real_t temp27 = temp23 + temp26;

 // temp28 = -0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)).
 real_t temp28 = temp15 * (context.stress_xz->readElem(t, x, y, z+1, __LINE__) - context.stress_xz->readElem(t, x, y, z-2, __LINE__));

 // temp29 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2))).
 real_t temp29 = temp27 + temp28;

 // temp30 = (delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))).
 real_t temp30 = temp12 * temp29;

 // temp31 = vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2))))).
 real_t temp31 = context.vel_x->readElem(t, x, y, z, __LINE__) + temp30;

 // temp32 = sponge(x, y, z).
 real_t temp32 = context.sponge->readElem(x, y, z, __LINE__);

 // temp33 = (vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))))) * sponge(x, y, z).
 real_t temp33 = temp31 * temp32;

 // temp34 = ((vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))))) * sponge(x, y, z)).
 real_t temp34 = temp33;

 // Save result to vel_x(t+1, x, y, z):
 context.vel_x->writeElem(temp34, t+1, x, y, z, __LINE__);

 // temp35 = rho(x+1, y, z).
 real_t temp35 = context.rho->readElem(x+1, y, z, __LINE__);

 // temp36 = rho(x, y, z) + rho(x+1, y, z).
 real_t temp36 = temp3 + temp35;

 // temp37 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1).
 real_t temp37 = temp36 + temp6;

 // temp38 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1).
 real_t temp38 = temp37 + context.rho->readElem(x+1, y, z-1, __LINE__);

 // temp39 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)).
 real_t temp39 = temp2 * temp38;

 // temp40 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25.
 real_t temp40 = temp39 * temp10;

 // temp41 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)).
 real_t temp41 = temp1 / temp40;

 // temp42 = (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z)).
 real_t temp42 = context.stress_xy->readElem(t, x+1, y, z, __LINE__) - temp18;

 // temp43 = 1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z)).
 real_t temp43 = temp13 * temp42;

 // temp44 = -0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z)).
 real_t temp44 = temp15 * (context.stress_xy->readElem(t, x+2, y, z, __LINE__) - context.stress_xy->readElem(t, x-1, y, z, __LINE__));

 // temp45 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))).
 real_t temp45 = temp43 + temp44;

 // temp46 = 1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z)).
 real_t temp46 = temp13 * (context.stress_yy->readElem(t, x, y+1, z, __LINE__) - context.stress_yy->readElem(t, x, y, z, __LINE__));

 // temp47 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))).
 real_t temp47 = temp45 + temp46;

 // temp48 = -0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z)).
 real_t temp48 = temp15 * (context.stress_yy->readElem(t, x, y+2, z, __LINE__) - context.stress_yy->readElem(t, x, y-1, z, __LINE__));

 // temp49 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))).
 real_t temp49 = temp47 + temp48;

 // temp50 = stress_yz(t, x, y, z).
 real_t temp50 = context.stress_yz->readElem(t, x, y, z, __LINE__);

 // temp51 = (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1)).
 real_t temp51 = temp50 - context.stress_yz->readElem(t, x, y, z-1, __LINE__);

 // temp52 = 1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1)).
 real_t temp52 = temp13 * temp51;

 // temp53 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))).
 real_t temp53 = temp49 + temp52;

 // temp54 = -0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)).
 real_t temp54 = temp15 * (context.stress_yz->readElem(t, x, y, z+1, __LINE__) - context.stress_yz->readElem(t, x, y, z-2, __LINE__));

 // temp55 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2))).
 real_t temp55 = temp53 + temp54;

 // temp56 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))).
 real_t temp56 = temp41 * temp55;

 // temp57 = vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2))))).
 real_t temp57 = context.vel_y->readElem(t, x, y, z, __LINE__) + temp56;

 // temp58 = (vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))))) * sponge(x, y, z).
 real_t temp58 = temp57 * temp32;

 // temp59 = ((vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))))) * sponge(x, y, z)).
 real_t temp59 = temp58;

 // Save result to vel_y(t+1, x, y, z):
 context.vel_y->writeElem(temp59, t+1, x, y, z, __LINE__);

 // temp60 = rho(x, y, z) + rho(x+1, y, z).
 real_t temp60 = temp3 + temp35;

 // temp61 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z).
 real_t temp61 = temp60 + temp4;

 // temp62 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z).
 real_t temp62 = temp61 + context.rho->readElem(x+1, y-1, z, __LINE__);

 // temp63 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)).
 real_t temp63 = temp2 * temp62;

 // temp64 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25.
 real_t temp64 = temp63 * temp10;

 // temp65 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)).
 real_t temp65 = temp1 / temp64;

 // temp66 = (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z)).
 real_t temp66 = context.stress_xz->readElem(t, x+1, y, z, __LINE__) - temp24;

 // temp67 = 1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z)).
 real_t temp67 = temp13 * temp66;

 // temp68 = -0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z)).
 real_t temp68 = temp15 * (context.stress_xz->readElem(t, x+2, y, z, __LINE__) - context.stress_xz->readElem(t, x-1, y, z, __LINE__));

 // temp69 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))).
 real_t temp69 = temp67 + temp68;

 // temp70 = (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z)).
 real_t temp70 = temp50 - context.stress_yz->readElem(t, x, y-1, z, __LINE__);

 // temp71 = 1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z)).
 real_t temp71 = temp13 * temp70;

 // temp72 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))).
 real_t temp72 = temp69 + temp71;

 // temp73 = -0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z)).
 real_t temp73 = temp15 * (context.stress_yz->readElem(t, x, y+1, z, __LINE__) - context.stress_yz->readElem(t, x, y-2, z, __LINE__));

 // temp74 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))).
 real_t temp74 = temp72 + temp73;

 // temp75 = 1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z)).
 real_t temp75 = temp13 * (context.stress_zz->readElem(t, x, y, z+1, __LINE__) - context.stress_zz->readElem(t, x, y, z, __LINE__));

 // temp76 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))).
 real_t temp76 = temp74 + temp75;

 // temp77 = -0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)).
 real_t temp77 = temp15 * (context.stress_zz->readElem(t, x, y, z+2, __LINE__) - context.stress_zz->readElem(t, x, y, z-1, __LINE__));

 // temp78 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1))).
 real_t temp78 = temp76 + temp77;

 // temp79 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))).
 real_t temp79 = temp65 * temp78;

 // temp80 = vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1))))).
 real_t temp80 = context.vel_z->readElem(t, x, y, z, __LINE__) + temp79;

 // temp81 = (vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))))) * sponge(x, y, z).
 real_t temp81 = temp80 * temp32;

 // temp82 = ((vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))))) * sponge(x, y, z)).
 real_t temp82 = temp81;

 // Save result to vel_z(t+1, x, y, z):
 context.vel_z->writeElem(temp82, t+1, x, y, z, __LINE__);
} // scalar calculation.

 // Calculate 16 result(s) relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 // SIMD calculations use 44 vector block(s) created from 38 aligned vector-block(s).
 // There are 1248 FP operation(s) per cluster.
 void calc_cluster(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {

 // Un-normalized indices.
 idx_t t = tv;
 idx_t x = xv * 4;
 idx_t y = yv * 4;
 idx_t z = zv * 1;

 // Read aligned vector block from vel_x at t, x, y, z.
 real_vec_t temp_vec1 = context.vel_x->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from rho at x, y, z.
 real_vec_t temp_vec2 = context.rho->readVecNorm(xv, yv, zv, __LINE__);

 // Read aligned vector block from rho at x, y-4, z.
 real_vec_t temp_vec3 = context.rho->readVecNorm(xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from rho at x, y-1, z.
 real_vec_t temp_vec4;
 // temp_vec4[0] = temp_vec3[12];  // for x, y-1, z;
 // temp_vec4[1] = temp_vec3[13];  // for x+1, y-1, z;
 // temp_vec4[2] = temp_vec3[14];  // for x+2, y-1, z;
 // temp_vec4[3] = temp_vec3[15];  // for x+3, y-1, z;
 // temp_vec4[4] = temp_vec2[0];  // for x, y, z;
 // temp_vec4[5] = temp_vec2[1];  // for x+1, y, z;
 // temp_vec4[6] = temp_vec2[2];  // for x+2, y, z;
 // temp_vec4[7] = temp_vec2[3];  // for x+3, y, z;
 // temp_vec4[8] = temp_vec2[4];  // for x, y+1, z;
 // temp_vec4[9] = temp_vec2[5];  // for x+1, y+1, z;
 // temp_vec4[10] = temp_vec2[6];  // for x+2, y+1, z;
 // temp_vec4[11] = temp_vec2[7];  // for x+3, y+1, z;
 // temp_vec4[12] = temp_vec2[8];  // for x, y+2, z;
 // temp_vec4[13] = temp_vec2[9];  // for x+1, y+2, z;
 // temp_vec4[14] = temp_vec2[10];  // for x+2, y+2, z;
 // temp_vec4[15] = temp_vec2[11];  // for x+3, y+2, z;
 // Get 12 element(s) from temp_vec2 and 4 from temp_vec3.
 real_vec_align<12>(temp_vec4, temp_vec2, temp_vec3);

 // Read aligned vector block from rho at x, y, z-1.
 real_vec_t temp_vec5 = context.rho->readVecNorm(xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from rho at x, y-4, z-1.
 real_vec_t temp_vec6 = context.rho->readVecNorm(xv, yv-(4/4), zv-(1/1), __LINE__);

 // Construct unaligned vector block from rho at x, y-1, z-1.
 real_vec_t temp_vec7;
 // temp_vec7[0] = temp_vec6[12];  // for x, y-1, z-1;
 // temp_vec7[1] = temp_vec6[13];  // for x+1, y-1, z-1;
 // temp_vec7[2] = temp_vec6[14];  // for x+2, y-1, z-1;
 // temp_vec7[3] = temp_vec6[15];  // for x+3, y-1, z-1;
 // temp_vec7[4] = temp_vec5[0];  // for x, y, z-1;
 // temp_vec7[5] = temp_vec5[1];  // for x+1, y, z-1;
 // temp_vec7[6] = temp_vec5[2];  // for x+2, y, z-1;
 // temp_vec7[7] = temp_vec5[3];  // for x+3, y, z-1;
 // temp_vec7[8] = temp_vec5[4];  // for x, y+1, z-1;
 // temp_vec7[9] = temp_vec5[5];  // for x+1, y+1, z-1;
 // temp_vec7[10] = temp_vec5[6];  // for x+2, y+1, z-1;
 // temp_vec7[11] = temp_vec5[7];  // for x+3, y+1, z-1;
 // temp_vec7[12] = temp_vec5[8];  // for x, y+2, z-1;
 // temp_vec7[13] = temp_vec5[9];  // for x+1, y+2, z-1;
 // temp_vec7[14] = temp_vec5[10];  // for x+2, y+2, z-1;
 // temp_vec7[15] = temp_vec5[11];  // for x+3, y+2, z-1;
 // Get 12 element(s) from temp_vec5 and 4 from temp_vec6.
 real_vec_align<12>(temp_vec7, temp_vec5, temp_vec6);

 // Read aligned vector block from stress_xx at t, x, y, z.
 real_vec_t temp_vec8 = context.stress_xx->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_xx at t, x-4, y, z.
 real_vec_t temp_vec9 = context.stress_xx->readVecNorm(tv, xv-(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from stress_xx at t, x-1, y, z.
 real_vec_t temp_vec10;
 // temp_vec10[0] = temp_vec9[3];  // for t, x-1, y, z;
 // temp_vec10[1] = temp_vec8[0];  // for t, x, y, z;
 // temp_vec10[2] = temp_vec8[1];  // for t, x+1, y, z;
 // temp_vec10[3] = temp_vec8[2];  // for t, x+2, y, z;
 // temp_vec10[4] = temp_vec9[7];  // for t, x-1, y+1, z;
 // temp_vec10[5] = temp_vec8[4];  // for t, x, y+1, z;
 // temp_vec10[6] = temp_vec8[5];  // for t, x+1, y+1, z;
 // temp_vec10[7] = temp_vec8[6];  // for t, x+2, y+1, z;
 // temp_vec10[8] = temp_vec9[11];  // for t, x-1, y+2, z;
 // temp_vec10[9] = temp_vec8[8];  // for t, x, y+2, z;
 // temp_vec10[10] = temp_vec8[9];  // for t, x+1, y+2, z;
 // temp_vec10[11] = temp_vec8[10];  // for t, x+2, y+2, z;
 // temp_vec10[12] = temp_vec9[15];  // for t, x-1, y+3, z;
 // temp_vec10[13] = temp_vec8[12];  // for t, x, y+3, z;
 // temp_vec10[14] = temp_vec8[13];  // for t, x+1, y+3, z;
 // temp_vec10[15] = temp_vec8[14];  // for t, x+2, y+3, z;
 const real_vec_t_data ctrl_data_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14 = { .ci = { 3, ctrl_sel_bit |0, ctrl_sel_bit |1, ctrl_sel_bit |2, 7, ctrl_sel_bit |4, ctrl_sel_bit |5, ctrl_sel_bit |6, 11, ctrl_sel_bit |8, ctrl_sel_bit |9, ctrl_sel_bit |10, 15, ctrl_sel_bit |12, ctrl_sel_bit |13, ctrl_sel_bit |14 } };
 const real_vec_t ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14(ctrl_data_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14);
 real_vec_permute2(temp_vec10, ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14, temp_vec9, temp_vec8);

 // Read aligned vector block from stress_xx at t, x+4, y, z.
 real_vec_t temp_vec11 = context.stress_xx->readVecNorm(tv, xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from stress_xx at t, x+1, y, z.
 real_vec_t temp_vec12;
 // temp_vec12[0] = temp_vec8[1];  // for t, x+1, y, z;
 // temp_vec12[1] = temp_vec8[2];  // for t, x+2, y, z;
 // temp_vec12[2] = temp_vec8[3];  // for t, x+3, y, z;
 // temp_vec12[3] = temp_vec11[0];  // for t, x+4, y, z;
 // temp_vec12[4] = temp_vec8[5];  // for t, x+1, y+1, z;
 // temp_vec12[5] = temp_vec8[6];  // for t, x+2, y+1, z;
 // temp_vec12[6] = temp_vec8[7];  // for t, x+3, y+1, z;
 // temp_vec12[7] = temp_vec11[4];  // for t, x+4, y+1, z;
 // temp_vec12[8] = temp_vec8[9];  // for t, x+1, y+2, z;
 // temp_vec12[9] = temp_vec8[10];  // for t, x+2, y+2, z;
 // temp_vec12[10] = temp_vec8[11];  // for t, x+3, y+2, z;
 // temp_vec12[11] = temp_vec11[8];  // for t, x+4, y+2, z;
 // temp_vec12[12] = temp_vec8[13];  // for t, x+1, y+3, z;
 // temp_vec12[13] = temp_vec8[14];  // for t, x+2, y+3, z;
 // temp_vec12[14] = temp_vec8[15];  // for t, x+3, y+3, z;
 // temp_vec12[15] = temp_vec11[12];  // for t, x+4, y+3, z;
 const real_vec_t_data ctrl_data_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12 = { .ci = { 1, 2, 3, ctrl_sel_bit |0, 5, 6, 7, ctrl_sel_bit |4, 9, 10, 11, ctrl_sel_bit |8, 13, 14, 15, ctrl_sel_bit |12 } };
 const real_vec_t ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12(ctrl_data_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12);
 real_vec_permute2(temp_vec12, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec8, temp_vec11);

 // Construct unaligned vector block from stress_xx at t, x-2, y, z.
 real_vec_t temp_vec13;
 // temp_vec13[0] = temp_vec9[2];  // for t, x-2, y, z;
 // temp_vec13[1] = temp_vec9[3];  // for t, x-1, y, z;
 // temp_vec13[2] = temp_vec8[0];  // for t, x, y, z;
 // temp_vec13[3] = temp_vec8[1];  // for t, x+1, y, z;
 // temp_vec13[4] = temp_vec9[6];  // for t, x-2, y+1, z;
 // temp_vec13[5] = temp_vec9[7];  // for t, x-1, y+1, z;
 // temp_vec13[6] = temp_vec8[4];  // for t, x, y+1, z;
 // temp_vec13[7] = temp_vec8[5];  // for t, x+1, y+1, z;
 // temp_vec13[8] = temp_vec9[10];  // for t, x-2, y+2, z;
 // temp_vec13[9] = temp_vec9[11];  // for t, x-1, y+2, z;
 // temp_vec13[10] = temp_vec8[8];  // for t, x, y+2, z;
 // temp_vec13[11] = temp_vec8[9];  // for t, x+1, y+2, z;
 // temp_vec13[12] = temp_vec9[14];  // for t, x-2, y+3, z;
 // temp_vec13[13] = temp_vec9[15];  // for t, x-1, y+3, z;
 // temp_vec13[14] = temp_vec8[12];  // for t, x, y+3, z;
 // temp_vec13[15] = temp_vec8[13];  // for t, x+1, y+3, z;
 const real_vec_t_data ctrl_data_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13 = { .ci = { 2, 3, ctrl_sel_bit |0, ctrl_sel_bit |1, 6, 7, ctrl_sel_bit |4, ctrl_sel_bit |5, 10, 11, ctrl_sel_bit |8, ctrl_sel_bit |9, 14, 15, ctrl_sel_bit |12, ctrl_sel_bit |13 } };
 const real_vec_t ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13(ctrl_data_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13);
 real_vec_permute2(temp_vec13, ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13, temp_vec9, temp_vec8);

 // Read aligned vector block from stress_xy at t, x, y, z.
 real_vec_t temp_vec14 = context.stress_xy->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_xy at t, x, y-4, z.
 real_vec_t temp_vec15 = context.stress_xy->readVecNorm(tv, xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from stress_xy at t, x, y-1, z.
 real_vec_t temp_vec16;
 // temp_vec16[0] = temp_vec15[12];  // for t, x, y-1, z;
 // temp_vec16[1] = temp_vec15[13];  // for t, x+1, y-1, z;
 // temp_vec16[2] = temp_vec15[14];  // for t, x+2, y-1, z;
 // temp_vec16[3] = temp_vec15[15];  // for t, x+3, y-1, z;
 // temp_vec16[4] = temp_vec14[0];  // for t, x, y, z;
 // temp_vec16[5] = temp_vec14[1];  // for t, x+1, y, z;
 // temp_vec16[6] = temp_vec14[2];  // for t, x+2, y, z;
 // temp_vec16[7] = temp_vec14[3];  // for t, x+3, y, z;
 // temp_vec16[8] = temp_vec14[4];  // for t, x, y+1, z;
 // temp_vec16[9] = temp_vec14[5];  // for t, x+1, y+1, z;
 // temp_vec16[10] = temp_vec14[6];  // for t, x+2, y+1, z;
 // temp_vec16[11] = temp_vec14[7];  // for t, x+3, y+1, z;
 // temp_vec16[12] = temp_vec14[8];  // for t, x, y+2, z;
 // temp_vec16[13] = temp_vec14[9];  // for t, x+1, y+2, z;
 // temp_vec16[14] = temp_vec14[10];  // for t, x+2, y+2, z;
 // temp_vec16[15] = temp_vec14[11];  // for t, x+3, y+2, z;
 // Get 12 element(s) from temp_vec14 and 4 from temp_vec15.
 real_vec_align<12>(temp_vec16, temp_vec14, temp_vec15);

 // Read aligned vector block from stress_xy at t, x, y+4, z.
 real_vec_t temp_vec17 = context.stress_xy->readVecNorm(tv, xv, yv+(4/4), zv, __LINE__);

 // Construct unaligned vector block from stress_xy at t, x, y+1, z.
 real_vec_t temp_vec18;
 // temp_vec18[0] = temp_vec14[4];  // for t, x, y+1, z;
 // temp_vec18[1] = temp_vec14[5];  // for t, x+1, y+1, z;
 // temp_vec18[2] = temp_vec14[6];  // for t, x+2, y+1, z;
 // temp_vec18[3] = temp_vec14[7];  // for t, x+3, y+1, z;
 // temp_vec18[4] = temp_vec14[8];  // for t, x, y+2, z;
 // temp_vec18[5] = temp_vec14[9];  // for t, x+1, y+2, z;
 // temp_vec18[6] = temp_vec14[10];  // for t, x+2, y+2, z;
 // temp_vec18[7] = temp_vec14[11];  // for t, x+3, y+2, z;
 // temp_vec18[8] = temp_vec14[12];  // for t, x, y+3, z;
 // temp_vec18[9] = temp_vec14[13];  // for t, x+1, y+3, z;
 // temp_vec18[10] = temp_vec14[14];  // for t, x+2, y+3, z;
 // temp_vec18[11] = temp_vec14[15];  // for t, x+3, y+3, z;
 // temp_vec18[12] = temp_vec17[0];  // for t, x, y+4, z;
 // temp_vec18[13] = temp_vec17[1];  // for t, x+1, y+4, z;
 // temp_vec18[14] = temp_vec17[2];  // for t, x+2, y+4, z;
 // temp_vec18[15] = temp_vec17[3];  // for t, x+3, y+4, z;
 // Get 4 element(s) from temp_vec17 and 12 from temp_vec14.
 real_vec_align<4>(temp_vec18, temp_vec17, temp_vec14);

 // Construct unaligned vector block from stress_xy at t, x, y-2, z.
 real_vec_t temp_vec19;
 // temp_vec19[0] = temp_vec15[8];  // for t, x, y-2, z;
 // temp_vec19[1] = temp_vec15[9];  // for t, x+1, y-2, z;
 // temp_vec19[2] = temp_vec15[10];  // for t, x+2, y-2, z;
 // temp_vec19[3] = temp_vec15[11];  // for t, x+3, y-2, z;
 // temp_vec19[4] = temp_vec15[12];  // for t, x, y-1, z;
 // temp_vec19[5] = temp_vec15[13];  // for t, x+1, y-1, z;
 // temp_vec19[6] = temp_vec15[14];  // for t, x+2, y-1, z;
 // temp_vec19[7] = temp_vec15[15];  // for t, x+3, y-1, z;
 // temp_vec19[8] = temp_vec14[0];  // for t, x, y, z;
 // temp_vec19[9] = temp_vec14[1];  // for t, x+1, y, z;
 // temp_vec19[10] = temp_vec14[2];  // for t, x+2, y, z;
 // temp_vec19[11] = temp_vec14[3];  // for t, x+3, y, z;
 // temp_vec19[12] = temp_vec14[4];  // for t, x, y+1, z;
 // temp_vec19[13] = temp_vec14[5];  // for t, x+1, y+1, z;
 // temp_vec19[14] = temp_vec14[6];  // for t, x+2, y+1, z;
 // temp_vec19[15] = temp_vec14[7];  // for t, x+3, y+1, z;
 // Get 8 element(s) from temp_vec14 and 8 from temp_vec15.
 real_vec_align<8>(temp_vec19, temp_vec14, temp_vec15);

 // Read aligned vector block from stress_xz at t, x, y, z.
 real_vec_t temp_vec20 = context.stress_xz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_xz at t, x, y, z-1.
 real_vec_t temp_vec21 = context.stress_xz->readVecNorm(tv, xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from stress_xz at t, x, y, z+1.
 real_vec_t temp_vec22 = context.stress_xz->readVecNorm(tv, xv, yv, zv+(1/1), __LINE__);

 // Read aligned vector block from stress_xz at t, x, y, z-2.
 real_vec_t temp_vec23 = context.stress_xz->readVecNorm(tv, xv, yv, zv-(2/1), __LINE__);

 // Read aligned vector block from sponge at x, y, z.
 real_vec_t temp_vec24 = context.sponge->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec25 = delta_t().
 real_vec_t temp_vec25 = (*context.delta_t)();

 // temp_vec26 = h().
 real_vec_t temp_vec26 = (*context.h)();

 // temp_vec27 = rho(x, y, z).
 real_vec_t temp_vec27 = temp_vec2;

 // temp_vec28 = rho(x, y-1, z).
 real_vec_t temp_vec28 = temp_vec4;

 // temp_vec29 = rho(x, y, z) + rho(x, y-1, z).
 real_vec_t temp_vec29 = temp_vec27 + temp_vec28;

 // temp_vec30 = rho(x, y, z-1).
 real_vec_t temp_vec30 = temp_vec5;

 // temp_vec31 = rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1).
 real_vec_t temp_vec31 = temp_vec29 + temp_vec30;

 // temp_vec32 = rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1).
 real_vec_t temp_vec32 = temp_vec31 + temp_vec7;

 // temp_vec33 = h() * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)).
 real_vec_t temp_vec33 = temp_vec26 * temp_vec32;

 // temp_vec34 = 0.25.
 real_vec_t temp_vec34 = 2.50000000000000000e-01;

 // temp_vec35 = h() * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25.
 real_vec_t temp_vec35 = temp_vec33 * temp_vec34;

 // temp_vec36 = (delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)).
 real_vec_t temp_vec36 = temp_vec25 / temp_vec35;

 // temp_vec37 = 1.125.
 real_vec_t temp_vec37 = 1.12500000000000000e+00;

 // temp_vec38 = 1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z)).
 real_vec_t temp_vec38 = temp_vec37 * (temp_vec8 - temp_vec10);

 // temp_vec39 = -0.0416667.
 real_vec_t temp_vec39 = -4.16666666666666644e-02;

 // temp_vec40 = -0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z)).
 real_vec_t temp_vec40 = temp_vec39 * (temp_vec12 - temp_vec13);

 // temp_vec41 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))).
 real_vec_t temp_vec41 = temp_vec38 + temp_vec40;

 // temp_vec42 = stress_xy(t, x, y, z).
 real_vec_t temp_vec42 = temp_vec14;

 // temp_vec43 = (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z)).
 real_vec_t temp_vec43 = temp_vec42 - temp_vec16;

 // temp_vec44 = 1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z)).
 real_vec_t temp_vec44 = temp_vec37 * temp_vec43;

 // temp_vec45 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))).
 real_vec_t temp_vec45 = temp_vec41 + temp_vec44;

 // temp_vec46 = -0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z)).
 real_vec_t temp_vec46 = temp_vec39 * (temp_vec18 - temp_vec19);

 // temp_vec47 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))).
 real_vec_t temp_vec47 = temp_vec45 + temp_vec46;

 // temp_vec48 = stress_xz(t, x, y, z).
 real_vec_t temp_vec48 = temp_vec20;

 // temp_vec49 = (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1)).
 real_vec_t temp_vec49 = temp_vec48 - temp_vec21;

 // temp_vec50 = 1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1)).
 real_vec_t temp_vec50 = temp_vec37 * temp_vec49;

 // temp_vec51 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))).
 real_vec_t temp_vec51 = temp_vec47 + temp_vec50;

 // temp_vec52 = -0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)).
 real_vec_t temp_vec52 = temp_vec39 * (temp_vec22 - temp_vec23);

 // temp_vec53 = (1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2))).
 real_vec_t temp_vec53 = temp_vec51 + temp_vec52;

 // temp_vec54 = (delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))).
 real_vec_t temp_vec54 = temp_vec36 * temp_vec53;

 // temp_vec55 = vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2))))).
 real_vec_t temp_vec55 = temp_vec1 + temp_vec54;

 // temp_vec56 = sponge(x, y, z).
 real_vec_t temp_vec56 = temp_vec24;

 // temp_vec57 = (vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))))) * sponge(x, y, z).
 real_vec_t temp_vec57 = temp_vec55 * temp_vec56;

 // temp_vec58 = ((vel_x(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x, y-1, z) + rho(x, y, z-1) + rho(x, y-1, z-1)) * 0.25)) * ((1.125 * (stress_xx(t, x, y, z) - stress_xx(t, x-1, y, z))) + (-0.0416667 * (stress_xx(t, x+1, y, z) - stress_xx(t, x-2, y, z))) + (1.125 * (stress_xy(t, x, y, z) - stress_xy(t, x, y-1, z))) + (-0.0416667 * (stress_xy(t, x, y+1, z) - stress_xy(t, x, y-2, z))) + (1.125 * (stress_xz(t, x, y, z) - stress_xz(t, x, y, z-1))) + (-0.0416667 * (stress_xz(t, x, y, z+1) - stress_xz(t, x, y, z-2)))))) * sponge(x, y, z)).
 real_vec_t temp_vec58 = temp_vec57;

 // Save result to vel_x(t+1, x, y, z):
 
 // Write aligned vector block to vel_x at t+1, x, y, z.
context.vel_x->writeVecNorm(temp_vec58, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from vel_y at t, x, y, z.
 real_vec_t temp_vec59 = context.vel_y->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from rho at x+4, y, z.
 real_vec_t temp_vec60 = context.rho->readVecNorm(xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from rho at x+1, y, z.
 real_vec_t temp_vec61;
 // temp_vec61[0] = temp_vec2[1];  // for x+1, y, z;
 // temp_vec61[1] = temp_vec2[2];  // for x+2, y, z;
 // temp_vec61[2] = temp_vec2[3];  // for x+3, y, z;
 // temp_vec61[3] = temp_vec60[0];  // for x+4, y, z;
 // temp_vec61[4] = temp_vec2[5];  // for x+1, y+1, z;
 // temp_vec61[5] = temp_vec2[6];  // for x+2, y+1, z;
 // temp_vec61[6] = temp_vec2[7];  // for x+3, y+1, z;
 // temp_vec61[7] = temp_vec60[4];  // for x+4, y+1, z;
 // temp_vec61[8] = temp_vec2[9];  // for x+1, y+2, z;
 // temp_vec61[9] = temp_vec2[10];  // for x+2, y+2, z;
 // temp_vec61[10] = temp_vec2[11];  // for x+3, y+2, z;
 // temp_vec61[11] = temp_vec60[8];  // for x+4, y+2, z;
 // temp_vec61[12] = temp_vec2[13];  // for x+1, y+3, z;
 // temp_vec61[13] = temp_vec2[14];  // for x+2, y+3, z;
 // temp_vec61[14] = temp_vec2[15];  // for x+3, y+3, z;
 // temp_vec61[15] = temp_vec60[12];  // for x+4, y+3, z;
 real_vec_permute2(temp_vec61, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec2, temp_vec60);

 // Read aligned vector block from rho at x+4, y, z-1.
 real_vec_t temp_vec62 = context.rho->readVecNorm(xv+(4/4), yv, zv-(1/1), __LINE__);

 // Construct unaligned vector block from rho at x+1, y, z-1.
 real_vec_t temp_vec63;
 // temp_vec63[0] = temp_vec5[1];  // for x+1, y, z-1;
 // temp_vec63[1] = temp_vec5[2];  // for x+2, y, z-1;
 // temp_vec63[2] = temp_vec5[3];  // for x+3, y, z-1;
 // temp_vec63[3] = temp_vec62[0];  // for x+4, y, z-1;
 // temp_vec63[4] = temp_vec5[5];  // for x+1, y+1, z-1;
 // temp_vec63[5] = temp_vec5[6];  // for x+2, y+1, z-1;
 // temp_vec63[6] = temp_vec5[7];  // for x+3, y+1, z-1;
 // temp_vec63[7] = temp_vec62[4];  // for x+4, y+1, z-1;
 // temp_vec63[8] = temp_vec5[9];  // for x+1, y+2, z-1;
 // temp_vec63[9] = temp_vec5[10];  // for x+2, y+2, z-1;
 // temp_vec63[10] = temp_vec5[11];  // for x+3, y+2, z-1;
 // temp_vec63[11] = temp_vec62[8];  // for x+4, y+2, z-1;
 // temp_vec63[12] = temp_vec5[13];  // for x+1, y+3, z-1;
 // temp_vec63[13] = temp_vec5[14];  // for x+2, y+3, z-1;
 // temp_vec63[14] = temp_vec5[15];  // for x+3, y+3, z-1;
 // temp_vec63[15] = temp_vec62[12];  // for x+4, y+3, z-1;
 real_vec_permute2(temp_vec63, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec5, temp_vec62);

 // Read aligned vector block from stress_xy at t, x+4, y, z.
 real_vec_t temp_vec64 = context.stress_xy->readVecNorm(tv, xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from stress_xy at t, x+1, y, z.
 real_vec_t temp_vec65;
 // temp_vec65[0] = temp_vec14[1];  // for t, x+1, y, z;
 // temp_vec65[1] = temp_vec14[2];  // for t, x+2, y, z;
 // temp_vec65[2] = temp_vec14[3];  // for t, x+3, y, z;
 // temp_vec65[3] = temp_vec64[0];  // for t, x+4, y, z;
 // temp_vec65[4] = temp_vec14[5];  // for t, x+1, y+1, z;
 // temp_vec65[5] = temp_vec14[6];  // for t, x+2, y+1, z;
 // temp_vec65[6] = temp_vec14[7];  // for t, x+3, y+1, z;
 // temp_vec65[7] = temp_vec64[4];  // for t, x+4, y+1, z;
 // temp_vec65[8] = temp_vec14[9];  // for t, x+1, y+2, z;
 // temp_vec65[9] = temp_vec14[10];  // for t, x+2, y+2, z;
 // temp_vec65[10] = temp_vec14[11];  // for t, x+3, y+2, z;
 // temp_vec65[11] = temp_vec64[8];  // for t, x+4, y+2, z;
 // temp_vec65[12] = temp_vec14[13];  // for t, x+1, y+3, z;
 // temp_vec65[13] = temp_vec14[14];  // for t, x+2, y+3, z;
 // temp_vec65[14] = temp_vec14[15];  // for t, x+3, y+3, z;
 // temp_vec65[15] = temp_vec64[12];  // for t, x+4, y+3, z;
 real_vec_permute2(temp_vec65, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec14, temp_vec64);

 // Construct unaligned vector block from stress_xy at t, x+2, y, z.
 real_vec_t temp_vec66;
 // temp_vec66[0] = temp_vec14[2];  // for t, x+2, y, z;
 // temp_vec66[1] = temp_vec14[3];  // for t, x+3, y, z;
 // temp_vec66[2] = temp_vec64[0];  // for t, x+4, y, z;
 // temp_vec66[3] = temp_vec64[1];  // for t, x+5, y, z;
 // temp_vec66[4] = temp_vec14[6];  // for t, x+2, y+1, z;
 // temp_vec66[5] = temp_vec14[7];  // for t, x+3, y+1, z;
 // temp_vec66[6] = temp_vec64[4];  // for t, x+4, y+1, z;
 // temp_vec66[7] = temp_vec64[5];  // for t, x+5, y+1, z;
 // temp_vec66[8] = temp_vec14[10];  // for t, x+2, y+2, z;
 // temp_vec66[9] = temp_vec14[11];  // for t, x+3, y+2, z;
 // temp_vec66[10] = temp_vec64[8];  // for t, x+4, y+2, z;
 // temp_vec66[11] = temp_vec64[9];  // for t, x+5, y+2, z;
 // temp_vec66[12] = temp_vec14[14];  // for t, x+2, y+3, z;
 // temp_vec66[13] = temp_vec14[15];  // for t, x+3, y+3, z;
 // temp_vec66[14] = temp_vec64[12];  // for t, x+4, y+3, z;
 // temp_vec66[15] = temp_vec64[13];  // for t, x+5, y+3, z;
 real_vec_permute2(temp_vec66, ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13, temp_vec14, temp_vec64);

 // Read aligned vector block from stress_xy at t, x-4, y, z.
 real_vec_t temp_vec67 = context.stress_xy->readVecNorm(tv, xv-(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from stress_xy at t, x-1, y, z.
 real_vec_t temp_vec68;
 // temp_vec68[0] = temp_vec67[3];  // for t, x-1, y, z;
 // temp_vec68[1] = temp_vec14[0];  // for t, x, y, z;
 // temp_vec68[2] = temp_vec14[1];  // for t, x+1, y, z;
 // temp_vec68[3] = temp_vec14[2];  // for t, x+2, y, z;
 // temp_vec68[4] = temp_vec67[7];  // for t, x-1, y+1, z;
 // temp_vec68[5] = temp_vec14[4];  // for t, x, y+1, z;
 // temp_vec68[6] = temp_vec14[5];  // for t, x+1, y+1, z;
 // temp_vec68[7] = temp_vec14[6];  // for t, x+2, y+1, z;
 // temp_vec68[8] = temp_vec67[11];  // for t, x-1, y+2, z;
 // temp_vec68[9] = temp_vec14[8];  // for t, x, y+2, z;
 // temp_vec68[10] = temp_vec14[9];  // for t, x+1, y+2, z;
 // temp_vec68[11] = temp_vec14[10];  // for t, x+2, y+2, z;
 // temp_vec68[12] = temp_vec67[15];  // for t, x-1, y+3, z;
 // temp_vec68[13] = temp_vec14[12];  // for t, x, y+3, z;
 // temp_vec68[14] = temp_vec14[13];  // for t, x+1, y+3, z;
 // temp_vec68[15] = temp_vec14[14];  // for t, x+2, y+3, z;
 real_vec_permute2(temp_vec68, ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14, temp_vec67, temp_vec14);

 // Read aligned vector block from stress_yy at t, x, y, z.
 real_vec_t temp_vec69 = context.stress_yy->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_yy at t, x, y+4, z.
 real_vec_t temp_vec70 = context.stress_yy->readVecNorm(tv, xv, yv+(4/4), zv, __LINE__);

 // Construct unaligned vector block from stress_yy at t, x, y+1, z.
 real_vec_t temp_vec71;
 // temp_vec71[0] = temp_vec69[4];  // for t, x, y+1, z;
 // temp_vec71[1] = temp_vec69[5];  // for t, x+1, y+1, z;
 // temp_vec71[2] = temp_vec69[6];  // for t, x+2, y+1, z;
 // temp_vec71[3] = temp_vec69[7];  // for t, x+3, y+1, z;
 // temp_vec71[4] = temp_vec69[8];  // for t, x, y+2, z;
 // temp_vec71[5] = temp_vec69[9];  // for t, x+1, y+2, z;
 // temp_vec71[6] = temp_vec69[10];  // for t, x+2, y+2, z;
 // temp_vec71[7] = temp_vec69[11];  // for t, x+3, y+2, z;
 // temp_vec71[8] = temp_vec69[12];  // for t, x, y+3, z;
 // temp_vec71[9] = temp_vec69[13];  // for t, x+1, y+3, z;
 // temp_vec71[10] = temp_vec69[14];  // for t, x+2, y+3, z;
 // temp_vec71[11] = temp_vec69[15];  // for t, x+3, y+3, z;
 // temp_vec71[12] = temp_vec70[0];  // for t, x, y+4, z;
 // temp_vec71[13] = temp_vec70[1];  // for t, x+1, y+4, z;
 // temp_vec71[14] = temp_vec70[2];  // for t, x+2, y+4, z;
 // temp_vec71[15] = temp_vec70[3];  // for t, x+3, y+4, z;
 // Get 4 element(s) from temp_vec70 and 12 from temp_vec69.
 real_vec_align<4>(temp_vec71, temp_vec70, temp_vec69);

 // Construct unaligned vector block from stress_yy at t, x, y+2, z.
 real_vec_t temp_vec72;
 // temp_vec72[0] = temp_vec69[8];  // for t, x, y+2, z;
 // temp_vec72[1] = temp_vec69[9];  // for t, x+1, y+2, z;
 // temp_vec72[2] = temp_vec69[10];  // for t, x+2, y+2, z;
 // temp_vec72[3] = temp_vec69[11];  // for t, x+3, y+2, z;
 // temp_vec72[4] = temp_vec69[12];  // for t, x, y+3, z;
 // temp_vec72[5] = temp_vec69[13];  // for t, x+1, y+3, z;
 // temp_vec72[6] = temp_vec69[14];  // for t, x+2, y+3, z;
 // temp_vec72[7] = temp_vec69[15];  // for t, x+3, y+3, z;
 // temp_vec72[8] = temp_vec70[0];  // for t, x, y+4, z;
 // temp_vec72[9] = temp_vec70[1];  // for t, x+1, y+4, z;
 // temp_vec72[10] = temp_vec70[2];  // for t, x+2, y+4, z;
 // temp_vec72[11] = temp_vec70[3];  // for t, x+3, y+4, z;
 // temp_vec72[12] = temp_vec70[4];  // for t, x, y+5, z;
 // temp_vec72[13] = temp_vec70[5];  // for t, x+1, y+5, z;
 // temp_vec72[14] = temp_vec70[6];  // for t, x+2, y+5, z;
 // temp_vec72[15] = temp_vec70[7];  // for t, x+3, y+5, z;
 // Get 8 element(s) from temp_vec70 and 8 from temp_vec69.
 real_vec_align<8>(temp_vec72, temp_vec70, temp_vec69);

 // Read aligned vector block from stress_yy at t, x, y-4, z.
 real_vec_t temp_vec73 = context.stress_yy->readVecNorm(tv, xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from stress_yy at t, x, y-1, z.
 real_vec_t temp_vec74;
 // temp_vec74[0] = temp_vec73[12];  // for t, x, y-1, z;
 // temp_vec74[1] = temp_vec73[13];  // for t, x+1, y-1, z;
 // temp_vec74[2] = temp_vec73[14];  // for t, x+2, y-1, z;
 // temp_vec74[3] = temp_vec73[15];  // for t, x+3, y-1, z;
 // temp_vec74[4] = temp_vec69[0];  // for t, x, y, z;
 // temp_vec74[5] = temp_vec69[1];  // for t, x+1, y, z;
 // temp_vec74[6] = temp_vec69[2];  // for t, x+2, y, z;
 // temp_vec74[7] = temp_vec69[3];  // for t, x+3, y, z;
 // temp_vec74[8] = temp_vec69[4];  // for t, x, y+1, z;
 // temp_vec74[9] = temp_vec69[5];  // for t, x+1, y+1, z;
 // temp_vec74[10] = temp_vec69[6];  // for t, x+2, y+1, z;
 // temp_vec74[11] = temp_vec69[7];  // for t, x+3, y+1, z;
 // temp_vec74[12] = temp_vec69[8];  // for t, x, y+2, z;
 // temp_vec74[13] = temp_vec69[9];  // for t, x+1, y+2, z;
 // temp_vec74[14] = temp_vec69[10];  // for t, x+2, y+2, z;
 // temp_vec74[15] = temp_vec69[11];  // for t, x+3, y+2, z;
 // Get 12 element(s) from temp_vec69 and 4 from temp_vec73.
 real_vec_align<12>(temp_vec74, temp_vec69, temp_vec73);

 // Read aligned vector block from stress_yz at t, x, y, z.
 real_vec_t temp_vec75 = context.stress_yz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_yz at t, x, y, z-1.
 real_vec_t temp_vec76 = context.stress_yz->readVecNorm(tv, xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from stress_yz at t, x, y, z+1.
 real_vec_t temp_vec77 = context.stress_yz->readVecNorm(tv, xv, yv, zv+(1/1), __LINE__);

 // Read aligned vector block from stress_yz at t, x, y, z-2.
 real_vec_t temp_vec78 = context.stress_yz->readVecNorm(tv, xv, yv, zv-(2/1), __LINE__);

 // temp_vec79 = rho(x+1, y, z).
 real_vec_t temp_vec79 = temp_vec61;

 // temp_vec80 = rho(x, y, z) + rho(x+1, y, z).
 real_vec_t temp_vec80 = temp_vec27 + temp_vec79;

 // temp_vec81 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1).
 real_vec_t temp_vec81 = temp_vec80 + temp_vec30;

 // temp_vec82 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1).
 real_vec_t temp_vec82 = temp_vec81 + temp_vec63;

 // temp_vec83 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)).
 real_vec_t temp_vec83 = temp_vec26 * temp_vec82;

 // temp_vec84 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25.
 real_vec_t temp_vec84 = temp_vec83 * temp_vec34;

 // temp_vec85 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)).
 real_vec_t temp_vec85 = temp_vec25 / temp_vec84;

 // temp_vec86 = (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z)).
 real_vec_t temp_vec86 = temp_vec65 - temp_vec42;

 // temp_vec87 = 1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z)).
 real_vec_t temp_vec87 = temp_vec37 * temp_vec86;

 // temp_vec88 = -0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z)).
 real_vec_t temp_vec88 = temp_vec39 * (temp_vec66 - temp_vec68);

 // temp_vec89 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))).
 real_vec_t temp_vec89 = temp_vec87 + temp_vec88;

 // temp_vec90 = 1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z)).
 real_vec_t temp_vec90 = temp_vec37 * (temp_vec71 - temp_vec69);

 // temp_vec91 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))).
 real_vec_t temp_vec91 = temp_vec89 + temp_vec90;

 // temp_vec92 = -0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z)).
 real_vec_t temp_vec92 = temp_vec39 * (temp_vec72 - temp_vec74);

 // temp_vec93 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))).
 real_vec_t temp_vec93 = temp_vec91 + temp_vec92;

 // temp_vec94 = stress_yz(t, x, y, z).
 real_vec_t temp_vec94 = temp_vec75;

 // temp_vec95 = (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1)).
 real_vec_t temp_vec95 = temp_vec94 - temp_vec76;

 // temp_vec96 = 1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1)).
 real_vec_t temp_vec96 = temp_vec37 * temp_vec95;

 // temp_vec97 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))).
 real_vec_t temp_vec97 = temp_vec93 + temp_vec96;

 // temp_vec98 = -0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)).
 real_vec_t temp_vec98 = temp_vec39 * (temp_vec77 - temp_vec78);

 // temp_vec99 = (1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2))).
 real_vec_t temp_vec99 = temp_vec97 + temp_vec98;

 // temp_vec100 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))).
 real_vec_t temp_vec100 = temp_vec85 * temp_vec99;

 // temp_vec101 = vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2))))).
 real_vec_t temp_vec101 = temp_vec59 + temp_vec100;

 // temp_vec102 = (vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))))) * sponge(x, y, z).
 real_vec_t temp_vec102 = temp_vec101 * temp_vec56;

 // temp_vec103 = ((vel_y(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y, z-1) + rho(x+1, y, z-1)) * 0.25)) * ((1.125 * (stress_xy(t, x+1, y, z) - stress_xy(t, x, y, z))) + (-0.0416667 * (stress_xy(t, x+2, y, z) - stress_xy(t, x-1, y, z))) + (1.125 * (stress_yy(t, x, y+1, z) - stress_yy(t, x, y, z))) + (-0.0416667 * (stress_yy(t, x, y+2, z) - stress_yy(t, x, y-1, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y, z-1))) + (-0.0416667 * (stress_yz(t, x, y, z+1) - stress_yz(t, x, y, z-2)))))) * sponge(x, y, z)).
 real_vec_t temp_vec103 = temp_vec102;

 // Save result to vel_y(t+1, x, y, z):
 
 // Write aligned vector block to vel_y at t+1, x, y, z.
context.vel_y->writeVecNorm(temp_vec103, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from vel_z at t, x, y, z.
 real_vec_t temp_vec104 = context.vel_z->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from rho at x+4, y-4, z.
 real_vec_t temp_vec105 = context.rho->readVecNorm(xv+(4/4), yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from rho at x+1, y-1, z.
 real_vec_t temp_vec106;
 // temp_vec106[0] = temp_vec3[13];  // for x+1, y-1, z;
 // temp_vec106[1] = temp_vec3[14];  // for x+2, y-1, z;
 // temp_vec106[2] = temp_vec3[15];  // for x+3, y-1, z;
 // temp_vec106[3] = temp_vec105[12];  // for x+4, y-1, z;
 // temp_vec106[4] = temp_vec2[1];  // for x+1, y, z;
 // temp_vec106[5] = temp_vec2[2];  // for x+2, y, z;
 // temp_vec106[6] = temp_vec2[3];  // for x+3, y, z;
 // temp_vec106[7] = temp_vec60[0];  // for x+4, y, z;
 // temp_vec106[8] = temp_vec2[5];  // for x+1, y+1, z;
 // temp_vec106[9] = temp_vec2[6];  // for x+2, y+1, z;
 // temp_vec106[10] = temp_vec2[7];  // for x+3, y+1, z;
 // temp_vec106[11] = temp_vec60[4];  // for x+4, y+1, z;
 // temp_vec106[12] = temp_vec2[9];  // for x+1, y+2, z;
 // temp_vec106[13] = temp_vec2[10];  // for x+2, y+2, z;
 // temp_vec106[14] = temp_vec2[11];  // for x+3, y+2, z;
 // temp_vec106[15] = temp_vec60[8];  // for x+4, y+2, z;
 // Get 9 element(s) from temp_vec2 and 3 from temp_vec3.
 real_vec_align<13>(temp_vec106, temp_vec2, temp_vec3);
 // Get 3 element(s) from temp_vec60 and 1 from temp_vec105.
 real_vec_align_masked<9>(temp_vec106, temp_vec60, temp_vec105, 0x8888);

 // Read aligned vector block from stress_xz at t, x+4, y, z.
 real_vec_t temp_vec107 = context.stress_xz->readVecNorm(tv, xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from stress_xz at t, x+1, y, z.
 real_vec_t temp_vec108;
 // temp_vec108[0] = temp_vec20[1];  // for t, x+1, y, z;
 // temp_vec108[1] = temp_vec20[2];  // for t, x+2, y, z;
 // temp_vec108[2] = temp_vec20[3];  // for t, x+3, y, z;
 // temp_vec108[3] = temp_vec107[0];  // for t, x+4, y, z;
 // temp_vec108[4] = temp_vec20[5];  // for t, x+1, y+1, z;
 // temp_vec108[5] = temp_vec20[6];  // for t, x+2, y+1, z;
 // temp_vec108[6] = temp_vec20[7];  // for t, x+3, y+1, z;
 // temp_vec108[7] = temp_vec107[4];  // for t, x+4, y+1, z;
 // temp_vec108[8] = temp_vec20[9];  // for t, x+1, y+2, z;
 // temp_vec108[9] = temp_vec20[10];  // for t, x+2, y+2, z;
 // temp_vec108[10] = temp_vec20[11];  // for t, x+3, y+2, z;
 // temp_vec108[11] = temp_vec107[8];  // for t, x+4, y+2, z;
 // temp_vec108[12] = temp_vec20[13];  // for t, x+1, y+3, z;
 // temp_vec108[13] = temp_vec20[14];  // for t, x+2, y+3, z;
 // temp_vec108[14] = temp_vec20[15];  // for t, x+3, y+3, z;
 // temp_vec108[15] = temp_vec107[12];  // for t, x+4, y+3, z;
 real_vec_permute2(temp_vec108, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec20, temp_vec107);

 // Construct unaligned vector block from stress_xz at t, x+2, y, z.
 real_vec_t temp_vec109;
 // temp_vec109[0] = temp_vec20[2];  // for t, x+2, y, z;
 // temp_vec109[1] = temp_vec20[3];  // for t, x+3, y, z;
 // temp_vec109[2] = temp_vec107[0];  // for t, x+4, y, z;
 // temp_vec109[3] = temp_vec107[1];  // for t, x+5, y, z;
 // temp_vec109[4] = temp_vec20[6];  // for t, x+2, y+1, z;
 // temp_vec109[5] = temp_vec20[7];  // for t, x+3, y+1, z;
 // temp_vec109[6] = temp_vec107[4];  // for t, x+4, y+1, z;
 // temp_vec109[7] = temp_vec107[5];  // for t, x+5, y+1, z;
 // temp_vec109[8] = temp_vec20[10];  // for t, x+2, y+2, z;
 // temp_vec109[9] = temp_vec20[11];  // for t, x+3, y+2, z;
 // temp_vec109[10] = temp_vec107[8];  // for t, x+4, y+2, z;
 // temp_vec109[11] = temp_vec107[9];  // for t, x+5, y+2, z;
 // temp_vec109[12] = temp_vec20[14];  // for t, x+2, y+3, z;
 // temp_vec109[13] = temp_vec20[15];  // for t, x+3, y+3, z;
 // temp_vec109[14] = temp_vec107[12];  // for t, x+4, y+3, z;
 // temp_vec109[15] = temp_vec107[13];  // for t, x+5, y+3, z;
 real_vec_permute2(temp_vec109, ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13, temp_vec20, temp_vec107);

 // Read aligned vector block from stress_xz at t, x-4, y, z.
 real_vec_t temp_vec110 = context.stress_xz->readVecNorm(tv, xv-(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from stress_xz at t, x-1, y, z.
 real_vec_t temp_vec111;
 // temp_vec111[0] = temp_vec110[3];  // for t, x-1, y, z;
 // temp_vec111[1] = temp_vec20[0];  // for t, x, y, z;
 // temp_vec111[2] = temp_vec20[1];  // for t, x+1, y, z;
 // temp_vec111[3] = temp_vec20[2];  // for t, x+2, y, z;
 // temp_vec111[4] = temp_vec110[7];  // for t, x-1, y+1, z;
 // temp_vec111[5] = temp_vec20[4];  // for t, x, y+1, z;
 // temp_vec111[6] = temp_vec20[5];  // for t, x+1, y+1, z;
 // temp_vec111[7] = temp_vec20[6];  // for t, x+2, y+1, z;
 // temp_vec111[8] = temp_vec110[11];  // for t, x-1, y+2, z;
 // temp_vec111[9] = temp_vec20[8];  // for t, x, y+2, z;
 // temp_vec111[10] = temp_vec20[9];  // for t, x+1, y+2, z;
 // temp_vec111[11] = temp_vec20[10];  // for t, x+2, y+2, z;
 // temp_vec111[12] = temp_vec110[15];  // for t, x-1, y+3, z;
 // temp_vec111[13] = temp_vec20[12];  // for t, x, y+3, z;
 // temp_vec111[14] = temp_vec20[13];  // for t, x+1, y+3, z;
 // temp_vec111[15] = temp_vec20[14];  // for t, x+2, y+3, z;
 real_vec_permute2(temp_vec111, ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14, temp_vec110, temp_vec20);

 // Read aligned vector block from stress_yz at t, x, y-4, z.
 real_vec_t temp_vec112 = context.stress_yz->readVecNorm(tv, xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from stress_yz at t, x, y-1, z.
 real_vec_t temp_vec113;
 // temp_vec113[0] = temp_vec112[12];  // for t, x, y-1, z;
 // temp_vec113[1] = temp_vec112[13];  // for t, x+1, y-1, z;
 // temp_vec113[2] = temp_vec112[14];  // for t, x+2, y-1, z;
 // temp_vec113[3] = temp_vec112[15];  // for t, x+3, y-1, z;
 // temp_vec113[4] = temp_vec75[0];  // for t, x, y, z;
 // temp_vec113[5] = temp_vec75[1];  // for t, x+1, y, z;
 // temp_vec113[6] = temp_vec75[2];  // for t, x+2, y, z;
 // temp_vec113[7] = temp_vec75[3];  // for t, x+3, y, z;
 // temp_vec113[8] = temp_vec75[4];  // for t, x, y+1, z;
 // temp_vec113[9] = temp_vec75[5];  // for t, x+1, y+1, z;
 // temp_vec113[10] = temp_vec75[6];  // for t, x+2, y+1, z;
 // temp_vec113[11] = temp_vec75[7];  // for t, x+3, y+1, z;
 // temp_vec113[12] = temp_vec75[8];  // for t, x, y+2, z;
 // temp_vec113[13] = temp_vec75[9];  // for t, x+1, y+2, z;
 // temp_vec113[14] = temp_vec75[10];  // for t, x+2, y+2, z;
 // temp_vec113[15] = temp_vec75[11];  // for t, x+3, y+2, z;
 // Get 12 element(s) from temp_vec75 and 4 from temp_vec112.
 real_vec_align<12>(temp_vec113, temp_vec75, temp_vec112);

 // Read aligned vector block from stress_yz at t, x, y+4, z.
 real_vec_t temp_vec114 = context.stress_yz->readVecNorm(tv, xv, yv+(4/4), zv, __LINE__);

 // Construct unaligned vector block from stress_yz at t, x, y+1, z.
 real_vec_t temp_vec115;
 // temp_vec115[0] = temp_vec75[4];  // for t, x, y+1, z;
 // temp_vec115[1] = temp_vec75[5];  // for t, x+1, y+1, z;
 // temp_vec115[2] = temp_vec75[6];  // for t, x+2, y+1, z;
 // temp_vec115[3] = temp_vec75[7];  // for t, x+3, y+1, z;
 // temp_vec115[4] = temp_vec75[8];  // for t, x, y+2, z;
 // temp_vec115[5] = temp_vec75[9];  // for t, x+1, y+2, z;
 // temp_vec115[6] = temp_vec75[10];  // for t, x+2, y+2, z;
 // temp_vec115[7] = temp_vec75[11];  // for t, x+3, y+2, z;
 // temp_vec115[8] = temp_vec75[12];  // for t, x, y+3, z;
 // temp_vec115[9] = temp_vec75[13];  // for t, x+1, y+3, z;
 // temp_vec115[10] = temp_vec75[14];  // for t, x+2, y+3, z;
 // temp_vec115[11] = temp_vec75[15];  // for t, x+3, y+3, z;
 // temp_vec115[12] = temp_vec114[0];  // for t, x, y+4, z;
 // temp_vec115[13] = temp_vec114[1];  // for t, x+1, y+4, z;
 // temp_vec115[14] = temp_vec114[2];  // for t, x+2, y+4, z;
 // temp_vec115[15] = temp_vec114[3];  // for t, x+3, y+4, z;
 // Get 4 element(s) from temp_vec114 and 12 from temp_vec75.
 real_vec_align<4>(temp_vec115, temp_vec114, temp_vec75);

 // Construct unaligned vector block from stress_yz at t, x, y-2, z.
 real_vec_t temp_vec116;
 // temp_vec116[0] = temp_vec112[8];  // for t, x, y-2, z;
 // temp_vec116[1] = temp_vec112[9];  // for t, x+1, y-2, z;
 // temp_vec116[2] = temp_vec112[10];  // for t, x+2, y-2, z;
 // temp_vec116[3] = temp_vec112[11];  // for t, x+3, y-2, z;
 // temp_vec116[4] = temp_vec112[12];  // for t, x, y-1, z;
 // temp_vec116[5] = temp_vec112[13];  // for t, x+1, y-1, z;
 // temp_vec116[6] = temp_vec112[14];  // for t, x+2, y-1, z;
 // temp_vec116[7] = temp_vec112[15];  // for t, x+3, y-1, z;
 // temp_vec116[8] = temp_vec75[0];  // for t, x, y, z;
 // temp_vec116[9] = temp_vec75[1];  // for t, x+1, y, z;
 // temp_vec116[10] = temp_vec75[2];  // for t, x+2, y, z;
 // temp_vec116[11] = temp_vec75[3];  // for t, x+3, y, z;
 // temp_vec116[12] = temp_vec75[4];  // for t, x, y+1, z;
 // temp_vec116[13] = temp_vec75[5];  // for t, x+1, y+1, z;
 // temp_vec116[14] = temp_vec75[6];  // for t, x+2, y+1, z;
 // temp_vec116[15] = temp_vec75[7];  // for t, x+3, y+1, z;
 // Get 8 element(s) from temp_vec75 and 8 from temp_vec112.
 real_vec_align<8>(temp_vec116, temp_vec75, temp_vec112);

 // Read aligned vector block from stress_zz at t, x, y, z+1.
 real_vec_t temp_vec117 = context.stress_zz->readVecNorm(tv, xv, yv, zv+(1/1), __LINE__);

 // Read aligned vector block from stress_zz at t, x, y, z.
 real_vec_t temp_vec118 = context.stress_zz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_zz at t, x, y, z+2.
 real_vec_t temp_vec119 = context.stress_zz->readVecNorm(tv, xv, yv, zv+(2/1), __LINE__);

 // Read aligned vector block from stress_zz at t, x, y, z-1.
 real_vec_t temp_vec120 = context.stress_zz->readVecNorm(tv, xv, yv, zv-(1/1), __LINE__);

 // temp_vec121 = rho(x, y, z) + rho(x+1, y, z).
 real_vec_t temp_vec121 = temp_vec27 + temp_vec79;

 // temp_vec122 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z).
 real_vec_t temp_vec122 = temp_vec121 + temp_vec28;

 // temp_vec123 = rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z).
 real_vec_t temp_vec123 = temp_vec122 + temp_vec106;

 // temp_vec124 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)).
 real_vec_t temp_vec124 = temp_vec26 * temp_vec123;

 // temp_vec125 = h() * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25.
 real_vec_t temp_vec125 = temp_vec124 * temp_vec34;

 // temp_vec126 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)).
 real_vec_t temp_vec126 = temp_vec25 / temp_vec125;

 // temp_vec127 = (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z)).
 real_vec_t temp_vec127 = temp_vec108 - temp_vec48;

 // temp_vec128 = 1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z)).
 real_vec_t temp_vec128 = temp_vec37 * temp_vec127;

 // temp_vec129 = -0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z)).
 real_vec_t temp_vec129 = temp_vec39 * (temp_vec109 - temp_vec111);

 // temp_vec130 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))).
 real_vec_t temp_vec130 = temp_vec128 + temp_vec129;

 // temp_vec131 = (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z)).
 real_vec_t temp_vec131 = temp_vec94 - temp_vec113;

 // temp_vec132 = 1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z)).
 real_vec_t temp_vec132 = temp_vec37 * temp_vec131;

 // temp_vec133 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))).
 real_vec_t temp_vec133 = temp_vec130 + temp_vec132;

 // temp_vec134 = -0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z)).
 real_vec_t temp_vec134 = temp_vec39 * (temp_vec115 - temp_vec116);

 // temp_vec135 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))).
 real_vec_t temp_vec135 = temp_vec133 + temp_vec134;

 // temp_vec136 = 1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z)).
 real_vec_t temp_vec136 = temp_vec37 * (temp_vec117 - temp_vec118);

 // temp_vec137 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))).
 real_vec_t temp_vec137 = temp_vec135 + temp_vec136;

 // temp_vec138 = -0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)).
 real_vec_t temp_vec138 = temp_vec39 * (temp_vec119 - temp_vec120);

 // temp_vec139 = (1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1))).
 real_vec_t temp_vec139 = temp_vec137 + temp_vec138;

 // temp_vec140 = (delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))).
 real_vec_t temp_vec140 = temp_vec126 * temp_vec139;

 // temp_vec141 = vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1))))).
 real_vec_t temp_vec141 = temp_vec104 + temp_vec140;

 // temp_vec142 = (vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))))) * sponge(x, y, z).
 real_vec_t temp_vec142 = temp_vec141 * temp_vec56;

 // temp_vec143 = ((vel_z(t, x, y, z) + ((delta_t / (h * (rho(x, y, z) + rho(x+1, y, z) + rho(x, y-1, z) + rho(x+1, y-1, z)) * 0.25)) * ((1.125 * (stress_xz(t, x+1, y, z) - stress_xz(t, x, y, z))) + (-0.0416667 * (stress_xz(t, x+2, y, z) - stress_xz(t, x-1, y, z))) + (1.125 * (stress_yz(t, x, y, z) - stress_yz(t, x, y-1, z))) + (-0.0416667 * (stress_yz(t, x, y+1, z) - stress_yz(t, x, y-2, z))) + (1.125 * (stress_zz(t, x, y, z+1) - stress_zz(t, x, y, z))) + (-0.0416667 * (stress_zz(t, x, y, z+2) - stress_zz(t, x, y, z-1)))))) * sponge(x, y, z)).
 real_vec_t temp_vec143 = temp_vec142;

 // Save result to vel_z(t+1, x, y, z):
 
 // Write aligned vector block to vel_z at t+1, x, y, z.
context.vel_z->writeVecNorm(temp_vec143, tv+(1/1), xv, yv, zv, __LINE__);
;
} // vector calculation.

 // Prefetches cache line(s) for entire stencil  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned rho at x, y-4, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x, y-4, z.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x, y, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y-4, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x-4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x+4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x-4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y-4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y+4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x+4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x-4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z-2.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z-1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z+1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x+4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y-4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y+4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y-4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z-2.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z-1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z+1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y+4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z-1.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+1.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+2.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t, x, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t, x, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t, x, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}

 // Prefetches cache line(s) for leading edge of stencil in '+x' direction  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster_x(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned rho at x, y-4, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y-4, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x+4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y-4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y+4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x+4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z-2.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z-1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z+1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x+4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y-4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y+4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y-4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z-2.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z-1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z+1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y+4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z-1.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+1.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+2.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t, x, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t, x, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t, x, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}

 // Prefetches cache line(s) for leading edge of stencil in '+y' direction  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster_y(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned rho at x, y, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z-1.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x-4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x+4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x-4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y+4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x+4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x-4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z-2.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z-1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z+1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x+4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y+4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z-2.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z-1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z+1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y+4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z-1.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+1.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+2.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t, x, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t, x, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t, x, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}

 // Prefetches cache line(s) for leading edge of stencil in '+z' direction  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster_z(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned rho at x, y-4, z.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y-4, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned rho at x+4, y, z.
 p = (const char*)context.rho->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x-4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x+4, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x-4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y-4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y+4, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x+4, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x-4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z+1.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x+4, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y-4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y+4, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y-4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z+1.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y+4, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z+2.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t, x, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t, x, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t, x, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}
};

 ////// Stencil equation 'stress' //////

struct Stencil_stress {
 std::string name = "stress";

 // 192 FP operation(s) per point:
 // stress_xx(t+1, x, y, z) = ((stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)))) * sponge(x, y, z)).
 // stress_yy(t+1, x, y, z) = ((stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)))) * sponge(x, y, z)).
 // stress_zz(t+1, x, y, z) = ((stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)))) * sponge(x, y, z)).
 // stress_xy(t+1, x, y, z) = ((stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)))) * sponge(x, y, z)).
 // stress_xz(t+1, x, y, z) = ((stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)))) * sponge(x, y, z)).
 // stress_yz(t+1, x, y, z) = ((stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)))) * sponge(x, y, z)).
 // stress_mem_xx(t+1, x, y, z) = ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 // stress_mem_yy(t+1, x, y, z) = ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 // stress_mem_zz(t+1, x, y, z) = ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 // stress_mem_xy(t+1, x, y, z) = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))).
 // stress_mem_xz(t+1, x, y, z) = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))).
 // stress_mem_yz(t+1, x, y, z) = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))).
 const int scalar_fp_ops = 192;

 // All grids updated by this equation.
 std::vector<RealVecGridBase*> eqGridPtrs;
 void init(StencilContext_awp& context) {
  eqGridPtrs.clear();
  eqGridPtrs.push_back(context.stress_xx);
  eqGridPtrs.push_back(context.stress_yy);
  eqGridPtrs.push_back(context.stress_zz);
  eqGridPtrs.push_back(context.stress_xy);
  eqGridPtrs.push_back(context.stress_xz);
  eqGridPtrs.push_back(context.stress_yz);
  eqGridPtrs.push_back(context.stress_mem_xx);
  eqGridPtrs.push_back(context.stress_mem_yy);
  eqGridPtrs.push_back(context.stress_mem_zz);
  eqGridPtrs.push_back(context.stress_mem_xy);
  eqGridPtrs.push_back(context.stress_mem_xz);
  eqGridPtrs.push_back(context.stress_mem_yz);
 }

 // Calculate one scalar result relative to indices t, x, y, z.
 void calc_scalar(StencilContext_awp& context, idx_t t, idx_t x, idx_t y, idx_t z) {

 // temp1 = delta_t().
 real_t temp1 = (*context.delta_t)();

 // temp2 = h().
 real_t temp2 = (*context.h)();

 // temp3 = (delta_t / h).
 real_t temp3 = temp1 / temp2;

 // temp4 = 2.
 real_t temp4 = 2.00000000000000000e+00;

 // temp5 = 8.
 real_t temp5 = 8.00000000000000000e+00;

 // temp6 = mu(x, y, z).
 real_t temp6 = context.mu->readElem(x, y, z, __LINE__);

 // temp7 = mu(x+1, y, z).
 real_t temp7 = context.mu->readElem(x+1, y, z, __LINE__);

 // temp8 = mu(x, y, z) + mu(x+1, y, z).
 real_t temp8 = temp6 + temp7;

 // temp9 = mu(x, y-1, z).
 real_t temp9 = context.mu->readElem(x, y-1, z, __LINE__);

 // temp10 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z).
 real_t temp10 = temp8 + temp9;

 // temp11 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z).
 real_t temp11 = temp10 + context.mu->readElem(x+1, y-1, z, __LINE__);

 // temp12 = mu(x, y, z-1).
 real_t temp12 = context.mu->readElem(x, y, z-1, __LINE__);

 // temp13 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1).
 real_t temp13 = temp11 + temp12;

 // temp14 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1).
 real_t temp14 = temp13 + context.mu->readElem(x+1, y, z-1, __LINE__);

 // temp15 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1).
 real_t temp15 = temp14 + context.mu->readElem(x, y-1, z-1, __LINE__);

 // temp16 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1).
 real_t temp16 = temp15 + context.mu->readElem(x+1, y-1, z-1, __LINE__);

 // temp17 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_t temp17 = temp5 / temp16;

 // temp18 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_t temp18 = temp4 * temp17;

 // temp19 = 1.125.
 real_t temp19 = 1.12500000000000000e+00;

 // temp20 = vel_x(t+1, x, y, z).
 real_t temp20 = context.vel_x->readElem(t+1, x, y, z, __LINE__);

 // temp21 = (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z)).
 real_t temp21 = context.vel_x->readElem(t+1, x+1, y, z, __LINE__) - temp20;

 // temp22 = 1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z)).
 real_t temp22 = temp19 * temp21;

 // temp23 = -0.0416667.
 real_t temp23 = -4.16666666666666644e-02;

 // temp24 = -0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z)).
 real_t temp24 = temp23 * (context.vel_x->readElem(t+1, x+2, y, z, __LINE__) - context.vel_x->readElem(t+1, x-1, y, z, __LINE__));

 // temp25 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_t temp25 = temp22 + temp24;

 // temp26 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z)))).
 real_t temp26 = temp18 * temp25;

 // temp27 = (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))).
 real_t temp27 = temp5 / (context.lambda->readElem(x, y, z, __LINE__) + context.lambda->readElem(x+1, y, z, __LINE__) + context.lambda->readElem(x, y-1, z, __LINE__) + context.lambda->readElem(x+1, y-1, z, __LINE__) + context.lambda->readElem(x, y, z-1, __LINE__) + context.lambda->readElem(x+1, y, z-1, __LINE__) + context.lambda->readElem(x, y-1, z-1, __LINE__) + context.lambda->readElem(x+1, y-1, z-1, __LINE__));

 // temp28 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_t temp28 = temp22 + temp24;

 // temp29 = vel_y(t+1, x, y, z).
 real_t temp29 = context.vel_y->readElem(t+1, x, y, z, __LINE__);

 // temp30 = (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z)).
 real_t temp30 = temp29 - context.vel_y->readElem(t+1, x, y-1, z, __LINE__);

 // temp31 = 1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z)).
 real_t temp31 = temp19 * temp30;

 // temp32 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))).
 real_t temp32 = temp28 + temp31;

 // temp33 = -0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z)).
 real_t temp33 = temp23 * (context.vel_y->readElem(t+1, x, y+1, z, __LINE__) - context.vel_y->readElem(t+1, x, y-2, z, __LINE__));

 // temp34 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_t temp34 = temp32 + temp33;

 // temp35 = vel_z(t+1, x, y, z).
 real_t temp35 = context.vel_z->readElem(t+1, x, y, z, __LINE__);

 // temp36 = (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1)).
 real_t temp36 = temp35 - context.vel_z->readElem(t+1, x, y, z-1, __LINE__);

 // temp37 = 1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1)).
 real_t temp37 = temp19 * temp36;

 // temp38 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))).
 real_t temp38 = temp34 + temp37;

 // temp39 = -0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)).
 real_t temp39 = temp23 * (context.vel_z->readElem(t+1, x, y, z+1, __LINE__) - context.vel_z->readElem(t+1, x, y, z-2, __LINE__));

 // temp40 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_t temp40 = temp38 + temp39;

 // temp41 = (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_t temp41 = temp27 * temp40;

 // temp42 = (2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))).
 real_t temp42 = temp26 + temp41;

 // temp43 = (delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp43 = temp3 * temp42;

 // temp44 = stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp44 = context.stress_xx->readElem(t, x, y, z, __LINE__) + temp43;

 // temp45 = tau2(x, y, z).
 real_t temp45 = context.tau2->readElem(x, y, z, __LINE__);

 // temp46 = stress_mem_xx(t, x, y, z).
 real_t temp46 = context.stress_mem_xx->readElem(t, x, y, z, __LINE__);

 // temp47 = tau2(x, y, z) * stress_mem_xx(t, x, y, z).
 real_t temp47 = temp45 * temp46;

 // temp48 = 1.
 real_t temp48 = 1.00000000000000000e+00;

 // temp49 = (1 / h).
 real_t temp49 = temp48 / temp2;

 // temp50 = (1 - tau2(x, y, z)).
 real_t temp50 = temp48 - temp45;

 // temp51 = (1 / h) * (1 - tau2(x, y, z)).
 real_t temp51 = temp49 * temp50;

 // temp52 = weight(x, y, z).
 real_t temp52 = context.weight->readElem(x, y, z, __LINE__);

 // temp53 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_t temp53 = temp51 * temp52;

 // temp54 = anelastic_as_diag(x, y, z).
 real_t temp54 = context.anelastic_as_diag->readElem(x, y, z, __LINE__);

 // temp55 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z).
 real_t temp55 = temp17 * temp54;

 // temp56 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_t temp56 = temp31 + temp33;

 // temp57 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))).
 real_t temp57 = temp56 + temp37;

 // temp58 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_t temp58 = temp57 + temp39;

 // temp59 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_t temp59 = temp55 * temp58;

 // temp60 = 0.5.
 real_t temp60 = 5.00000000000000000e-01;

 // temp61 = 0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))).
 real_t temp61 = temp60 * temp27;

 // temp62 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1)))).
 real_t temp62 = temp17 + temp61;

 // temp63 = ((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z).
 real_t temp63 = temp62 * context.anelastic_ap->readElem(x, y, z, __LINE__);

 // temp64 = ((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_t temp64 = temp63 * temp40;

 // temp65 = (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp65 = temp59 - temp64;

 // temp66 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp66 = temp53 * temp65;

 // temp67 = (tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp67 = temp47 + temp66;

 // temp68 = (tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z).
 real_t temp68 = temp67 + temp46;

 // temp69 = delta_t() * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)).
 real_t temp69 = temp1 * temp68;

 // temp70 = stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z))).
 real_t temp70 = temp44 + temp69;

 // temp71 = sponge(x, y, z).
 real_t temp71 = context.sponge->readElem(x, y, z, __LINE__);

 // temp72 = (stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)))) * sponge(x, y, z).
 real_t temp72 = temp70 * temp71;

 // temp73 = ((stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)))) * sponge(x, y, z)).
 real_t temp73 = temp72;

 // Save result to stress_xx(t+1, x, y, z):
 context.stress_xx->writeElem(temp73, t+1, x, y, z, __LINE__);

 // temp74 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_t temp74 = temp4 * temp17;

 // temp75 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_t temp75 = temp31 + temp33;

 // temp76 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z)))).
 real_t temp76 = temp74 * temp75;

 // temp77 = (2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))).
 real_t temp77 = temp76 + temp41;

 // temp78 = (delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp78 = temp3 * temp77;

 // temp79 = stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp79 = context.stress_yy->readElem(t, x, y, z, __LINE__) + temp78;

 // temp80 = stress_mem_yy(t, x, y, z).
 real_t temp80 = context.stress_mem_yy->readElem(t, x, y, z, __LINE__);

 // temp81 = tau2(x, y, z) * stress_mem_yy(t, x, y, z).
 real_t temp81 = temp45 * temp80;

 // temp82 = (1 / h) * (1 - tau2(x, y, z)).
 real_t temp82 = temp49 * temp50;

 // temp83 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_t temp83 = temp82 * temp52;

 // temp84 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z).
 real_t temp84 = temp17 * temp54;

 // temp85 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_t temp85 = temp22 + temp24;

 // temp86 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))).
 real_t temp86 = temp85 + temp37;

 // temp87 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_t temp87 = temp86 + temp39;

 // temp88 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_t temp88 = temp84 * temp87;

 // temp89 = (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp89 = temp88 - temp64;

 // temp90 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp90 = temp83 * temp89;

 // temp91 = (tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp91 = temp81 + temp90;

 // temp92 = (tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z).
 real_t temp92 = temp91 + temp80;

 // temp93 = delta_t() * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)).
 real_t temp93 = temp1 * temp92;

 // temp94 = stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z))).
 real_t temp94 = temp79 + temp93;

 // temp95 = (stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)))) * sponge(x, y, z).
 real_t temp95 = temp94 * temp71;

 // temp96 = ((stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)))) * sponge(x, y, z)).
 real_t temp96 = temp95;

 // Save result to stress_yy(t+1, x, y, z):
 context.stress_yy->writeElem(temp96, t+1, x, y, z, __LINE__);

 // temp97 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_t temp97 = temp4 * temp17;

 // temp98 = (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_t temp98 = temp37 + temp39;

 // temp99 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_t temp99 = temp97 * temp98;

 // temp100 = (2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))).
 real_t temp100 = temp99 + temp41;

 // temp101 = (delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp101 = temp3 * temp100;

 // temp102 = stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp102 = context.stress_zz->readElem(t, x, y, z, __LINE__) + temp101;

 // temp103 = stress_mem_zz(t, x, y, z).
 real_t temp103 = context.stress_mem_zz->readElem(t, x, y, z, __LINE__);

 // temp104 = tau2(x, y, z) * stress_mem_zz(t, x, y, z).
 real_t temp104 = temp45 * temp103;

 // temp105 = (1 / h) * (1 - tau2(x, y, z)).
 real_t temp105 = temp49 * temp50;

 // temp106 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_t temp106 = temp105 * temp52;

 // temp107 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z).
 real_t temp107 = temp17 * temp54;

 // temp108 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_t temp108 = temp22 + temp24;

 // temp109 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))).
 real_t temp109 = temp108 + temp31;

 // temp110 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_t temp110 = temp109 + temp33;

 // temp111 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z)))).
 real_t temp111 = temp107 * temp110;

 // temp112 = (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp112 = temp111 - temp64;

 // temp113 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_t temp113 = temp106 * temp112;

 // temp114 = (tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp114 = temp104 + temp113;

 // temp115 = (tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z).
 real_t temp115 = temp114 + temp103;

 // temp116 = delta_t() * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)).
 real_t temp116 = temp1 * temp115;

 // temp117 = stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z))).
 real_t temp117 = temp102 + temp116;

 // temp118 = (stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)))) * sponge(x, y, z).
 real_t temp118 = temp117 * temp71;

 // temp119 = ((stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)))) * sponge(x, y, z)).
 real_t temp119 = temp118;

 // Save result to stress_zz(t+1, x, y, z):
 context.stress_zz->writeElem(temp119, t+1, x, y, z, __LINE__);

 // temp120 = mu(x, y, z) + mu(x, y, z-1).
 real_t temp120 = temp6 + temp12;

 // temp121 = (2 / (mu(x, y, z) + mu(x, y, z-1))).
 real_t temp121 = temp4 / temp120;

 // temp122 = (2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t().
 real_t temp122 = temp121 * temp1;

 // temp123 = (((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h).
 real_t temp123 = temp122 / temp2;

 // temp124 = (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z)).
 real_t temp124 = context.vel_x->readElem(t+1, x, y+1, z, __LINE__) - temp20;

 // temp125 = 1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z)).
 real_t temp125 = temp19 * temp124;

 // temp126 = -0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z)).
 real_t temp126 = temp23 * (context.vel_x->readElem(t+1, x, y+2, z, __LINE__) - context.vel_x->readElem(t+1, x, y-1, z, __LINE__));

 // temp127 = (1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))).
 real_t temp127 = temp125 + temp126;

 // temp128 = (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z)).
 real_t temp128 = temp29 - context.vel_y->readElem(t+1, x-1, y, z, __LINE__);

 // temp129 = 1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z)).
 real_t temp129 = temp19 * temp128;

 // temp130 = (1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))).
 real_t temp130 = temp127 + temp129;

 // temp131 = -0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)).
 real_t temp131 = temp23 * (context.vel_y->readElem(t+1, x+1, y, z, __LINE__) - context.vel_y->readElem(t+1, x-2, y, z, __LINE__));

 // temp132 = (1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))).
 real_t temp132 = temp130 + temp131;

 // temp133 = (((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))).
 real_t temp133 = temp123 * temp132;

 // temp134 = stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))).
 real_t temp134 = context.stress_xy->readElem(t, x, y, z, __LINE__) + temp133;

 // temp135 = stress_mem_xy(t, x, y, z).
 real_t temp135 = context.stress_mem_xy->readElem(t, x, y, z, __LINE__);

 // temp136 = tau2(x, y, z) * stress_mem_xy(t, x, y, z).
 real_t temp136 = temp45 * temp135;

 // temp137 = (0.5 / h).
 real_t temp137 = temp60 / temp2;

 // temp138 = (0.5 / h) * (1 - tau2(x, y, z)).
 real_t temp138 = temp137 * temp50;

 // temp139 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_t temp139 = temp138 * temp52;

 // temp140 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))).
 real_t temp140 = temp139 * temp121;

 // temp141 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z).
 real_t temp141 = temp140 * context.anelastic_xy->readElem(x, y, z, __LINE__);

 // temp142 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))).
 real_t temp142 = temp141 * temp132;

 // temp143 = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))).
 real_t temp143 = temp136 - temp142;

 // temp144 = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z).
 real_t temp144 = temp143 + temp135;

 // temp145 = delta_t() * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)).
 real_t temp145 = temp1 * temp144;

 // temp146 = stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z))).
 real_t temp146 = temp134 + temp145;

 // temp147 = (stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)))) * sponge(x, y, z).
 real_t temp147 = temp146 * temp71;

 // temp148 = ((stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)))) * sponge(x, y, z)).
 real_t temp148 = temp147;

 // Save result to stress_xy(t+1, x, y, z):
 context.stress_xy->writeElem(temp148, t+1, x, y, z, __LINE__);

 // temp149 = mu(x, y, z) + mu(x, y-1, z).
 real_t temp149 = temp6 + temp9;

 // temp150 = (2 / (mu(x, y, z) + mu(x, y-1, z))).
 real_t temp150 = temp4 / temp149;

 // temp151 = (2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t().
 real_t temp151 = temp150 * temp1;

 // temp152 = (((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h).
 real_t temp152 = temp151 / temp2;

 // temp153 = (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z)).
 real_t temp153 = context.vel_x->readElem(t+1, x, y, z+1, __LINE__) - temp20;

 // temp154 = 1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z)).
 real_t temp154 = temp19 * temp153;

 // temp155 = -0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1)).
 real_t temp155 = temp23 * (context.vel_x->readElem(t+1, x, y, z+2, __LINE__) - context.vel_x->readElem(t+1, x, y, z-1, __LINE__));

 // temp156 = (1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))).
 real_t temp156 = temp154 + temp155;

 // temp157 = (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z)).
 real_t temp157 = temp35 - context.vel_z->readElem(t+1, x-1, y, z, __LINE__);

 // temp158 = 1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z)).
 real_t temp158 = temp19 * temp157;

 // temp159 = (1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))).
 real_t temp159 = temp156 + temp158;

 // temp160 = -0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)).
 real_t temp160 = temp23 * (context.vel_z->readElem(t+1, x+1, y, z, __LINE__) - context.vel_z->readElem(t+1, x-2, y, z, __LINE__));

 // temp161 = (1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))).
 real_t temp161 = temp159 + temp160;

 // temp162 = (((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))).
 real_t temp162 = temp152 * temp161;

 // temp163 = stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))).
 real_t temp163 = context.stress_xz->readElem(t, x, y, z, __LINE__) + temp162;

 // temp164 = stress_mem_xz(t, x, y, z).
 real_t temp164 = context.stress_mem_xz->readElem(t, x, y, z, __LINE__);

 // temp165 = tau2(x, y, z) * stress_mem_xz(t, x, y, z).
 real_t temp165 = temp45 * temp164;

 // temp166 = (0.5 / h) * (1 - tau2(x, y, z)).
 real_t temp166 = temp137 * temp50;

 // temp167 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_t temp167 = temp166 * temp52;

 // temp168 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))).
 real_t temp168 = temp167 * temp150;

 // temp169 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z).
 real_t temp169 = temp168 * context.anelastic_xz->readElem(x, y, z, __LINE__);

 // temp170 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))).
 real_t temp170 = temp169 * temp161;

 // temp171 = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))).
 real_t temp171 = temp165 - temp170;

 // temp172 = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z).
 real_t temp172 = temp171 + temp164;

 // temp173 = delta_t() * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)).
 real_t temp173 = temp1 * temp172;

 // temp174 = stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z))).
 real_t temp174 = temp163 + temp173;

 // temp175 = (stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)))) * sponge(x, y, z).
 real_t temp175 = temp174 * temp71;

 // temp176 = ((stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)))) * sponge(x, y, z)).
 real_t temp176 = temp175;

 // Save result to stress_xz(t+1, x, y, z):
 context.stress_xz->writeElem(temp176, t+1, x, y, z, __LINE__);

 // temp177 = mu(x, y, z) + mu(x+1, y, z).
 real_t temp177 = temp6 + temp7;

 // temp178 = (2 / (mu(x, y, z) + mu(x+1, y, z))).
 real_t temp178 = temp4 / temp177;

 // temp179 = (2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t().
 real_t temp179 = temp178 * temp1;

 // temp180 = (((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h).
 real_t temp180 = temp179 / temp2;

 // temp181 = (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z)).
 real_t temp181 = context.vel_y->readElem(t+1, x, y, z+1, __LINE__) - temp29;

 // temp182 = 1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z)).
 real_t temp182 = temp19 * temp181;

 // temp183 = -0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1)).
 real_t temp183 = temp23 * (context.vel_y->readElem(t+1, x, y, z+2, __LINE__) - context.vel_y->readElem(t+1, x, y, z-1, __LINE__));

 // temp184 = (1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))).
 real_t temp184 = temp182 + temp183;

 // temp185 = (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z)).
 real_t temp185 = context.vel_z->readElem(t+1, x, y+1, z, __LINE__) - temp35;

 // temp186 = 1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z)).
 real_t temp186 = temp19 * temp185;

 // temp187 = (1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))).
 real_t temp187 = temp184 + temp186;

 // temp188 = -0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)).
 real_t temp188 = temp23 * (context.vel_z->readElem(t+1, x, y+2, z, __LINE__) - context.vel_z->readElem(t+1, x, y-1, z, __LINE__));

 // temp189 = (1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))).
 real_t temp189 = temp187 + temp188;

 // temp190 = (((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))).
 real_t temp190 = temp180 * temp189;

 // temp191 = stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))).
 real_t temp191 = context.stress_yz->readElem(t, x, y, z, __LINE__) + temp190;

 // temp192 = stress_mem_yz(t, x, y, z).
 real_t temp192 = context.stress_mem_yz->readElem(t, x, y, z, __LINE__);

 // temp193 = tau2(x, y, z) * stress_mem_yz(t, x, y, z).
 real_t temp193 = temp45 * temp192;

 // temp194 = (0.5 / h) * (1 - tau2(x, y, z)).
 real_t temp194 = temp137 * temp50;

 // temp195 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_t temp195 = temp194 * temp52;

 // temp196 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))).
 real_t temp196 = temp195 * temp178;

 // temp197 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z).
 real_t temp197 = temp196 * context.anelastic_yz->readElem(x, y, z, __LINE__);

 // temp198 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))).
 real_t temp198 = temp197 * temp189;

 // temp199 = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))).
 real_t temp199 = temp193 - temp198;

 // temp200 = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z).
 real_t temp200 = temp199 + temp192;

 // temp201 = delta_t() * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)).
 real_t temp201 = temp1 * temp200;

 // temp202 = stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z))).
 real_t temp202 = temp191 + temp201;

 // temp203 = (stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)))) * sponge(x, y, z).
 real_t temp203 = temp202 * temp71;

 // temp204 = ((stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)))) * sponge(x, y, z)).
 real_t temp204 = temp203;

 // Save result to stress_yz(t+1, x, y, z):
 context.stress_yz->writeElem(temp204, t+1, x, y, z, __LINE__);

 // temp205 = (tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp205 = temp47 + temp66;

 // temp206 = ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 real_t temp206 = temp205;

 // Save result to stress_mem_xx(t+1, x, y, z):
 context.stress_mem_xx->writeElem(temp206, t+1, x, y, z, __LINE__);

 // temp207 = (tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp207 = temp81 + temp90;

 // temp208 = ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 real_t temp208 = temp207;

 // Save result to stress_mem_yy(t+1, x, y, z):
 context.stress_mem_yy->writeElem(temp208, t+1, x, y, z, __LINE__);

 // temp209 = (tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_t temp209 = temp104 + temp113;

 // temp210 = ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 real_t temp210 = temp209;

 // Save result to stress_mem_zz(t+1, x, y, z):
 context.stress_mem_zz->writeElem(temp210, t+1, x, y, z, __LINE__);

 // temp211 = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))).
 real_t temp211 = temp143;

 // Save result to stress_mem_xy(t+1, x, y, z):
 context.stress_mem_xy->writeElem(temp211, t+1, x, y, z, __LINE__);

 // temp212 = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))).
 real_t temp212 = temp171;

 // Save result to stress_mem_xz(t+1, x, y, z):
 context.stress_mem_xz->writeElem(temp212, t+1, x, y, z, __LINE__);

 // temp213 = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))).
 real_t temp213 = temp199;

 // Save result to stress_mem_yz(t+1, x, y, z):
 context.stress_mem_yz->writeElem(temp213, t+1, x, y, z, __LINE__);
} // scalar calculation.

 // Calculate 16 result(s) relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 // SIMD calculations use 66 vector block(s) created from 60 aligned vector-block(s).
 // There are 3072 FP operation(s) per cluster.
 void calc_cluster(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {

 // Un-normalized indices.
 idx_t t = tv;
 idx_t x = xv * 4;
 idx_t y = yv * 4;
 idx_t z = zv * 1;

 // Read aligned vector block from stress_xx at t, x, y, z.
 real_vec_t temp_vec1 = context.stress_xx->readVecNorm(tv, xv, yv, zv, __LINE__);

 // temp_vec2 = delta_t().
 real_vec_t temp_vec2 = (*context.delta_t)();

 // temp_vec3 = h().
 real_vec_t temp_vec3 = (*context.h)();

 // temp_vec4 = (delta_t / h).
 real_vec_t temp_vec4 = temp_vec2 / temp_vec3;

 // Read aligned vector block from mu at x, y, z.
 real_vec_t temp_vec5 = context.mu->readVecNorm(xv, yv, zv, __LINE__);

 // Read aligned vector block from mu at x+4, y, z.
 real_vec_t temp_vec6 = context.mu->readVecNorm(xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from mu at x+1, y, z.
 real_vec_t temp_vec7;
 // temp_vec7[0] = temp_vec5[1];  // for x+1, y, z;
 // temp_vec7[1] = temp_vec5[2];  // for x+2, y, z;
 // temp_vec7[2] = temp_vec5[3];  // for x+3, y, z;
 // temp_vec7[3] = temp_vec6[0];  // for x+4, y, z;
 // temp_vec7[4] = temp_vec5[5];  // for x+1, y+1, z;
 // temp_vec7[5] = temp_vec5[6];  // for x+2, y+1, z;
 // temp_vec7[6] = temp_vec5[7];  // for x+3, y+1, z;
 // temp_vec7[7] = temp_vec6[4];  // for x+4, y+1, z;
 // temp_vec7[8] = temp_vec5[9];  // for x+1, y+2, z;
 // temp_vec7[9] = temp_vec5[10];  // for x+2, y+2, z;
 // temp_vec7[10] = temp_vec5[11];  // for x+3, y+2, z;
 // temp_vec7[11] = temp_vec6[8];  // for x+4, y+2, z;
 // temp_vec7[12] = temp_vec5[13];  // for x+1, y+3, z;
 // temp_vec7[13] = temp_vec5[14];  // for x+2, y+3, z;
 // temp_vec7[14] = temp_vec5[15];  // for x+3, y+3, z;
 // temp_vec7[15] = temp_vec6[12];  // for x+4, y+3, z;
 const real_vec_t_data ctrl_data_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12 = { .ci = { 1, 2, 3, ctrl_sel_bit |0, 5, 6, 7, ctrl_sel_bit |4, 9, 10, 11, ctrl_sel_bit |8, 13, 14, 15, ctrl_sel_bit |12 } };
 const real_vec_t ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12(ctrl_data_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12);
 real_vec_permute2(temp_vec7, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec5, temp_vec6);

 // Read aligned vector block from mu at x, y-4, z.
 real_vec_t temp_vec8 = context.mu->readVecNorm(xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from mu at x, y-1, z.
 real_vec_t temp_vec9;
 // temp_vec9[0] = temp_vec8[12];  // for x, y-1, z;
 // temp_vec9[1] = temp_vec8[13];  // for x+1, y-1, z;
 // temp_vec9[2] = temp_vec8[14];  // for x+2, y-1, z;
 // temp_vec9[3] = temp_vec8[15];  // for x+3, y-1, z;
 // temp_vec9[4] = temp_vec5[0];  // for x, y, z;
 // temp_vec9[5] = temp_vec5[1];  // for x+1, y, z;
 // temp_vec9[6] = temp_vec5[2];  // for x+2, y, z;
 // temp_vec9[7] = temp_vec5[3];  // for x+3, y, z;
 // temp_vec9[8] = temp_vec5[4];  // for x, y+1, z;
 // temp_vec9[9] = temp_vec5[5];  // for x+1, y+1, z;
 // temp_vec9[10] = temp_vec5[6];  // for x+2, y+1, z;
 // temp_vec9[11] = temp_vec5[7];  // for x+3, y+1, z;
 // temp_vec9[12] = temp_vec5[8];  // for x, y+2, z;
 // temp_vec9[13] = temp_vec5[9];  // for x+1, y+2, z;
 // temp_vec9[14] = temp_vec5[10];  // for x+2, y+2, z;
 // temp_vec9[15] = temp_vec5[11];  // for x+3, y+2, z;
 // Get 12 element(s) from temp_vec5 and 4 from temp_vec8.
 real_vec_align<12>(temp_vec9, temp_vec5, temp_vec8);

 // Read aligned vector block from mu at x+4, y-4, z.
 real_vec_t temp_vec10 = context.mu->readVecNorm(xv+(4/4), yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from mu at x+1, y-1, z.
 real_vec_t temp_vec11;
 // temp_vec11[0] = temp_vec8[13];  // for x+1, y-1, z;
 // temp_vec11[1] = temp_vec8[14];  // for x+2, y-1, z;
 // temp_vec11[2] = temp_vec8[15];  // for x+3, y-1, z;
 // temp_vec11[3] = temp_vec10[12];  // for x+4, y-1, z;
 // temp_vec11[4] = temp_vec5[1];  // for x+1, y, z;
 // temp_vec11[5] = temp_vec5[2];  // for x+2, y, z;
 // temp_vec11[6] = temp_vec5[3];  // for x+3, y, z;
 // temp_vec11[7] = temp_vec6[0];  // for x+4, y, z;
 // temp_vec11[8] = temp_vec5[5];  // for x+1, y+1, z;
 // temp_vec11[9] = temp_vec5[6];  // for x+2, y+1, z;
 // temp_vec11[10] = temp_vec5[7];  // for x+3, y+1, z;
 // temp_vec11[11] = temp_vec6[4];  // for x+4, y+1, z;
 // temp_vec11[12] = temp_vec5[9];  // for x+1, y+2, z;
 // temp_vec11[13] = temp_vec5[10];  // for x+2, y+2, z;
 // temp_vec11[14] = temp_vec5[11];  // for x+3, y+2, z;
 // temp_vec11[15] = temp_vec6[8];  // for x+4, y+2, z;
 // Get 9 element(s) from temp_vec5 and 3 from temp_vec8.
 real_vec_align<13>(temp_vec11, temp_vec5, temp_vec8);
 // Get 3 element(s) from temp_vec6 and 1 from temp_vec10.
 real_vec_align_masked<9>(temp_vec11, temp_vec6, temp_vec10, 0x8888);

 // Read aligned vector block from mu at x, y, z-1.
 real_vec_t temp_vec12 = context.mu->readVecNorm(xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from mu at x+4, y, z-1.
 real_vec_t temp_vec13 = context.mu->readVecNorm(xv+(4/4), yv, zv-(1/1), __LINE__);

 // Construct unaligned vector block from mu at x+1, y, z-1.
 real_vec_t temp_vec14;
 // temp_vec14[0] = temp_vec12[1];  // for x+1, y, z-1;
 // temp_vec14[1] = temp_vec12[2];  // for x+2, y, z-1;
 // temp_vec14[2] = temp_vec12[3];  // for x+3, y, z-1;
 // temp_vec14[3] = temp_vec13[0];  // for x+4, y, z-1;
 // temp_vec14[4] = temp_vec12[5];  // for x+1, y+1, z-1;
 // temp_vec14[5] = temp_vec12[6];  // for x+2, y+1, z-1;
 // temp_vec14[6] = temp_vec12[7];  // for x+3, y+1, z-1;
 // temp_vec14[7] = temp_vec13[4];  // for x+4, y+1, z-1;
 // temp_vec14[8] = temp_vec12[9];  // for x+1, y+2, z-1;
 // temp_vec14[9] = temp_vec12[10];  // for x+2, y+2, z-1;
 // temp_vec14[10] = temp_vec12[11];  // for x+3, y+2, z-1;
 // temp_vec14[11] = temp_vec13[8];  // for x+4, y+2, z-1;
 // temp_vec14[12] = temp_vec12[13];  // for x+1, y+3, z-1;
 // temp_vec14[13] = temp_vec12[14];  // for x+2, y+3, z-1;
 // temp_vec14[14] = temp_vec12[15];  // for x+3, y+3, z-1;
 // temp_vec14[15] = temp_vec13[12];  // for x+4, y+3, z-1;
 real_vec_permute2(temp_vec14, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec12, temp_vec13);

 // Read aligned vector block from mu at x, y-4, z-1.
 real_vec_t temp_vec15 = context.mu->readVecNorm(xv, yv-(4/4), zv-(1/1), __LINE__);

 // Construct unaligned vector block from mu at x, y-1, z-1.
 real_vec_t temp_vec16;
 // temp_vec16[0] = temp_vec15[12];  // for x, y-1, z-1;
 // temp_vec16[1] = temp_vec15[13];  // for x+1, y-1, z-1;
 // temp_vec16[2] = temp_vec15[14];  // for x+2, y-1, z-1;
 // temp_vec16[3] = temp_vec15[15];  // for x+3, y-1, z-1;
 // temp_vec16[4] = temp_vec12[0];  // for x, y, z-1;
 // temp_vec16[5] = temp_vec12[1];  // for x+1, y, z-1;
 // temp_vec16[6] = temp_vec12[2];  // for x+2, y, z-1;
 // temp_vec16[7] = temp_vec12[3];  // for x+3, y, z-1;
 // temp_vec16[8] = temp_vec12[4];  // for x, y+1, z-1;
 // temp_vec16[9] = temp_vec12[5];  // for x+1, y+1, z-1;
 // temp_vec16[10] = temp_vec12[6];  // for x+2, y+1, z-1;
 // temp_vec16[11] = temp_vec12[7];  // for x+3, y+1, z-1;
 // temp_vec16[12] = temp_vec12[8];  // for x, y+2, z-1;
 // temp_vec16[13] = temp_vec12[9];  // for x+1, y+2, z-1;
 // temp_vec16[14] = temp_vec12[10];  // for x+2, y+2, z-1;
 // temp_vec16[15] = temp_vec12[11];  // for x+3, y+2, z-1;
 // Get 12 element(s) from temp_vec12 and 4 from temp_vec15.
 real_vec_align<12>(temp_vec16, temp_vec12, temp_vec15);

 // Read aligned vector block from mu at x+4, y-4, z-1.
 real_vec_t temp_vec17 = context.mu->readVecNorm(xv+(4/4), yv-(4/4), zv-(1/1), __LINE__);

 // Construct unaligned vector block from mu at x+1, y-1, z-1.
 real_vec_t temp_vec18;
 // temp_vec18[0] = temp_vec15[13];  // for x+1, y-1, z-1;
 // temp_vec18[1] = temp_vec15[14];  // for x+2, y-1, z-1;
 // temp_vec18[2] = temp_vec15[15];  // for x+3, y-1, z-1;
 // temp_vec18[3] = temp_vec17[12];  // for x+4, y-1, z-1;
 // temp_vec18[4] = temp_vec12[1];  // for x+1, y, z-1;
 // temp_vec18[5] = temp_vec12[2];  // for x+2, y, z-1;
 // temp_vec18[6] = temp_vec12[3];  // for x+3, y, z-1;
 // temp_vec18[7] = temp_vec13[0];  // for x+4, y, z-1;
 // temp_vec18[8] = temp_vec12[5];  // for x+1, y+1, z-1;
 // temp_vec18[9] = temp_vec12[6];  // for x+2, y+1, z-1;
 // temp_vec18[10] = temp_vec12[7];  // for x+3, y+1, z-1;
 // temp_vec18[11] = temp_vec13[4];  // for x+4, y+1, z-1;
 // temp_vec18[12] = temp_vec12[9];  // for x+1, y+2, z-1;
 // temp_vec18[13] = temp_vec12[10];  // for x+2, y+2, z-1;
 // temp_vec18[14] = temp_vec12[11];  // for x+3, y+2, z-1;
 // temp_vec18[15] = temp_vec13[8];  // for x+4, y+2, z-1;
 // Get 9 element(s) from temp_vec12 and 3 from temp_vec15.
 real_vec_align<13>(temp_vec18, temp_vec12, temp_vec15);
 // Get 3 element(s) from temp_vec13 and 1 from temp_vec17.
 real_vec_align_masked<9>(temp_vec18, temp_vec13, temp_vec17, 0x8888);

 // Read aligned vector block from vel_x at t+1, x, y, z.
 real_vec_t temp_vec19 = context.vel_x->readVecNorm(tv+(1/1), xv, yv, zv, __LINE__);

 // Read aligned vector block from vel_x at t+1, x+4, y, z.
 real_vec_t temp_vec20 = context.vel_x->readVecNorm(tv+(1/1), xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from vel_x at t+1, x+1, y, z.
 real_vec_t temp_vec21;
 // temp_vec21[0] = temp_vec19[1];  // for t+1, x+1, y, z;
 // temp_vec21[1] = temp_vec19[2];  // for t+1, x+2, y, z;
 // temp_vec21[2] = temp_vec19[3];  // for t+1, x+3, y, z;
 // temp_vec21[3] = temp_vec20[0];  // for t+1, x+4, y, z;
 // temp_vec21[4] = temp_vec19[5];  // for t+1, x+1, y+1, z;
 // temp_vec21[5] = temp_vec19[6];  // for t+1, x+2, y+1, z;
 // temp_vec21[6] = temp_vec19[7];  // for t+1, x+3, y+1, z;
 // temp_vec21[7] = temp_vec20[4];  // for t+1, x+4, y+1, z;
 // temp_vec21[8] = temp_vec19[9];  // for t+1, x+1, y+2, z;
 // temp_vec21[9] = temp_vec19[10];  // for t+1, x+2, y+2, z;
 // temp_vec21[10] = temp_vec19[11];  // for t+1, x+3, y+2, z;
 // temp_vec21[11] = temp_vec20[8];  // for t+1, x+4, y+2, z;
 // temp_vec21[12] = temp_vec19[13];  // for t+1, x+1, y+3, z;
 // temp_vec21[13] = temp_vec19[14];  // for t+1, x+2, y+3, z;
 // temp_vec21[14] = temp_vec19[15];  // for t+1, x+3, y+3, z;
 // temp_vec21[15] = temp_vec20[12];  // for t+1, x+4, y+3, z;
 real_vec_permute2(temp_vec21, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec19, temp_vec20);

 // Construct unaligned vector block from vel_x at t+1, x+2, y, z.
 real_vec_t temp_vec22;
 // temp_vec22[0] = temp_vec19[2];  // for t+1, x+2, y, z;
 // temp_vec22[1] = temp_vec19[3];  // for t+1, x+3, y, z;
 // temp_vec22[2] = temp_vec20[0];  // for t+1, x+4, y, z;
 // temp_vec22[3] = temp_vec20[1];  // for t+1, x+5, y, z;
 // temp_vec22[4] = temp_vec19[6];  // for t+1, x+2, y+1, z;
 // temp_vec22[5] = temp_vec19[7];  // for t+1, x+3, y+1, z;
 // temp_vec22[6] = temp_vec20[4];  // for t+1, x+4, y+1, z;
 // temp_vec22[7] = temp_vec20[5];  // for t+1, x+5, y+1, z;
 // temp_vec22[8] = temp_vec19[10];  // for t+1, x+2, y+2, z;
 // temp_vec22[9] = temp_vec19[11];  // for t+1, x+3, y+2, z;
 // temp_vec22[10] = temp_vec20[8];  // for t+1, x+4, y+2, z;
 // temp_vec22[11] = temp_vec20[9];  // for t+1, x+5, y+2, z;
 // temp_vec22[12] = temp_vec19[14];  // for t+1, x+2, y+3, z;
 // temp_vec22[13] = temp_vec19[15];  // for t+1, x+3, y+3, z;
 // temp_vec22[14] = temp_vec20[12];  // for t+1, x+4, y+3, z;
 // temp_vec22[15] = temp_vec20[13];  // for t+1, x+5, y+3, z;
 const real_vec_t_data ctrl_data_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13 = { .ci = { 2, 3, ctrl_sel_bit |0, ctrl_sel_bit |1, 6, 7, ctrl_sel_bit |4, ctrl_sel_bit |5, 10, 11, ctrl_sel_bit |8, ctrl_sel_bit |9, 14, 15, ctrl_sel_bit |12, ctrl_sel_bit |13 } };
 const real_vec_t ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13(ctrl_data_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13);
 real_vec_permute2(temp_vec22, ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13, temp_vec19, temp_vec20);

 // Read aligned vector block from vel_x at t+1, x-4, y, z.
 real_vec_t temp_vec23 = context.vel_x->readVecNorm(tv+(1/1), xv-(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from vel_x at t+1, x-1, y, z.
 real_vec_t temp_vec24;
 // temp_vec24[0] = temp_vec23[3];  // for t+1, x-1, y, z;
 // temp_vec24[1] = temp_vec19[0];  // for t+1, x, y, z;
 // temp_vec24[2] = temp_vec19[1];  // for t+1, x+1, y, z;
 // temp_vec24[3] = temp_vec19[2];  // for t+1, x+2, y, z;
 // temp_vec24[4] = temp_vec23[7];  // for t+1, x-1, y+1, z;
 // temp_vec24[5] = temp_vec19[4];  // for t+1, x, y+1, z;
 // temp_vec24[6] = temp_vec19[5];  // for t+1, x+1, y+1, z;
 // temp_vec24[7] = temp_vec19[6];  // for t+1, x+2, y+1, z;
 // temp_vec24[8] = temp_vec23[11];  // for t+1, x-1, y+2, z;
 // temp_vec24[9] = temp_vec19[8];  // for t+1, x, y+2, z;
 // temp_vec24[10] = temp_vec19[9];  // for t+1, x+1, y+2, z;
 // temp_vec24[11] = temp_vec19[10];  // for t+1, x+2, y+2, z;
 // temp_vec24[12] = temp_vec23[15];  // for t+1, x-1, y+3, z;
 // temp_vec24[13] = temp_vec19[12];  // for t+1, x, y+3, z;
 // temp_vec24[14] = temp_vec19[13];  // for t+1, x+1, y+3, z;
 // temp_vec24[15] = temp_vec19[14];  // for t+1, x+2, y+3, z;
 const real_vec_t_data ctrl_data_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14 = { .ci = { 3, ctrl_sel_bit |0, ctrl_sel_bit |1, ctrl_sel_bit |2, 7, ctrl_sel_bit |4, ctrl_sel_bit |5, ctrl_sel_bit |6, 11, ctrl_sel_bit |8, ctrl_sel_bit |9, ctrl_sel_bit |10, 15, ctrl_sel_bit |12, ctrl_sel_bit |13, ctrl_sel_bit |14 } };
 const real_vec_t ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14(ctrl_data_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14);
 real_vec_permute2(temp_vec24, ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14, temp_vec23, temp_vec19);

 // temp_vec25 = 2.
 real_vec_t temp_vec25 = 2.00000000000000000e+00;

 // temp_vec26 = 8.
 real_vec_t temp_vec26 = 8.00000000000000000e+00;

 // temp_vec27 = mu(x, y, z).
 real_vec_t temp_vec27 = temp_vec5;

 // temp_vec28 = mu(x+1, y, z).
 real_vec_t temp_vec28 = temp_vec7;

 // temp_vec29 = mu(x, y, z) + mu(x+1, y, z).
 real_vec_t temp_vec29 = temp_vec27 + temp_vec28;

 // temp_vec30 = mu(x, y-1, z).
 real_vec_t temp_vec30 = temp_vec9;

 // temp_vec31 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z).
 real_vec_t temp_vec31 = temp_vec29 + temp_vec30;

 // temp_vec32 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z).
 real_vec_t temp_vec32 = temp_vec31 + temp_vec11;

 // temp_vec33 = mu(x, y, z-1).
 real_vec_t temp_vec33 = temp_vec12;

 // temp_vec34 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1).
 real_vec_t temp_vec34 = temp_vec32 + temp_vec33;

 // temp_vec35 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1).
 real_vec_t temp_vec35 = temp_vec34 + temp_vec14;

 // temp_vec36 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1).
 real_vec_t temp_vec36 = temp_vec35 + temp_vec16;

 // temp_vec37 = mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1).
 real_vec_t temp_vec37 = temp_vec36 + temp_vec18;

 // temp_vec38 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_vec_t temp_vec38 = temp_vec26 / temp_vec37;

 // temp_vec39 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_vec_t temp_vec39 = temp_vec25 * temp_vec38;

 // temp_vec40 = 1.125.
 real_vec_t temp_vec40 = 1.12500000000000000e+00;

 // temp_vec41 = vel_x(t+1, x, y, z).
 real_vec_t temp_vec41 = temp_vec19;

 // temp_vec42 = (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z)).
 real_vec_t temp_vec42 = temp_vec21 - temp_vec41;

 // temp_vec43 = 1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z)).
 real_vec_t temp_vec43 = temp_vec40 * temp_vec42;

 // temp_vec44 = -0.0416667.
 real_vec_t temp_vec44 = -4.16666666666666644e-02;

 // temp_vec45 = -0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z)).
 real_vec_t temp_vec45 = temp_vec44 * (temp_vec22 - temp_vec24);

 // temp_vec46 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_vec_t temp_vec46 = temp_vec43 + temp_vec45;

 // temp_vec47 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z)))).
 real_vec_t temp_vec47 = temp_vec39 * temp_vec46;

 // Read aligned vector block from lambda at x, y, z.
 real_vec_t temp_vec48 = context.lambda->readVecNorm(xv, yv, zv, __LINE__);

 // Read aligned vector block from lambda at x+4, y, z.
 real_vec_t temp_vec49 = context.lambda->readVecNorm(xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from lambda at x+1, y, z.
 real_vec_t temp_vec50;
 // temp_vec50[0] = temp_vec48[1];  // for x+1, y, z;
 // temp_vec50[1] = temp_vec48[2];  // for x+2, y, z;
 // temp_vec50[2] = temp_vec48[3];  // for x+3, y, z;
 // temp_vec50[3] = temp_vec49[0];  // for x+4, y, z;
 // temp_vec50[4] = temp_vec48[5];  // for x+1, y+1, z;
 // temp_vec50[5] = temp_vec48[6];  // for x+2, y+1, z;
 // temp_vec50[6] = temp_vec48[7];  // for x+3, y+1, z;
 // temp_vec50[7] = temp_vec49[4];  // for x+4, y+1, z;
 // temp_vec50[8] = temp_vec48[9];  // for x+1, y+2, z;
 // temp_vec50[9] = temp_vec48[10];  // for x+2, y+2, z;
 // temp_vec50[10] = temp_vec48[11];  // for x+3, y+2, z;
 // temp_vec50[11] = temp_vec49[8];  // for x+4, y+2, z;
 // temp_vec50[12] = temp_vec48[13];  // for x+1, y+3, z;
 // temp_vec50[13] = temp_vec48[14];  // for x+2, y+3, z;
 // temp_vec50[14] = temp_vec48[15];  // for x+3, y+3, z;
 // temp_vec50[15] = temp_vec49[12];  // for x+4, y+3, z;
 real_vec_permute2(temp_vec50, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec48, temp_vec49);

 // Read aligned vector block from lambda at x, y-4, z.
 real_vec_t temp_vec51 = context.lambda->readVecNorm(xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from lambda at x, y-1, z.
 real_vec_t temp_vec52;
 // temp_vec52[0] = temp_vec51[12];  // for x, y-1, z;
 // temp_vec52[1] = temp_vec51[13];  // for x+1, y-1, z;
 // temp_vec52[2] = temp_vec51[14];  // for x+2, y-1, z;
 // temp_vec52[3] = temp_vec51[15];  // for x+3, y-1, z;
 // temp_vec52[4] = temp_vec48[0];  // for x, y, z;
 // temp_vec52[5] = temp_vec48[1];  // for x+1, y, z;
 // temp_vec52[6] = temp_vec48[2];  // for x+2, y, z;
 // temp_vec52[7] = temp_vec48[3];  // for x+3, y, z;
 // temp_vec52[8] = temp_vec48[4];  // for x, y+1, z;
 // temp_vec52[9] = temp_vec48[5];  // for x+1, y+1, z;
 // temp_vec52[10] = temp_vec48[6];  // for x+2, y+1, z;
 // temp_vec52[11] = temp_vec48[7];  // for x+3, y+1, z;
 // temp_vec52[12] = temp_vec48[8];  // for x, y+2, z;
 // temp_vec52[13] = temp_vec48[9];  // for x+1, y+2, z;
 // temp_vec52[14] = temp_vec48[10];  // for x+2, y+2, z;
 // temp_vec52[15] = temp_vec48[11];  // for x+3, y+2, z;
 // Get 12 element(s) from temp_vec48 and 4 from temp_vec51.
 real_vec_align<12>(temp_vec52, temp_vec48, temp_vec51);

 // Read aligned vector block from lambda at x+4, y-4, z.
 real_vec_t temp_vec53 = context.lambda->readVecNorm(xv+(4/4), yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from lambda at x+1, y-1, z.
 real_vec_t temp_vec54;
 // temp_vec54[0] = temp_vec51[13];  // for x+1, y-1, z;
 // temp_vec54[1] = temp_vec51[14];  // for x+2, y-1, z;
 // temp_vec54[2] = temp_vec51[15];  // for x+3, y-1, z;
 // temp_vec54[3] = temp_vec53[12];  // for x+4, y-1, z;
 // temp_vec54[4] = temp_vec48[1];  // for x+1, y, z;
 // temp_vec54[5] = temp_vec48[2];  // for x+2, y, z;
 // temp_vec54[6] = temp_vec48[3];  // for x+3, y, z;
 // temp_vec54[7] = temp_vec49[0];  // for x+4, y, z;
 // temp_vec54[8] = temp_vec48[5];  // for x+1, y+1, z;
 // temp_vec54[9] = temp_vec48[6];  // for x+2, y+1, z;
 // temp_vec54[10] = temp_vec48[7];  // for x+3, y+1, z;
 // temp_vec54[11] = temp_vec49[4];  // for x+4, y+1, z;
 // temp_vec54[12] = temp_vec48[9];  // for x+1, y+2, z;
 // temp_vec54[13] = temp_vec48[10];  // for x+2, y+2, z;
 // temp_vec54[14] = temp_vec48[11];  // for x+3, y+2, z;
 // temp_vec54[15] = temp_vec49[8];  // for x+4, y+2, z;
 // Get 9 element(s) from temp_vec48 and 3 from temp_vec51.
 real_vec_align<13>(temp_vec54, temp_vec48, temp_vec51);
 // Get 3 element(s) from temp_vec49 and 1 from temp_vec53.
 real_vec_align_masked<9>(temp_vec54, temp_vec49, temp_vec53, 0x8888);

 // Read aligned vector block from lambda at x, y, z-1.
 real_vec_t temp_vec55 = context.lambda->readVecNorm(xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from lambda at x+4, y, z-1.
 real_vec_t temp_vec56 = context.lambda->readVecNorm(xv+(4/4), yv, zv-(1/1), __LINE__);

 // Construct unaligned vector block from lambda at x+1, y, z-1.
 real_vec_t temp_vec57;
 // temp_vec57[0] = temp_vec55[1];  // for x+1, y, z-1;
 // temp_vec57[1] = temp_vec55[2];  // for x+2, y, z-1;
 // temp_vec57[2] = temp_vec55[3];  // for x+3, y, z-1;
 // temp_vec57[3] = temp_vec56[0];  // for x+4, y, z-1;
 // temp_vec57[4] = temp_vec55[5];  // for x+1, y+1, z-1;
 // temp_vec57[5] = temp_vec55[6];  // for x+2, y+1, z-1;
 // temp_vec57[6] = temp_vec55[7];  // for x+3, y+1, z-1;
 // temp_vec57[7] = temp_vec56[4];  // for x+4, y+1, z-1;
 // temp_vec57[8] = temp_vec55[9];  // for x+1, y+2, z-1;
 // temp_vec57[9] = temp_vec55[10];  // for x+2, y+2, z-1;
 // temp_vec57[10] = temp_vec55[11];  // for x+3, y+2, z-1;
 // temp_vec57[11] = temp_vec56[8];  // for x+4, y+2, z-1;
 // temp_vec57[12] = temp_vec55[13];  // for x+1, y+3, z-1;
 // temp_vec57[13] = temp_vec55[14];  // for x+2, y+3, z-1;
 // temp_vec57[14] = temp_vec55[15];  // for x+3, y+3, z-1;
 // temp_vec57[15] = temp_vec56[12];  // for x+4, y+3, z-1;
 real_vec_permute2(temp_vec57, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec55, temp_vec56);

 // Read aligned vector block from lambda at x, y-4, z-1.
 real_vec_t temp_vec58 = context.lambda->readVecNorm(xv, yv-(4/4), zv-(1/1), __LINE__);

 // Construct unaligned vector block from lambda at x, y-1, z-1.
 real_vec_t temp_vec59;
 // temp_vec59[0] = temp_vec58[12];  // for x, y-1, z-1;
 // temp_vec59[1] = temp_vec58[13];  // for x+1, y-1, z-1;
 // temp_vec59[2] = temp_vec58[14];  // for x+2, y-1, z-1;
 // temp_vec59[3] = temp_vec58[15];  // for x+3, y-1, z-1;
 // temp_vec59[4] = temp_vec55[0];  // for x, y, z-1;
 // temp_vec59[5] = temp_vec55[1];  // for x+1, y, z-1;
 // temp_vec59[6] = temp_vec55[2];  // for x+2, y, z-1;
 // temp_vec59[7] = temp_vec55[3];  // for x+3, y, z-1;
 // temp_vec59[8] = temp_vec55[4];  // for x, y+1, z-1;
 // temp_vec59[9] = temp_vec55[5];  // for x+1, y+1, z-1;
 // temp_vec59[10] = temp_vec55[6];  // for x+2, y+1, z-1;
 // temp_vec59[11] = temp_vec55[7];  // for x+3, y+1, z-1;
 // temp_vec59[12] = temp_vec55[8];  // for x, y+2, z-1;
 // temp_vec59[13] = temp_vec55[9];  // for x+1, y+2, z-1;
 // temp_vec59[14] = temp_vec55[10];  // for x+2, y+2, z-1;
 // temp_vec59[15] = temp_vec55[11];  // for x+3, y+2, z-1;
 // Get 12 element(s) from temp_vec55 and 4 from temp_vec58.
 real_vec_align<12>(temp_vec59, temp_vec55, temp_vec58);

 // Read aligned vector block from lambda at x+4, y-4, z-1.
 real_vec_t temp_vec60 = context.lambda->readVecNorm(xv+(4/4), yv-(4/4), zv-(1/1), __LINE__);

 // Construct unaligned vector block from lambda at x+1, y-1, z-1.
 real_vec_t temp_vec61;
 // temp_vec61[0] = temp_vec58[13];  // for x+1, y-1, z-1;
 // temp_vec61[1] = temp_vec58[14];  // for x+2, y-1, z-1;
 // temp_vec61[2] = temp_vec58[15];  // for x+3, y-1, z-1;
 // temp_vec61[3] = temp_vec60[12];  // for x+4, y-1, z-1;
 // temp_vec61[4] = temp_vec55[1];  // for x+1, y, z-1;
 // temp_vec61[5] = temp_vec55[2];  // for x+2, y, z-1;
 // temp_vec61[6] = temp_vec55[3];  // for x+3, y, z-1;
 // temp_vec61[7] = temp_vec56[0];  // for x+4, y, z-1;
 // temp_vec61[8] = temp_vec55[5];  // for x+1, y+1, z-1;
 // temp_vec61[9] = temp_vec55[6];  // for x+2, y+1, z-1;
 // temp_vec61[10] = temp_vec55[7];  // for x+3, y+1, z-1;
 // temp_vec61[11] = temp_vec56[4];  // for x+4, y+1, z-1;
 // temp_vec61[12] = temp_vec55[9];  // for x+1, y+2, z-1;
 // temp_vec61[13] = temp_vec55[10];  // for x+2, y+2, z-1;
 // temp_vec61[14] = temp_vec55[11];  // for x+3, y+2, z-1;
 // temp_vec61[15] = temp_vec56[8];  // for x+4, y+2, z-1;
 // Get 9 element(s) from temp_vec55 and 3 from temp_vec58.
 real_vec_align<13>(temp_vec61, temp_vec55, temp_vec58);
 // Get 3 element(s) from temp_vec56 and 1 from temp_vec60.
 real_vec_align_masked<9>(temp_vec61, temp_vec56, temp_vec60, 0x8888);

 // Read aligned vector block from vel_y at t+1, x, y, z.
 real_vec_t temp_vec62 = context.vel_y->readVecNorm(tv+(1/1), xv, yv, zv, __LINE__);

 // Read aligned vector block from vel_y at t+1, x, y-4, z.
 real_vec_t temp_vec63 = context.vel_y->readVecNorm(tv+(1/1), xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from vel_y at t+1, x, y-1, z.
 real_vec_t temp_vec64;
 // temp_vec64[0] = temp_vec63[12];  // for t+1, x, y-1, z;
 // temp_vec64[1] = temp_vec63[13];  // for t+1, x+1, y-1, z;
 // temp_vec64[2] = temp_vec63[14];  // for t+1, x+2, y-1, z;
 // temp_vec64[3] = temp_vec63[15];  // for t+1, x+3, y-1, z;
 // temp_vec64[4] = temp_vec62[0];  // for t+1, x, y, z;
 // temp_vec64[5] = temp_vec62[1];  // for t+1, x+1, y, z;
 // temp_vec64[6] = temp_vec62[2];  // for t+1, x+2, y, z;
 // temp_vec64[7] = temp_vec62[3];  // for t+1, x+3, y, z;
 // temp_vec64[8] = temp_vec62[4];  // for t+1, x, y+1, z;
 // temp_vec64[9] = temp_vec62[5];  // for t+1, x+1, y+1, z;
 // temp_vec64[10] = temp_vec62[6];  // for t+1, x+2, y+1, z;
 // temp_vec64[11] = temp_vec62[7];  // for t+1, x+3, y+1, z;
 // temp_vec64[12] = temp_vec62[8];  // for t+1, x, y+2, z;
 // temp_vec64[13] = temp_vec62[9];  // for t+1, x+1, y+2, z;
 // temp_vec64[14] = temp_vec62[10];  // for t+1, x+2, y+2, z;
 // temp_vec64[15] = temp_vec62[11];  // for t+1, x+3, y+2, z;
 // Get 12 element(s) from temp_vec62 and 4 from temp_vec63.
 real_vec_align<12>(temp_vec64, temp_vec62, temp_vec63);

 // Read aligned vector block from vel_y at t+1, x, y+4, z.
 real_vec_t temp_vec65 = context.vel_y->readVecNorm(tv+(1/1), xv, yv+(4/4), zv, __LINE__);

 // Construct unaligned vector block from vel_y at t+1, x, y+1, z.
 real_vec_t temp_vec66;
 // temp_vec66[0] = temp_vec62[4];  // for t+1, x, y+1, z;
 // temp_vec66[1] = temp_vec62[5];  // for t+1, x+1, y+1, z;
 // temp_vec66[2] = temp_vec62[6];  // for t+1, x+2, y+1, z;
 // temp_vec66[3] = temp_vec62[7];  // for t+1, x+3, y+1, z;
 // temp_vec66[4] = temp_vec62[8];  // for t+1, x, y+2, z;
 // temp_vec66[5] = temp_vec62[9];  // for t+1, x+1, y+2, z;
 // temp_vec66[6] = temp_vec62[10];  // for t+1, x+2, y+2, z;
 // temp_vec66[7] = temp_vec62[11];  // for t+1, x+3, y+2, z;
 // temp_vec66[8] = temp_vec62[12];  // for t+1, x, y+3, z;
 // temp_vec66[9] = temp_vec62[13];  // for t+1, x+1, y+3, z;
 // temp_vec66[10] = temp_vec62[14];  // for t+1, x+2, y+3, z;
 // temp_vec66[11] = temp_vec62[15];  // for t+1, x+3, y+3, z;
 // temp_vec66[12] = temp_vec65[0];  // for t+1, x, y+4, z;
 // temp_vec66[13] = temp_vec65[1];  // for t+1, x+1, y+4, z;
 // temp_vec66[14] = temp_vec65[2];  // for t+1, x+2, y+4, z;
 // temp_vec66[15] = temp_vec65[3];  // for t+1, x+3, y+4, z;
 // Get 4 element(s) from temp_vec65 and 12 from temp_vec62.
 real_vec_align<4>(temp_vec66, temp_vec65, temp_vec62);

 // Construct unaligned vector block from vel_y at t+1, x, y-2, z.
 real_vec_t temp_vec67;
 // temp_vec67[0] = temp_vec63[8];  // for t+1, x, y-2, z;
 // temp_vec67[1] = temp_vec63[9];  // for t+1, x+1, y-2, z;
 // temp_vec67[2] = temp_vec63[10];  // for t+1, x+2, y-2, z;
 // temp_vec67[3] = temp_vec63[11];  // for t+1, x+3, y-2, z;
 // temp_vec67[4] = temp_vec63[12];  // for t+1, x, y-1, z;
 // temp_vec67[5] = temp_vec63[13];  // for t+1, x+1, y-1, z;
 // temp_vec67[6] = temp_vec63[14];  // for t+1, x+2, y-1, z;
 // temp_vec67[7] = temp_vec63[15];  // for t+1, x+3, y-1, z;
 // temp_vec67[8] = temp_vec62[0];  // for t+1, x, y, z;
 // temp_vec67[9] = temp_vec62[1];  // for t+1, x+1, y, z;
 // temp_vec67[10] = temp_vec62[2];  // for t+1, x+2, y, z;
 // temp_vec67[11] = temp_vec62[3];  // for t+1, x+3, y, z;
 // temp_vec67[12] = temp_vec62[4];  // for t+1, x, y+1, z;
 // temp_vec67[13] = temp_vec62[5];  // for t+1, x+1, y+1, z;
 // temp_vec67[14] = temp_vec62[6];  // for t+1, x+2, y+1, z;
 // temp_vec67[15] = temp_vec62[7];  // for t+1, x+3, y+1, z;
 // Get 8 element(s) from temp_vec62 and 8 from temp_vec63.
 real_vec_align<8>(temp_vec67, temp_vec62, temp_vec63);

 // Read aligned vector block from vel_z at t+1, x, y, z.
 real_vec_t temp_vec68 = context.vel_z->readVecNorm(tv+(1/1), xv, yv, zv, __LINE__);

 // Read aligned vector block from vel_z at t+1, x, y, z-1.
 real_vec_t temp_vec69 = context.vel_z->readVecNorm(tv+(1/1), xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from vel_z at t+1, x, y, z+1.
 real_vec_t temp_vec70 = context.vel_z->readVecNorm(tv+(1/1), xv, yv, zv+(1/1), __LINE__);

 // Read aligned vector block from vel_z at t+1, x, y, z-2.
 real_vec_t temp_vec71 = context.vel_z->readVecNorm(tv+(1/1), xv, yv, zv-(2/1), __LINE__);

 // temp_vec72 = (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))).
 real_vec_t temp_vec72 = temp_vec26 / (temp_vec48 + temp_vec50 + temp_vec52 + temp_vec54 + temp_vec55 + temp_vec57 + temp_vec59 + temp_vec61);

 // temp_vec73 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_vec_t temp_vec73 = temp_vec43 + temp_vec45;

 // temp_vec74 = vel_y(t+1, x, y, z).
 real_vec_t temp_vec74 = temp_vec62;

 // temp_vec75 = (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z)).
 real_vec_t temp_vec75 = temp_vec74 - temp_vec64;

 // temp_vec76 = 1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z)).
 real_vec_t temp_vec76 = temp_vec40 * temp_vec75;

 // temp_vec77 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))).
 real_vec_t temp_vec77 = temp_vec73 + temp_vec76;

 // temp_vec78 = -0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z)).
 real_vec_t temp_vec78 = temp_vec44 * (temp_vec66 - temp_vec67);

 // temp_vec79 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_vec_t temp_vec79 = temp_vec77 + temp_vec78;

 // temp_vec80 = vel_z(t+1, x, y, z).
 real_vec_t temp_vec80 = temp_vec68;

 // temp_vec81 = (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1)).
 real_vec_t temp_vec81 = temp_vec80 - temp_vec69;

 // temp_vec82 = 1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1)).
 real_vec_t temp_vec82 = temp_vec40 * temp_vec81;

 // temp_vec83 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))).
 real_vec_t temp_vec83 = temp_vec79 + temp_vec82;

 // temp_vec84 = -0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)).
 real_vec_t temp_vec84 = temp_vec44 * (temp_vec70 - temp_vec71);

 // temp_vec85 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_vec_t temp_vec85 = temp_vec83 + temp_vec84;

 // temp_vec86 = (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_vec_t temp_vec86 = temp_vec72 * temp_vec85;

 // temp_vec87 = (2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))).
 real_vec_t temp_vec87 = temp_vec47 + temp_vec86;

 // temp_vec88 = (delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec88 = temp_vec4 * temp_vec87;

 // temp_vec89 = stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec89 = temp_vec1 + temp_vec88;

 // Read aligned vector block from tau2 at x, y, z.
 real_vec_t temp_vec90 = context.tau2->readVecNorm(xv, yv, zv, __LINE__);

 // Read aligned vector block from stress_mem_xx at t, x, y, z.
 real_vec_t temp_vec91 = context.stress_mem_xx->readVecNorm(tv, xv, yv, zv, __LINE__);

 // temp_vec92 = tau2(x, y, z).
 real_vec_t temp_vec92 = temp_vec90;

 // temp_vec93 = stress_mem_xx(t, x, y, z).
 real_vec_t temp_vec93 = temp_vec91;

 // temp_vec94 = tau2(x, y, z) * stress_mem_xx(t, x, y, z).
 real_vec_t temp_vec94 = temp_vec92 * temp_vec93;

 // temp_vec95 = 1.
 real_vec_t temp_vec95 = 1.00000000000000000e+00;

 // temp_vec96 = (1 / h).
 real_vec_t temp_vec96 = temp_vec95 / temp_vec3;

 // temp_vec97 = (1 - tau2(x, y, z)).
 real_vec_t temp_vec97 = temp_vec95 - temp_vec92;

 // temp_vec98 = (1 / h) * (1 - tau2(x, y, z)).
 real_vec_t temp_vec98 = temp_vec96 * temp_vec97;

 // Read aligned vector block from weight at x, y, z.
 real_vec_t temp_vec99 = context.weight->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec100 = weight(x, y, z).
 real_vec_t temp_vec100 = temp_vec99;

 // temp_vec101 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_vec_t temp_vec101 = temp_vec98 * temp_vec100;

 // Read aligned vector block from anelastic_as_diag at x, y, z.
 real_vec_t temp_vec102 = context.anelastic_as_diag->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec103 = anelastic_as_diag(x, y, z).
 real_vec_t temp_vec103 = temp_vec102;

 // temp_vec104 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z).
 real_vec_t temp_vec104 = temp_vec38 * temp_vec103;

 // temp_vec105 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_vec_t temp_vec105 = temp_vec76 + temp_vec78;

 // temp_vec106 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))).
 real_vec_t temp_vec106 = temp_vec105 + temp_vec82;

 // temp_vec107 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_vec_t temp_vec107 = temp_vec106 + temp_vec84;

 // temp_vec108 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_vec_t temp_vec108 = temp_vec104 * temp_vec107;

 // temp_vec109 = 0.5.
 real_vec_t temp_vec109 = 5.00000000000000000e-01;

 // temp_vec110 = 0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))).
 real_vec_t temp_vec110 = temp_vec109 * temp_vec72;

 // temp_vec111 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1)))).
 real_vec_t temp_vec111 = temp_vec38 + temp_vec110;

 // Read aligned vector block from anelastic_ap at x, y, z.
 real_vec_t temp_vec112 = context.anelastic_ap->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec113 = ((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z).
 real_vec_t temp_vec113 = temp_vec111 * temp_vec112;

 // temp_vec114 = ((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_vec_t temp_vec114 = temp_vec113 * temp_vec85;

 // temp_vec115 = (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec115 = temp_vec108 - temp_vec114;

 // temp_vec116 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec116 = temp_vec101 * temp_vec115;

 // temp_vec117 = (tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec117 = temp_vec94 + temp_vec116;

 // temp_vec118 = (tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z).
 real_vec_t temp_vec118 = temp_vec117 + temp_vec93;

 // temp_vec119 = delta_t() * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)).
 real_vec_t temp_vec119 = temp_vec2 * temp_vec118;

 // temp_vec120 = stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z))).
 real_vec_t temp_vec120 = temp_vec89 + temp_vec119;

 // Read aligned vector block from sponge at x, y, z.
 real_vec_t temp_vec121 = context.sponge->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec122 = sponge(x, y, z).
 real_vec_t temp_vec122 = temp_vec121;

 // temp_vec123 = (stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)))) * sponge(x, y, z).
 real_vec_t temp_vec123 = temp_vec120 * temp_vec122;

 // temp_vec124 = ((stress_xx(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_xx(t, x, y, z)))) * sponge(x, y, z)).
 real_vec_t temp_vec124 = temp_vec123;

 // Save result to stress_xx(t+1, x, y, z):
 
 // Write aligned vector block to stress_xx at t+1, x, y, z.
context.stress_xx->writeVecNorm(temp_vec124, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from stress_yy at t, x, y, z.
 real_vec_t temp_vec125 = context.stress_yy->readVecNorm(tv, xv, yv, zv, __LINE__);

 // temp_vec126 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_vec_t temp_vec126 = temp_vec25 * temp_vec38;

 // temp_vec127 = (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_vec_t temp_vec127 = temp_vec76 + temp_vec78;

 // temp_vec128 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z)))).
 real_vec_t temp_vec128 = temp_vec126 * temp_vec127;

 // temp_vec129 = (2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))).
 real_vec_t temp_vec129 = temp_vec128 + temp_vec86;

 // temp_vec130 = (delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec130 = temp_vec4 * temp_vec129;

 // temp_vec131 = stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec131 = temp_vec125 + temp_vec130;

 // Read aligned vector block from stress_mem_yy at t, x, y, z.
 real_vec_t temp_vec132 = context.stress_mem_yy->readVecNorm(tv, xv, yv, zv, __LINE__);

 // temp_vec133 = stress_mem_yy(t, x, y, z).
 real_vec_t temp_vec133 = temp_vec132;

 // temp_vec134 = tau2(x, y, z) * stress_mem_yy(t, x, y, z).
 real_vec_t temp_vec134 = temp_vec92 * temp_vec133;

 // temp_vec135 = (1 / h) * (1 - tau2(x, y, z)).
 real_vec_t temp_vec135 = temp_vec96 * temp_vec97;

 // temp_vec136 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_vec_t temp_vec136 = temp_vec135 * temp_vec100;

 // temp_vec137 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z).
 real_vec_t temp_vec137 = temp_vec38 * temp_vec103;

 // temp_vec138 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_vec_t temp_vec138 = temp_vec43 + temp_vec45;

 // temp_vec139 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))).
 real_vec_t temp_vec139 = temp_vec138 + temp_vec82;

 // temp_vec140 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_vec_t temp_vec140 = temp_vec139 + temp_vec84;

 // temp_vec141 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_vec_t temp_vec141 = temp_vec137 * temp_vec140;

 // temp_vec142 = (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec142 = temp_vec141 - temp_vec114;

 // temp_vec143 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec143 = temp_vec136 * temp_vec142;

 // temp_vec144 = (tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec144 = temp_vec134 + temp_vec143;

 // temp_vec145 = (tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z).
 real_vec_t temp_vec145 = temp_vec144 + temp_vec133;

 // temp_vec146 = delta_t() * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)).
 real_vec_t temp_vec146 = temp_vec2 * temp_vec145;

 // temp_vec147 = stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z))).
 real_vec_t temp_vec147 = temp_vec131 + temp_vec146;

 // temp_vec148 = (stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)))) * sponge(x, y, z).
 real_vec_t temp_vec148 = temp_vec147 * temp_vec122;

 // temp_vec149 = ((stress_yy(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_yy(t, x, y, z)))) * sponge(x, y, z)).
 real_vec_t temp_vec149 = temp_vec148;

 // Save result to stress_yy(t+1, x, y, z):
 
 // Write aligned vector block to stress_yy at t+1, x, y, z.
context.stress_yy->writeVecNorm(temp_vec149, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from stress_zz at t, x, y, z.
 real_vec_t temp_vec150 = context.stress_zz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // temp_vec151 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))).
 real_vec_t temp_vec151 = temp_vec25 * temp_vec38;

 // temp_vec152 = (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))).
 real_vec_t temp_vec152 = temp_vec82 + temp_vec84;

 // temp_vec153 = 2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))).
 real_vec_t temp_vec153 = temp_vec151 * temp_vec152;

 // temp_vec154 = (2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))).
 real_vec_t temp_vec154 = temp_vec153 + temp_vec86;

 // temp_vec155 = (delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec155 = temp_vec4 * temp_vec154;

 // temp_vec156 = stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec156 = temp_vec150 + temp_vec155;

 // Read aligned vector block from stress_mem_zz at t, x, y, z.
 real_vec_t temp_vec157 = context.stress_mem_zz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // temp_vec158 = stress_mem_zz(t, x, y, z).
 real_vec_t temp_vec158 = temp_vec157;

 // temp_vec159 = tau2(x, y, z) * stress_mem_zz(t, x, y, z).
 real_vec_t temp_vec159 = temp_vec92 * temp_vec158;

 // temp_vec160 = (1 / h) * (1 - tau2(x, y, z)).
 real_vec_t temp_vec160 = temp_vec96 * temp_vec97;

 // temp_vec161 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_vec_t temp_vec161 = temp_vec160 * temp_vec100;

 // temp_vec162 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z).
 real_vec_t temp_vec162 = temp_vec38 * temp_vec103;

 // temp_vec163 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))).
 real_vec_t temp_vec163 = temp_vec43 + temp_vec45;

 // temp_vec164 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))).
 real_vec_t temp_vec164 = temp_vec163 + temp_vec76;

 // temp_vec165 = (1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))).
 real_vec_t temp_vec165 = temp_vec164 + temp_vec78;

 // temp_vec166 = (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z)))).
 real_vec_t temp_vec166 = temp_vec162 * temp_vec165;

 // temp_vec167 = (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec167 = temp_vec166 - temp_vec114;

 // temp_vec168 = (1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))).
 real_vec_t temp_vec168 = temp_vec161 * temp_vec167;

 // temp_vec169 = (tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec169 = temp_vec159 + temp_vec168;

 // temp_vec170 = (tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z).
 real_vec_t temp_vec170 = temp_vec169 + temp_vec158;

 // temp_vec171 = delta_t() * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)).
 real_vec_t temp_vec171 = temp_vec2 * temp_vec170;

 // temp_vec172 = stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z))).
 real_vec_t temp_vec172 = temp_vec156 + temp_vec171;

 // temp_vec173 = (stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)))) * sponge(x, y, z).
 real_vec_t temp_vec173 = temp_vec172 * temp_vec122;

 // temp_vec174 = ((stress_zz(t, x, y, z) + ((delta_t / h) * ((2 * (8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * ((1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) + ((8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + (delta_t * ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))) + stress_mem_zz(t, x, y, z)))) * sponge(x, y, z)).
 real_vec_t temp_vec174 = temp_vec173;

 // Save result to stress_zz(t+1, x, y, z):
 
 // Write aligned vector block to stress_zz at t+1, x, y, z.
context.stress_zz->writeVecNorm(temp_vec174, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from stress_xy at t, x, y, z.
 real_vec_t temp_vec175 = context.stress_xy->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from vel_x at t+1, x, y+4, z.
 real_vec_t temp_vec176 = context.vel_x->readVecNorm(tv+(1/1), xv, yv+(4/4), zv, __LINE__);

 // Construct unaligned vector block from vel_x at t+1, x, y+1, z.
 real_vec_t temp_vec177;
 // temp_vec177[0] = temp_vec19[4];  // for t+1, x, y+1, z;
 // temp_vec177[1] = temp_vec19[5];  // for t+1, x+1, y+1, z;
 // temp_vec177[2] = temp_vec19[6];  // for t+1, x+2, y+1, z;
 // temp_vec177[3] = temp_vec19[7];  // for t+1, x+3, y+1, z;
 // temp_vec177[4] = temp_vec19[8];  // for t+1, x, y+2, z;
 // temp_vec177[5] = temp_vec19[9];  // for t+1, x+1, y+2, z;
 // temp_vec177[6] = temp_vec19[10];  // for t+1, x+2, y+2, z;
 // temp_vec177[7] = temp_vec19[11];  // for t+1, x+3, y+2, z;
 // temp_vec177[8] = temp_vec19[12];  // for t+1, x, y+3, z;
 // temp_vec177[9] = temp_vec19[13];  // for t+1, x+1, y+3, z;
 // temp_vec177[10] = temp_vec19[14];  // for t+1, x+2, y+3, z;
 // temp_vec177[11] = temp_vec19[15];  // for t+1, x+3, y+3, z;
 // temp_vec177[12] = temp_vec176[0];  // for t+1, x, y+4, z;
 // temp_vec177[13] = temp_vec176[1];  // for t+1, x+1, y+4, z;
 // temp_vec177[14] = temp_vec176[2];  // for t+1, x+2, y+4, z;
 // temp_vec177[15] = temp_vec176[3];  // for t+1, x+3, y+4, z;
 // Get 4 element(s) from temp_vec176 and 12 from temp_vec19.
 real_vec_align<4>(temp_vec177, temp_vec176, temp_vec19);

 // Construct unaligned vector block from vel_x at t+1, x, y+2, z.
 real_vec_t temp_vec178;
 // temp_vec178[0] = temp_vec19[8];  // for t+1, x, y+2, z;
 // temp_vec178[1] = temp_vec19[9];  // for t+1, x+1, y+2, z;
 // temp_vec178[2] = temp_vec19[10];  // for t+1, x+2, y+2, z;
 // temp_vec178[3] = temp_vec19[11];  // for t+1, x+3, y+2, z;
 // temp_vec178[4] = temp_vec19[12];  // for t+1, x, y+3, z;
 // temp_vec178[5] = temp_vec19[13];  // for t+1, x+1, y+3, z;
 // temp_vec178[6] = temp_vec19[14];  // for t+1, x+2, y+3, z;
 // temp_vec178[7] = temp_vec19[15];  // for t+1, x+3, y+3, z;
 // temp_vec178[8] = temp_vec176[0];  // for t+1, x, y+4, z;
 // temp_vec178[9] = temp_vec176[1];  // for t+1, x+1, y+4, z;
 // temp_vec178[10] = temp_vec176[2];  // for t+1, x+2, y+4, z;
 // temp_vec178[11] = temp_vec176[3];  // for t+1, x+3, y+4, z;
 // temp_vec178[12] = temp_vec176[4];  // for t+1, x, y+5, z;
 // temp_vec178[13] = temp_vec176[5];  // for t+1, x+1, y+5, z;
 // temp_vec178[14] = temp_vec176[6];  // for t+1, x+2, y+5, z;
 // temp_vec178[15] = temp_vec176[7];  // for t+1, x+3, y+5, z;
 // Get 8 element(s) from temp_vec176 and 8 from temp_vec19.
 real_vec_align<8>(temp_vec178, temp_vec176, temp_vec19);

 // Read aligned vector block from vel_x at t+1, x, y-4, z.
 real_vec_t temp_vec179 = context.vel_x->readVecNorm(tv+(1/1), xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from vel_x at t+1, x, y-1, z.
 real_vec_t temp_vec180;
 // temp_vec180[0] = temp_vec179[12];  // for t+1, x, y-1, z;
 // temp_vec180[1] = temp_vec179[13];  // for t+1, x+1, y-1, z;
 // temp_vec180[2] = temp_vec179[14];  // for t+1, x+2, y-1, z;
 // temp_vec180[3] = temp_vec179[15];  // for t+1, x+3, y-1, z;
 // temp_vec180[4] = temp_vec19[0];  // for t+1, x, y, z;
 // temp_vec180[5] = temp_vec19[1];  // for t+1, x+1, y, z;
 // temp_vec180[6] = temp_vec19[2];  // for t+1, x+2, y, z;
 // temp_vec180[7] = temp_vec19[3];  // for t+1, x+3, y, z;
 // temp_vec180[8] = temp_vec19[4];  // for t+1, x, y+1, z;
 // temp_vec180[9] = temp_vec19[5];  // for t+1, x+1, y+1, z;
 // temp_vec180[10] = temp_vec19[6];  // for t+1, x+2, y+1, z;
 // temp_vec180[11] = temp_vec19[7];  // for t+1, x+3, y+1, z;
 // temp_vec180[12] = temp_vec19[8];  // for t+1, x, y+2, z;
 // temp_vec180[13] = temp_vec19[9];  // for t+1, x+1, y+2, z;
 // temp_vec180[14] = temp_vec19[10];  // for t+1, x+2, y+2, z;
 // temp_vec180[15] = temp_vec19[11];  // for t+1, x+3, y+2, z;
 // Get 12 element(s) from temp_vec19 and 4 from temp_vec179.
 real_vec_align<12>(temp_vec180, temp_vec19, temp_vec179);

 // Read aligned vector block from vel_y at t+1, x-4, y, z.
 real_vec_t temp_vec181 = context.vel_y->readVecNorm(tv+(1/1), xv-(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from vel_y at t+1, x-1, y, z.
 real_vec_t temp_vec182;
 // temp_vec182[0] = temp_vec181[3];  // for t+1, x-1, y, z;
 // temp_vec182[1] = temp_vec62[0];  // for t+1, x, y, z;
 // temp_vec182[2] = temp_vec62[1];  // for t+1, x+1, y, z;
 // temp_vec182[3] = temp_vec62[2];  // for t+1, x+2, y, z;
 // temp_vec182[4] = temp_vec181[7];  // for t+1, x-1, y+1, z;
 // temp_vec182[5] = temp_vec62[4];  // for t+1, x, y+1, z;
 // temp_vec182[6] = temp_vec62[5];  // for t+1, x+1, y+1, z;
 // temp_vec182[7] = temp_vec62[6];  // for t+1, x+2, y+1, z;
 // temp_vec182[8] = temp_vec181[11];  // for t+1, x-1, y+2, z;
 // temp_vec182[9] = temp_vec62[8];  // for t+1, x, y+2, z;
 // temp_vec182[10] = temp_vec62[9];  // for t+1, x+1, y+2, z;
 // temp_vec182[11] = temp_vec62[10];  // for t+1, x+2, y+2, z;
 // temp_vec182[12] = temp_vec181[15];  // for t+1, x-1, y+3, z;
 // temp_vec182[13] = temp_vec62[12];  // for t+1, x, y+3, z;
 // temp_vec182[14] = temp_vec62[13];  // for t+1, x+1, y+3, z;
 // temp_vec182[15] = temp_vec62[14];  // for t+1, x+2, y+3, z;
 real_vec_permute2(temp_vec182, ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14, temp_vec181, temp_vec62);

 // Read aligned vector block from vel_y at t+1, x+4, y, z.
 real_vec_t temp_vec183 = context.vel_y->readVecNorm(tv+(1/1), xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from vel_y at t+1, x+1, y, z.
 real_vec_t temp_vec184;
 // temp_vec184[0] = temp_vec62[1];  // for t+1, x+1, y, z;
 // temp_vec184[1] = temp_vec62[2];  // for t+1, x+2, y, z;
 // temp_vec184[2] = temp_vec62[3];  // for t+1, x+3, y, z;
 // temp_vec184[3] = temp_vec183[0];  // for t+1, x+4, y, z;
 // temp_vec184[4] = temp_vec62[5];  // for t+1, x+1, y+1, z;
 // temp_vec184[5] = temp_vec62[6];  // for t+1, x+2, y+1, z;
 // temp_vec184[6] = temp_vec62[7];  // for t+1, x+3, y+1, z;
 // temp_vec184[7] = temp_vec183[4];  // for t+1, x+4, y+1, z;
 // temp_vec184[8] = temp_vec62[9];  // for t+1, x+1, y+2, z;
 // temp_vec184[9] = temp_vec62[10];  // for t+1, x+2, y+2, z;
 // temp_vec184[10] = temp_vec62[11];  // for t+1, x+3, y+2, z;
 // temp_vec184[11] = temp_vec183[8];  // for t+1, x+4, y+2, z;
 // temp_vec184[12] = temp_vec62[13];  // for t+1, x+1, y+3, z;
 // temp_vec184[13] = temp_vec62[14];  // for t+1, x+2, y+3, z;
 // temp_vec184[14] = temp_vec62[15];  // for t+1, x+3, y+3, z;
 // temp_vec184[15] = temp_vec183[12];  // for t+1, x+4, y+3, z;
 real_vec_permute2(temp_vec184, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec62, temp_vec183);

 // Construct unaligned vector block from vel_y at t+1, x-2, y, z.
 real_vec_t temp_vec185;
 // temp_vec185[0] = temp_vec181[2];  // for t+1, x-2, y, z;
 // temp_vec185[1] = temp_vec181[3];  // for t+1, x-1, y, z;
 // temp_vec185[2] = temp_vec62[0];  // for t+1, x, y, z;
 // temp_vec185[3] = temp_vec62[1];  // for t+1, x+1, y, z;
 // temp_vec185[4] = temp_vec181[6];  // for t+1, x-2, y+1, z;
 // temp_vec185[5] = temp_vec181[7];  // for t+1, x-1, y+1, z;
 // temp_vec185[6] = temp_vec62[4];  // for t+1, x, y+1, z;
 // temp_vec185[7] = temp_vec62[5];  // for t+1, x+1, y+1, z;
 // temp_vec185[8] = temp_vec181[10];  // for t+1, x-2, y+2, z;
 // temp_vec185[9] = temp_vec181[11];  // for t+1, x-1, y+2, z;
 // temp_vec185[10] = temp_vec62[8];  // for t+1, x, y+2, z;
 // temp_vec185[11] = temp_vec62[9];  // for t+1, x+1, y+2, z;
 // temp_vec185[12] = temp_vec181[14];  // for t+1, x-2, y+3, z;
 // temp_vec185[13] = temp_vec181[15];  // for t+1, x-1, y+3, z;
 // temp_vec185[14] = temp_vec62[12];  // for t+1, x, y+3, z;
 // temp_vec185[15] = temp_vec62[13];  // for t+1, x+1, y+3, z;
 real_vec_permute2(temp_vec185, ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13, temp_vec181, temp_vec62);

 // Read aligned vector block from stress_mem_xy at t, x, y, z.
 real_vec_t temp_vec186 = context.stress_mem_xy->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from anelastic_xy at x, y, z.
 real_vec_t temp_vec187 = context.anelastic_xy->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec188 = mu(x, y, z) + mu(x, y, z-1).
 real_vec_t temp_vec188 = temp_vec27 + temp_vec33;

 // temp_vec189 = (2 / (mu(x, y, z) + mu(x, y, z-1))).
 real_vec_t temp_vec189 = temp_vec25 / temp_vec188;

 // temp_vec190 = (2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t().
 real_vec_t temp_vec190 = temp_vec189 * temp_vec2;

 // temp_vec191 = (((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h).
 real_vec_t temp_vec191 = temp_vec190 / temp_vec3;

 // temp_vec192 = (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z)).
 real_vec_t temp_vec192 = temp_vec177 - temp_vec41;

 // temp_vec193 = 1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z)).
 real_vec_t temp_vec193 = temp_vec40 * temp_vec192;

 // temp_vec194 = -0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z)).
 real_vec_t temp_vec194 = temp_vec44 * (temp_vec178 - temp_vec180);

 // temp_vec195 = (1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))).
 real_vec_t temp_vec195 = temp_vec193 + temp_vec194;

 // temp_vec196 = (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z)).
 real_vec_t temp_vec196 = temp_vec74 - temp_vec182;

 // temp_vec197 = 1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z)).
 real_vec_t temp_vec197 = temp_vec40 * temp_vec196;

 // temp_vec198 = (1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))).
 real_vec_t temp_vec198 = temp_vec195 + temp_vec197;

 // temp_vec199 = -0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)).
 real_vec_t temp_vec199 = temp_vec44 * (temp_vec184 - temp_vec185);

 // temp_vec200 = (1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))).
 real_vec_t temp_vec200 = temp_vec198 + temp_vec199;

 // temp_vec201 = (((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))).
 real_vec_t temp_vec201 = temp_vec191 * temp_vec200;

 // temp_vec202 = stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))).
 real_vec_t temp_vec202 = temp_vec175 + temp_vec201;

 // temp_vec203 = stress_mem_xy(t, x, y, z).
 real_vec_t temp_vec203 = temp_vec186;

 // temp_vec204 = tau2(x, y, z) * stress_mem_xy(t, x, y, z).
 real_vec_t temp_vec204 = temp_vec92 * temp_vec203;

 // temp_vec205 = (0.5 / h).
 real_vec_t temp_vec205 = temp_vec109 / temp_vec3;

 // temp_vec206 = (0.5 / h) * (1 - tau2(x, y, z)).
 real_vec_t temp_vec206 = temp_vec205 * temp_vec97;

 // temp_vec207 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_vec_t temp_vec207 = temp_vec206 * temp_vec100;

 // temp_vec208 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))).
 real_vec_t temp_vec208 = temp_vec207 * temp_vec189;

 // temp_vec209 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z).
 real_vec_t temp_vec209 = temp_vec208 * temp_vec187;

 // temp_vec210 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))).
 real_vec_t temp_vec210 = temp_vec209 * temp_vec200;

 // temp_vec211 = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))).
 real_vec_t temp_vec211 = temp_vec204 - temp_vec210;

 // temp_vec212 = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z).
 real_vec_t temp_vec212 = temp_vec211 + temp_vec203;

 // temp_vec213 = delta_t() * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)).
 real_vec_t temp_vec213 = temp_vec2 * temp_vec212;

 // temp_vec214 = stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z))).
 real_vec_t temp_vec214 = temp_vec202 + temp_vec213;

 // temp_vec215 = (stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)))) * sponge(x, y, z).
 real_vec_t temp_vec215 = temp_vec214 * temp_vec122;

 // temp_vec216 = ((stress_xy(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y, z-1))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))) + stress_mem_xy(t, x, y, z)))) * sponge(x, y, z)).
 real_vec_t temp_vec216 = temp_vec215;

 // Save result to stress_xy(t+1, x, y, z):
 
 // Write aligned vector block to stress_xy at t+1, x, y, z.
context.stress_xy->writeVecNorm(temp_vec216, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from stress_xz at t, x, y, z.
 real_vec_t temp_vec217 = context.stress_xz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from vel_x at t+1, x, y, z+1.
 real_vec_t temp_vec218 = context.vel_x->readVecNorm(tv+(1/1), xv, yv, zv+(1/1), __LINE__);

 // Read aligned vector block from vel_x at t+1, x, y, z+2.
 real_vec_t temp_vec219 = context.vel_x->readVecNorm(tv+(1/1), xv, yv, zv+(2/1), __LINE__);

 // Read aligned vector block from vel_x at t+1, x, y, z-1.
 real_vec_t temp_vec220 = context.vel_x->readVecNorm(tv+(1/1), xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from vel_z at t+1, x-4, y, z.
 real_vec_t temp_vec221 = context.vel_z->readVecNorm(tv+(1/1), xv-(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from vel_z at t+1, x-1, y, z.
 real_vec_t temp_vec222;
 // temp_vec222[0] = temp_vec221[3];  // for t+1, x-1, y, z;
 // temp_vec222[1] = temp_vec68[0];  // for t+1, x, y, z;
 // temp_vec222[2] = temp_vec68[1];  // for t+1, x+1, y, z;
 // temp_vec222[3] = temp_vec68[2];  // for t+1, x+2, y, z;
 // temp_vec222[4] = temp_vec221[7];  // for t+1, x-1, y+1, z;
 // temp_vec222[5] = temp_vec68[4];  // for t+1, x, y+1, z;
 // temp_vec222[6] = temp_vec68[5];  // for t+1, x+1, y+1, z;
 // temp_vec222[7] = temp_vec68[6];  // for t+1, x+2, y+1, z;
 // temp_vec222[8] = temp_vec221[11];  // for t+1, x-1, y+2, z;
 // temp_vec222[9] = temp_vec68[8];  // for t+1, x, y+2, z;
 // temp_vec222[10] = temp_vec68[9];  // for t+1, x+1, y+2, z;
 // temp_vec222[11] = temp_vec68[10];  // for t+1, x+2, y+2, z;
 // temp_vec222[12] = temp_vec221[15];  // for t+1, x-1, y+3, z;
 // temp_vec222[13] = temp_vec68[12];  // for t+1, x, y+3, z;
 // temp_vec222[14] = temp_vec68[13];  // for t+1, x+1, y+3, z;
 // temp_vec222[15] = temp_vec68[14];  // for t+1, x+2, y+3, z;
 real_vec_permute2(temp_vec222, ctrl_A3_B0_B1_B2_A7_B4_B5_B6_A11_B8_B9_B10_A15_B12_B13_B14, temp_vec221, temp_vec68);

 // Read aligned vector block from vel_z at t+1, x+4, y, z.
 real_vec_t temp_vec223 = context.vel_z->readVecNorm(tv+(1/1), xv+(4/4), yv, zv, __LINE__);

 // Construct unaligned vector block from vel_z at t+1, x+1, y, z.
 real_vec_t temp_vec224;
 // temp_vec224[0] = temp_vec68[1];  // for t+1, x+1, y, z;
 // temp_vec224[1] = temp_vec68[2];  // for t+1, x+2, y, z;
 // temp_vec224[2] = temp_vec68[3];  // for t+1, x+3, y, z;
 // temp_vec224[3] = temp_vec223[0];  // for t+1, x+4, y, z;
 // temp_vec224[4] = temp_vec68[5];  // for t+1, x+1, y+1, z;
 // temp_vec224[5] = temp_vec68[6];  // for t+1, x+2, y+1, z;
 // temp_vec224[6] = temp_vec68[7];  // for t+1, x+3, y+1, z;
 // temp_vec224[7] = temp_vec223[4];  // for t+1, x+4, y+1, z;
 // temp_vec224[8] = temp_vec68[9];  // for t+1, x+1, y+2, z;
 // temp_vec224[9] = temp_vec68[10];  // for t+1, x+2, y+2, z;
 // temp_vec224[10] = temp_vec68[11];  // for t+1, x+3, y+2, z;
 // temp_vec224[11] = temp_vec223[8];  // for t+1, x+4, y+2, z;
 // temp_vec224[12] = temp_vec68[13];  // for t+1, x+1, y+3, z;
 // temp_vec224[13] = temp_vec68[14];  // for t+1, x+2, y+3, z;
 // temp_vec224[14] = temp_vec68[15];  // for t+1, x+3, y+3, z;
 // temp_vec224[15] = temp_vec223[12];  // for t+1, x+4, y+3, z;
 real_vec_permute2(temp_vec224, ctrl_A1_A2_A3_B0_A5_A6_A7_B4_A9_A10_A11_B8_A13_A14_A15_B12, temp_vec68, temp_vec223);

 // Construct unaligned vector block from vel_z at t+1, x-2, y, z.
 real_vec_t temp_vec225;
 // temp_vec225[0] = temp_vec221[2];  // for t+1, x-2, y, z;
 // temp_vec225[1] = temp_vec221[3];  // for t+1, x-1, y, z;
 // temp_vec225[2] = temp_vec68[0];  // for t+1, x, y, z;
 // temp_vec225[3] = temp_vec68[1];  // for t+1, x+1, y, z;
 // temp_vec225[4] = temp_vec221[6];  // for t+1, x-2, y+1, z;
 // temp_vec225[5] = temp_vec221[7];  // for t+1, x-1, y+1, z;
 // temp_vec225[6] = temp_vec68[4];  // for t+1, x, y+1, z;
 // temp_vec225[7] = temp_vec68[5];  // for t+1, x+1, y+1, z;
 // temp_vec225[8] = temp_vec221[10];  // for t+1, x-2, y+2, z;
 // temp_vec225[9] = temp_vec221[11];  // for t+1, x-1, y+2, z;
 // temp_vec225[10] = temp_vec68[8];  // for t+1, x, y+2, z;
 // temp_vec225[11] = temp_vec68[9];  // for t+1, x+1, y+2, z;
 // temp_vec225[12] = temp_vec221[14];  // for t+1, x-2, y+3, z;
 // temp_vec225[13] = temp_vec221[15];  // for t+1, x-1, y+3, z;
 // temp_vec225[14] = temp_vec68[12];  // for t+1, x, y+3, z;
 // temp_vec225[15] = temp_vec68[13];  // for t+1, x+1, y+3, z;
 real_vec_permute2(temp_vec225, ctrl_A2_A3_B0_B1_A6_A7_B4_B5_A10_A11_B8_B9_A14_A15_B12_B13, temp_vec221, temp_vec68);

 // Read aligned vector block from stress_mem_xz at t, x, y, z.
 real_vec_t temp_vec226 = context.stress_mem_xz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from anelastic_xz at x, y, z.
 real_vec_t temp_vec227 = context.anelastic_xz->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec228 = mu(x, y, z) + mu(x, y-1, z).
 real_vec_t temp_vec228 = temp_vec27 + temp_vec30;

 // temp_vec229 = (2 / (mu(x, y, z) + mu(x, y-1, z))).
 real_vec_t temp_vec229 = temp_vec25 / temp_vec228;

 // temp_vec230 = (2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t().
 real_vec_t temp_vec230 = temp_vec229 * temp_vec2;

 // temp_vec231 = (((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h).
 real_vec_t temp_vec231 = temp_vec230 / temp_vec3;

 // temp_vec232 = (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z)).
 real_vec_t temp_vec232 = temp_vec218 - temp_vec41;

 // temp_vec233 = 1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z)).
 real_vec_t temp_vec233 = temp_vec40 * temp_vec232;

 // temp_vec234 = -0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1)).
 real_vec_t temp_vec234 = temp_vec44 * (temp_vec219 - temp_vec220);

 // temp_vec235 = (1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))).
 real_vec_t temp_vec235 = temp_vec233 + temp_vec234;

 // temp_vec236 = (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z)).
 real_vec_t temp_vec236 = temp_vec80 - temp_vec222;

 // temp_vec237 = 1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z)).
 real_vec_t temp_vec237 = temp_vec40 * temp_vec236;

 // temp_vec238 = (1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))).
 real_vec_t temp_vec238 = temp_vec235 + temp_vec237;

 // temp_vec239 = -0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)).
 real_vec_t temp_vec239 = temp_vec44 * (temp_vec224 - temp_vec225);

 // temp_vec240 = (1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))).
 real_vec_t temp_vec240 = temp_vec238 + temp_vec239;

 // temp_vec241 = (((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))).
 real_vec_t temp_vec241 = temp_vec231 * temp_vec240;

 // temp_vec242 = stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))).
 real_vec_t temp_vec242 = temp_vec217 + temp_vec241;

 // temp_vec243 = stress_mem_xz(t, x, y, z).
 real_vec_t temp_vec243 = temp_vec226;

 // temp_vec244 = tau2(x, y, z) * stress_mem_xz(t, x, y, z).
 real_vec_t temp_vec244 = temp_vec92 * temp_vec243;

 // temp_vec245 = (0.5 / h) * (1 - tau2(x, y, z)).
 real_vec_t temp_vec245 = temp_vec205 * temp_vec97;

 // temp_vec246 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_vec_t temp_vec246 = temp_vec245 * temp_vec100;

 // temp_vec247 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))).
 real_vec_t temp_vec247 = temp_vec246 * temp_vec229;

 // temp_vec248 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z).
 real_vec_t temp_vec248 = temp_vec247 * temp_vec227;

 // temp_vec249 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))).
 real_vec_t temp_vec249 = temp_vec248 * temp_vec240;

 // temp_vec250 = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))).
 real_vec_t temp_vec250 = temp_vec244 - temp_vec249;

 // temp_vec251 = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z).
 real_vec_t temp_vec251 = temp_vec250 + temp_vec243;

 // temp_vec252 = delta_t() * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)).
 real_vec_t temp_vec252 = temp_vec2 * temp_vec251;

 // temp_vec253 = stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z))).
 real_vec_t temp_vec253 = temp_vec242 + temp_vec252;

 // temp_vec254 = (stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)))) * sponge(x, y, z).
 real_vec_t temp_vec254 = temp_vec253 * temp_vec122;

 // temp_vec255 = ((stress_xz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x, y-1, z))) * delta_t) / h) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))) + stress_mem_xz(t, x, y, z)))) * sponge(x, y, z)).
 real_vec_t temp_vec255 = temp_vec254;

 // Save result to stress_xz(t+1, x, y, z):
 
 // Write aligned vector block to stress_xz at t+1, x, y, z.
context.stress_xz->writeVecNorm(temp_vec255, tv+(1/1), xv, yv, zv, __LINE__);
;

 // Read aligned vector block from stress_yz at t, x, y, z.
 real_vec_t temp_vec256 = context.stress_yz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from vel_y at t+1, x, y, z+1.
 real_vec_t temp_vec257 = context.vel_y->readVecNorm(tv+(1/1), xv, yv, zv+(1/1), __LINE__);

 // Read aligned vector block from vel_y at t+1, x, y, z+2.
 real_vec_t temp_vec258 = context.vel_y->readVecNorm(tv+(1/1), xv, yv, zv+(2/1), __LINE__);

 // Read aligned vector block from vel_y at t+1, x, y, z-1.
 real_vec_t temp_vec259 = context.vel_y->readVecNorm(tv+(1/1), xv, yv, zv-(1/1), __LINE__);

 // Read aligned vector block from vel_z at t+1, x, y+4, z.
 real_vec_t temp_vec260 = context.vel_z->readVecNorm(tv+(1/1), xv, yv+(4/4), zv, __LINE__);

 // Construct unaligned vector block from vel_z at t+1, x, y+1, z.
 real_vec_t temp_vec261;
 // temp_vec261[0] = temp_vec68[4];  // for t+1, x, y+1, z;
 // temp_vec261[1] = temp_vec68[5];  // for t+1, x+1, y+1, z;
 // temp_vec261[2] = temp_vec68[6];  // for t+1, x+2, y+1, z;
 // temp_vec261[3] = temp_vec68[7];  // for t+1, x+3, y+1, z;
 // temp_vec261[4] = temp_vec68[8];  // for t+1, x, y+2, z;
 // temp_vec261[5] = temp_vec68[9];  // for t+1, x+1, y+2, z;
 // temp_vec261[6] = temp_vec68[10];  // for t+1, x+2, y+2, z;
 // temp_vec261[7] = temp_vec68[11];  // for t+1, x+3, y+2, z;
 // temp_vec261[8] = temp_vec68[12];  // for t+1, x, y+3, z;
 // temp_vec261[9] = temp_vec68[13];  // for t+1, x+1, y+3, z;
 // temp_vec261[10] = temp_vec68[14];  // for t+1, x+2, y+3, z;
 // temp_vec261[11] = temp_vec68[15];  // for t+1, x+3, y+3, z;
 // temp_vec261[12] = temp_vec260[0];  // for t+1, x, y+4, z;
 // temp_vec261[13] = temp_vec260[1];  // for t+1, x+1, y+4, z;
 // temp_vec261[14] = temp_vec260[2];  // for t+1, x+2, y+4, z;
 // temp_vec261[15] = temp_vec260[3];  // for t+1, x+3, y+4, z;
 // Get 4 element(s) from temp_vec260 and 12 from temp_vec68.
 real_vec_align<4>(temp_vec261, temp_vec260, temp_vec68);

 // Construct unaligned vector block from vel_z at t+1, x, y+2, z.
 real_vec_t temp_vec262;
 // temp_vec262[0] = temp_vec68[8];  // for t+1, x, y+2, z;
 // temp_vec262[1] = temp_vec68[9];  // for t+1, x+1, y+2, z;
 // temp_vec262[2] = temp_vec68[10];  // for t+1, x+2, y+2, z;
 // temp_vec262[3] = temp_vec68[11];  // for t+1, x+3, y+2, z;
 // temp_vec262[4] = temp_vec68[12];  // for t+1, x, y+3, z;
 // temp_vec262[5] = temp_vec68[13];  // for t+1, x+1, y+3, z;
 // temp_vec262[6] = temp_vec68[14];  // for t+1, x+2, y+3, z;
 // temp_vec262[7] = temp_vec68[15];  // for t+1, x+3, y+3, z;
 // temp_vec262[8] = temp_vec260[0];  // for t+1, x, y+4, z;
 // temp_vec262[9] = temp_vec260[1];  // for t+1, x+1, y+4, z;
 // temp_vec262[10] = temp_vec260[2];  // for t+1, x+2, y+4, z;
 // temp_vec262[11] = temp_vec260[3];  // for t+1, x+3, y+4, z;
 // temp_vec262[12] = temp_vec260[4];  // for t+1, x, y+5, z;
 // temp_vec262[13] = temp_vec260[5];  // for t+1, x+1, y+5, z;
 // temp_vec262[14] = temp_vec260[6];  // for t+1, x+2, y+5, z;
 // temp_vec262[15] = temp_vec260[7];  // for t+1, x+3, y+5, z;
 // Get 8 element(s) from temp_vec260 and 8 from temp_vec68.
 real_vec_align<8>(temp_vec262, temp_vec260, temp_vec68);

 // Read aligned vector block from vel_z at t+1, x, y-4, z.
 real_vec_t temp_vec263 = context.vel_z->readVecNorm(tv+(1/1), xv, yv-(4/4), zv, __LINE__);

 // Construct unaligned vector block from vel_z at t+1, x, y-1, z.
 real_vec_t temp_vec264;
 // temp_vec264[0] = temp_vec263[12];  // for t+1, x, y-1, z;
 // temp_vec264[1] = temp_vec263[13];  // for t+1, x+1, y-1, z;
 // temp_vec264[2] = temp_vec263[14];  // for t+1, x+2, y-1, z;
 // temp_vec264[3] = temp_vec263[15];  // for t+1, x+3, y-1, z;
 // temp_vec264[4] = temp_vec68[0];  // for t+1, x, y, z;
 // temp_vec264[5] = temp_vec68[1];  // for t+1, x+1, y, z;
 // temp_vec264[6] = temp_vec68[2];  // for t+1, x+2, y, z;
 // temp_vec264[7] = temp_vec68[3];  // for t+1, x+3, y, z;
 // temp_vec264[8] = temp_vec68[4];  // for t+1, x, y+1, z;
 // temp_vec264[9] = temp_vec68[5];  // for t+1, x+1, y+1, z;
 // temp_vec264[10] = temp_vec68[6];  // for t+1, x+2, y+1, z;
 // temp_vec264[11] = temp_vec68[7];  // for t+1, x+3, y+1, z;
 // temp_vec264[12] = temp_vec68[8];  // for t+1, x, y+2, z;
 // temp_vec264[13] = temp_vec68[9];  // for t+1, x+1, y+2, z;
 // temp_vec264[14] = temp_vec68[10];  // for t+1, x+2, y+2, z;
 // temp_vec264[15] = temp_vec68[11];  // for t+1, x+3, y+2, z;
 // Get 12 element(s) from temp_vec68 and 4 from temp_vec263.
 real_vec_align<12>(temp_vec264, temp_vec68, temp_vec263);

 // Read aligned vector block from stress_mem_yz at t, x, y, z.
 real_vec_t temp_vec265 = context.stress_mem_yz->readVecNorm(tv, xv, yv, zv, __LINE__);

 // Read aligned vector block from anelastic_yz at x, y, z.
 real_vec_t temp_vec266 = context.anelastic_yz->readVecNorm(xv, yv, zv, __LINE__);

 // temp_vec267 = mu(x, y, z) + mu(x+1, y, z).
 real_vec_t temp_vec267 = temp_vec27 + temp_vec28;

 // temp_vec268 = (2 / (mu(x, y, z) + mu(x+1, y, z))).
 real_vec_t temp_vec268 = temp_vec25 / temp_vec267;

 // temp_vec269 = (2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t().
 real_vec_t temp_vec269 = temp_vec268 * temp_vec2;

 // temp_vec270 = (((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h).
 real_vec_t temp_vec270 = temp_vec269 / temp_vec3;

 // temp_vec271 = (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z)).
 real_vec_t temp_vec271 = temp_vec257 - temp_vec74;

 // temp_vec272 = 1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z)).
 real_vec_t temp_vec272 = temp_vec40 * temp_vec271;

 // temp_vec273 = -0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1)).
 real_vec_t temp_vec273 = temp_vec44 * (temp_vec258 - temp_vec259);

 // temp_vec274 = (1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))).
 real_vec_t temp_vec274 = temp_vec272 + temp_vec273;

 // temp_vec275 = (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z)).
 real_vec_t temp_vec275 = temp_vec261 - temp_vec80;

 // temp_vec276 = 1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z)).
 real_vec_t temp_vec276 = temp_vec40 * temp_vec275;

 // temp_vec277 = (1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))).
 real_vec_t temp_vec277 = temp_vec274 + temp_vec276;

 // temp_vec278 = -0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)).
 real_vec_t temp_vec278 = temp_vec44 * (temp_vec262 - temp_vec264);

 // temp_vec279 = (1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))).
 real_vec_t temp_vec279 = temp_vec277 + temp_vec278;

 // temp_vec280 = (((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))).
 real_vec_t temp_vec280 = temp_vec270 * temp_vec279;

 // temp_vec281 = stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))).
 real_vec_t temp_vec281 = temp_vec256 + temp_vec280;

 // temp_vec282 = stress_mem_yz(t, x, y, z).
 real_vec_t temp_vec282 = temp_vec265;

 // temp_vec283 = tau2(x, y, z) * stress_mem_yz(t, x, y, z).
 real_vec_t temp_vec283 = temp_vec92 * temp_vec282;

 // temp_vec284 = (0.5 / h) * (1 - tau2(x, y, z)).
 real_vec_t temp_vec284 = temp_vec205 * temp_vec97;

 // temp_vec285 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z).
 real_vec_t temp_vec285 = temp_vec284 * temp_vec100;

 // temp_vec286 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))).
 real_vec_t temp_vec286 = temp_vec285 * temp_vec268;

 // temp_vec287 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z).
 real_vec_t temp_vec287 = temp_vec286 * temp_vec266;

 // temp_vec288 = (0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))).
 real_vec_t temp_vec288 = temp_vec287 * temp_vec279;

 // temp_vec289 = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))).
 real_vec_t temp_vec289 = temp_vec283 - temp_vec288;

 // temp_vec290 = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z).
 real_vec_t temp_vec290 = temp_vec289 + temp_vec282;

 // temp_vec291 = delta_t() * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)).
 real_vec_t temp_vec291 = temp_vec2 * temp_vec290;

 // temp_vec292 = stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z))).
 real_vec_t temp_vec292 = temp_vec281 + temp_vec291;

 // temp_vec293 = (stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)))) * sponge(x, y, z).
 real_vec_t temp_vec293 = temp_vec292 * temp_vec122;

 // temp_vec294 = ((stress_yz(t, x, y, z) + ((((2 / (mu(x, y, z) + mu(x+1, y, z))) * delta_t) / h) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z))))) + (delta_t * (((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))) + stress_mem_yz(t, x, y, z)))) * sponge(x, y, z)).
 real_vec_t temp_vec294 = temp_vec293;

 // Save result to stress_yz(t+1, x, y, z):
 
 // Write aligned vector block to stress_yz at t+1, x, y, z.
context.stress_yz->writeVecNorm(temp_vec294, tv+(1/1), xv, yv, zv, __LINE__);
;

 // temp_vec295 = (tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec295 = temp_vec94 + temp_vec116;

 // temp_vec296 = ((tau2(x, y, z) * stress_mem_xx(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 real_vec_t temp_vec296 = temp_vec295;

 // Save result to stress_mem_xx(t+1, x, y, z):
 
 // Write aligned vector block to stress_mem_xx at t+1, x, y, z.
context.stress_mem_xx->writeVecNorm(temp_vec296, tv+(1/1), xv, yv, zv, __LINE__);
;

 // temp_vec297 = (tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec297 = temp_vec134 + temp_vec143;

 // temp_vec298 = ((tau2(x, y, z) * stress_mem_yy(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 real_vec_t temp_vec298 = temp_vec297;

 // Save result to stress_mem_yy(t+1, x, y, z):
 
 // Write aligned vector block to stress_mem_yy at t+1, x, y, z.
context.stress_mem_yy->writeVecNorm(temp_vec298, tv+(1/1), xv, yv, zv, __LINE__);
;

 // temp_vec299 = (tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2))))))).
 real_vec_t temp_vec299 = temp_vec159 + temp_vec168;

 // temp_vec300 = ((tau2(x, y, z) * stress_mem_zz(t, x, y, z)) + ((1 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) * anelastic_as_diag(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))))) - (((8 / (mu(x, y, z) + mu(x+1, y, z) + mu(x, y-1, z) + mu(x+1, y-1, z) + mu(x, y, z-1) + mu(x+1, y, z-1) + mu(x, y-1, z-1) + mu(x+1, y-1, z-1))) + (0.5 * (8 / (lambda(x, y, z) + lambda(x+1, y, z) + lambda(x, y-1, z) + lambda(x+1, y-1, z) + lambda(x, y, z-1) + lambda(x+1, y, z-1) + lambda(x, y-1, z-1) + lambda(x+1, y-1, z-1))))) * anelastic_ap(x, y, z) * ((1.125 * (vel_x(t+1, x+1, y, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x+2, y, z) - vel_x(t+1, x-1, y, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x, y-1, z))) + (-0.0416667 * (vel_y(t+1, x, y+1, z) - vel_y(t+1, x, y-2, z))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x, y, z-1))) + (-0.0416667 * (vel_z(t+1, x, y, z+1) - vel_z(t+1, x, y, z-2)))))))).
 real_vec_t temp_vec300 = temp_vec299;

 // Save result to stress_mem_zz(t+1, x, y, z):
 
 // Write aligned vector block to stress_mem_zz at t+1, x, y, z.
context.stress_mem_zz->writeVecNorm(temp_vec300, tv+(1/1), xv, yv, zv, __LINE__);
;

 // temp_vec301 = ((tau2(x, y, z) * stress_mem_xy(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y, z-1))) * anelastic_xy(x, y, z) * ((1.125 * (vel_x(t+1, x, y+1, z) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y+2, z) - vel_x(t+1, x, y-1, z))) + (1.125 * (vel_y(t+1, x, y, z) - vel_y(t+1, x-1, y, z))) + (-0.0416667 * (vel_y(t+1, x+1, y, z) - vel_y(t+1, x-2, y, z)))))).
 real_vec_t temp_vec301 = temp_vec211;

 // Save result to stress_mem_xy(t+1, x, y, z):
 
 // Write aligned vector block to stress_mem_xy at t+1, x, y, z.
context.stress_mem_xy->writeVecNorm(temp_vec301, tv+(1/1), xv, yv, zv, __LINE__);
;

 // temp_vec302 = ((tau2(x, y, z) * stress_mem_xz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x, y-1, z))) * anelastic_xz(x, y, z) * ((1.125 * (vel_x(t+1, x, y, z+1) - vel_x(t+1, x, y, z))) + (-0.0416667 * (vel_x(t+1, x, y, z+2) - vel_x(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y, z) - vel_z(t+1, x-1, y, z))) + (-0.0416667 * (vel_z(t+1, x+1, y, z) - vel_z(t+1, x-2, y, z)))))).
 real_vec_t temp_vec302 = temp_vec250;

 // Save result to stress_mem_xz(t+1, x, y, z):
 
 // Write aligned vector block to stress_mem_xz at t+1, x, y, z.
context.stress_mem_xz->writeVecNorm(temp_vec302, tv+(1/1), xv, yv, zv, __LINE__);
;

 // temp_vec303 = ((tau2(x, y, z) * stress_mem_yz(t, x, y, z)) - ((0.5 / h) * (1 - tau2(x, y, z)) * weight(x, y, z) * (2 / (mu(x, y, z) + mu(x+1, y, z))) * anelastic_yz(x, y, z) * ((1.125 * (vel_y(t+1, x, y, z+1) - vel_y(t+1, x, y, z))) + (-0.0416667 * (vel_y(t+1, x, y, z+2) - vel_y(t+1, x, y, z-1))) + (1.125 * (vel_z(t+1, x, y+1, z) - vel_z(t+1, x, y, z))) + (-0.0416667 * (vel_z(t+1, x, y+2, z) - vel_z(t+1, x, y-1, z)))))).
 real_vec_t temp_vec303 = temp_vec289;

 // Save result to stress_mem_yz(t+1, x, y, z):
 
 // Write aligned vector block to stress_mem_yz at t+1, x, y, z.
context.stress_mem_yz->writeVecNorm(temp_vec303, tv+(1/1), xv, yv, zv, __LINE__);
;
} // vector calculation.

 // Prefetches cache line(s) for entire stencil  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned anelastic_ap at x, y, z.
 p = (const char*)context.anelastic_ap->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_as_diag at x, y, z.
 p = (const char*)context.anelastic_as_diag->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xy at x, y, z.
 p = (const char*)context.anelastic_xy->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xz at x, y, z.
 p = (const char*)context.anelastic_xz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_yz at x, y, z.
 p = (const char*)context.anelastic_yz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y-4, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y-4, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y-4, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y-4, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y-4, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y-4, z.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y-4, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y-4, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xx at t, x, y, z.
 p = (const char*)context.stress_mem_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xy at t, x, y, z.
 p = (const char*)context.stress_mem_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xz at t, x, y, z.
 p = (const char*)context.stress_mem_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yy at t, x, y, z.
 p = (const char*)context.stress_mem_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yz at t, x, y, z.
 p = (const char*)context.stress_mem_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_zz at t, x, y, z.
 p = (const char*)context.stress_mem_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned tau2 at x, y, z.
 p = (const char*)context.tau2->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x-4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y-4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z-1.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+1.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+2.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y+4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x+4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x-4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y-4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z-1.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+1.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+2.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y+4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x+4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x-4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y-4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z-2.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z-1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z+1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y+4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x+4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned weight at x, y, z.
 p = (const char*)context.weight->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}

 // Prefetches cache line(s) for leading edge of stencil in '+x' direction  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster_x(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned anelastic_ap at x, y, z.
 p = (const char*)context.anelastic_ap->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_as_diag at x, y, z.
 p = (const char*)context.anelastic_as_diag->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xy at x, y, z.
 p = (const char*)context.anelastic_xy->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xz at x, y, z.
 p = (const char*)context.anelastic_xz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_yz at x, y, z.
 p = (const char*)context.anelastic_yz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y-4, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y-4, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y-4, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv-(4/4), zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y-4, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xx at t, x, y, z.
 p = (const char*)context.stress_mem_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xy at t, x, y, z.
 p = (const char*)context.stress_mem_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xz at t, x, y, z.
 p = (const char*)context.stress_mem_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yy at t, x, y, z.
 p = (const char*)context.stress_mem_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yz at t, x, y, z.
 p = (const char*)context.stress_mem_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_zz at t, x, y, z.
 p = (const char*)context.stress_mem_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned tau2 at x, y, z.
 p = (const char*)context.tau2->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y-4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z-1.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+1.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+2.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y+4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x+4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y-4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z-1.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+1.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+2.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y+4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x+4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y-4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z-2.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z-1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z+1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y+4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x+4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned weight at x, y, z.
 p = (const char*)context.weight->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}

 // Prefetches cache line(s) for leading edge of stencil in '+y' direction  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster_y(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned anelastic_ap at x, y, z.
 p = (const char*)context.anelastic_ap->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_as_diag at x, y, z.
 p = (const char*)context.anelastic_as_diag->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xy at x, y, z.
 p = (const char*)context.anelastic_xy->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xz at x, y, z.
 p = (const char*)context.anelastic_xz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_yz at x, y, z.
 p = (const char*)context.anelastic_yz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z-1.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z-1.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xx at t, x, y, z.
 p = (const char*)context.stress_mem_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xy at t, x, y, z.
 p = (const char*)context.stress_mem_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xz at t, x, y, z.
 p = (const char*)context.stress_mem_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yy at t, x, y, z.
 p = (const char*)context.stress_mem_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yz at t, x, y, z.
 p = (const char*)context.stress_mem_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_zz at t, x, y, z.
 p = (const char*)context.stress_mem_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned tau2 at x, y, z.
 p = (const char*)context.tau2->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x-4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z-1.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+1.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+2.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y+4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x+4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x-4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z-1.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+1.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+2.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y+4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x+4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x-4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z-2.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv-(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z-1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv-(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z+1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y+4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x+4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned weight at x, y, z.
 p = (const char*)context.weight->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}

 // Prefetches cache line(s) for leading edge of stencil in '+z' direction  relative to indices t, x, y, z in a 'x=1 * y=1 * z=1' cluster of 'x=4 * y=4 * z=1' vector(s).
 // Indices must be normalized, i.e., already divided by VLEN_*.
 template<int level> void prefetch_cluster_z(StencilContext_awp& context, idx_t tv, idx_t xv, idx_t yv, idx_t zv) {
 const char* p = 0;

 // Aligned anelastic_ap at x, y, z.
 p = (const char*)context.anelastic_ap->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_as_diag at x, y, z.
 p = (const char*)context.anelastic_as_diag->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xy at x, y, z.
 p = (const char*)context.anelastic_xy->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_xz at x, y, z.
 p = (const char*)context.anelastic_xz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned anelastic_yz at x, y, z.
 p = (const char*)context.anelastic_yz->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y-4, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y-4, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned lambda at x+4, y, z.
 p = (const char*)context.lambda->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y-4, z.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y-4, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned mu at x+4, y, z.
 p = (const char*)context.mu->getVecPtrNorm(xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned sponge at x, y, z.
 p = (const char*)context.sponge->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xx at t, x, y, z.
 p = (const char*)context.stress_mem_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xy at t, x, y, z.
 p = (const char*)context.stress_mem_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_xz at t, x, y, z.
 p = (const char*)context.stress_mem_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yy at t, x, y, z.
 p = (const char*)context.stress_mem_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_yz at t, x, y, z.
 p = (const char*)context.stress_mem_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_mem_zz at t, x, y, z.
 p = (const char*)context.stress_mem_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xx at t, x, y, z.
 p = (const char*)context.stress_xx->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xy at t, x, y, z.
 p = (const char*)context.stress_xy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_xz at t, x, y, z.
 p = (const char*)context.stress_xz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yy at t, x, y, z.
 p = (const char*)context.stress_yy->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_yz at t, x, y, z.
 p = (const char*)context.stress_yz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned stress_zz at t, x, y, z.
 p = (const char*)context.stress_zz->getVecPtrNorm(tv, xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned tau2 at x, y, z.
 p = (const char*)context.tau2->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x-4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y-4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y, z+2.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x, y+4, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_x at t+1, x+4, y, z.
 p = (const char*)context.vel_x->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x-4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y-4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y, z+2.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv, zv+(2/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x, y+4, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_y at t+1, x+4, y, z.
 p = (const char*)context.vel_y->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x-4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv-(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y-4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv-(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y, z+1.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv, zv+(1/1), false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x, y+4, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv, yv+(4/4), zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned vel_z at t+1, x+4, y, z.
 p = (const char*)context.vel_z->getVecPtrNorm(tv+(1/1), xv+(4/4), yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);

 // Aligned weight at x, y, z.
 p = (const char*)context.weight->getVecPtrNorm(xv, yv, zv, false);
 MCP(p, level, __LINE__);
 _mm_prefetch(p, level);
}
};

 ////// Overall stencil-equations class //////
template <typename ContextClass>
struct StencilEquations_awp : public StencilEquations {

 // Stencils.
 StencilTemplate<Stencil_velocity,StencilContext_awp> stencil_velocity;
 StencilTemplate<Stencil_stress,StencilContext_awp> stencil_stress;

 StencilEquations_awp() {
name = "awp";
  stencils.push_back(&stencil_velocity);
  stencils.push_back(&stencil_stress);
 }
};
} // namespace yask.
