
#ifndef _EXTERN_VAR_H
#define _EXTERN_VAR_H

#include <complex.h>
#include <omp.h>
#include <stdio.h>
#include <cuda.h>
// #include <cuda_runtime.h>
// #include <device_launch_parameters.h>
//#include<conio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>


#define DOUBLEPRECISION
//#define FlOATPRECISION

//#define DOUBLECOMPLEX

#if defined(DOUBLEPRECISION)
	typedef double comp;
	typedef double real;
	typedef float real2;
#elif defined(FlOATPRECISION)
	typedef float comp;
	typedef float real;
#elif defined(DOUBLECOMPLEX)
	typedef double complex comp;
	typedef double complex real;
	typedef float real2;
#endif

extern FILE *Source;

extern void READ_DATA_FILE(void);
extern void SOURCE_SETUP(void);
extern void SOURCE_IN();

extern int StaticField;

extern int PBC_CTW;
extern real  Nonlocalend;
extern int Nonlocalend_int;

extern int Diverge_Gradient,Laplacian;
extern int material, WL_or_freq,Test_offset;
extern int NONLOCAL;
extern int Tend_inc,NUM_freq_inc,Freq_start,Spect_loc,centerx,centery;
extern int trigger;
extern comp *Hx_source, *Hy_source, *Hz_source;
extern comp *Ex_source, *Ey_source, *Ez_source;

extern int * mat_matrix,* mat_matrixX,* mat_matrixY,* mat_matrixZ, first_medium,first_medium_max,nano_sphere;
extern int * mat_matrixdev,* mat_matrixXdev,* mat_matrixYdev,* mat_matrixZdev;

extern real  nano_sphere_radius;

extern real  BandWidth;
extern real MAX_AMP;
extern int num_trials,min_trials,max_trials,TE_TM,trials;
extern int inf_disp_slab;
//Polarization type
extern int TMz,TEz;
extern int Snap_in;

//number of cells in simulation domain
extern int NCELLX, NCELLY, NCELLZ;
//incident plane position
extern int inc_plane,source_end;
extern real  f_0;
//number of TFSF cells in each dimension
extern int NtfsfX, NtfsfY,NtfsfZ;
extern int cpml_N_X, cpml_F_X,cpml_N_Y,cpml_F_Y,cpml_N_Z,cpml_F_Z;
extern int cpml_x_lim,cpml_y_lim,cpml_z_lim;

//number of CPML cells in each domain
extern int NcpmlX,NcpmlY,NcpmlZ;
//E and H fields
extern comp *ex,*ey,*ez,*hx,*hy,*hz;
extern comp  *exdev,*eydev,*ezdev,*hxdev,*hydev,*hzdev;
extern comp *ex_n_1dev,*ey_n_1dev,*ez_n_1dev;
extern comp *ex_ndev,*ey_ndev,*ez_ndev;

extern comp *Dx,*Dy,*Dz;
extern comp *ex_n,*ey_n,*ez_n;
extern comp *ex_n_1,*ey_n_1,*ez_n_1;
extern comp *hxPrev,*hyPrev,*hzPrev;
extern comp *hxPrevdev,*hyPrevdev,*hzPrevdev;
extern comp *NDx,*NDy,*NDz;
extern comp *NDxdev,*NDydev,*NDzdev;
extern comp *NDx_prev,*NDy_prev,*NDz_prev;
extern comp *NDx_prevdev,*NDy_prevdev,*NDz_prevdev;

extern comp *Hx_w_FT,*Hy_w_FT,*Hz_w_FT,*Ex_w_FT,*Ey_w_FT,*Ez_w_FT;
extern comp *Hx_t_FT,*Hy_t_FT,*Hz_t_FT,*Ex_t_FT,*Ey_t_FT,*Ez_t_FT;

//update coefficients
extern real  *Cexe,*Cexh,*Ceye,*Ceyh,*Ceze,*Cezh;
extern real  *Chxe,*Chxh,*Chye,*Chyh,*Chze,*Chzh;
extern comp *Chxedev,*Chxhdev,*Chyedev,*Chyhdev,*Chzedev,*Chzhdev;
extern real  *Cexedev,*Cexhdev,*Ceyedev,*Ceyhdev,*Cezedev,*Cezhdev;

extern real  *sigma_e, *sigma_m;
extern real  *eps,*mu;
//CPML Stuff
extern comp *psi_Ex_y_N,*psi_Ex_z_N,*psi_Ey_z_N,*psi_Ey_x_N,*psi_Ez_x_N,*psi_Ez_y_N;
extern comp *psi_Hx_y_N,*psi_Hx_z_N,*psi_Hy_z_N,*psi_Hy_x_N,*psi_Hz_x_N,*psi_Hz_y_N;
extern comp *psi_Ex_y_F,*psi_Ex_z_F,*psi_Ey_z_F,*psi_Ey_x_F,*psi_Ez_x_F,*psi_Ez_y_F;
extern comp *psi_Hx_y_F,*psi_Hx_z_F,*psi_Hy_z_F,*psi_Hy_x_F,*psi_Hz_x_F,*psi_Hz_y_F;

extern comp  *psi_Ex_y_Ndev,*psi_Ex_z_Ndev,*psi_Ey_z_Ndev,*psi_Ey_x_Ndev,*psi_Ez_x_Ndev,*psi_Ez_y_Ndev;
extern comp  *psi_Hx_y_Ndev,*psi_Hx_z_Ndev,*psi_Hy_z_Ndev,*psi_Hy_x_Ndev,*psi_Hz_x_Ndev,*psi_Hz_y_Ndev;
extern comp  *psi_Ex_y_Fdev,*psi_Ex_z_Fdev,*psi_Ey_z_Fdev,*psi_Ey_x_Fdev,*psi_Ez_x_Fdev,*psi_Ez_y_Fdev;
extern comp  *psi_Hx_y_Fdev,*psi_Hx_z_Fdev,*psi_Hy_z_Fdev,*psi_Hy_x_Fdev,*psi_Hz_x_Fdev,*psi_Hz_y_Fdev;
extern real  *kedx,*kedy,*kedz;
extern real  *khdx,*khdy,*khdz;
extern real  *be_x_N,*be_y_N,*be_z_N;
extern real  *bh_x_N,*bh_y_N,*bh_z_N;
extern real  *ce_x_N,*ce_y_N,*ce_z_N;
extern real  *ch_x_N,*ch_y_N,*ch_z_N;
extern real  *be_x_F,*be_y_F,*be_z_F;
extern real  *bh_x_F,*bh_y_F,*bh_z_F;
extern real  *ce_x_F,*ce_y_F,*ce_z_F;
extern real  *ch_x_F,*ch_y_F,*ch_z_F;

extern real  *kedxdev,*kedydev,*kedzdev;
extern real  *khdxdev,*khdydev,*khdzdev;
extern real  *be_x_Ndev,*be_y_Ndev,*be_z_Ndev;
extern real  *bh_x_Ndev,*bh_y_Ndev,*bh_z_Ndev;
extern real  *ce_x_Ndev,*ce_y_Ndev,*ce_z_Ndev;
extern real  *ch_x_Ndev,*ch_y_Ndev,*ch_z_Ndev;
extern real  *be_x_Fdev,*be_y_Fdev,*be_z_Fdev;
extern real  *bh_x_Fdev,*bh_y_Fdev,*bh_z_Fdev;
extern real  *ce_x_Fdev,*ce_y_Fdev,*ce_z_Fdev;
extern real  *ch_x_Fdev,*ch_y_Fdev,*ch_z_Fdev;


extern real  *sigma_e_x, *sigma_e_y,*sigma_e_z;
extern real  *sigma_h_x, *sigma_h_y,*sigma_h_z;
extern real  cpml_exp;
extern real  max_stretch_factor_x,max_stretch_factor_y,max_stretch_factor_z;
extern real  max_sigma_cpml_x,max_sigma_cpml_y,max_sigma_cpml_z;
extern real  max_alpha_x, max_alpha_y, max_alpha_z;
extern real  exp_alpha_x, exp_alpha_y, exp_alpha_z;

//temporal and spacial step sizes
extern real  dt,dx,dy,dz,d_1D;
//Time Stepping parameters
extern int t,Tend;
//Snapshot parameters
extern int x_skip,y_skip,z_skip,t_skip,snapshot_count;

//source fields and parameters;
extern real  inc_phi, inc_theta, polar_psi;
extern comp *e_inc, *h_inc;
extern comp *e_incdev, *h_incdev;

extern comp *hx_inc, *hy_inc, *hz_inc;
extern comp *ex_inc, *ey_inc, *ez_inc;
extern int inc_Length,m0,i_0,j_0,k_0;
extern real  delay, width, width_real ;
//TFSF buffers
extern real  e1,e2;
extern int inc_Length;

//constants
extern real  z0,ep0,mu0,c0,pi,me,e0;

//Periodic Boundary Conditions
extern real  k_x,k_y,k_z,k_rho,period_x,period_y,period_z;
extern int Periodic_XY,Periodic_YZ,Periodic_XZ;


//FUNCTIONS
extern real * MALLOC1D(real  *ARRAY, int i);
extern real * MALLOC1D_double (real  *ARRAY, int i);
extern comp* MALLOC1D_Complex(comp *ARRAY, int i);
extern double complex* MALLOC1D_Complex2(double complex *ARRAY, int i);
extern real2* MALLOC1D_Real2(real2 *ARRAY, int i);


extern real * MALLOC2D(real  *ARRAY, int i, int j);
extern real * MALLOC2D_double (real  *ARRAY, int i, int j);
extern comp* MALLOC2D_Complex(comp *ARRAY, int i, int j);


extern real * MALLOC3D(real  *ARRAY, int i, int j,int k);
extern real * MALLOC3D_double (real  *ARRAY, int i, int j,int k);
extern comp* MALLOC3D_Complex(comp *ARRAY, int i, int j,int k);
extern double complex* MALLOC3D_Complex2(double complex *ARRAY, int i, int j,int k);
extern real2* MALLOC3D_Real2(real2 *ARRAY, int i, int j,int k);


extern real * MALLOC4D(real  *ARRAY, int i, int j,int k, int size4);
extern real * MALLOC4D_double (real  *ARRAY, int i, int j,int k, int size4);
extern comp* MALLOC4D_Complex(comp *ARRAY, int i, int j,int k, int size4);

extern void FREE4D(real  *ARRAY,int i,int j,int size4);
extern void FREE4D_double (real  *ARRAY,int i,int j,int size4);
extern void FREE4D_Complex(comp *ARRAY,int i, int j,int size4);

extern void FREE3D(real  *ARRAY,int i,int j);
extern void FREE3D_double (real  *ARRAY,int i,int j);
extern void FREE3D_Complex(comp *ARRAY,int i, int j);
extern void FREE3D_Complex2(double complex *ARRAY,int i, int j);


extern void FREE2D(real  *ARRAY,int i);
extern void FREE2D_double (real  *ARRAY,int i);
extern void FREE2D_Complex(comp *ARRAY,int i);

extern void FREE1D(real  *ARRAY);
extern void FREE1D_double (real  *ARRAY);
extern void FREE1D_Complex(comp *ARRAY);
extern void FREE1D_Complex2(double complex *ARRAY);

extern comp Calc_DIV_GRADz(int i, int j, int k);
extern comp Calc_DIV_GRADx(int i, int j, int k);
extern comp Calc_DIV_GRADy(int i, int j, int k);
extern comp Calc_DIV_GRADz2(int i, int j, int k);
extern comp Calc_DIV_GRADx2(int i, int j, int k);
extern comp Calc_DIV_GRADy2(int i, int j, int k);

// extern void UPDATE_ex();
extern __global__ void UPDATE_ex(real  *ex,real  *ex_n,real  *ex_n_1,real  *hy,real  *hz,real  *Cexe,real  *Cexh,real  *kedy,real  *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real  *psi_Ex_z_N,
real  *psi_Ex_z_F,real  *psi_Ex_y_N,real  *psi_Ex_y_F,real  *Px_cp,real  *Px_cp_n,real  *Px_cp_n_1,real  *Px_d,real  *Px_d_n,real  *Px_d_n_1,real  *C_1_cp,real  *C_2_cp,real  *C_3_cp,real  *C_4_cp,real  *C_5_cp,real  *d_1_d,
real  *d_2_d,real  *d_3_d,real  *d_4_d,real  *d_5_d,real  *d_NL,real  C_E,real  z0,int N_CP_poles,int N_drude_poles,real  *ce_z_N,real  *ce_z_F,real  *be_z_N,real  *be_z_F,real  *ce_y_N,real  *ce_y_F,real  *be_y_N,real  *be_y_F,
real  dx,real  dy,real  dz,real  dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real  C_E_1,real  C_E_2,int Periodic_XY);
//extern void UPDATE_ey();
extern __global__ void UPDATE_ey(real  *ex,real  *ex_n,real  *ex_n_1,real  *hy,real  *hz,real  *Cexe,real  *Cexh,real  *kedy,real  *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real  *psi_Ex_z_N,
real  *psi_Ex_z_F,real  *psi_Ex_y_N,real  *psi_Ex_y_F,real  *Px_cp,real  *Px_cp_n,real  *Px_cp_n_1,real  *Px_d,real  *Px_d_n,real  *Px_d_n_1,real  *C_1_cp,real  *C_2_cp,real  *C_3_cp,real  *C_4_cp,real  *C_5_cp,real  *d_1_d,
real  *d_2_d,real  *d_3_d,real  *d_4_d,real  *d_5_d,real  *d_NL,real  C_E,real  z0,int N_CP_poles,int N_drude_poles,real  *ce_z_N,real  *ce_z_F,real  *be_z_N,real  *be_z_F,real  *ce_y_N,real  *ce_y_F,real  *be_y_N,real  *be_y_F,
real  dx,real  dy,real  dz,real  dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real  C_E_1,real  C_E_2,int Periodic_XY);


extern __global__ void UPDATE_ez(real  *ex,real  *ex_n,real  *ex_n_1,real  *hy,real  *hz,real  *Cexe,real  *Cexh,real  *kedy,real  *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real  *psi_Ex_z_N,
real  *psi_Ex_z_F,real  *psi_Ex_y_N,real  *psi_Ex_y_F,real  *Px_cp,real  *Px_cp_n,real  *Px_cp_n_1,real  *Px_d,real  *Px_d_n,real  *Px_d_n_1,real  *C_1_cp,real  *C_2_cp,real  *C_3_cp,real  *C_4_cp,real  *C_5_cp,real  *d_1_d,
real  *d_2_d,real  *d_3_d,real  *d_4_d,real  *d_5_d,real  *d_NL,real  C_E,real  z0,int N_CP_poles,int N_drude_poles,real  *ce_z_N,real  *ce_z_F,real  *be_z_N,real  *be_z_F,real  *ce_y_N,real  *ce_y_F,real  *be_y_N,real  *be_y_F,
real  dx,real  dy,real  dz,real  dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real  C_E_1,real  C_E_2,int Periodic_XY);

// extern void UPDATE_ez();
extern void UPDATE_dx(void);
extern void UPDATE_dy(void);
extern void UPDATE_dz(void);
extern void UPDATE_PML_E_X_F(void);
extern void UPDATE_PML_E_Y_F(void);
extern void UPDATE_PML_E_Z_F(void);
extern void UPDATE_PML_E_X_N(void);
extern void UPDATE_PML_E_Y_N(void);
extern void UPDATE_PML_E_Z_N(void);
extern void UPDATE_Px(void);
extern void UPDATE_Py(void);
extern void UPDATE_Pz(void);
extern void PBC_E_X(void);
extern void PBC_E_Y(void);
extern void PBC_E_Z(void);
extern void CP_D_ex(int i, int j, int k, comp Curl_H,comp Div_Grad);
extern void CP_D_ey(int i, int j, int k, comp Curl_H,comp Div_Grad);
extern void CP_D_ez(int i, int j, int k, comp Curl_H,comp Div_Grad);

// extern void UPDATE_hx();

  extern __global__ void UPDATE_hx(real  *hx,real  *hxPrev,real  *ex,real  *ey,real  *Chxh,real  *Chxe,real  *psi_Hx_z_N,real  *psi_Hx_z_F,real  *psi_Hx_y_N,real  *psi_Hx_y_F,real  *khdy,real
    *khdz,real  *bh_z_N,real  *bh_z_F,real  *ch_z_N,real  *ch_z_F,real  *bh_y_N,real  *bh_y_F,real  *ch_y_N,real  *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real  dx,real  dy,real  dz,real  dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics);
// extern void UPDATE_hy();
extern __global__ void UPDATE_hy(real  *hx,real  *hxPrev,real  *ex,real  *ey,real  *Chxh,real  *Chxe,real  *psi_Hx_z_N,real  *psi_Hx_z_F,real  *psi_Hx_y_N,real  *psi_Hx_y_F,real  *khdy,real
	*khdz,real  *bh_z_N,real  *bh_z_F,real  *ch_z_N,real  *ch_z_F,real  *bh_y_N,real  *bh_y_F,real  *ch_y_N,real  *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real  dx,real  dy,real  dz,real  dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics);

// extern void UPDATE_hz();
extern __global__ void UPDATE_hz(real  *hx,real  *hxPrev,real  *ex,real  *ey,real  *Chxh,real  *Chxe,real  *psi_Hx_z_N,real  *psi_Hx_z_F,real  *psi_Hx_y_N,real  *psi_Hx_y_F,real  *khdy,real
	*khdz,real  *bh_z_N,real  *bh_z_F,real  *ch_z_N,real  *ch_z_F,real  *bh_y_N,real  *bh_y_F,real  *ch_y_N,real  *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real  dx,real  dy,real  dz,real  dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics);


extern void UPDATE_PML_H_X_F(void);
extern void UPDATE_PML_H_Y_F(void);
extern void UPDATE_PML_H_Z_F(void);
extern void UPDATE_PML_H_X_N(void);
extern void UPDATE_PML_H_Y_N(void);
extern void UPDATE_PML_H_Z_N(void);
extern void PBC_B_X(void);
extern void PBC_B_Y(void);
extern void PBC_B_Z(void);

extern void SETUP_CPML_X(void);
extern void SETUP_CPML_Y(void);
extern void SETUP_CPML_Z(void);

//TFSF corrections
extern void SETUP_TFSF(void);
extern void TFSF_CORRECT(void);
extern void SOURCE_IN_E(void);
extern void SOURCE_IN_B(void);

extern void FT_source_calc();

extern int* MALLOC3D_int(int *grid, int sizeX,int sizeY,int sizeZ);
extern void MATERIAL_MATRIX(void);

extern void CORRECT_Y(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real  *e_inc,real  *h_inc,real  *ex,real  *ey,real  *ez,real  *hx,real  *hy,real  *hz,real  inc_theta,
  real  inc_phi,real  polar_psi,real  polar_theta,real  dx,real  dy,real  dz,real  dt,int i0,int j0,int k0,real  d_1D,int m0,real  *Cexe,real  *Ceye,real  *Ceze,real  *Cexh,real  *Ceyh,real  *Cezh,real  *Chxe,real  *Chye,real  *Chze,real  *Chxh,real  *Chyh,real  *Chzh);
extern void CORRECT_Z(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real  *e_inc,real  *h_inc,real  *ex,real  *ey,real  *ez,real  *hx,real  *hy,real  *hz,real  inc_theta,
  real  inc_phi,real  polar_psi,real  polar_theta,real  dx,real  dy,real  dz,real  dt,int i0,int j0,int k0,real  d_1D,int m0,real  *Cexe,real  *Ceye,real  *Ceze,real  *Cexh,real  *Ceyh,real  *Cezh,real  *Chxe,real  *Chye,real  *Chze,real  *Chxh,real  *Chyh,real  *Chzh);
extern void CORRECT_X(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real  *e_inc,real  *h_inc,real  *ex,real  *ey,real  *ez,real  *hx,real  *hy,real  *hz,real  inc_theta,
  real  inc_phi,real  polar_psi,real  polar_theta,real  dx,real  dy,real  dz,real  dt,int i0,int j0,int k0,real  d_1D,int m0,real  *Cexe,real  *Ceye,real  *Ceze,real  *Cexh,real  *Ceyh,real  *Cezh,real  *Chxe,real  *Chye,real  *Chze,real  *Chxh,real  *Chyh,real  *Chzh);


	extern __global__  void UPDATE_e_inc(real * e_inc,real * h_inc,int inc_length,real  d_1D,real  c0,real  dt,real  z0,real  ep0,int t,real delay,real  width,real  pi,real  f_0,int m0);
extern __global__  void UPDATE_h_inc(real * e_inc,real * h_inc,int inc_length,real  d_1D,real  c0,real  dt,real  z0,real  ep0,int t,real delay,real  width,real  pi,real  f_0,real  mu0);
extern void SETUP_t_inc(void);


//Drude and Criticsl Point multi-pole
extern void SETUP_Drude_CP(void);
extern int N_drude_poles,N_CP_poles,N_lorentz_poles;
extern real  *w_D,*gamma_d,*OMEGA_cp,*A_cp,*phi_cp,*GAMMA_cp,eps_inf,*d_eps_L,*omg_L,*delta_L;
extern real  *a_0_cp,*a_1_cp,*b_0_cp,*b_1_cp,*b_2_cp;
extern real  *C_1_cp,*C_2_cp,*C_3_cp,*C_4_cp,*C_5_cp,*C_cp;
extern real  C_E,C_E_1,C_E_2,C1_NL,C2_NL;
extern real  *d_1_d,*d_2_d,*d_3_d,*d_4_d,*d_5_d,*d_d, *d_NL;
extern real  *alpha_L,*psi_L,*eta_L,*alpha_HD1,*alpha_HD2,*psi_HD,*eta_HD;

extern comp *Px_cp,*Px_cp_n,*Px_cp_n_1,*Py_cp,*Py_cp_n,*Py_cp_n_1,*Pz_cp,*Pz_cp_n,*Pz_cp_n_1;
extern comp *Px_d,*Px_d_n,*Px_d_n_1,*Py_d,*Py_d_n,*Py_d_n_1,*Pz_d,*Pz_d_n,*Pz_d_n_1;
extern comp *Px_NL,*Px_NL_n,*Px_NL_n_1,*Py_NL,*Py_NL_n,*Py_NL_n_1,*Pz_NL,*Pz_NL_n,*Pz_NL_n_1;
extern comp *Jx_NL,*Jx_NL_n,*Jx_NL_n_1,*Jy_NL,*Jy_NL_n,*Jy_NL_n_1,*Jz_NL,*Jz_NL_n,*Jz_NL_n_1;
extern comp *Jx_Lo,*Jx_Lo_n,*Jx_Lo_n_1,*Jy_Lo,*Jy_Lo_n,*Jy_Lo_n_1,*Jz_Lo,*Jz_Lo_n,*Jz_Lo_n_1;
extern int dispersive_slab;

extern real TimeFactor;


extern comp *Px_cpdev,*Px_cp_ndev,*Px_cp_n_1dev,*Py_cpdev, *Py_cp_ndev, *Py_cp_n_1dev,*Pz_cpdev,*Pz_cp_ndev,*Pz_cp_n_1dev;
extern comp *Px_ddev,*Px_d_ndev,*Px_d_n_1dev,*Py_ddev,*Py_d_ndev,*Py_d_n_1dev,*Pz_ddev,*Pz_d_ndev,*Pz_d_n_1dev,*Pz_d_n_2dev,*Py_d_n_2dev,*Px_d_n_2dev;
extern real  *C_1_cpdev,*C_2_cpdev,*C_3_cpdev,*C_4_cpdev,*C_5_cpdev,*C_cpdev;
extern real  *d_1_ddev,*d_2_ddev,*d_3_ddev,*d_4_ddev,*d_5_ddev,*d_ddev, *d_NLdev;

extern void Fourier_Transform(void);
extern void Reflected_Calculation(void);
extern comp PULSE(int t);
extern void Reflectance_XZ(void);
extern double complex e_reflected,*E_reflected,*E_transmitted, *Ex_Reflected, *Hx_Reflected,*Ey_Reflected, *Hy_Reflected;
extern double complex e_incident,*E_incident, *E_Incident, *H_Incident;
extern double complex  *Ex_Transmitted, *Hx_Transmitted,*Ey_Transmitted, *Hy_Transmitted;
extern real  *t_inc;
extern comp e_total;
extern int NUM_freq;
extern real  dw, *freq, f_min;

extern int Scattering,Absorption;
extern int XSTARTAbs,XSTARTSca, XENDAbs,XENDSca,YSTARTAbs,YSTARTSca,YENDAbs,YENDSca,ZSTARTAbs,ZSTARTSca,ZENDAbs,ZENDSca;
extern int ZNEARAbs,ZFARAbs,ZNEARSca,ZFARSca,YNEARAbs,YFARAbs,YNEARSca,YFARSca,XNEARAbs,XFARAbs,XNEARSca,XFARSca;

extern real2 *ExTransformNearZScaRe;
extern real2 *EyTransformNearZScaRe;
extern real2 *EzTransformNearZScaRe;
extern real2 *HxTransformNearZScaRe;
extern real2 *HyTransformNearZScaRe;
extern real2 *HzTransformNearZScaRe;

extern real2 *ExTransformFarZScaRe;
extern real2 *EyTransformFarZScaRe;
extern real2 *EzTransformFarZScaRe;
extern real2 *HxTransformFarZScaRe;
extern real2 *HyTransformFarZScaRe;
extern real2 *HzTransformFarZScaRe;

extern real2 *ExTransformNearYScaRe;
extern real2 *EyTransformNearYScaRe;
extern real2 *EzTransformNearYScaRe;
extern real2 *HxTransformNearYScaRe;
extern real2 *HyTransformNearYScaRe;
extern real2 *HzTransformNearYScaRe;

extern real2 *ExTransformFarYScaRe;
extern real2 *EyTransformFarYScaRe;
extern real2 *EzTransformFarYScaRe;
extern real2 *HxTransformFarYScaRe;
extern real2 *HyTransformFarYScaRe;
extern real2 *HzTransformFarYScaRe;

extern real2 *ExTransformNearXScaRe;
extern real2 *EyTransformNearXScaRe;
extern real2 *EzTransformNearXScaRe;
extern real2 *HxTransformNearXScaRe;
extern real2 *HyTransformNearXScaRe;
extern real2 *HzTransformNearXScaRe;

extern real2 *ExTransformFarXScaRe;
extern real2 *EyTransformFarXScaRe;
extern real2 *EzTransformFarXScaRe;
extern real2 *HxTransformFarXScaRe;
extern real2 *HyTransformFarXScaRe;
extern real2 *HzTransformFarXScaRe;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
extern real2 *ExTransformNearZAbsRe;
extern real2 *EyTransformNearZAbsRe;
extern real2 *EzTransformNearZAbsRe;
extern real2 *HxTransformNearZAbsRe;
extern real2 *HyTransformNearZAbsRe;
extern real2 *HzTransformNearZAbsRe;

extern real2 *ExTransformFarZAbsRe;
extern real2 *EyTransformFarZAbsRe;
extern real2 *EzTransformFarZAbsRe;
extern real2 *HxTransformFarZAbsRe;
extern real2 *HyTransformFarZAbsRe;
extern real2 *HzTransformFarZAbsRe;

extern real2 *ExTransformNearYAbsRe;
extern real2 *EyTransformNearYAbsRe;
extern real2 *EzTransformNearYAbsRe;
extern real2 *HxTransformNearYAbsRe;
extern real2 *HyTransformNearYAbsRe;
extern real2 *HzTransformNearYAbsRe;

extern real2 *ExTransformFarYAbsRe;
extern real2 *EyTransformFarYAbsRe;
extern real2 *EzTransformFarYAbsRe;
extern real2 *HxTransformFarYAbsRe;
extern real2 *HyTransformFarYAbsRe;
extern real2 *HzTransformFarYAbsRe;

extern real2 *ExTransformNearXAbsRe;
extern real2 *EyTransformNearXAbsRe;
extern real2 *EzTransformNearXAbsRe;
extern real2 *HxTransformNearXAbsRe;
extern real2 *HyTransformNearXAbsRe;
extern real2 *HzTransformNearXAbsRe;

extern real2 *ExTransformFarXAbsRe;
extern real2 *EyTransformFarXAbsRe;
extern real2 *EzTransformFarXAbsRe;
extern real2 *HxTransformFarXAbsRe;
extern real2 *HyTransformFarXAbsRe;
extern real2 *HzTransformFarXAbsRe;

extern real2 *ExTransformNearZScaIm;
extern real2 *EyTransformNearZScaIm;
extern real2 *EzTransformNearZScaIm;
extern real2 *HxTransformNearZScaIm;
extern real2 *HyTransformNearZScaIm;
extern real2 *HzTransformNearZScaIm;

extern real2 *ExTransformFarZScaIm;
extern real2 *EyTransformFarZScaIm;
extern real2 *EzTransformFarZScaIm;
extern real2 *HxTransformFarZScaIm;
extern real2 *HyTransformFarZScaIm;
extern real2 *HzTransformFarZScaIm;

extern real2 *ExTransformNearYScaIm;
extern real2 *EyTransformNearYScaIm;
extern real2 *EzTransformNearYScaIm;
extern real2 *HxTransformNearYScaIm;
extern real2 *HyTransformNearYScaIm;
extern real2 *HzTransformNearYScaIm;

extern real2 *ExTransformFarYScaIm;
extern real2 *EyTransformFarYScaIm;
extern real2 *EzTransformFarYScaIm;
extern real2 *HxTransformFarYScaIm;
extern real2 *HyTransformFarYScaIm;
extern real2 *HzTransformFarYScaIm;

extern real2 *ExTransformNearXScaIm;
extern real2 *EyTransformNearXScaIm;
extern real2 *EzTransformNearXScaIm;
extern real2 *HxTransformNearXScaIm;
extern real2 *HyTransformNearXScaIm;
extern real2 *HzTransformNearXScaIm;

extern real2 *ExTransformFarXScaIm;
extern real2 *EyTransformFarXScaIm;
extern real2 *EzTransformFarXScaIm;
extern real2 *HxTransformFarXScaIm;
extern real2 *HyTransformFarXScaIm;
extern real2 *HzTransformFarXScaIm;



extern real2 *ExTransformNearZAbsIm;
extern real2 *EyTransformNearZAbsIm;
extern real2 *EzTransformNearZAbsIm;
extern real2 *HxTransformNearZAbsIm;
extern real2 *HyTransformNearZAbsIm;
extern real2 *HzTransformNearZAbsIm;

extern real2 *ExTransformFarZAbsIm;
extern real2 *EyTransformFarZAbsIm;
extern real2 *EzTransformFarZAbsIm;
extern real2 *HxTransformFarZAbsIm;
extern real2 *HyTransformFarZAbsIm;
extern real2 *HzTransformFarZAbsIm;

extern real2 *ExTransformNearYAbsIm;
extern real2 *EyTransformNearYAbsIm;
extern real2 *EzTransformNearYAbsIm;
extern real2 *HxTransformNearYAbsIm;
extern real2 *HyTransformNearYAbsIm;
extern real2 *HzTransformNearYAbsIm;

extern real2 *ExTransformFarYAbsIm;
extern real2 *EyTransformFarYAbsIm;
extern real2 *EzTransformFarYAbsIm;
extern real2 *HxTransformFarYAbsIm;
extern real2 *HyTransformFarYAbsIm;
extern real2 *HzTransformFarYAbsIm;

extern real2 *ExTransformNearXAbsIm;
extern real2 *EyTransformNearXAbsIm;
extern real2 *EzTransformNearXAbsIm;
extern real2 *HxTransformNearXAbsIm;
extern real2 *HyTransformNearXAbsIm;
extern real2 *HzTransformNearXAbsIm;

extern real2 *ExTransformFarXAbsIm;
extern real2 *EyTransformFarXAbsIm;
extern real2 *EzTransformFarXAbsIm;
extern real2 *HxTransformFarXAbsIm;
extern real2 *HyTransformFarXAbsIm;
extern real2 *HzTransformFarXAbsIm;







extern __global__ void ScattAbs(real  *ex,real  *ey,real  *ez,real  *hx,real  *hy,real  *hz,int NUM_freq,int t,real  dt,real * freq,real  pi,int XSTARTAbs,int XENDAbs,int YSTARTAbs,int YENDAbs,int ZSTARTAbs,int ZENDAbs,int XSTARTSca,int XENDSca,int YSTARTSca,int YENDSca,int ZSTARTSca,int ZENDSca,int XNEARAbs,int XFARAbs,int YNEARAbs,int YFARAbs,int ZNEARAbs,int ZFARAbs,
int XNEARSca,int XFARSca,int YNEARSca,int YFARSca,int ZNEARSca,int ZFARSca,real2  *exTransformNearZAbsRe,real2  *exTransformNearZAbsIm,real2  *eyTransformNearZAbsRe,real2  *eyTransformNearZAbsIm,real2  *hxTransformNearZAbsRe,real2  *hxTransformNearZAbsIm,real2  *hyTransformNearZAbsRe,real2  *hyTransformNearZAbsIm,
real2  *exTransformFarZAbsRe,real2  *exTransformFarZAbsIm,real2  *eyTransformFarZAbsRe,real2  *eyTransformFarZAbsIm,real2  *hxTransformFarZAbsRe,real2  *hxTransformFarZAbsIm,real2  *hyTransformFarZAbsRe,real2  *hyTransformFarZAbsIm,
real2  *exTransformNearYAbsRe,real2  *exTransformNearYAbsIm,real2  *ezTransformNearYAbsRe,real2  *ezTransformNearYAbsIm,real2  *hxTransformNearYAbsRe,real2  *hxTransformNearYAbsIm,real2  *hzTransformNearYAbsRe,real2  *hzTransformNearYAbsIm,
real2  *exTransformFarYAbsRe,real2  *exTransformFarYAbsIm,real2  *ezTransformFarYAbsRe,real2  *ezTransformFarYAbsIm,real2  *hxTransformFarYAbsRe,real2  *hxTransformFarYAbsIm,real2  *hzTransformFarYAbsRe,real2  *hzTransformFarYAbsIm,
real2  *eyTransformNearXAbsRe,real2  *eyTransformNearXAbsIm,real2  *ezTransformNearXAbsRe,real2  *ezTransformNearXAbsIm,real2  *hyTransformNearXAbsRe,real2  *hyTransformNearXAbsIm,real2  *hzTransformNearXAbsRe,real2  *hzTransformNearXAbsIm,
real2  *eyTransformFarXAbsRe,real2  *eyTransformFarXAbsIm,real2  *ezTransformFarXAbsRe,real2  *ezTransformFarXAbsIm,real2  *hyTransformFarXAbsRe,real2  *hyTransformFarXAbsIm,real2  *hzTransformFarXAbsRe,real2  *hzTransformFarXAbsIm,
real2  *exTransformNearZScaRe,real2  *exTransformNearZScaIm,real2  *eyTransformNearZScaRe,real2  *eyTransformNearZScaIm,real2  *hxTransformNearZScaRe,real2  *hxTransformNearZScaIm,real2  *hyTransformNearZScaRe,real2  *hyTransformNearZScaIm,
real2  *exTransformFarZScaRe,real2  *exTransformFarZScaIm,real2  *eyTransformFarZScaRe,real2  *eyTransformFarZScaIm,real2  *hxTransformFarZScaRe,real2  *hxTransformFarZScaIm,real2  *hyTransformFarZScaRe,real2  *hyTransformFarZScaIm,
real2  *exTransformNearYScaRe,real2  *exTransformNearYScaIm,real2  *ezTransformNearYScaRe,real2  *ezTransformNearYScaIm,real2  *hxTransformNearYScaRe,real2  *hxTransformNearYScaIm,real2  *hzTransformNearYScaRe,real2  *hzTransformNearYScaIm,
real2  *exTransformFarYScaRe,real2  *exTransformFarYScaIm,real2  *ezTransformFarYScaRe,real2  *ezTransformFarYScaIm,real2  *hxTransformFarYScaRe,real2  *hxTransformFarYScaIm,real2  *hzTransformFarYScaRe,real2  *hzTransformFarYScaIm,
real2  *eyTransformNearXScaRe,real2  *eyTransformNearXScaIm,real2  *ezTransformNearXScaRe,real2  *ezTransformNearXScaIm,real2  *hyTransformNearXScaRe,real2  *hyTransformNearXScaIm,real2  *hzTransformNearXScaRe,real2  *hzTransformNearXScaIm,
real2  *eyTransformFarXScaRe,real2  *eyTransformFarXScaIm,real2  *ezTransformFarXScaRe,real2  *ezTransformFarXScaIm,real2  *hyTransformFarXScaRe,real2  *hyTransformFarXScaIm,real2  *hzTransformFarXScaRe,real2  *hzTransformFarXScaIm,int NCELLX,int NCELLY,int NCELLZ);










extern real2 *ExTransformNearZScaRedev;
extern real2 *EyTransformNearZScaRedev;
extern real2 *EzTransformNearZScaRedev;
extern real2 *HxTransformNearZScaRedev;
extern real2 *HyTransformNearZScaRedev;
extern real2 *HzTransformNearZScaRedev;

extern real2 *ExTransformFarZScaRedev;
extern real2 *EyTransformFarZScaRedev;
extern real2 *EzTransformFarZScaRedev;
extern real2 *HxTransformFarZScaRedev;
extern real2 *HyTransformFarZScaRedev;
extern real2 *HzTransformFarZScaRedev;

extern real2 *ExTransformNearYScaRedev;
extern real2 *EyTransformNearYScaRedev;
extern real2 *EzTransformNearYScaRedev;
extern real2 *HxTransformNearYScaRedev;
extern real2 *HyTransformNearYScaRedev;
extern real2 *HzTransformNearYScaRedev;

extern real2 *ExTransformFarYScaRedev;
extern real2 *EyTransformFarYScaRedev;
extern real2 *EzTransformFarYScaRedev;
extern real2 *HxTransformFarYScaRedev;
extern real2 *HyTransformFarYScaRedev;
extern real2 *HzTransformFarYScaRedev;

extern real2 *ExTransformNearXScaRedev;
extern real2 *EyTransformNearXScaRedev;
extern real2 *EzTransformNearXScaRedev;
extern real2 *HxTransformNearXScaRedev;
extern real2 *HyTransformNearXScaRedev;
extern real2 *HzTransformNearXScaRedev;

extern real2 *ExTransformFarXScaRedev;
extern real2 *EyTransformFarXScaRedev;
extern real2 *EzTransformFarXScaRedev;
extern real2 *HxTransformFarXScaRedev;
extern real2 *HyTransformFarXScaRedev;
extern real2 *HzTransformFarXScaRedev;

////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
extern real2 *ExTransformNearZAbsRedev;
extern real2 *EyTransformNearZAbsRedev;
extern real2 *EzTransformNearZAbsRedev;
extern real2 *HxTransformNearZAbsRedev;
extern real2 *HyTransformNearZAbsRedev;
extern real2 *HzTransformNearZAbsRedev;

extern real2 *ExTransformFarZAbsRedev;
extern real2 *EyTransformFarZAbsRedev;
extern real2 *EzTransformFarZAbsRedev;
extern real2 *HxTransformFarZAbsRedev;
extern real2 *HyTransformFarZAbsRedev;
extern real2 *HzTransformFarZAbsRedev;

extern real2 *ExTransformNearYAbsRedev;
extern real2 *EyTransformNearYAbsRedev;
extern real2 *EzTransformNearYAbsRedev;
extern real2 *HxTransformNearYAbsRedev;
extern real2 *HyTransformNearYAbsRedev;
extern real2 *HzTransformNearYAbsRedev;

extern real2 *ExTransformFarYAbsRedev;
extern real2 *EyTransformFarYAbsRedev;
extern real2 *EzTransformFarYAbsRedev;
extern real2 *HxTransformFarYAbsRedev;
extern real2 *HyTransformFarYAbsRedev;
extern real2 *HzTransformFarYAbsRedev;

extern real2 *ExTransformNearXAbsRedev;
extern real2 *EyTransformNearXAbsRedev;
extern real2 *EzTransformNearXAbsRedev;
extern real2 *HxTransformNearXAbsRedev;
extern real2 *HyTransformNearXAbsRedev;
extern real2 *HzTransformNearXAbsRedev;

extern real2 *ExTransformFarXAbsRedev;
extern real2 *EyTransformFarXAbsRedev;
extern real2 *EzTransformFarXAbsRedev;
extern real2 *HxTransformFarXAbsRedev;
extern real2 *HyTransformFarXAbsRedev;
extern real2 *HzTransformFarXAbsRedev;

extern real2 *ExTransformNearZScaImdev;
extern real2 *EyTransformNearZScaImdev;
extern real2 *EzTransformNearZScaImdev;
extern real2 *HxTransformNearZScaImdev;
extern real2 *HyTransformNearZScaImdev;
extern real2 *HzTransformNearZScaImdev;

extern real2 *ExTransformFarZScaImdev;
extern real2 *EyTransformFarZScaImdev;
extern real2 *EzTransformFarZScaImdev;
extern real2 *HxTransformFarZScaImdev;
extern real2 *HyTransformFarZScaImdev;
extern real2 *HzTransformFarZScaImdev;

extern real2 *ExTransformNearYScaImdev;
extern real2 *EyTransformNearYScaImdev;
extern real2 *EzTransformNearYScaImdev;
extern real2 *HxTransformNearYScaImdev;
extern real2 *HyTransformNearYScaImdev;
extern real2 *HzTransformNearYScaImdev;

extern real2 *ExTransformFarYScaImdev;
extern real2 *EyTransformFarYScaImdev;
extern real2 *EzTransformFarYScaImdev;
extern real2 *HxTransformFarYScaImdev;
extern real2 *HyTransformFarYScaImdev;
extern real2 *HzTransformFarYScaImdev;

extern real2 *ExTransformNearXScaImdev;
extern real2 *EyTransformNearXScaImdev;
extern real2 *EzTransformNearXScaImdev;
extern real2 *HxTransformNearXScaImdev;
extern real2 *HyTransformNearXScaImdev;
extern real2 *HzTransformNearXScaImdev;

extern real2 *ExTransformFarXScaImdev;
extern real2 *EyTransformFarXScaImdev;
extern real2 *EzTransformFarXScaImdev;
extern real2 *HxTransformFarXScaImdev;
extern real2 *HyTransformFarXScaImdev;
extern real2 *HzTransformFarXScaImdev;



extern real2 *ExTransformNearZAbsImdev;
extern real2 *EyTransformNearZAbsImdev;
extern real2 *EzTransformNearZAbsImdev;
extern real2 *HxTransformNearZAbsImdev;
extern real2 *HyTransformNearZAbsImdev;
extern real2 *HzTransformNearZAbsImdev;

extern real2 *ExTransformFarZAbsImdev;
extern real2 *EyTransformFarZAbsImdev;
extern real2 *EzTransformFarZAbsImdev;
extern real2 *HxTransformFarZAbsImdev;
extern real2 *HyTransformFarZAbsImdev;
extern real2 *HzTransformFarZAbsImdev;

extern real2 *ExTransformNearYAbsImdev;
extern real2 *EyTransformNearYAbsImdev;
extern real2 *EzTransformNearYAbsImdev;
extern real2 *HxTransformNearYAbsImdev;
extern real2 *HyTransformNearYAbsImdev;
extern real2 *HzTransformNearYAbsImdev;

extern real2 *ExTransformFarYAbsImdev;
extern real2 *EyTransformFarYAbsImdev;
extern real2 *EzTransformFarYAbsImdev;
extern real2 *HxTransformFarYAbsImdev;
extern real2 *HyTransformFarYAbsImdev;
extern real2 *HzTransformFarYAbsImdev;

extern real2 *ExTransformNearXAbsImdev;
extern real2 *EyTransformNearXAbsImdev;
extern real2 *EzTransformNearXAbsImdev;
extern real2 *HxTransformNearXAbsImdev;
extern real2 *HyTransformNearXAbsImdev;
extern real2 *HzTransformNearXAbsImdev;

extern real2 *ExTransformFarXAbsImdev;
extern real2 *EyTransformFarXAbsImdev;
extern real2 *EzTransformFarXAbsImdev;
extern real2 *HxTransformFarXAbsImdev;
extern real2 *HyTransformFarXAbsImdev;
extern real2 *HzTransformFarXAbsImdev;





extern void UpdateHydroPx(void);
extern void UpdateHydroPy(void);
extern void UpdateHydroPz(void);

extern void CalculateAbsScatt(void);
extern int Hydrodynamics;
extern int WithMagField;
extern int WithConvection;
extern __host__ __device__   int ThreeDMap(int i,int j, int k,int X,int Y);
extern __host__ __device__   int FourDMap(int i,int j, int k, int n,int X,int Y,int Z);
extern __host__ __device__   int TwoDMap(int i,int j,int X);

extern __device__ int ThreeDMapD(int i,int j, int k,int X,int Y);
extern __device__ int FourDMapD(int i,int j, int k, int n,int X,int Y,int Z);
extern __device__ int TwoDMapD(int i,int j,int X);

extern __global__ void UpdateHydroPx(real *ex,real *ex_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
	real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me);
	extern __global__ void UpdateHydroPy(real *ex,real *ex_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
		real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me);
	extern	__global__ void UpdateHydroPz(real *ex,real *ex_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
			real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me);

// int Hydrodynamics;
extern real N_EQ;

#endif // _EXTERN_VAR_H
