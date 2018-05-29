#ifndef _GLOBAL_VAR_H
#define _GLOBAL_VAR_H

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

FILE *Source;

int StaticField;

int PBC_CTW;
int Laplacian,Diverge_Gradient;
void READ_DATA_FILE(void);
void SOURCE_SETUP(void);

int material,WL_or_freq;
int trigger;
real Nonlocalend;
int Nonlocalend_int;

int Snap_in,Test_offset;
int Tend_inc,NUM_freq_inc,Freq_start,Spect_loc;

comp *Hx_source, *Hy_source, *Hz_source;
comp *Ex_source, *Ey_source, *Ez_source;

int * mat_matrix,* mat_matrixX,* mat_matrixY,* mat_matrixZ, first_medium,first_medium_max,nano_sphere;
int * mat_matrixdev,* mat_matrixXdev,* mat_matrixYdev,* mat_matrixZdev;
real nano_sphere_radius;

real BandWidth;
int NONLOCAL;
int num_trials,min_trials,max_trials,TE_TM,trials;
int inf_disp_slab,centerx,centery;
//polarization type;
int TMz,TEz;
//number of cells in simulation domain
int NCELLX, NCELLY, NCELLZ;
//number of TFSF cells in each dimension
int NtfsfX, NtfsfY,NtfsfZ;
//number of CPML cells in each domain
int NcpmlX,NcpmlY,NcpmlZ;
int cpml_N_X, cpml_F_X,cpml_N_Y,cpml_F_Y,cpml_N_Z,cpml_F_Z;
int cpml_x_lim,cpml_y_lim,cpml_z_lim;
//incident plane position
int inc_plane,source_end;
real f_0;
//E and H fields
comp *ex,*ey,*ez,*hx,*hy,*hz;
comp  *exdev,*eydev,*ezdev,*hxdev,*hydev,*hzdev;

comp *Dx,*Dy,*Dz;
comp *ex_n,*ey_n,*ez_n;
comp *ex_ndev,*ey_ndev,*ez_ndev;

comp *ex_n_1,*ey_n_1,*ez_n_1;
comp *ex_n_1dev,*ey_n_1dev,*ez_n_1dev;

comp *hxPrev,*hyPrev,*hzPrev;
comp *hxPrevdev,*hyPrevdev,*hzPrevdev;
comp *NDx,*NDy,*NDz;
comp *NDxdev,*NDydev,*NDzdev;
comp *NDx_prev,*NDy_prev,*NDz_prev;
comp *NDx_prevdev,*NDy_prevdev,*NDz_prevdev;



comp *Hx_w_FT,*Hy_w_FT,*Hz_w_FT,*Ex_w_FT,*Ey_w_FT,*Ez_w_FT;
comp *Hx_t_FT,*Hy_t_FT,*Hz_t_FT,*Ex_t_FT,*Ey_t_FT,*Ez_t_FT;


//update coefficients
real *Cexe,*Cexh,*Ceye,*Ceyh,*Ceze,*Cezh;
real *Chxe,*Chxh,*Chye,*Chyh,*Chze,*Chzh;
real *Chxedev,*Chxhdev,*Chyedev,*Chyhdev,*Chzedev,*Chzhdev;
real *Cexedev,*Cexhdev,*Ceyedev,*Ceyhdev,*Cezedev,*Cezhdev;


real *sigma_e, *sigma_m;
real *eps,*mu;
//CPML stuff
comp *psi_Ex_y_N,*psi_Ex_z_N,*psi_Ey_z_N,*psi_Ey_x_N,*psi_Ez_x_N,*psi_Ez_y_N;
comp *psi_Hx_y_N,*psi_Hx_z_N,*psi_Hy_z_N,*psi_Hy_x_N,*psi_Hz_x_N,*psi_Hz_y_N;
comp *psi_Ex_y_F,*psi_Ex_z_F,*psi_Ey_z_F,*psi_Ey_x_F,*psi_Ez_x_F,*psi_Ez_y_F;
comp *psi_Hx_y_F,*psi_Hx_z_F,*psi_Hy_z_F,*psi_Hy_x_F,*psi_Hz_x_F,*psi_Hz_y_F;

comp *psi_Ex_y_Ndev,*psi_Ex_z_Ndev,*psi_Ey_z_Ndev,*psi_Ey_x_Ndev,*psi_Ez_x_Ndev,*psi_Ez_y_Ndev;
comp *psi_Hx_y_Ndev,*psi_Hx_z_Ndev,*psi_Hy_z_Ndev,*psi_Hy_x_Ndev,*psi_Hz_x_Ndev,*psi_Hz_y_Ndev;
comp *psi_Ex_y_Fdev,*psi_Ex_z_Fdev,*psi_Ey_z_Fdev,*psi_Ey_x_Fdev,*psi_Ez_x_Fdev,*psi_Ez_y_Fdev;
comp *psi_Hx_y_Fdev,*psi_Hx_z_Fdev,*psi_Hy_z_Fdev,*psi_Hy_x_Fdev,*psi_Hz_x_Fdev,*psi_Hz_y_Fdev;

real *kedx,*kedy,*kedz;
real *khdx,*khdy,*khdz;
real *be_x_N,*be_y_N,*be_z_N;
real *bh_x_N,*bh_y_N,*bh_z_N;
real *ce_x_N,*ce_y_N,*ce_z_N;
real *ch_x_N,*ch_y_N,*ch_z_N;
real *be_x_F,*be_y_F,*be_z_F;
real *bh_x_F,*bh_y_F,*bh_z_F;
real *ce_x_F,*ce_y_F,*ce_z_F;
real *ch_x_F,*ch_y_F,*ch_z_F;

real *kedxdev,*kedydev,*kedzdev;
real *khdxdev,*khdydev,*khdzdev;
real *be_x_Ndev,*be_y_Ndev,*be_z_Ndev;
real *bh_x_Ndev,*bh_y_Ndev,*bh_z_Ndev;
real *ce_x_Ndev,*ce_y_Ndev,*ce_z_Ndev;
real *ch_x_Ndev,*ch_y_Ndev,*ch_z_Ndev;
real *be_x_Fdev,*be_y_Fdev,*be_z_Fdev;
real *bh_x_Fdev,*bh_y_Fdev,*bh_z_Fdev;
real *ce_x_Fdev,*ce_y_Fdev,*ce_z_Fdev;
real *ch_x_Fdev,*ch_y_Fdev,*ch_z_Fdev;

real *sigma_e_x, *sigma_e_y,*sigma_e_z;
real *sigma_h_x, *sigma_h_y,*sigma_h_z;
real cpml_exp;
real max_stretch_factor_x,max_stretch_factor_y,max_stretch_factor_z;
real max_sigma_cpml_x,max_sigma_cpml_y,max_sigma_cpml_z;
real max_alpha_x, max_alpha_y, max_alpha_z;
real exp_alpha_x, exp_alpha_y, exp_alpha_z;
real MAX_AMP;

//temporal and spacial step sizes
real dt,dx,dy,dz,d_1D;

//time stepping variables
int t,Tend;

//Snapshot variables
int x_skip,y_skip,z_skip,t_skip,snapshot_count;

//source fields and parameters;
real inc_phi, inc_theta, polar_psi; //incident angles and polarization
comp *e_inc, *h_inc;
comp *e_incdev, *h_incdev;

comp *hx_inc, *hy_inc, *hz_inc;
comp *ex_inc, *ey_inc, *ez_inc;
int inc_Length,m0,i_0,j_0,k_0;
//buffers for auxilliary field
real e1,e2;
real delay, width, width_real;

//repeat of constants
real z0,ep0,mu0,c0,pi,me,e0;

//Periodic Boundary Conditions
real  k_x,k_y,k_z,k_rho,period_x,period_y,period_z;
int Periodic_XY,Periodic_YZ,Periodic_XZ;

//Functions
real* MALLOC1D(real *ARRAY, int i);
real* MALLOC1D_double(real *ARRAY, int i);
comp* MALLOC1D_Complex(comp *ARRAY, int i);
double complex* MALLOC1D_Complex2(double complex *ARRAY, int i);
real2* MALLOC1D_Real2(real2 *ARRAY, int i);


real* MALLOC2D(real *ARRAY, int i, int j);
real* MALLOC2D_double(real *ARRAY, int i, int j);
comp* MALLOC2D_Complex(comp *ARRAY, int i, int j);


real* MALLOC3D(real *ARRAY, int i, int j,int k);
real* MALLOC3D_double(real *ARRAY, int i, int j,int k);
comp* MALLOC3D_Complex(comp *ARRAY, int i, int j,int k);
real2* MALLOC3D_Real2(real2 *ARRAY, int i, int j,int k);
double complex* MALLOC3D_Complex2(double complex *ARRAY, int i, int j,int k);


real* MALLOC4D(real *ARRAY, int i, int j,int k, int size4);
real* MALLOC4D_double(real *ARRAY, int i, int j,int k, int size4);
comp* MALLOC4D_Complex(comp *ARRAY, int i, int j,int k, int size4);

void FREE4D(real *ARRAY,int i,int j,int size4);
void FREE4D_double(real *ARRAY,int i,int j,int size4);
void FREE4D_Complex(comp *ARRAY,int i, int j,int size4);

void FREE3D(real *ARRAY,int i,int j);
void FREE3D_double(real *ARRAY,int i,int j);
void FREE3D_Complex(comp *ARRAY,int i, int j);
void FREE3D_Complex2(double complex *ARRAY,int i, int j);


void FREE2D(real *ARRAY,int i);
void FREE2D_double(real *ARRAY,int i);
void FREE2D_Complex(comp *ARRAY,int i);

void FREE1D(real *ARRAY);
void FREE1D_double(real *ARRAY);
void FREE1D_Complex(comp *ARRAY);
void FREE1D_Complex2(double complex *ARRAY);


real* ZERO_VECTORS3D(real *ARRAY, int i, int j, int k);
real* ZERO_VECTORS3D_double(real *ARRAY, int i, int j, int k);
comp* ZERO_VECTORS3D_Complex(comp *ARRAY, int i, int j, int k);
real2* ZERO_VECTORS3D_Real2(real2 *ARRAY, int i, int j, int k);

double complex* ZERO_VECTORS3D_Complex2(double complex *ARRAY, int i, int j, int k);

real* ZERO_VECTORS1D(real *ARRAY);
real* ZERO_VECTORS1D_double(real *ARRAY);
comp* ZERO_VECTORS1D_Complex(comp *ARRAY,int length);
double complex* ZERO_VECTORS1D_Complex2(double complex *ARRAY,int length);

comp* ZERO_VECTORS2D_Complex(comp *ARRAY, int SizeX, int SizeY);

void ALLOCATE_MEM(void);
void FREE_MEM(void);

void SOURCE_IN_B();
void SOURCE_IN_E();
void SOURCE_IN();

comp Calc_DIV_GRADz(int i, int j, int k);
comp Calc_DIV_GRADx(int i, int j, int k);
comp Calc_DIV_GRADy(int i, int j, int k);
comp Calc_DIV_GRADz2(int i, int j, int k);
comp Calc_DIV_GRADx2(int i, int j, int k);
comp Calc_DIV_GRADy2(int i, int j, int k);

void FT_source_calc();

void SETUP_CONST();

void SETUP_CPML(void);
void SETUP_CPML_X(void);
void SETUP_CPML_Y(void);
void SETUP_CPML_Z(void);

void SETUP_SNAPSHOT(void);
void SNAPSHOT_2D();
void SNAPSHOT_2D_N();
void SNAPSHOT_1D();

real * DEF_UPDATE_COEFF_EonH(real *ARRAY);
real * DEF_UPDATE_COEFF_HonE(real *ARRAY);
real * DEF_UPDATE_COEFF_EonE(real *ARRAY);
real * DEF_UPDATE_COEFF_HonH(real *ARRAY);

void DEF_EPS(void);
void DEF_MU(void);
void DEF_SIGMA_E(void);
void DEF_SIGMA_M(void);

void UPDATE_E(void);
// void UPDATE_ex();
__global__ void UPDATE_ex(real *ex,real *ex_n,real *ex_n_1,real *hy,real *hz,real *Cexe,real *Cexh,real *kedy,real *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real *psi_Ex_z_N,
real *psi_Ex_z_F,real *psi_Ex_y_N,real *psi_Ex_y_F,real *Px_cp,real *Px_cp_n,real *Px_cp_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,
real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY,  real *NDx,real *NDy,real *NDz,real *NDx_prev,real *NDy_prev,real *NDz_prev,real e0,real N_EQ);
//void UPDATE_ey();
__global__ void UPDATE_ey(real *ex,real *ex_n,real *ex_n_1,real *hy,real *hz,real *Cexe,real *Cexh,real *kedy,real *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real *psi_Ex_z_N,
real *psi_Ex_z_F,real *psi_Ex_y_N,real *psi_Ex_y_F,real *Px_cp,real *Px_cp_n,real *Px_cp_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,
real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY,  real *NDx,real *NDy,real *NDz,real *NDx_prev,real *NDy_prev,real *NDz_prev,real e0,real N_EQ);

__global__ void UPDATE_ez(real *ex,real *ex_n,real *ex_n_1,real *hy,real *hz,real *Cexe,real *Cexh,real *kedy,real *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real *psi_Ex_z_N,
real *psi_Ex_z_F,real *psi_Ex_y_N,real *psi_Ex_y_F,real *Px_cp,real *Px_cp_n,real *Px_cp_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,
real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY,  real *NDx,real *NDy,real *NDz,real *NDx_prev,real *NDy_prev,real *NDz_prev,real e0,real N_EQ);


// void UPDATE_ez();
void UPDATE_PML_E_X_F(void);
void UPDATE_PML_E_Y_F(void);
void UPDATE_PML_E_Z_F(void);
void UPDATE_PML_E_X_N(void);
void UPDATE_PML_E_Y_N(void);
void UPDATE_PML_E_Z_N(void);
void UPDATE_Px(void);
void UPDATE_Py(void);
void UPDATE_Pz(void);
void PBC_E_X(void);
void PBC_E_Y(void);
void PBC_E_Z(void);
void CP_D_ex(int i, int j, int k, comp Curl_H,comp Div_Grad);
void CP_D_ey(int i, int j, int k, comp Curl_H,comp Div_Grad);
void CP_D_ez(int i, int j, int k, comp Curl_H,comp Div_Grad);


void UPDATE_B(void);
// void UPDATE_hx();

  __global__ void UPDATE_hx(real *hx,real *hxPrev,real *ex,real *ey,real *Chxh,real *Chxe,real *psi_Hx_z_N,real *psi_Hx_z_F,real *psi_Hx_y_N,real *psi_Hx_y_F,real *khdy,real
    *khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics);
// void UPDATE_hy();

__global__ void UPDATE_hy(real *hx,real *hxPrev,real *ex,real *ey,real *Chxh,real *Chxe,real *psi_Hx_z_N,real *psi_Hx_z_F,real *psi_Hx_y_N,real *psi_Hx_y_F,real *khdy,real
	*khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics);

// void UPDATE_hz();
__global__ void UPDATE_hz(real *hx,real *hxPrev,real *ex,real *ey,real *Chxh,real *Chxe,real *psi_Hx_z_N,real *psi_Hx_z_F,real *psi_Hx_y_N,real *psi_Hx_y_F,real *khdy,real
	*khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics);


void UPDATE_PML_H_X_F(void);
void UPDATE_PML_H_Y_F(void);
void UPDATE_PML_H_Z_F(void);
void UPDATE_PML_H_X_N(void);
void UPDATE_PML_H_Y_N(void);
void UPDATE_PML_H_Z_N(void);
void PBC_B_X(void);
void PBC_B_Y(void);
void PBC_B_Z(void);

void Fourier_Transform(void);
void Reflected_Calculation(void);
comp PULSE(int t);

__host__ __device__ int ThreeDMap(int i,int j, int k,int X,int Y);
__host__ __device__ int FourDMap(int i,int j, int k, int n,int X,int Y,int Z);
__host__ __device__ int TwoDMap(int i,int j,int X);

 __device__ int ThreeDMapD(int i,int j, int k,int X,int Y);
 __device__ int FourDMapD(int i,int j, int k, int n,int X,int Y,int Z);
 __device__ int TwoDMapD(int i,int j,int X);


double complex e_reflected,*E_reflected,*E_transmitted, *Ex_Reflected, *Hx_Reflected,*Ey_Reflected, *Hy_Reflected;
double complex e_incident,*E_incident, *E_Incident, *H_Incident;
double complex  *Ex_Transmitted, *Hx_Transmitted,*Ey_Transmitted, *Hy_Transmitted;


real *t_inc;
comp e_total;
real f_min;
real dw;
real *freq;
int NUM_freq;


//TFSF Functions
void SETUP_TFSF(void);
void TFSF_CORRECT(void);

 void CORRECT_Y(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real *e_inc,real *h_inc,real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,real inc_theta,
  real inc_phi,real polar_psi,real polar_theta,real dx,real dy,real dz,real dt,int i0,int j0,int k0,real d_1D,int m0,real *Cexe,real *Ceye,real *Ceze,real *Cexh,real *Ceyh,real *Cezh,real *Chxe,real *Chye,real *Chze,real *Chxh,real *Chyh,real *Chzh);
 void CORRECT_Z(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real *e_inc,real *h_inc,real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,real inc_theta,
  real inc_phi,real polar_psi,real polar_theta,real dx,real dy,real dz,real dt,int i0,int j0,int k0,real d_1D,int m0,real *Cexe,real *Ceye,real *Ceze,real *Cexh,real *Ceyh,real *Cezh,real *Chxe,real *Chye,real *Chze,real *Chxh,real *Chyh,real *Chzh);
 void CORRECT_X(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real *e_inc,real *h_inc,real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,real inc_theta,
  real inc_phi,real polar_psi,real polar_theta,real dx,real dy,real dz,real dt,int i0,int j0,int k0,real d_1D,int m0,real *Cexe,real *Ceye,real *Ceze,real *Cexh,real *Ceyh,real *Cezh,real *Chxe,real *Chye,real *Chze,real *Chxh,real *Chyh,real *Chzh);



	__global__  void UPDATE_e_inc(real* e_inc,real* h_inc,int inc_length,real d_1D,real c0,real dt,real z0,real ep0,int t,real delay,real width,real pi,real f_0,int m0);
__global__  void UPDATE_h_inc(real* e_inc,real* h_inc,int inc_length,real d_1D,real c0,real dt,real z0,real ep0,int t,real delay,real width,real pi,real f_0,real mu0);

void DIELECTRIC_SLAB(void);
int* MALLOC3D_int(int *grid, int sizeX,int sizeY,int sizeZ);
void MATERIAL_MATRIX(void);

real TimeFactor;

void Reflectance_XZ(void);

void SETUP_t_inc(void);

//drude and Critical point multi-pole
void SETUP_Drude_CP(void);
int N_drude_poles,N_CP_poles,N_lorentz_poles;
real *w_D,*gamma_d,*OMEGA_cp,*A_cp,*phi_cp,*GAMMA_cp,eps_inf,*d_eps_L,*omg_L,*delta_L;
real *a_0_cp,*a_1_cp,*b_0_cp,*b_1_cp,*b_2_cp;
real *C_1_cp,*C_2_cp,*C_3_cp,*C_4_cp,*C_5_cp,*C_cp;
real *d_1_d,*d_2_d,*d_3_d,*d_4_d,*d_5_d,*d_d, *d_NL;
real *C_1_cpdev,*C_2_cpdev,*C_3_cpdev,*C_4_cpdev,*C_5_cpdev,*C_cpdev;
real *d_1_ddev,*d_2_ddev,*d_3_ddev,*d_4_ddev,*d_5_ddev,*d_ddev, *d_NLdev;
real *alpha_L,*psi_L,*eta_L,*alpha_HD1,*alpha_HD2,*psi_HD,*eta_HD;
real C_E,C_E_1,C_E_2,C1_NL,C2_NL;
comp *Px_cp,*Px_cp_n,*Px_cp_n_1,*Py_cp, *Py_cp_n, *Py_cp_n_1,*Pz_cp,*Pz_cp_n,*Pz_cp_n_1;
comp *Px_d,*Px_d_n,*Px_d_n_1,*Py_d,*Py_d_n,*Py_d_n_1,*Pz_d,*Pz_d_n,*Pz_d_n_1;
comp *Px_cpdev,*Px_cp_ndev,*Px_cp_n_1dev,*Py_cpdev, *Py_cp_ndev, *Py_cp_n_1dev,*Pz_cpdev,*Pz_cp_ndev,*Pz_cp_n_1dev;
comp *Px_ddev,*Px_d_ndev,*Px_d_n_1dev,*Py_ddev,*Py_d_ndev,*Py_d_n_1dev,*Pz_ddev,*Pz_d_ndev,*Pz_d_n_1dev,*Pz_d_n_2dev,*Py_d_n_2dev,*Px_d_n_2dev;
comp *Px_NL,*Px_NL_n,*Px_NL_n_1,*Py_NL,*Py_NL_n,*Py_NL_n_1,*Pz_NL,*Pz_NL_n,*Pz_NL_n_1;
comp *Jx_NL,*Jx_NL_n,*Jx_NL_n_1,*Jy_NL,*Jy_NL_n,*Jy_NL_n_1,*Jz_NL,*Jz_NL_n,*Jz_NL_n_1;
comp *Jx_Lo,*Jx_Lo_n,*Jx_Lo_n_1,*Jy_Lo,*Jy_Lo_n,*Jy_Lo_n_1,*Jz_Lo,*Jz_Lo_n,*Jz_Lo_n_1;

int dispersive_slab;



int Scattering,Absorption;
int XSTARTAbs,XSTARTSca, XENDAbs,XENDSca,YSTARTAbs,YSTARTSca,YENDAbs,YENDSca,ZSTARTAbs,ZSTARTSca,ZENDAbs,ZENDSca;
int ZNEARAbs,ZFARAbs,ZNEARSca,ZFARSca,YNEARAbs,YFARAbs,YNEARSca,YFARSca,XNEARAbs,XFARAbs,XNEARSca,XFARSca;

real2 *ExTransformNearZScaRe;
real2 *EyTransformNearZScaRe;
real2 *EzTransformNearZScaRe;
real2 *HxTransformNearZScaRe;
real2 *HyTransformNearZScaRe;
real2 *HzTransformNearZScaRe;

real2 *ExTransformFarZScaRe;
real2 *EyTransformFarZScaRe;
real2 *EzTransformFarZScaRe;
real2 *HxTransformFarZScaRe;
real2 *HyTransformFarZScaRe;
real2 *HzTransformFarZScaRe;

real2 *ExTransformNearYScaRe;
real2 *EyTransformNearYScaRe;
real2 *EzTransformNearYScaRe;
real2 *HxTransformNearYScaRe;
real2 *HyTransformNearYScaRe;
real2 *HzTransformNearYScaRe;

real2 *ExTransformFarYScaRe;
real2 *EyTransformFarYScaRe;
real2 *EzTransformFarYScaRe;
real2 *HxTransformFarYScaRe;
real2 *HyTransformFarYScaRe;
real2 *HzTransformFarYScaRe;

real2 *ExTransformNearXScaRe;
real2 *EyTransformNearXScaRe;
real2 *EzTransformNearXScaRe;
real2 *HxTransformNearXScaRe;
real2 *HyTransformNearXScaRe;
real2 *HzTransformNearXScaRe;

real2 *ExTransformFarXScaRe;
real2 *EyTransformFarXScaRe;
real2 *EzTransformFarXScaRe;
real2 *HxTransformFarXScaRe;
real2 *HyTransformFarXScaRe;
real2 *HzTransformFarXScaRe;



real2 *ExTransformNearZScaIm;
real2 *EyTransformNearZScaIm;
real2 *EzTransformNearZScaIm;
real2 *HxTransformNearZScaIm;
real2 *HyTransformNearZScaIm;
real2 *HzTransformNearZScaIm;

real2 *ExTransformFarZScaIm;
real2 *EyTransformFarZScaIm;
real2 *EzTransformFarZScaIm;
real2 *HxTransformFarZScaIm;
real2 *HyTransformFarZScaIm;
real2 *HzTransformFarZScaIm;

real2 *ExTransformNearYScaIm;
real2 *EyTransformNearYScaIm;
real2 *EzTransformNearYScaIm;
real2 *HxTransformNearYScaIm;
real2 *HyTransformNearYScaIm;
real2 *HzTransformNearYScaIm;

real2 *ExTransformFarYScaIm;
real2 *EyTransformFarYScaIm;
real2 *EzTransformFarYScaIm;
real2 *HxTransformFarYScaIm;
real2 *HyTransformFarYScaIm;
real2 *HzTransformFarYScaIm;

real2 *ExTransformNearXScaIm;
real2 *EyTransformNearXScaIm;
real2 *EzTransformNearXScaIm;
real2 *HxTransformNearXScaIm;
real2 *HyTransformNearXScaIm;
real2 *HzTransformNearXScaIm;

real2 *ExTransformFarXScaIm;
real2 *EyTransformFarXScaIm;
real2 *EzTransformFarXScaIm;
real2 *HxTransformFarXScaIm;
real2 *HyTransformFarXScaIm;
real2 *HzTransformFarXScaIm;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
real2 *ExTransformNearZAbsRe;
real2 *EyTransformNearZAbsRe;
real2 *EzTransformNearZAbsRe;
real2 *HxTransformNearZAbsRe;
real2 *HyTransformNearZAbsRe;
real2 *HzTransformNearZAbsRe;

real2 *ExTransformFarZAbsRe;
real2 *EyTransformFarZAbsRe;
real2 *EzTransformFarZAbsRe;
real2 *HxTransformFarZAbsRe;
real2 *HyTransformFarZAbsRe;
real2 *HzTransformFarZAbsRe;

real2 *ExTransformNearYAbsRe;
real2 *EyTransformNearYAbsRe;
real2 *EzTransformNearYAbsRe;
real2 *HxTransformNearYAbsRe;
real2 *HyTransformNearYAbsRe;
real2 *HzTransformNearYAbsRe;

real2 *ExTransformFarYAbsRe;
real2 *EyTransformFarYAbsRe;
real2 *EzTransformFarYAbsRe;
real2 *HxTransformFarYAbsRe;
real2 *HyTransformFarYAbsRe;
real2 *HzTransformFarYAbsRe;

real2 *ExTransformNearXAbsRe;
real2 *EyTransformNearXAbsRe;
real2 *EzTransformNearXAbsRe;
real2 *HxTransformNearXAbsRe;
real2 *HyTransformNearXAbsRe;
real2 *HzTransformNearXAbsRe;

real2 *ExTransformFarXAbsRe;
real2 *EyTransformFarXAbsRe;
real2 *EzTransformFarXAbsRe;
real2 *HxTransformFarXAbsRe;
real2 *HyTransformFarXAbsRe;
real2 *HzTransformFarXAbsRe;


real2 *ExTransformNearZAbsIm;
real2 *EyTransformNearZAbsIm;
real2 *EzTransformNearZAbsIm;
real2 *HxTransformNearZAbsIm;
real2 *HyTransformNearZAbsIm;
real2 *HzTransformNearZAbsIm;

real2 *ExTransformFarZAbsIm;
real2 *EyTransformFarZAbsIm;
real2 *EzTransformFarZAbsIm;
real2 *HxTransformFarZAbsIm;
real2 *HyTransformFarZAbsIm;
real2 *HzTransformFarZAbsIm;

real2 *ExTransformNearYAbsIm;
real2 *EyTransformNearYAbsIm;
real2 *EzTransformNearYAbsIm;
real2 *HxTransformNearYAbsIm;
real2 *HyTransformNearYAbsIm;
real2 *HzTransformNearYAbsIm;

real2 *ExTransformFarYAbsIm;
real2 *EyTransformFarYAbsIm;
real2 *EzTransformFarYAbsIm;
real2 *HxTransformFarYAbsIm;
real2 *HyTransformFarYAbsIm;
real2 *HzTransformFarYAbsIm;

real2 *ExTransformNearXAbsIm;
real2 *EyTransformNearXAbsIm;
real2 *EzTransformNearXAbsIm;
real2 *HxTransformNearXAbsIm;
real2 *HyTransformNearXAbsIm;
real2 *HzTransformNearXAbsIm;

real2 *ExTransformFarXAbsIm;
real2 *EyTransformFarXAbsIm;
real2 *EzTransformFarXAbsIm;
real2 *HxTransformFarXAbsIm;
real2 *HyTransformFarXAbsIm;
real2 *HzTransformFarXAbsIm;




__global__ void ScattAbs(real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,int NUM_freq,int t,real dt,real* freq,real pi,int XSTARTAbs,int XENDAbs,int YSTARTAbs,int YENDAbs,int ZSTARTAbs,int ZENDAbs,int XSTARTSca,int XENDSca,int YSTARTSca,int YENDSca,int ZSTARTSca,int ZENDSca,int XNEARAbs,int XFARAbs,int YNEARAbs,int YFARAbs,int ZNEARAbs,int ZFARAbs,
int XNEARSca,int XFARSca,int YNEARSca,int YFARSca,int ZNEARSca,int ZFARSca,real2 *exTransformNearZAbsRe,real2 *exTransformNearZAbsIm,real2 *eyTransformNearZAbsRe,real2 *eyTransformNearZAbsIm,real2 *hxTransformNearZAbsRe,real2 *hxTransformNearZAbsIm,real2 *hyTransformNearZAbsRe,real2 *hyTransformNearZAbsIm,
real2 *exTransformFarZAbsRe,real2 *exTransformFarZAbsIm,real2 *eyTransformFarZAbsRe,real2 *eyTransformFarZAbsIm,real2 *hxTransformFarZAbsRe,real2 *hxTransformFarZAbsIm,real2 *hyTransformFarZAbsRe,real2 *hyTransformFarZAbsIm,
real2 *exTransformNearYAbsRe,real2 *exTransformNearYAbsIm,real2 *ezTransformNearYAbsRe,real2 *ezTransformNearYAbsIm,real2 *hxTransformNearYAbsRe,real2 *hxTransformNearYAbsIm,real2 *hzTransformNearYAbsRe,real2 *hzTransformNearYAbsIm,
real2 *exTransformFarYAbsRe,real2 *exTransformFarYAbsIm,real2 *ezTransformFarYAbsRe,real2 *ezTransformFarYAbsIm,real2 *hxTransformFarYAbsRe,real2 *hxTransformFarYAbsIm,real2 *hzTransformFarYAbsRe,real2 *hzTransformFarYAbsIm,
real2 *eyTransformNearXAbsRe,real2 *eyTransformNearXAbsIm,real2 *ezTransformNearXAbsRe,real2 *ezTransformNearXAbsIm,real2 *hyTransformNearXAbsRe,real2 *hyTransformNearXAbsIm,real2 *hzTransformNearXAbsRe,real2 *hzTransformNearXAbsIm,
real2 *eyTransformFarXAbsRe,real2 *eyTransformFarXAbsIm,real2 *ezTransformFarXAbsRe,real2 *ezTransformFarXAbsIm,real2 *hyTransformFarXAbsRe,real2 *hyTransformFarXAbsIm,real2 *hzTransformFarXAbsRe,real2 *hzTransformFarXAbsIm,
real2 *exTransformNearZScaRe,real2 *exTransformNearZScaIm,real2 *eyTransformNearZScaRe,real2 *eyTransformNearZScaIm,real2 *hxTransformNearZScaRe,real2 *hxTransformNearZScaIm,real2 *hyTransformNearZScaRe,real2 *hyTransformNearZScaIm,
real2 *exTransformFarZScaRe,real2 *exTransformFarZScaIm,real2 *eyTransformFarZScaRe,real2 *eyTransformFarZScaIm,real2 *hxTransformFarZScaRe,real2 *hxTransformFarZScaIm,real2 *hyTransformFarZScaRe,real2 *hyTransformFarZScaIm,
real2 *exTransformNearYScaRe,real2 *exTransformNearYScaIm,real2 *ezTransformNearYScaRe,real2 *ezTransformNearYScaIm,real2 *hxTransformNearYScaRe,real2 *hxTransformNearYScaIm,real2 *hzTransformNearYScaRe,real2 *hzTransformNearYScaIm,
real2 *exTransformFarYScaRe,real2 *exTransformFarYScaIm,real2 *ezTransformFarYScaRe,real2 *ezTransformFarYScaIm,real2 *hxTransformFarYScaRe,real2 *hxTransformFarYScaIm,real2 *hzTransformFarYScaRe,real2 *hzTransformFarYScaIm,
real2 *eyTransformNearXScaRe,real2 *eyTransformNearXScaIm,real2 *ezTransformNearXScaRe,real2 *ezTransformNearXScaIm,real2 *hyTransformNearXScaRe,real2 *hyTransformNearXScaIm,real2 *hzTransformNearXScaRe,real2 *hzTransformNearXScaIm,
real2 *eyTransformFarXScaRe,real2 *eyTransformFarXScaIm,real2 *ezTransformFarXScaRe,real2 *ezTransformFarXScaIm,real2 *hyTransformFarXScaRe,real2 *hyTransformFarXScaIm,real2 *hzTransformFarXScaRe,real2 *hzTransformFarXScaIm,int NCELLX,int NCELLY,int NCELLZ);






real2 *ExTransformNearZScaRedev;
real2 *EyTransformNearZScaRedev;
real2 *EzTransformNearZScaRedev;
real2 *HxTransformNearZScaRedev;
real2 *HyTransformNearZScaRedev;
real2 *HzTransformNearZScaRedev;

real2 *ExTransformFarZScaRedev;
real2 *EyTransformFarZScaRedev;
real2 *EzTransformFarZScaRedev;
real2 *HxTransformFarZScaRedev;
real2 *HyTransformFarZScaRedev;
real2 *HzTransformFarZScaRedev;

real2 *ExTransformNearYScaRedev;
real2 *EyTransformNearYScaRedev;
real2 *EzTransformNearYScaRedev;
real2 *HxTransformNearYScaRedev;
real2 *HyTransformNearYScaRedev;
real2 *HzTransformNearYScaRedev;

real2 *ExTransformFarYScaRedev;
real2 *EyTransformFarYScaRedev;
real2 *EzTransformFarYScaRedev;
real2 *HxTransformFarYScaRedev;
real2 *HyTransformFarYScaRedev;
real2 *HzTransformFarYScaRedev;

real2 *ExTransformNearXScaRedev;
real2 *EyTransformNearXScaRedev;
real2 *EzTransformNearXScaRedev;
real2 *HxTransformNearXScaRedev;
real2 *HyTransformNearXScaRedev;
real2 *HzTransformNearXScaRedev;

real2 *ExTransformFarXScaRedev;
real2 *EyTransformFarXScaRedev;
real2 *EzTransformFarXScaRedev;
real2 *HxTransformFarXScaRedev;
real2 *HyTransformFarXScaRedev;
real2 *HzTransformFarXScaRedev;



real2 *ExTransformNearZScaImdev;
real2 *EyTransformNearZScaImdev;
real2 *EzTransformNearZScaImdev;
real2 *HxTransformNearZScaImdev;
real2 *HyTransformNearZScaImdev;
real2 *HzTransformNearZScaImdev;

real2 *ExTransformFarZScaImdev;
real2 *EyTransformFarZScaImdev;
real2 *EzTransformFarZScaImdev;
real2 *HxTransformFarZScaImdev;
real2 *HyTransformFarZScaImdev;
real2 *HzTransformFarZScaImdev;

real2 *ExTransformNearYScaImdev;
real2 *EyTransformNearYScaImdev;
real2 *EzTransformNearYScaImdev;
real2 *HxTransformNearYScaImdev;
real2 *HyTransformNearYScaImdev;
real2 *HzTransformNearYScaImdev;

real2 *ExTransformFarYScaImdev;
real2 *EyTransformFarYScaImdev;
real2 *EzTransformFarYScaImdev;
real2 *HxTransformFarYScaImdev;
real2 *HyTransformFarYScaImdev;
real2 *HzTransformFarYScaImdev;

real2 *ExTransformNearXScaImdev;
real2 *EyTransformNearXScaImdev;
real2 *EzTransformNearXScaImdev;
real2 *HxTransformNearXScaImdev;
real2 *HyTransformNearXScaImdev;
real2 *HzTransformNearXScaImdev;

real2 *ExTransformFarXScaImdev;
real2 *EyTransformFarXScaImdev;
real2 *EzTransformFarXScaImdev;
real2 *HxTransformFarXScaImdev;
real2 *HyTransformFarXScaImdev;
real2 *HzTransformFarXScaImdev;


////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////
real2 *ExTransformNearZAbsRedev;
real2 *EyTransformNearZAbsRedev;
real2 *EzTransformNearZAbsRedev;
real2 *HxTransformNearZAbsRedev;
real2 *HyTransformNearZAbsRedev;
real2 *HzTransformNearZAbsRedev;

real2 *ExTransformFarZAbsRedev;
real2 *EyTransformFarZAbsRedev;
real2 *EzTransformFarZAbsRedev;
real2 *HxTransformFarZAbsRedev;
real2 *HyTransformFarZAbsRedev;
real2 *HzTransformFarZAbsRedev;

real2 *ExTransformNearYAbsRedev;
real2 *EyTransformNearYAbsRedev;
real2 *EzTransformNearYAbsRedev;
real2 *HxTransformNearYAbsRedev;
real2 *HyTransformNearYAbsRedev;
real2 *HzTransformNearYAbsRedev;

real2 *ExTransformFarYAbsRedev;
real2 *EyTransformFarYAbsRedev;
real2 *EzTransformFarYAbsRedev;
real2 *HxTransformFarYAbsRedev;
real2 *HyTransformFarYAbsRedev;
real2 *HzTransformFarYAbsRedev;

real2 *ExTransformNearXAbsRedev;
real2 *EyTransformNearXAbsRedev;
real2 *EzTransformNearXAbsRedev;
real2 *HxTransformNearXAbsRedev;
real2 *HyTransformNearXAbsRedev;
real2 *HzTransformNearXAbsRedev;

real2 *ExTransformFarXAbsRedev;
real2 *EyTransformFarXAbsRedev;
real2 *EzTransformFarXAbsRedev;
real2 *HxTransformFarXAbsRedev;
real2 *HyTransformFarXAbsRedev;
real2 *HzTransformFarXAbsRedev;


real2 *ExTransformNearZAbsImdev;
real2 *EyTransformNearZAbsImdev;
real2 *EzTransformNearZAbsImdev;
real2 *HxTransformNearZAbsImdev;
real2 *HyTransformNearZAbsImdev;
real2 *HzTransformNearZAbsImdev;

real2 *ExTransformFarZAbsImdev;
real2 *EyTransformFarZAbsImdev;
real2 *EzTransformFarZAbsImdev;
real2 *HxTransformFarZAbsImdev;
real2 *HyTransformFarZAbsImdev;
real2 *HzTransformFarZAbsImdev;

real2 *ExTransformNearYAbsImdev;
real2 *EyTransformNearYAbsImdev;
real2 *EzTransformNearYAbsImdev;
real2 *HxTransformNearYAbsImdev;
real2 *HyTransformNearYAbsImdev;
real2 *HzTransformNearYAbsImdev;

real2 *ExTransformFarYAbsImdev;
real2 *EyTransformFarYAbsImdev;
real2 *EzTransformFarYAbsImdev;
real2 *HxTransformFarYAbsImdev;
real2 *HyTransformFarYAbsImdev;
real2 *HzTransformFarYAbsImdev;

real2 *ExTransformNearXAbsImdev;
real2 *EyTransformNearXAbsImdev;
real2 *EzTransformNearXAbsImdev;
real2 *HxTransformNearXAbsImdev;
real2 *HyTransformNearXAbsImdev;
real2 *HzTransformNearXAbsImdev;

real2 *ExTransformFarXAbsImdev;
real2 *EyTransformFarXAbsImdev;
real2 *EzTransformFarXAbsImdev;
real2 *HxTransformFarXAbsImdev;
real2 *HyTransformFarXAbsImdev;
real2 *HzTransformFarXAbsImdev;





















void CalculateAbsScatt(void);
__global__ void UpdateHydroPx(real *ex,real *ex_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
	real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me);
	__global__ void UpdateHydroPy(real *ex,real *ex_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
		real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me);
		__global__ void UpdateHydroPz(real *ex,real *ex_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
			real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me);

int Hydrodynamics;
// int Hydrodynamics;
int WithMagField;
int WithConvection;

real N_EQ;

#endif // _GLOBAL_VAR_H
