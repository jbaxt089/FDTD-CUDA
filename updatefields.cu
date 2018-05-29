#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "extern_var.h"
#include "cuda_profiler_api.h"
#include <cuda.h>
// #include <cuda_runtime.h>
// #include <device_launch_parameters.h>
//#include<conio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
#include <cudaProfiler.h>
//  __device__  int ThreeDMapD(int i,int j,int k,int SizeZ,int SizeY){
//   int num = k + SizeZ*j +SizeY*SizeZ*i;
//   return num;
// }
//
//
//  __device__  int FourDMapD(int i,int j,int k,int n,int SizeN,int SizeZ,int SizeY){
//   int num = n + SizeN*( k + SizeZ*j +SizeY*SizeZ*i);
//   return num;
// }
//
//  __device__  int TwoDMapD(int i,int j,int size){
//   int num = j + i*size;
//   return num;
// }


//update B-field
void UPDATE_B(){
	// if(TEz && polar_psi==0){
	// 	UPDATE_hx();
	// 	UPDATE_hz();
	// }
	// else if(TMz && polar_psi==0){
	// 	UPDATE_hy();
	// }
	//else{
  int Number;
	int threadsPerBlock = 256;
  Number = NCELLX * NCELLY *NCELLZ;
	int blocksPerGrid = Number/threadsPerBlock + 1;
  // cudaProfilerStart() ;
        //
		    // UPDATE_hx <<<blocksPerGrid, threadsPerBlock>>> (hx,ez,ey,Chxh,Chxe,psi_Hx_z_N,psi_Hx_z_F,psi_Hx_y_N,psi_Hx_y_F,khdy,khdz,bh_z_N,bh_z_F,ch_z_N,ch_z_F,bh_y_N,bh_y_F,ch_y_N,ch_y_F,NCELLX,NCELLY,NCELLZ,Periodic_XY,dx,dy,dz,dt,cpml_N_Z,cpml_F_Z,cpml_N_Y,cpml_F_Y,cpml_z_lim,cpml_y_lim,cpml_x_lim,NcpmlZ,NcpmlY);
        // UPDATE_hy <<<blocksPerGrid, threadsPerBlock>>> (hy,ez,ex,Chyh,Chye,psi_Hy_z_N,psi_Hy_z_F,psi_Hy_x_N,psi_Hy_x_F,khdx,khdz,bh_z_N,bh_z_F,ch_z_N,ch_z_F,bh_x_N,bh_x_F,ch_x_N,ch_x_F,NCELLX,NCELLY,NCELLZ,Periodic_XY,dx,dy,dz,dt,cpml_N_Z,cpml_F_Z,cpml_N_X,cpml_F_X,cpml_z_lim,cpml_y_lim,cpml_x_lim,NcpmlZ,NcpmlX);
        // UPDATE_hz <<<blocksPerGrid, threadsPerBlock>>> (hz,ey,ex,Chzh,Chze,psi_Hz_x_N,psi_Hz_x_F,psi_Hz_y_N,psi_Hz_y_F,khdx,khdy,bh_x_N,bh_x_F,ch_x_N,ch_x_F,bh_y_N,bh_y_F,ch_y_N,ch_y_F,NCELLX,NCELLY,NCELLZ,Periodic_XY,dx,dy,dz,dt,cpml_N_X,cpml_F_X,cpml_N_Y,cpml_F_Y,cpml_z_lim,cpml_y_lim,cpml_x_lim,NcpmlY,NcpmlX);

        cudaDeviceSynchronize();
    // cudaProfilerStop() ;

	//			UPDATE_hz();


//	}

}
//update E-field
void UPDATE_E(){

		int i,j,k,n;
comp hold;
int Number;
int threadsPerBlock = 256;
Number = NCELLX * NCELLY *NCELLZ;
int blocksPerGrid = Number/threadsPerBlock + 1;
     //
			// // UPDATE_ex();
     //  UPDATE_ex <<<blocksPerGrid, threadsPerBlock>>> (ex,ex_n,ex_n_1,hy,hz,Cexe,Cexh,kedy,kedz,mat_matrix,mat_matrixX,first_medium_max,psi_Ex_z_N,psi_Ex_z_F,psi_Ex_y_N,psi_Ex_y_F,Px_cp,Px_cp_n,Px_cp_n_1,Px_d,Px_d_n,Px_d_n_1,
     // C_1_cp,C_2_cp,C_3_cp,C_4_cp,C_5_cp,d_1_d,d_2_d,d_3_d,d_4_d,d_5_d,d_NL,C_E,z0,N_CP_poles,N_drude_poles,ce_z_N,ce_z_F,be_z_N,be_z_F,ce_y_N,ce_y_F,be_y_N,be_y_F,dx,dy,dz,dt,NCELLX,NCELLY,NCELLZ,
     // Hydrodynamics,cpml_x_lim,cpml_y_lim,cpml_z_lim,cpml_N_Y,cpml_F_Y,cpml_N_Z,cpml_F_Z,NcpmlY,NcpmlZ, C_E_1,C_E_2,Periodic_XY);
			// //UPDATE_ey();
     //  UPDATE_ey <<<blocksPerGrid, threadsPerBlock>>> (ey,ey_n,ey_n_1,hx,hz,Ceye,Ceyh,kedx,kedz,mat_matrix,mat_matrixY,first_medium_max,psi_Ey_z_N,psi_Ey_z_F,psi_Ey_x_N,psi_Ey_x_F,Py_cp,Py_cp_n,Py_cp_n_1,Py_d,Py_d_n,Py_d_n_1,
     // C_1_cp,C_2_cp,C_3_cp,C_4_cp,C_5_cp,d_1_d,d_2_d,d_3_d,d_4_d,d_5_d,d_NL,C_E,z0,N_CP_poles,N_drude_poles,ce_z_N,ce_z_F,be_z_N,be_z_F,ce_y_N,ce_x_F,be_x_N,be_x_F,dx,dy,dz,dt,NCELLX,NCELLY,NCELLZ,
     // Hydrodynamics,cpml_x_lim,cpml_y_lim,cpml_z_lim,cpml_N_X,cpml_F_X,cpml_N_Z,cpml_F_Z,NcpmlX,NcpmlZ, C_E_1,C_E_2,Periodic_XY);
     //
			// // UPDATE_ez();
     //  UPDATE_ez <<<blocksPerGrid, threadsPerBlock>>> (ez,ez_n,ez_n_1,hx,hy,Ceze,Cezh,kedx,kedy,mat_matrix,mat_matrixZ,first_medium_max,psi_Ez_y_N,psi_Ez_y_F,psi_Ez_x_N,psi_Ez_x_F,Pz_cp,Pz_cp_n,Pz_cp_n_1,Pz_d,Pz_d_n,Pz_d_n_1,
     // C_1_cp,C_2_cp,C_3_cp,C_4_cp,C_5_cp,d_1_d,d_2_d,d_3_d,d_4_d,d_5_d,d_NL,C_E,z0,N_CP_poles,N_drude_poles,ce_y_N,ce_y_F,be_y_N,be_y_F,ce_x_N,ce_x_F,be_x_N,be_x_F,dx,dy,dz,dt,NCELLX,NCELLY,NCELLZ,
     // Hydrodynamics,cpml_x_lim,cpml_y_lim,cpml_z_lim,cpml_N_X,cpml_F_X,cpml_N_Y,cpml_F_Y,NcpmlX,NcpmlY, C_E_1,C_E_2,Periodic_XY);
     //  cudaDeviceSynchronize();

// for(i=0;i<NCELLX;i++){
// 	for(j=0;j<NCELLY;j++){
// 	printf("%e\t",Px_d[i][j][dispersive_slab+20][0]);
// 	//printf("%e\t",hx[i][j][dispersive_slab+20]);
// 	}
// 	printf("\n");
// }

		////#pragma omp parallel for collapse(3) private(i,j,k,n) //// schedule(guided)
		for(i=0;i<NCELLX;i++){
			for(j=0;j<NCELLY;j++){
				for(k=0;k<NCELLZ;k++){
					if(mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] < 6 || mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] <6 || mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] <6 || mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] <6){

								for(n=0;n<N_drude_poles;n++){
									Px_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
									Px_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
									Py_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
									Py_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
									Pz_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
									Pz_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];

								}
								for(n=0;n<N_CP_poles;n++){
									Px_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
									Px_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
									Py_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
									Py_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
									Pz_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
									Pz_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
								}

								if(Hydrodynamics){
									hxPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									hyPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									hzPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									hold = NDx[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									NDx[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = NDx_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									NDx_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=hold;
									hold = NDy[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									NDy[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = NDy_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									NDy_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=hold;
									hold = NDz[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									NDz[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = NDz_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
									NDz_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=hold;
								}
							}
					}
			}
		}
	}

  //
  // __global__ void UPDATE_ex(real *ex,real *ex_n,real *ex_n_1,real *hy,real *hz,real *Cexe,real *Cexh,real *kedy,real *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real *psi_Ex_z_N,
  // real *psi_Ex_z_F,real *psi_Ex_y_N,real *psi_Ex_y_F,real *Px_cp,real *Px_cp_n,real *Px_cp_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
  // real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,
  // real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY){
  //     int i,j,k,n,k2,j2;
  //     comp Curl_H, Div_Grad=0.0,J_T,dummy_var;
  // 		comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2;
  // 			    comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
  // 					double INV_DX = 1.0/dx;
  // 					double INV_DY = 1.0/dy;
  // 					double INV_DZ = 1.0/dz;
  //
  //           int idx = blockDim.x * blockIdx.x + threadIdx.x;
  //
  //           i = idx / (NCELLZ*NCELLY);
  //           j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
  //           k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
  //
  //     if(Periodic_XY){
  // 			////#pragma omp parallel for collapse(3) private(Curl_H,i,j,j2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
  // 			// for(k=1;k<NCELLZ-1;k++){
  // 			// 	for(i=0;i<NCELLX;i++){
  // 		  //       for(j=0;j<NCELLY;j++){
  //               if(i<NCELLX && j<NCELLY && k>0 && k<(NCELLZ-1)){
  // 	                //for(k=1;k<NCELLZ-1;k++){
  // 	                	if(j==0){
  // 											#ifdef DOUBLECOMPLEX
  // 												Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)]*cexp(I*k_y*period_y))/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
  // 											#endif
  // 											#ifdef DOUBLEPRECISION
  // 												Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)])/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
  // 											#endif
  //
  // 	                	}
  // 	                	else{
  // 	                		Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
  // 	                	}
  //
  // 										if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){
  //
  //
  // 											//  Div_Grad = Calc_DIV_GRADx(i,j,k);
  //                       Div_Grad = 0.0;
  //
  //                       C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //
  //                       for(n=0;n<N_drude_poles;n++){
  //                           C_P_1+=(d_1_d[n]-1)*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                           C_P_3+=(d_2_d[n])*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                          C_P_NL += d_NL[n]*Div_Grad;
  //                       }
  //                       for(n=0;n<N_CP_poles;n++){
  //                           C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                           C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                       }
  //                       ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                       ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                       ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4 - C_P_NL);
  //
  //
  // 												//printf("%e\n",Div_Grad);
  // 												//Z-CPML
  // 												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 													//Near-Z-PML
  // 													psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
  // 												}
  // 												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 													k2 = k - cpml_F_Z ;
  // 													psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=(1/C_E)*dt*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 												}
  //
  // 											 	for(n=0;n<N_CP_poles;n++){
  // 											 							Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											 	}
  //
  // 												for(n=0;n<N_drude_poles;n++){
  // 													Px_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
  //
  // 												}
  //
  // 										}
  //
  // 										else{
  // 												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Cexe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
  // 												//Z-CPML
  // 												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 													//Near-Z-PML
  // 													psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												}
  // 												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 													k2 = k - cpml_F_Z ;
  // 													psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
  // 												}
  // 										}
  //                   }
  //
  //
  //
  // 	  //       }
  // 	  //   }
  //     // }
  // 	}
  //     //No PBCs
  //     else{
  // 			////#pragma omp target teams distribute parallel for collapse(3) schedule(static,1) private(Curl_H,i,j,j2,k,k2,n,dummy_var,J_T,Div_Grad)
  // 		//	//#pragma omp parallel for collapse(3) private(Curl_H,i,j,j2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
  // 	    // for(i=0;i<NCELLX-1;i++){
  // 	    //     for(j=1;j<NCELLY-1;j++){
  // 	    //             for(k=1;k<NCELLZ-1;k++){
  //                     if(i<(NCELLX-1) && j>0 && j<(NCELLY-1) && k>0 && k<(NCELLZ-1)){
  // 	                    Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
  //
  // 											if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){
  //
  //
  // 													if(Hydrodynamics == 0)
  // 													{
  // 														//Div_Grad = Calc_DIV_GRADx(i,j,k);
  //                             Div_Grad = 0.0;
  // 														// CP_D_ex(i,j,k,Curl_H,Div_Grad);
  //
  //                              C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //
  //                              for(n=0;n<N_drude_poles;n++){
  //                                  C_P_1+=(d_1_d[n]-1)*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                                  C_P_3+=(d_2_d[n])*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                                 C_P_NL += d_NL[n]*Div_Grad;
  //                              }
  //                              for(n=0;n<N_CP_poles;n++){
  //                                  C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                                  C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                              }
  //                              ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                              ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                              ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4 - C_P_NL);
  //
  //
  //
  //
  // 														for(n=0;n<N_CP_poles;n++){
  // 																				Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 														}
  //
  // 														for(n=0;n<N_drude_poles;n++){
  // 																				Px_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
  //
  // 														}
  // 													}
  // 													// else{
  //                           //
  // 													// 	Vx1 = Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)];
  // 													// 	Vx2 = Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)];
  // 													// 	Nx1 = NDx_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)];
  // 													// 	Nx2 = NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)];
  //                           //
  // 													// 	Vy1 = 0.5*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 													// 	Vy2 = 0.5*(Py_d_n[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 													// 	Ny1 = 0.5*(NDy_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 													// 	Ny2 = 0.5*(NDy_prev[ThreeDMapD(i+1,j-1,k,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);
  //                           //
  // 													// 	Vz1 = 0.5*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 													// 	Vz2 = 0.5*(Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 													// 	Nz1 = 0.5*(NDz_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 													// 	Nz2 = 0.5*(NDz_prev[ThreeDMapD(i+1,j,k-1,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);
  //                           //
  // 													// 	NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] - 2.0*dt*INV_DX*(0.5*(Nx1*Vx1-Nx2*Vx2) + (Ny1*Vy1-Ny2*Vy2) + (Nz1*Vz1-Nz2*Vz2) + (0.5*(Vx1-Vx2) + (Vy1-Vy2) + (Vz1-Vz2))*N_EQ);
  //                           //
  //                           //
  // 												  //   C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //                           //
  //                           //
  // 												  //   for(n=0;n<N_CP_poles;n++){
  // 												  //       C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  // 												  //       C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  // 												  //   }
  // 												  //   ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												  //   ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												  //   ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0 + dt*Px_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*e0*(N_EQ + NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]) -C_E_2*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_2-C_P_4);
  //                           //
  // 													// 	for(n=0;n<N_CP_poles;n++){
  // 													// 							Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 													// 	}
  //                           //
  //                           //
  // 													// }
  //
  // 												//	printf("%e\n",Div_Grad);
  // 													//Z-CPML
  // 													// if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 													// 	//Near-Z-PML
  // 													// 	psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													// 	ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 													// }
  // 													// if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 													// 	k2 = k - cpml_F_Z ;
  // 													// 	psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													// 	ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=(1/C_E)*dt*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 													// }
  //
  //
  //
  // 											}
  //
  // 	                    else{
  // 	                        ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Cexe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
  // 	                    }
  // 											//Z-CPML
  // 											if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 												//Near-Z-PML
  // 												psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
  // 											}
  // 											if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){ //Far Z PML
  // 												k2 = k - cpml_F_Z;
  // 												psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 												// if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
  // 												// 	ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=(1/C_E)*dt*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 												// }
  // 												// else{
  // 													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
  // 												//}
  // 											}
  // 											//Y PML
  // 											if(j<cpml_N_Y+1 && i<cpml_x_lim && k<cpml_z_lim){ //Near Y PML
  // 												psi_Ex_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]=be_y_N[j]*psi_Ex_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]+ce_y_N[j]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
  // 												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)];
  // 											}
  // 											if(j>=cpml_F_Y && i<cpml_x_lim && k<cpml_z_lim){ //Far Y PML
  // 												j2 = j - cpml_F_Y;
  // 												psi_Ex_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]=be_y_F[j2]*psi_Ex_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]+ce_y_F[j2]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
  // 												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)];
  // 											}
  // 									// 		if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
  // 									// 				for(n=0;n<N_CP_poles;n++){
  // 									// 										Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 									// 				}
  // 									// 				for(n=0;n<N_drude_poles;n++){
  // 									// 										Px_d[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=d_1_d[n]*Px_d_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+d_2_d[n]*Px_d_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+d_3_d[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 									// 				}
  // 	                // }
  //                 }
  // 	//         }
  // 	//     }
  // 	// }
  // }
  //
  //   return;
  // }
  //
  //
  //
  //
  //
  //
  // __global__ void UPDATE_ey(real *ey,real *ey_n,real *ey_n_1,real *hx,real *hz,real *Ceye,real *Ceyh,real *kedx,real *kedz,int *mat_matrix,int *mat_matrixY,int first_medium_max,real *psi_Ey_z_N,
  // real *psi_Ey_z_F,real *psi_Ey_x_N,real *psi_Ey_x_F,real *Py_cp,real *Py_cp_n,real *Py_cp_n_1,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
  // real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_x_N,real *ce_x_F,real *be_x_N,real *be_x_F,
  // real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_X,int cpml_F_X,int cpml_N_Z,int cpml_F_Z,int NcpmlX,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY){
  //     int i,j,k,i2,k2,n;
  //     comp Curl_H,Div_Grad=0.0,dummy_var,J_T;
  // 		comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
  // 		double INV_DX = 1.0/dx;
  // 		double INV_DY = 1.0/dy;
  // 		double INV_DZ = 1.0/dz;
  // 		comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2;
  //     int idx = blockDim.x * blockIdx.x + threadIdx.x;
  //
  //     i = idx / (NCELLZ*NCELLY);
  //     j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
  //     k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
  //
  //     if(Periodic_XY){
  //     	////#pragma omp parallel for collapse(3) private(Curl_H,i,j,i2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
  // 			// for(k=0;k<NCELLZ-1;k++){
  // 			//   for(i=0;i<NCELLX;i++){
  // 		  //       for(j=0;j<NCELLY;j++){
  //               if(i<NCELLX && j<NCELLY && k>0 && k<(NCELLZ-1)){
  // 	               // for(k=1;k<NCELLZ-1;k++){
  // 	                	if(i==0){
  // 											#ifdef DOUBLECOMPLEX
  // 	                		Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]*cexp(I*period_x*k_x))/kedx[i];
  // 											#endif
  // 											#ifdef DOUBLEPRECISION
  // 											Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)])/kedx[i];
  // 											#endif
  // 											  //printf("%d,%d,%d \t %f\t%f\t%f\n",i,j,k,creal(Curl_H),creal(hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]), creal(hz[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]));
  //
  // 	                	}
  //
  // 	                	else{
  // 	                		Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i];
  // 	                	}
  // 										if(mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] <6){
  //
  // 												//Div_Grad = Calc_DIV_GRADy(i,j,k);
  //                         Div_Grad = 0.0;
  // 											// printf("%e\n",d_NL[0]*Div_Grad);
  // 												// CP_D_ey(i,j,k,Curl_H,Div_Grad);
  //                         C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //                         //printf("here");
  //                         for(n=0;n<N_drude_poles;n++){
  //                             C_P_1+=(d_1_d[n]-1.0)*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                             C_P_3+=(d_2_d[n])*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                             C_P_NL += d_NL[n]*Div_Grad;
  //                         }
  //                         for(n=0;n<N_CP_poles;n++){
  //                             C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                             C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                         }
  //                         ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                         ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                         ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);
  //
  //
  // 											//	printf("%e\n",ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 												//Z-CPML
  // 												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 													//Here we are in the near Z-PML
  // 													psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
  // 												}
  // 												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 													//Here we are in the far Z-PML
  // 														k2 = k - cpml_F_Z ;
  // 														psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 															ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 												}
  //
  // 														for(n=0;n<N_CP_poles;n++){
  // 																				Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 														}
  //
  // 														//if(Hydrodynamics == 0){
  // 														for(n=0;n<N_drude_poles;n++){
  // 															Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
  //
  // 														}
  // 													//}
  //
  // 										}
  //
  // 										else{
  // 												ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
  // 												//Z-CPML
  // 												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 													//Here we are in the near Z-PML
  // 													psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
  // 												}
  // 												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 													//Here we are in the far Z-PML
  // 														k2 = k - cpml_F_Z ;
  // 														psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 														if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
  // 															ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 														}
  // 														else{
  // 															ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
  // 														}
  // 												}
  // 										}
  //
  //                   }
  //
  // 	      //       }
  // 	      //   }
  //    	    // }
  //     }
  //
  //     else{
  // 		//	//#pragma omp parallel for collapse(3) private(Curl_H,i,j,i2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
  // 	    // for(i=1;i<NCELLX-1;i++){
  // 	    //     for(j=0;j<NCELLY-1;j++){
  // 	    //             for(k=1;k<NCELLZ-1;k++){
  //                     if(i>0 && i<(NCELLX-1) && j<(NCELLY-1) && k>0 && k<(NCELLZ-1)){
  // 	                    Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i];
  // 											if(mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] <6){
  //
  // 												if(Hydrodynamics == 0){
  // 													//Div_Grad = Calc_DIV_GRADy(i,j,k);
  //                           Div_Grad= 0.0;
  // 												// printf("%e\n",d_NL[0]*Div_Grad);
  // 													// CP_D_ey(i,j,k,Curl_H,Div_Grad);
  //                           C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //                       		//printf("here");
  //                           for(n=0;n<N_drude_poles;n++){
  //                               C_P_1+=(d_1_d[n]-1.0)*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                               C_P_3+=(d_2_d[n])*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                       				C_P_NL += d_NL[n]*Div_Grad;
  //                           }
  //                           for(n=0;n<N_CP_poles;n++){
  //                               C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                               C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                           }
  //                           ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                           ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                           ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);
  //
  // 												//	printf("%e\n",ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 													//Z-CPML
  // 													// if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 													// 	//Here we are in the near Z-PML
  // 													// 	psi_Ey_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													// 	ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 													// }
  // 													// if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 													// 	//Here we are in the far Z-PML
  // 													// 		k2 = k - cpml_F_Z ;
  // 													// 		psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 													// 			ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 													// }
  //
  // 															for(n=0;n<N_CP_poles;n++){
  // 																					Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 															}
  //
  // 															for(n=0;n<N_drude_poles;n++){
  // 																Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
  //
  // 															}
  //
  // 													}
  // 											// 		else{
  //                       //
  // 											// 			Vy1 = Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)];
  // 											// 			Vy2 = Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)];
  // 											// 			Ny1 = NDy_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)];
  // 											// 			Ny2 = NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)];
  //                       //
  // 											// 			Vx1 = 0.5*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 											// 			Vx2 = 0.5*(Px_d_n[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 											// 			Nx1 = 0.5*(NDx_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 											// 			Nx2 = 0.5*(NDx_prev[ThreeDMapD(i-1,j+1,k,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);
  //                       //
  // 											// 			Vz1 = 0.5*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 											// 			Vz2 = 0.5*(Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 											// 			Nz1 = 0.5*(NDz_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 											// 			Nz2 = 0.5*(NDz_prev[ThreeDMapD(i,j+1,k-1,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);
  //                       //
  // 											// 			NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] - 2.0*dt*INV_DX*((Nx1*Vx1-Nx2*Vx2) + 0.5*(Ny1*Vy1-Ny2*Vy2) + (Nz1*Vz1-Nz2*Vz2) + ((Vx1-Vx2) + 0.5*(Vy1-Vy2) + (Vz1-Vz2))*N_EQ);
  //                       //
  //                       //
  // 											//     C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //                       //
  // 											//     for(n=0;n<N_CP_poles;n++){
  // 											//         C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  // 											//         C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  // 											//     }
  // 											//     ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											//     ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											//     ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0 + dt*Py_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*e0*(N_EQ + NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]) + C_E_1*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_2-C_P_4);
  // 											// 		for(n=0;n<N_CP_poles;n++){
  // 											// 								Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											// 		}
  //                       //
  // 											// }
  // 										}
  //
  //
  // 										else{
  // 												ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
  //
  // 										}
  //
  // 										if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
  // 											//Here we are in the near Z-PML
  // 											psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 											ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
  // 										}
  // 										if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
  // 												//Here we are in the far Z-PML
  // 												k2 = k - cpml_F_Z;
  // 												psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
  // 												// if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
  // 												// 	ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
  // 												// }
  // 												// else{
  // 													ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
  // 												//}
  // 										}
  // 										//X-CPML
  // 										if(i<cpml_N_X+1 && j<cpml_y_lim && k<cpml_z_lim){
  // 											//Here we are in the near-X-PML
  // 											psi_Ey_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_x_N[i]*psi_Ey_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_x_N[i]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
  // 											ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 										}
  // 										if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){
  // 											//Here we are in the far-X-PML
  // 											i2 = i - cpml_F_X;
  // 											psi_Ey_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=be_x_F[i2]*psi_Ey_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ce_x_F[i2]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
  // 											ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
  // 										}
  //
  //
  //                   }
  //
  // 								}
  // 	  //       }
  //     //
  // 	  //   }
  //     // }
  //
  //   return;
  // }
  //
  //
  // //void UPDATE_ez(void){
  //   __global__ void UPDATE_ez(real *ez,real *ez_n,real *ez_n_1,real *hx,real *hy,real *Ceze,real *Cezh,real *kedx,real *kedy,int *mat_matrix,int *mat_matrixZ,int first_medium_max,real *psi_Ez_y_N,
  //   real *psi_Ez_y_F,real *psi_Ez_x_N,real *psi_Ez_x_F,real *Pz_cp,real *Pz_cp_n,real *Pz_cp_n_1,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
  //   real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,real *ce_x_N,real *ce_x_F,real *be_x_N,real *be_x_F,
  //   real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_X,int cpml_F_X,int cpml_N_Y,int cpml_F_Y,int NcpmlX,int NcpmlY,real C_E_1,real C_E_2,int Periodic_XY){
  //
  //     int i,j,k,i2,j2,n;
  //     comp Curl_H,Div_Grad=0.0,dummy_var,J_T;
  // 		comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
  // 		comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2;
  //
  // 		double INV_DX = 1.0/dx;
  // 		double INV_DY = 1.0/dy;
  // 		double INV_DZ = 1.0/dz;
  //     int idx = blockDim.x * blockIdx.x + threadIdx.x;
  //
  //     i = idx / (NCELLZ*NCELLY);
  //     j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
  //     k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
  //
  //     if(Periodic_XY){
  // 	//	//#pragma omp parallel for collapse(3) private(Curl_H,i,j,k,i2,j2,n,dummy_var,J_T,Div_Grad) // schedule(static)
  // 		// for(k=0;k<NCELLZ-1;k++){
  // 		// for(i=0;i<NCELLX;i++){
  // 		//         for(j=0;j<NCELLY;j++){
  //               if(i<NCELLX,j<NCELLY,k<NCELLZ){
  // 		           //     for(k=0;k<NCELLZ-1;k++){
  // 											//printf("Thread %d, ready to work\n",omp_get_thread_num());
  //
  // 		                	if(i==0 || j==0){
  // 		                		if(i==0 && j==0){
  // 													#ifdef DOUBLECOMPLEX
  // 				        						Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]*cexp(I*k_x*period_x))/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)]*cexp(I*k_y*period_y))/kedy[j];
  // 														#endif
  // 														#ifdef DOUBLEPRECISION
  // 														Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)])/kedy[j];
  // 														#endif
  //
  //
  // 			               	    }
  // 			               	    else if(i==0){
  // 													#ifdef DOUBLECOMPLEX
  // 				                 	Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]*cexp(I*k_x*period_x))/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
  // 													#endif
  // 													#ifdef DOUBLEPRECISION
  // 													Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
  // 													#endif
  //
  //
  // 				                }
  // 				                else{
  // 													#ifdef DOUBLECOMPLEX
  // 				                	Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)]*cexp(I*k_y*period_y))/kedy[j];
  // 													#endif
  // 													#ifdef DOUBLEPRECISION
  // 													Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)])/kedy[j];
  // 													#endif
  //
  //
  // 				                }
  // 		        			}
  //
  // 		        			else{
  // 		        				Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
  //
  // 		        			}
  // 									if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]<6){
  //
  // 											// Div_Grad  = Calc_DIV_GRADz(i,j,k);
  //                       Div_Grad = 0.0;
  //
  // 											// CP_D_ez(i,j,k,Curl_H,Div_Grad);
  //                       C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //
  //                       for(n=0;n<N_drude_poles;n++){
  //                           C_P_1+=(d_1_d[n]-1)*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                           C_P_3+=(d_2_d[n])*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                           C_P_NL += d_NL[n]*Div_Grad;
  //                       }
  //                       for(n=0;n<N_CP_poles;n++){
  //                           C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                           C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                       }
  //                       ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                       ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                       ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);
  //
  //
  // 												for(n=0;n<N_CP_poles;n++){
  // 																		Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												}
  //
  // 												for(n=0;n<N_drude_poles;n++){
  // 													Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
  // 												}
  // 									}
  //
  // 									else{
  // 											ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
  // 									}
  //
  //                 }
  // 		      //     }
  // 		      //   }
  // 					// }
  //     }
  //
  //     else{
  // 	//		//#pragma omp parallel for collapse(3) private(Curl_H,i,j,k,i2,j2,n,dummy_var,J_T,Div_Grad) // schedule(static)
  // 	    // for(i=1;i<NCELLX-1;i++){
  // 	    //     for(j=1;j<NCELLY-1;j++){
  // 	    //             for(k=0;k<NCELLZ-1;k++){
  //                     if(i>0 && i<(NCELLX-1) && j>0 && j<(NCELLY-1) && k<(NCELLZ-1)){
  // 	                    Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
  // 											if(mat_matrixZ[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixZ[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]<6){
  //
  // 												if(Hydrodynamics == 0){
  // 			 									 //Div_Grad  = Calc_DIV_GRADz(i,j,k);
  //                          Div_Grad = 0.0;
  // 			 									// CP_D_ez(i,j,k,Curl_H,Div_Grad);
  //                         C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //
  //                         for(n=0;n<N_drude_poles;n++){
  //                             C_P_1+=(d_1_d[n]-1)*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                             C_P_3+=(d_2_d[n])*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  //                             C_P_NL += d_NL[n]*Div_Grad;
  //                         }
  //                         for(n=0;n<N_CP_poles;n++){
  //                             C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                             C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  //                         }
  //                         ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                         ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  //                         ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);
  //
  //
  // 			 										 for(n=0;n<N_CP_poles;n++){
  // 			 																 Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 			 										 }
  // 			 										 for(n=0;n<N_drude_poles;n++){
  // 			 											 Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
  // 			 										 }
  //
  // 											 }
  // 											 // else{
  //                        //
  // 												//  C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  //                        //
  // 												//  Vz1 = Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)];
  // 												//  Vz2 = Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)];
  // 												//  Nz1 = NDz_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)];
  // 												//  Nz2 = NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)];
  //                        //
  // 												//  Vx1 = 0.5*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 												//  Vx2 = 0.5*(Px_d_n[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 												//  Nx1 = 0.5*(NDx_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 												//  Nx2 = 0.5*(NDx_prev[ThreeDMapD(i-1,j,k+1,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);
  //                        //
  // 												//  Vy1 = 0.5*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 												//  Vy2 = 0.5*(Py_d_n[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 												//  Ny1 = 0.5*(NDy_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
  // 												//  Ny2 = 0.5*(NDy_prev[ThreeDMapD(i,j-1,k+1,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);
  //                        //
  // 												//  NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] - 2.0*dt*INV_DX*((Nx1*Vx1-Nx2*Vx2) + (Ny1*Vy1-Ny2*Vy2) + 0.5*(Vz1-Vz2) + ((Vx1-Vx2) + (Vy1-Vy2) + 0.5*(Vz1-Vz2))*N_EQ);
  //                        //
  //                        //
  // 												//  for(n=0;n<N_CP_poles;n++){
  // 												// 		 C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  // 												// 		 C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
  // 												//  }
  // 												//  ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												//  ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												//  ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0 + dt*Pz_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*e0*(NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + N_EQ)+C_E_1*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_2-C_P_4);
  //                        //
  // 												//  for(n=0;n<N_CP_poles;n++){
  // 												// 						 Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 												//  }
  // 											 // }
  // 			 							 }
  //
  //
  // 	                    else{
  // 	                        ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
  // 	                    }
  // 											//Y CPML
  // 											if(j<cpml_N_Y && i<cpml_x_lim && k<cpml_z_lim){ //Near Y PML
  // 												psi_Ez_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]=be_y_N[j]*psi_Ez_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]+ce_y_N[j]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
  // 												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)];
  // 											}
  // 											if(j>=cpml_F_Y && i<cpml_x_lim &&  k<cpml_z_lim){
  // 												j2 = j - cpml_F_Y;
  // 												psi_Ez_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]=be_y_F[j2]*psi_Ez_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]+ce_y_F[j2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
  // 												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)];
  // 											}
  // 											//X PML
  // 											if(i<cpml_N_X+1 && j<cpml_y_lim && k<cpml_z_lim){//Near X-PML
  // 												psi_Ez_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_x_N[i]*psi_Ez_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_x_N[i]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
  // 												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											}
  // 											if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){//far X-PML
  // 												i2 = i - cpml_F_X;
  // 												psi_Ez_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=be_x_F[i2]*psi_Ez_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ce_x_F[i2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
  // 												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
  // 											}
  // 											// if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
  // 											// 		 for(n=0;n<N_CP_poles;n++){
  // 											// 				Pz_cp[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											// 				}
  // 											// 		for(n=0;n<N_drude_poles;n++){
  // 											// 				Pz_d[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Pz_d_n[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Pz_d_n_1[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  // 											// 				}
  // 											// }
  //
  // 	                }
  //       //           }
  // 	    //     }
  //       //
  // 	    // }
  //     }
  //
  //   return;
  // }
  //





void CP_D_ex(int i,int j, int k, comp Curl_H,comp Div_Grad){
      int n;
      comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
      C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

      for(n=0;n<N_drude_poles;n++){
          C_P_1+=(d_1_d[n]-1)*Px_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
          C_P_3+=(d_2_d[n])*Px_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  				C_P_NL += d_NL[n]*Div_Grad;
      }
      for(n=0;n<N_CP_poles;n++){
          C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
          C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
      }
      ex_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
      ex_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
      ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ex_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-C_E_2*ex_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4 - C_P_NL);
  }

  void CP_D_ey(int i,int j, int k, comp Curl_H,comp Div_Grad ){
      int n;
      comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
      C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
  		//printf("here");
      for(n=0;n<N_drude_poles;n++){
          C_P_1+=(d_1_d[n]-1.0)*Py_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
          C_P_3+=(d_2_d[n])*Py_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  				C_P_NL += d_NL[n]*Div_Grad;
      }
      for(n=0;n<N_CP_poles;n++){
          C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
          C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
      }
      ey_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
      ey_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
      ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ey_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);

  }

  void CP_D_ez(int i,int j, int k, comp Curl_H, comp Div_Grad){
      int n;
      comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
      C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

      for(n=0;n<N_drude_poles;n++){
          C_P_1+=(d_1_d[n]-1)*Pz_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
          C_P_3+=(d_2_d[n])*Pz_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
  				C_P_NL += d_NL[n]*Div_Grad;
      }
      for(n=0;n<N_CP_poles;n++){
          C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
          C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
      }
      ez_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
      ez_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
      ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ez_n[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);

  }

//
//
//   __global__ void UPDATE_hx(real *hx,real *ez,real *ey,real *Chxh,real *Chxe,real *psi_Hx_z_N,real *psi_Hx_z_F,real *psi_Hx_y_N,real *psi_Hx_y_F,real *khdy,real
//     *khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY){
// //void UPDATE_hx(void){
// // cudaProfilerStart();
//
//       int i,j,k,j2,k2;
//       comp Curl_E;
//       int idx = blockDim.x * blockIdx.x + threadIdx.x;
//
//       i = idx / (NCELLZ*NCELLY);
//       j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
//       k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
//
//       if(Periodic_XY){
//   	////#pragma omp parallel for collapse(3) private(i,j,k,Curl_E,j2,k2) // schedule(static)
//   	// for(k=0;k<NCELLZ;k++){
//   	// 	for(i=0;i<NCELLX;i++){
//   	//         for(j=0;j<NCELLY;j++){
//
//               if(i<NCELLX && j<NCELLY && k<NCELLZ){
//   	           //     for(k=0;k<NCELLZ-1;k++){
//   									 //if(i==1) printf("%d %d %d\n",i,j,k);
//   		                	if(j==NCELLY-1){
//   												#ifdef DOUBLECOMPLEX
//   		                				Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]*cexp(-I*k_y*period_y)-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
//   												#endif
//   												#ifdef DOUBLEPRECISION
//   		                				Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
//   												#endif
//   		                   	    hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chxh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//
//   		                	}
//   		                	 else{
//   														Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
//   		                   	    hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chxh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   		                	}
//   										 //Z-CPML
//   										 if(k<cpml_N_Z && i<cpml_x_lim && j<cpml_y_lim){
//   											 	//Near Z-PML
//   											 		psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   													hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
//   										 }
//   										 if(k>=cpml_F_Z && j<cpml_y_lim && i<cpml_x_lim){
//   											 //Far Z-PML
//   											 		k2 = k - cpml_F_Z;
//   													psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   													hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
//   										 }
//                      }
//   	    //             }
//   	    //     }
//   	    // }
//       }
//
//       else{
//   			//  //#pragma omp target device(0) MapD(Chxe[:NCELLX-1][:NCELLY-1][:NCELLZ-1],Chxh[:NCELLX-1][:NCELLY-1][:NCELLZ-1],ez[:NCELLX-1][:NCELLY-1][:NCELLZ-1],ey[:NCELLX-1][:NCELLY-1][:NCELLZ-1],khdy[:NCELLY-1],khdz[:NCELLZ-1],bh_z_N[:NcpmlZ-1],bh_z_F[:NcpmlZ-1],ch_z_N[:NcpmlZ-1],ch_z_F[:NcpmlZ-1],bh_y_N[:NcpmlY-1],bh_y_F[:NcpmlY-1],ch_y_N[:NcpmlY-1],ch_y_F[:NcpmlY-1]) 		MapD(tofrom:hx[:NCELLX-1][:NCELLY-1][:NCELLZ-1],psi_Hx_z_N[:NCELLX-1][:NCELLY-1][:cpml_N_Z-1],psi_Hx_z_F[:NCELLX-1][:NCELLY-1][:cpml_N_Z-1],psi_Hx_y_N[:NCELLX-1][:cpml_N_Y-1][:NCELLZ-1],psi_Hx_y_F[:NCELLX-1][:cpml_N_Y-1][:NCELLZ-1])
//   			//  {
//   			// //#pragma omp parallel for collapse(3) private(i,j,k,Curl_E,j2,k2) // schedule(static)
//   	    // for(i=0;i<NCELLX;i++){
//   	    //     for(j=0;j<NCELLY-1;j++){
//   	    //             for(k=0;k<NCELLZ-1;k++){
//                       if(i<NCELLX && j<(NCELLY-1) && k<(NCELLZ-1)){
//
//   	                    Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
//   	                    hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chxh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   											//Z-CPML
//   											if(k<cpml_N_Z && i<cpml_x_lim && j<cpml_y_lim){
//   												 //Near Z-PML
//   													 psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   													 hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
//   											}
//   											if(k>=cpml_F_Z && j<cpml_y_lim && i<cpml_x_lim){
//   												//Far Z-PML
//   													 k2 = k - cpml_F_Z;
//   													 psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   													 hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
//   											}
//   										   //Y- PML
//   											if(j<cpml_N_Y && i<cpml_x_lim && j<cpml_y_lim){
//   													 psi_Hx_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]=bh_y_N[j]*psi_Hx_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]+ch_y_N[j]*(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
//   													 hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)];
//   											}
//   											if(j>=cpml_F_Y && i<cpml_x_lim && k<cpml_z_lim){
//   													j2 = j - cpml_F_Y;
//   													psi_Hx_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]=bh_y_F[j2]*psi_Hx_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]+ch_y_F[j2]*(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
//   													hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)];
//   											}
//                       }
//   	  //               }
//   	  //       }
//   	  //  // }
//   		// }
//     }
//     // cudaProfilerStop();
//
//     return;
//   }
//
//   __global__   void UPDATE_hy(real *hy,real *ez,real *ex,real *Chyh,real *Chye,real *psi_Hy_z_N,real *psi_Hy_z_F,real *psi_Hy_x_N,real *psi_Hy_x_F,real *khdx,real
//     *khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_x_N,real *bh_x_F,real *ch_x_N,real *ch_x_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_X,int cpml_F_X,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlX){
//
//       int i,j,k,n,i2,k2;
//       comp Curl_E;
//       int idx = blockDim.x * blockIdx.x + threadIdx.x;
//
//       i = idx / (NCELLZ*NCELLY);
//       j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
//       k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
//       if(Periodic_XY){
//   		////#pragma omp parallel for collapse(3) private(Curl_E,i,i2,j,k,k2) // schedule(static)
//   		// for(k=0;k<NCELLZ;k++){
//   		// for(i=0;i<NCELLX;i++){
//   	  //       for(j=0;j<NCELLY;j++){
//               if(i<NCELLX && j<NCELLY && k<NCELLZ){
//   	             //   for(k=0;k<NCELLZ-1;k++){
//   	                	if(i==NCELLX-1){
//   											#ifdef DOUBLECOMPLEX
//   	                		 Curl_E=(ez[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]*cexp(-I*k_x*period_x)-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
//   											 #endif
//   											 #ifdef DOUBLEPRECISION
//   											 Curl_E=(ez[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
//   											 #endif
//   	                    	 hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   	                	}
//   	                	else{
//   	                		 Curl_E=(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
//   	                    	 hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   	                }
//   										//Z-PML
//   										if(k<cpml_N_Z && j<cpml_y_lim && k<cpml_z_lim){
//   											psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   											hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
//   										}
//   										if(k>=cpml_F_Z && j<cpml_y_lim && k<cpml_z_lim){
//   											k2 = k - cpml_F_Z;
//   											psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   											hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
//   										}
//                     }
//   	    //             }
//   	    //     }
//   	    // }
//       }
//
//       else{
//   			////#pragma omp parallel for collapse(3) private(Curl_E,i,i2,j,k,k2) // schedule(static)
//   	    // for(i=0;i<NCELLX-1;i++){
//   	    //     for(j=0;j<NCELLY;j++){
//   	    //             for(k=0;k<NCELLZ-1;k++){
//                       if(i<(NCELLX-1) && j<NCELLY && k<(NCELLZ-1)){
//
//   	                    Curl_E=(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
//   	                    hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//
//   											//Z-PML
//   											if(k<cpml_N_Z && j<cpml_y_lim && k<cpml_z_lim){
//   												psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
//   											}
//   											if(k>=cpml_F_Z && j<cpml_y_lim && k<cpml_z_lim){
//   												k2 = k - cpml_F_Z;
//   												psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
//   												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
//   											}
//   											//X-PML
//   											if(i<cpml_N_X && j<cpml_y_lim && k<cpml_z_lim){
//   												psi_Hy_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=bh_x_N[i]*psi_Hy_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ch_x_N[i]*(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
//   												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
//   											}
//   											if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){
//   												i2 = i - cpml_F_X;
//   												psi_Hy_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=bh_x_F[i2]*psi_Hy_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ch_x_F[i2]*(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
//   												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
//   											}
//                       }
//   	    //             }
//   	    //     }
//   	    // }
//       }
//     return;
//   }
//
//   //void UPDATE_hz(void){
//   __global__   void UPDATE_hz(real *hz,real *ey,real *ex,real *Chzh,real *Chze,real *psi_Hz_x_N,real *psi_Hz_x_F,real *psi_Hz_y_N,real *psi_Hz_y_F,real *khdx,real
//     *khdy,real *bh_x_N,real *bh_x_F,real *ch_x_N,real *ch_x_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_X,int cpml_F_X,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlY,int NcpmlX){
//
//       int i,j,k,j2,i2;
//       comp Curl_E;
//       int idx = blockDim.x * blockIdx.x + threadIdx.x;
//
//       i = idx / (NCELLZ*NCELLY);
//       j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
//       k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
//       if(Periodic_XY){
//     //  //#pragma omp parallel for collapse(3) private(Curl_E,i,j,k,j2,i2) // schedule(static)
//   		// for(k=0;k<NCELLZ;k++){
//   	  //   for(i=0;i<NCELLX;i++){
//   	  //       for(j=0;j<NCELLY;j++){
//           if(i<NCELLX && j<NCELLY && k<NCELLZ){
//   	          //      for(k=0;k<NCELLZ;k++){
//   	                	if(i==NCELLX-1 || j== NCELLY-1){
//   	                		if(i==NCELLX-1 && j==NCELLY-1){
//   	                			//printf("%d,%d,%d\n",i,j,k);
//   												#ifdef DOUBLECOMPLEX
//   	                			 Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]*cexp(-I*k_y*period_y)-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]*cexp(-I*k_x*period_x)-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   												 #endif
//   												 #ifdef DOUBLEPRECISION
//   													Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   													#endif
//   	                   			 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   	                		}
//   	                		else if(i==NCELLX-1){
//   												#ifdef DOUBLECOMPLEX
//   	                			 Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]*cexp(-I*k_x*period_x)-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   												 #endif
//   												 #ifdef DOUBLEPRECISION
//   												 Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   												 #endif
//   	                    		 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   	                		}
//   	                		else{
//   												#ifdef DOUBLECOMPLEX
//   	                			 Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]*cexp(-I*k_y*period_y)-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   												 #endif
//   												 #ifdef DOUBLEPRECISION
//   												 Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   												 #endif
//   	                    		 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   	                		}
//   	                	}
//   	                	else{
//   							 Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   	                    	 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   	               	}
//                   }
//   	    //             }
//   	    //     }
//   	    // }
//       }
//
//       else{
//   		//	//#pragma omp parallel for collapse(3) private(Curl_E,i,j,k,j2,i2) // schedule(static)
//   	    // for(i=0;i<NCELLX-1;i++){
//   	    //     for(j=0;j<NCELLY-1;j++){
//   	    //             for(k=0;k<NCELLZ;k++){
//                       if(i<(NCELLX-1) && j<(NCELLY-1) && k<NCELLZ){
//
//   	                    Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
//   	                    hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
//   											//X-PML
//   											if(i<cpml_N_X && j<cpml_y_lim && k<cpml_z_lim){
//   												psi_Hz_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=bh_x_N[i]*psi_Hz_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ch_x_N[i]*(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
//   												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
//   											}
//   											if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){
//   												i2 = i - cpml_F_X;
//   												psi_Hz_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=bh_x_F[i2]*psi_Hz_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ch_x_F[i2]*(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
//   												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
//   											}
//   											//Y-PML
//   											if(j<cpml_N_Y && i<cpml_x_lim && k<cpml_z_lim){
//   												psi_Hz_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]=bh_y_N[j]*psi_Hz_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]+ch_y_N[j]*(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
//   												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)];
//   											}
//   											if(j>=cpml_F_Y && i<cpml_x_lim && k<cpml_z_lim){
//   												j2 = j - cpml_F_Y;
//   												psi_Hz_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]=bh_y_F[j2]*psi_Hz_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]+ch_y_F[j2]*(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
//   												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)];
//   											}
//                       }
//   	    //             }
//   	    //     }
//   	    // }
//       }
//     return;
//   }























  // void UpdateHydroPx(void){
  // 	int i,j,k;
  // 	real Vx1,Vx2,Vx3,Vy1,Vy2,Vy3,Vz1,Vz2,Vz3,Hx1,Hz1,Hy1,Hx2,Hz2,Hy2,Ex1,Ey1,Ez1,VdotGrad,VdotGrad2,DivV,VcrossH,VcrossH2,Pressure,ND1,ND2,ND3,Grad_Div,Grad_Div2;
  // 	real INV_DX,INV_DY,INV_DZ;
  // 	INV_DX = 1.0/dx;
  // 	INV_DY = 1.0/dy;
  // 	INV_DZ = 1.0/dz;
  // 	Grad_Div = 0.0;
  // 	Grad_Div2 =0.0;
  // 	////#pragma omp parallel for collapse(3)  // schedule(static)
  // 	for(i=0;i<NCELLX-1;i++){
  // 			for(j=1;j<NCELLY-1;j++){
  // 							for(k=1;k<NCELLZ-1;k++){
  // 								if(mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)]> first_medium && mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] < 6){
  // 									ND1 = N_EQ + NDx_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
  //
  // 									Grad_Div = Calc_DIV_GRADx(i,j,k);
  // 									Grad_Div2 = Calc_DIV_GRADx2(i,j,k);
  //
  //
  // 									Vx1 = 0.5 * (Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									Vx2 = 0.5 * (Px_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  //
  // 									Vy1 = 0.25 * (Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									Vy2 = 0.25 * (Py_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									//Vy1 = 0.5 * (Vy1 + Vy2);
  // 									Vz1 = 0.25 * (Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									Vz2 = 0.25 * (Pz_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									//Vz1 = 0.5 * (Vz1 + Vz2);
  // 									Hy2 = 0.5 * (hyPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hyPrev[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]);
  // 									Hz2 = 0.5 * (hzPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hzPrev[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]);
  // 									Hy1 = 0.5 * (hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hy[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]);
  // 									Hz1 = 0.5 * (hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hz[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]);
  //
  // 									if(WithConvection) {
  // 										VdotGrad = 0.5*(Vx1*(Px_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
  // 										VdotGrad2 = 0.5*(Vx2*(Px_d_n_1[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Px_d_n_1[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Px_d_n_1[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
  // 										VdotGrad = (VdotGrad - VdotGrad2)/dt;
  // 									}
  // 									else VdotGrad = 0.0;
  //
  // 									if(WithMagField){
  // 										VcrossH = Vy1*Hz1 - Vz1*Hy1;
  // 										VcrossH2 = Vy2*Hz2 - Vz2*Hy2;
  // 										VcrossH = (VcrossH - VcrossH2)/dt;
  // 									}
  // 									else VcrossH = 0.0;
  //
  // 										 Px_d[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] = d_1_d[0]*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_2_d[0]*Px_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] +d_3_d[0]*ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + d_4_d[0]*ex_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + d_NL[0]*(Grad_Div + Grad_Div2/N_EQ)/pow(ND1,1.0/3.0) + d_5_d[0]*(VdotGrad + VcrossH);
  //
  // 								}
  //
  // 							}
  // 						}
  // 					}
  // }
  //
  // void UpdateHydroPy(void){
  // 	int i,j,k;
  // 	real Vx1,Vx2,Vx3,Vy1,Vy2,Vy3,Vz1,Vz2,Vz3,Hx1,Hz1,Hy1,Hx2,Hz2,Hy2,Ex1,Ey1,Ez1,VdotGrad,VdotGrad2,DivV,VcrossH,VcrossH2,Pressure,ND1,ND2,ND3,Grad_Div,Grad_Div2;
  // 	real INV_DX,INV_DY,INV_DZ;
  // 	INV_DX = 1.0/dx;
  // 	INV_DY = 1.0/dy;
  // 	INV_DZ = 1.0/dz;
  // 	Grad_Div = 0.0;
  // 	Grad_Div2 =0.0;
  // 	////#pragma omp parallel for collapse(3) // schedule(static)
  // 	for(i=1;i<NCELLX-1;i++){
  // 			for(j=0;j<NCELLY-1;j++){
  // 							for(k=1;k<NCELLZ-1;k++){
  // 								if(mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)]> first_medium && mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] < 6){
  //
  // 									Grad_Div = Calc_DIV_GRADy(i,j,k);
  // 									Grad_Div2 = Calc_DIV_GRADy2(i,j,k);
  //
  // 									ND1 = N_EQ + NDy_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
  //
  // 									Vy1 = 0.5 * (Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									Vy2 = 0.5 * (Py_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  //
  // 									Vx1 = 0.25 * (Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									Vx2 = 0.25 * (Px_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 								 //	Vx1 = 0.5 * (Vx1 + Vx2);
  // 									Vz1 = 0.25 * (Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									Vz2 = 0.25 * (Pz_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 									//Vz1 = 0.5 * (Vz1 + Vz2);
  // 									Hx1 = 0.5 * (hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hx[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]);
  // 									Hz1 = 0.5 * (hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hz[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]);
  // 									Hx2 = 0.5 * (hxPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hxPrev[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]);
  // 									Hz2 = 0.5 * (hzPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hzPrev[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]);
  //
  //
  // 									 if(WithMagField){
  // 										 VcrossH = Vz1*Hx1 - Vx1*Hz1;
  // 										 VcrossH2 = Vz2*Hx2 - Vx2*Hz2;
  // 										 VcrossH = (VcrossH - VcrossH2)/dt;
  // 									 }
  // 									 else VcrossH = 0.0;
  // 									 if(WithConvection) {
  // 										 VdotGrad = 0.5*(Vx1*(Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Py_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
  // 										 VdotGrad2 = 0.5*(Vx2*(Py_d_n_1[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Py_d_n_1[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Py_d_n_1[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
  // 										 VdotGrad = (VdotGrad - VdotGrad2)/dt;
  // 									 }
  // 									 else VdotGrad = 0.0;
  //
  // 								 Py_d[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] = d_1_d[0]*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_2_d[0]*Py_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_3_d[0]*ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + d_4_d[0]*ey_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + d_NL[0]*(Grad_Div + Grad_Div2/N_EQ)/pow(ND1,1.0/3.0) + d_5_d[0]*(VdotGrad + VcrossH);
  //
  // 								}
  //
  //
  // 							}
  // 						}
  // 					}
  // }
  //
  // void UpdateHydroPz(void){
  // 	int i,j,k;
  // 	real Vx1,Vx2,Vx3,Vy1,Vy2,Vy3,Vz1,Vz2,Vz3,Hx1,Hz1,Hy1,Hx2,Hz2,Hy2,Ex1,Ey1,Ez1,VdotGrad,VdotGrad2,DivV,VcrossH,VcrossH2,Pressure,ND1,ND2,ND3,Grad_Div,Grad_Div2;
  // 	real INV_DX,INV_DY,INV_DZ;
  // 	INV_DX = 1.0/dx;
  // 	INV_DY = 1.0/dy;
  // 	INV_DZ = 1.0/dz;
  // 	Grad_Div = 0.0;
  // 	Grad_Div2 =0.0;
  // 	////#pragma omp parallel for collapse(3) // schedule(static)
  // 	for(i=1;i<NCELLX-1;i++){
  // 			for(j=1;j<NCELLY-1;j++){
  // 							for(k=0;k<NCELLZ-1;k++){
  // 								if(mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)]> first_medium && mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] < 6){
  //
  // 									Grad_Div = Calc_DIV_GRADz(i,j,k);
  // 									Grad_Div2 = Calc_DIV_GRADz2(i,j,k);
  //
  // 									ND1 = N_EQ + NDz_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
  //
  // 								 Vz1 = 0.5 * (Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 								 Vz2 = 0.5 * (Pz_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  //
  // 								 Vx1 = 0.25 * (Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 								 Vx2 = 0.25 * (Px_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 								 //Vx1 = 0.5 * (Vx1 + Vx2);
  // 								 Vy1 = 0.25 * (Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 								 Vy2 = 0.25 * (Py_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  // 								 //	Vy1 = 0.5 * (Vy1 + Vy2);
  // 								 Hy1 = 0.5 * (hy[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)] + hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]);
  // 								 Hx1 = 0.5 * (hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hx[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]);
  // 								 Hy2 = 0.5 * (hyPrev[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)] + hyPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]);
  // 								 Hx2 = 0.5 * (hxPrev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + hxPrev[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]);
  //
  // 								 if(WithMagField){
  // 								 	VcrossH = Vx1*Hy1 - Vy1*Hx1;
  // 								 	VcrossH2 = Vx2*Hy2 - Vy2*Hx2;
  // 								 	VcrossH = (VcrossH - VcrossH2)/dt;
  // 								 }
  // 								 else VcrossH = 0.0;
  //
  // 								 if(WithConvection){
  // 								 	VdotGrad = 0.5*(Vx1*(Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
  // 								 	VdotGrad2 = 0.5*(Vx2*(Pz_d_n_1[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Pz_d_n_1[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Pz_d_n_1[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
  // 								 	VdotGrad = (VdotGrad - VdotGrad2)/dt;
  // 								 }
  // 								 else VdotGrad = 0.0;
  //
  // 								 Pz_d[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] = d_1_d[0]*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_2_d[0]*Pz_d_n_1[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_3_d[0]*ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + d_4_d[0]*ez_n_1[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + d_NL[0]*(Grad_Div + Grad_Div2/N_EQ)/pow(ND1,1.0/3.0) + d_5_d[0]*(VdotGrad + VcrossH);
  //
  // 								}
  //
  //
  // 							}
  // 						}
  // 					}
  // }
  //
  //
  //













  comp Calc_DIV_GRADx2(int i,int j, int k){
  	comp Div_Grad;
  	real INV_DX = 1.0/dx;
  	real INV_DY = 1.0/dy;
  	real INV_DZ = 1.0/dz;
  	int n;

  	Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i+1,j,k,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)])
  						 + INV_DX*INV_DY*(Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i+1,j,k,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i+1,j-1,k,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)])
  						 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i+1,j,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i+1,j,k-1,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]);
  }

  comp Calc_DIV_GRADx(int i,int j,int k){
  	comp Div_Grad;
  	real INV_DX = 1.0/dx;
  	real INV_DY = 1.0/dy;
  	real INV_DZ = 1.0/dz;
  	int n;


  	for(n=0;n<N_drude_poles;n++){
  	if(Diverge_Gradient){

  		if(i==0 && j==0){
  				Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  				 					 + INV_DX*INV_DY*(Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i+1,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}

  		else if(i==NCELLX-1 && j==0){
  				Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  				 					 + INV_DX*INV_DY*(Py_d_n[FourDMap(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(0,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(0,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else if(i==NCELLX-1){
  				Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DY*(Py_d_n[FourDMap(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(0,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(0,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else if(i==0){
  				Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DY*(Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else if(j==0){
  				Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  				 					 + INV_DX*INV_DY*(Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i+1,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else{
  				Div_Grad = INV_DX*INV_DX*(Px_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  				 					 + INV_DX*INV_DY*(Py_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DZ*(Pz_d_n[FourDMap(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}

  	}

   }
  	return Div_Grad;
  }


  comp Calc_DIV_GRADy2(int i,int j,int k){
  	comp Div_Grad;
  	real INV_DX = 1.0/dx;
  	real INV_DY = 1.0/dy;
  	real INV_DZ = 1.0/dz;
  	int n;

  	Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j+1,k,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)])
  						 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i,j+1,k,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i-1,j+1,k,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)])
  						 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j+1,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j+1,k-1,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]);

  	return Div_Grad;
  }


  comp Calc_DIV_GRADy(int i,int j,int k){
  	comp Div_Grad;
  	real INV_DX = 1.0/dx;
  	real INV_DY = 1.0/dy;
  	real INV_DZ = 1.0/dz;
  	int n;
  	for(n=0;n<N_drude_poles;n++){
  	if(Diverge_Gradient){
  		if(i==0 && j==0){
  			Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(NCELLX-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}

  		else if(i==0 && j==NCELLY-1){
  			Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(NCELLX-1,0,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,0,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  			}
  			else if(j==NCELLY-1){
  				Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,0,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  									 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,0,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  			}
  		else if(j==0){
  			Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else if(i==0){
  			Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(NCELLX-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else{
  			Div_Grad = INV_DY*INV_DY*(Py_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DY*(Px_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DY*INV_DZ*(Pz_d_n[FourDMap(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}

  	}


  }
  	//printf("%e\n",Div_Grad);
  	return Div_Grad;
  }


  comp Calc_DIV_GRADz2(int i,int j,int k){
  	comp Div_Grad;
  	real INV_DX = 1.0/dx;
  	real INV_DY = 1.0/dy;
  	real INV_DZ = 1.0/dz;
  	int n;

  	Div_Grad = INV_DZ*INV_DZ*(Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k+1,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)])
  						 + INV_DZ*INV_DY*(Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j,k+1,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j-1,k+1,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)])
  						 + INV_DX*INV_DZ*(Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i,j,k+1,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i-1,j,k+1,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]);
  return Div_Grad;
  }

  comp Calc_DIV_GRADz(int i,int j,int k){
  	comp Div_Grad;
  	real INV_DX = 1.0/dx;
  	real INV_DY = 1.0/dy;
  	real INV_DZ = 1.0/dz;
  	int n;
  	for(n=0;n<N_drude_poles;n++){
  	if(Diverge_Gradient){

  		if(i==0 && j==0){
  			Div_Grad = INV_DZ*INV_DZ*(Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
  			 					 + INV_DZ*INV_DY*(Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,NCELLY-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DZ*(Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(NCELLX-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else if(i==0){
  			Div_Grad = INV_DZ*INV_DZ*(Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
  			 					 + INV_DZ*INV_DY*(Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DZ*(Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(NCELLX-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else if(j==0){
  			Div_Grad = INV_DZ*INV_DZ*(Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DZ*INV_DY*(Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,NCELLY-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DZ*(Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
  		}
  		else{
  			Div_Grad = INV_DZ*INV_DZ*(Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
  			 					 + INV_DZ*INV_DY*(Py_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
  								 + INV_DX*INV_DZ*(Px_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMap(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);

  								//  if(t==173 && k==102){
  								// 	 printf("%e\t%e\t%e\n",Pz_d_n[FourDMap(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMap(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMap(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)],Py_d_n[i][j+1][k+1][0]-Py_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[i][j-1][k-1][0],Px_d_n[i+1][j][k+1][0]-Px_d_n[FourDMap(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMap(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[i-1][j][k-1][0]);
  								// 	 printf("%e\n",Py_d_n[i][j+1][k+1][0]-(Py_d_n[FourDMap(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMap(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)])+Py_d_n[i][j-1][k-1][0]);
  								 //
  								//  }
  								 //1.186946e-66
  		}
  	//	printf("%e\n",Div_Grad);
  	}

  }
  return Div_Grad;
  }
