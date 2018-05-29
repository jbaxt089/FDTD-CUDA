#include <stdio.h>
#include <stdlib.h>
#include "extern_var.h"
#include <math.h>
#include "define.h"
#include <cuda.h>



void READ_DATA_FILE(void){
  FILE *Input;
  Input=fopen("Input.txt", "r");
  char string[300];
  real lam_0;
//printf("HERE\n");
#ifndef FlOATPRECISION
  fscanf(Input,"%s", &string);
  fscanf(Input,"%d %s",&Tend, &string);
  fscanf(Input,"%lf %s",&dx, &string);
  fscanf(Input,"%lf %s",&dy, &string);
  fscanf(Input,"%lf %s",&dz, &string);
  fscanf(Input,"%d %s",&num_trials, &string);
  fscanf(Input,"%d %s",&min_trials, &string);
  fscanf(Input,"%d %s",&max_trials, &string);
  fscanf(Input,"%d %s",&TE_TM, &string);
  fscanf(Input,"%d %s",&NCELLX, &string);
  fscanf(Input,"%d %s",&NCELLY, &string);
  fscanf(Input,"%d %s",&NCELLZ, &string);
  fscanf(Input,"%d %s",&NcpmlX, &string);
  fscanf(Input,"%d %s",&NcpmlY, &string);
  fscanf(Input,"%d %s",&NcpmlZ, &string);
  fscanf(Input,"%d %s",&inf_disp_slab, &string);
  fscanf(Input,"%d %s",&nano_sphere, &string);
  fscanf(Input,"%lf %s",&nano_sphere_radius, &string);
  fscanf(Input,"%d %s",&material, &string);
  fscanf(Input,"%d %s",&first_medium, &string);
  fscanf(Input,"%d %s",&inc_plane, &string);
  fscanf(Input,"%d %s",&dispersive_slab, &string);
  fscanf(Input,"%d %s", &NUM_freq, &string);
  fscanf(Input,"%d %s", &WL_or_freq, &string);
  fscanf(Input,"%d %s", &NUM_freq_inc,&string);
  fscanf(Input,"%d %s", &Tend_inc, &string);
  fscanf(Input,"%d %s", &Snap_in, &string);
  fscanf(Input,"%d %s", &Freq_start, &string);
  fscanf(Input,"%d %s", &Test_offset, &string);
  fscanf(Input,"%d %s", &Spect_loc, &string);
  fscanf(Input,"%lf %s", &BandWidth, &string);
  fscanf(Input,"%d %s", &NONLOCAL, &string);
  fscanf(Input,"%d %s", &Diverge_Gradient, &string);
  fscanf(Input,"%d %s", &t_skip, &string);
  fscanf(Input,"%d %s", &x_skip, &string);
  fscanf(Input,"%d %s", &y_skip, &string);
  fscanf(Input,"%d %s", &z_skip, &string);
  fscanf(Input,"%lf %s", &Nonlocalend, &string);
  fscanf(Input,"%lf %s", &MAX_AMP, &string);
  fscanf(Input,"%lf %s", &TimeFactor, &string);
  fscanf(Input,"%lf %s", &lam_0, &string);
  fscanf(Input,"%d %s", &StaticField, &string);

  printf("%d\n",StaticField);

  f_0 = C0/lam_0;
  printf("Domain Size: %d,%d,%d\n", NCELLX,NCELLY,NCELLZ);
  printf("Sim Time: %d\n", Tend);
  printf("dx: %e\tdy: %e\tdz: %e\n",dx,dy,dz);
  if(inf_disp_slab && material == 2) printf("Infinite dispersive slab: material: Silver\n");
  if(inf_disp_slab && material == 3) printf("Infinite dispersive slab: material: Gold\n");
  if(!inf_disp_slab) printf("No Dispersive Slab\n");
  printf("Number of Trials: %d, Min Trial: %d, Max Trial %d\n", num_trials,min_trials,max_trials);
  printf("Freq Count = %d\n",NUM_freq);
  printf("Inc Freq Count = %d\n",NUM_freq_inc);
  printf("Inc Tend: %d\n",Tend_inc);
  if(WL_or_freq == 1) printf("Wavelength Plot\n");
  else printf("Frequency Plot\n");
  // printf("Sphere Radius: %e\n",nano_sphere_radius);
#else
fscanf(Input,"%s", &string);
fscanf(Input,"%d %s",&Tend, &string);
fscanf(Input,"%f %s",&dx, &string);
fscanf(Input,"%f %s",&dy, &string);
fscanf(Input,"%f %s",&dz, &string);
fscanf(Input,"%d %s",&num_trials, &string);
fscanf(Input,"%d %s",&min_trials, &string);
fscanf(Input,"%d %s",&max_trials, &string);
fscanf(Input,"%d %s",&TE_TM, &string);
fscanf(Input,"%d %s",&NCELLX, &string);
fscanf(Input,"%d %s",&NCELLY, &string);
fscanf(Input,"%d %s",&NCELLZ, &string);
fscanf(Input,"%d %s",&NcpmlX, &string);
fscanf(Input,"%d %s",&NcpmlY, &string);
fscanf(Input,"%d %s",&NcpmlZ, &string);
fscanf(Input,"%d %s",&inf_disp_slab, &string);
fscanf(Input,"%d %s",&nano_sphere, &string);
fscanf(Input,"%f %s",&nano_sphere_radius, &string);
fscanf(Input,"%d %s",&material, &string);
fscanf(Input,"%d %s",&first_medium, &string);
fscanf(Input,"%d %s",&inc_plane, &string);
fscanf(Input,"%d %s",&dispersive_slab, &string);
fscanf(Input,"%d %s", &NUM_freq, &string);
fscanf(Input,"%d %s", &WL_or_freq, &string);
fscanf(Input,"%d %s", &NUM_freq_inc,&string);
fscanf(Input,"%d %s", &Tend_inc, &string);
fscanf(Input,"%d %s", &Snap_in, &string);
fscanf(Input,"%d %s", &Freq_start, &string);
fscanf(Input,"%d %s", &Test_offset, &string);
fscanf(Input,"%d %s", &Spect_loc, &string);
fscanf(Input,"%f %s", &BandWidth, &string);
fscanf(Input,"%d %s", &NONLOCAL, &string);
fscanf(Input,"%d %s", &Diverge_Gradient, &string);
fscanf(Input,"%d %s", &t_skip, &string);
fscanf(Input,"%d %s", &x_skip, &string);
fscanf(Input,"%d %s", &y_skip, &string);
fscanf(Input,"%d %s", &z_skip, &string);
fscanf(Input,"%f %s", &Nonlocalend, &string);
fscanf(Input,"%lf %s", &MAX_AMP, &string);
fscanf(Input,"%lf %s", &TimeFactor, &string);
fscanf(Input,"%lf %s", &lam_0, &string);
fscanf(Input,"%d %s", &StaticField, &string);
#endif

  if(Diverge_Gradient == 1)
{
  Laplacian =0;
}
else if(Diverge_Gradient==0) Laplacian =1;
  //printf("HERE2\n");

  fclose(Input);

  //printf("HERE3\n");


}

void MATERIAL_MATRIX(void){
cudaError_t err;
  mat_matrix = MALLOC3D_int(mat_matrix, NCELLX,NCELLY,NCELLZ);
  mat_matrixX = MALLOC3D_int(mat_matrixX, NCELLX,NCELLY,NCELLZ);
  mat_matrixY = MALLOC3D_int(mat_matrixY, NCELLX,NCELLY,NCELLZ);
  mat_matrixZ = MALLOC3D_int(mat_matrixZ, NCELLX,NCELLY,NCELLZ);

err=  cudaMalloc(&mat_matrixdev,NCELLX*NCELLY*NCELLZ*sizeof(int));
err=  cudaMalloc(&mat_matrixXdev,NCELLX*NCELLY*NCELLZ*sizeof(int));
err=   cudaMalloc(&mat_matrixYdev,NCELLX*NCELLY*NCELLZ*sizeof(int));
err=  cudaMalloc(&mat_matrixZdev,NCELLX*NCELLY*NCELLZ*sizeof(int));
if( cudaSuccess != err)
{
    printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
    exit(-1);
}

  int i,j,k;

  printf("Material: %d\n",material);

  if(inf_disp_slab==1){
    if(PBC_CTW==1) dispersive_slab=inc_plane+dispersive_slab;
    if(PBC_CTW==0) dispersive_slab=dispersive_slab + NtfsfY;
    for(i=0; i<NCELLX; i++){
      for(j=0;j<NCELLY; j++){
        for(k=0; k<NCELLZ; k++){

            if(k>=dispersive_slab){
              mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
              mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
              mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
              mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
            }
            else{
              mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;


            }
        }
      }
    }
  }
  else if(inf_disp_slab == 2){
    if(PBC_CTW==1) dispersive_slab=inc_plane+dispersive_slab;
    if(PBC_CTW==0) dispersive_slab=dispersive_slab + NtfsfY;
    for(i=0; i<NCELLX; i++){
      for(j=0;j<NCELLY; j++){
        for(k=0; k<NCELLZ; k++){
            if(k>=dispersive_slab && k <(dispersive_slab + (int)(nano_sphere_radius/dz))){
              mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
              mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
              mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
              mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;

            }
            else{
              mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;

            }
        }
      }
    }
  }

  else{
    for(i=0; i<NCELLX; i++){
      for(j=0;j<NCELLY; j++){
        for(k=0; k<NCELLZ; k++){
              mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
              mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = first_medium;
        }
      }
    }
  }

  if(nano_sphere==1){
    int nano_sphere_radius_int,center_point;
    nano_sphere_radius_int = (int)(nano_sphere_radius/dx);
    double center_point_z;
    if(Periodic_XY)  center_point_z = dispersive_slab - nano_sphere_radius_int;
    else  center_point_z = NCELLZ/2;
    double center_point_x = NCELLX/2;
    double center_point_y = NCELLY/2;
    printf("%f\t%f\t%f\n",center_point_x,center_point_y,center_point_z);
    for(i=0; i<NCELLX; i++){
      for(j=0;j<NCELLY; j++){
        for(k=0; k<NCELLZ; k++){
          if(pow((i-center_point_x)*(i-center_point_x)+(j-center_point_y)*(j-center_point_y)+(k-center_point_z)*(k-center_point_z),0.5) < nano_sphere_radius_int){
            mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }
          if(pow((i+0.5-center_point_x)*(i+0.5-center_point_x)+(j-center_point_y)*(j-center_point_y)+(k-center_point_z)*(k-center_point_z),0.5) < nano_sphere_radius_int){
            mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }

          if(pow((i-center_point_x)*(i-center_point_x)+(j+0.5-center_point_y)*(j+0.5-center_point_y)+(k-center_point_z)*(k-center_point_z),0.5) < nano_sphere_radius_int){
            mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }

          if(pow((i-center_point_x)*(i-center_point_x)+(j-center_point_y)*(j-center_point_y)+(k+0.5-center_point_z)*(k+0.5-center_point_z),0.5) < nano_sphere_radius_int){
            mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }
        }
      }
    }
  }
  else if(nano_sphere==2){
    for(i=0; i<NCELLX; i++){
      for(j=0;j<NCELLY; j++){
        for(k=0; k<NCELLZ; k++){
          double center_point_z = dispersive_slab;
          double center_point_x = NCELLX/2;
          double center_point_y = NCELLY/2;
          if(i<center_point_x+2 && i>center_point_x-2 && j<center_point_y+2 && j>center_point_y-2 && k>center_point_z-20 ){
            mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }
        }
      }
    }
  }

  if(nano_sphere==3){
    int nano_sphere_radius_int,center_point;
    nano_sphere_radius_int = (int)(nano_sphere_radius/dx);
    double center_point_z;
    if(Periodic_XY)  center_point_z = dispersive_slab - nano_sphere_radius_int;
    else  center_point_z = NCELLZ/2;
    double center_point_x = NCELLX/2;
    double center_point_y = NCELLY/2;
    printf("%f\t%f\t%f\n",center_point_x,center_point_y,center_point_z);
    for(i=0; i<NCELLX; i++){
      for(j=0;j<NCELLY; j++){
        for(k=0; k<NCELLZ; k++){
          if((i-center_point_x)*(i-center_point_x) < nano_sphere_radius_int*nano_sphere_radius_int && (j-center_point_y)*(j-center_point_y) < nano_sphere_radius_int*nano_sphere_radius_int && (k-center_point_z)*(k-center_point_z) < nano_sphere_radius_int*nano_sphere_radius_int){
            mat_matrix[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }
          if((i+0.5-center_point_x)*(i+0.5-center_point_x) < nano_sphere_radius_int*nano_sphere_radius_int&& (j-center_point_y)*(j-center_point_y) < nano_sphere_radius_int*nano_sphere_radius_int && (k-center_point_z)*(k-center_point_z) < nano_sphere_radius_int*nano_sphere_radius_int){
            mat_matrixX[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }

          if((i-center_point_x)*(i-center_point_x)  < nano_sphere_radius_int*nano_sphere_radius_int && (j+0.5-center_point_y)*(j+0.5-center_point_y) < nano_sphere_radius_int*nano_sphere_radius_int && (k-center_point_z)*(k-center_point_z)< nano_sphere_radius_int*nano_sphere_radius_int){
            mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }

          if((i-center_point_x)*(i-center_point_x) < nano_sphere_radius_int*nano_sphere_radius_int && (j-center_point_y)*(j-center_point_y)  < nano_sphere_radius_int*nano_sphere_radius_int&& (k+0.5-center_point_z)*(k+0.5-center_point_z)  < nano_sphere_radius_int*nano_sphere_radius_int){
            mat_matrixZ[ThreeDMap(i,j,k,NCELLZ,NCELLY)] = material;
          }
        }
      }
    }
  }

  //Output MaterialPlot
  FILE *MaterialPlot;
  MaterialPlot = fopen("Material_XZ.txt","w");
  j = (int)NCELLY/2;
  for(i=0;i<NCELLX;i++){
      for(k=0;k<NCELLZ;k++){
        fprintf(MaterialPlot,"%d\t",mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)]);
      }
      fprintf(MaterialPlot,"\n");
    }
    fclose(MaterialPlot);

    MaterialPlot = fopen("Material_XY.txt","w");
    k = (int)NCELLZ/2;
    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
          fprintf(MaterialPlot,"%d\t",mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)]);
        }
        fprintf(MaterialPlot,"\n");
      }
      fclose(MaterialPlot);


      MaterialPlot = fopen("Material_YZ.txt","w");
      i = (int)NCELLX/2;

          for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
            fprintf(MaterialPlot,"%d\t",mat_matrixY[ThreeDMap(i,j,k,NCELLZ,NCELLY)]);
          }
          fprintf(MaterialPlot,"\n");
        }
        fclose(MaterialPlot);





    Nonlocalend_int = dispersive_slab + ceil(Nonlocalend/dz);
    printf("Nonlocal Barrier: %d\n",Nonlocalend_int);


    cudaMemcpy(mat_matrixdev,mat_matrix,NCELLX*NCELLY*NCELLZ*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(mat_matrixXdev,mat_matrixX,NCELLX*NCELLY*NCELLZ*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(mat_matrixYdev,mat_matrixY,NCELLX*NCELLY*NCELLZ*sizeof(int),cudaMemcpyHostToDevice);
    cudaMemcpy(mat_matrixZdev,mat_matrixZ,NCELLX*NCELLY*NCELLZ*sizeof(int),cudaMemcpyHostToDevice);

}



//Sets variables
void SETUP_CONST(){

    source_end=1230;
    inc_Length=NCELLZ*6;

    int i;
    //only if dx=dy=dz
    if(NCELLX==1 || NCELLY==1 || NCELLZ==1){
        dt=0.9*dy/(C0*pow(2,0.5));
    }
    else{

        dt=0.5*dy/C0/(real)TimeFactor;
         printf("dt = %e\n",dt);
    }

    //constants
    z0=Z0;
  //  z0 = 1.0;
    ep0=EP0;
    mu0=MU0;
    c0=C0;
    pi=PI;
    me = ME;
    e0 = E0;
    //CPML constants
    // NcpmlX=0;
    // NcpmlY=0;
    // NcpmlZ=100;
    //add the CPML cells to the regular domain
    NCELLX+=2*NcpmlX;
    NCELLY+=2*NcpmlY;
    NCELLZ+=2*NcpmlZ;

    cpml_N_X=NcpmlX;
    cpml_F_X=NCELLX-NcpmlX;
    cpml_N_Y=NcpmlY;
    cpml_F_Y=NCELLY-NcpmlY;
    cpml_N_Z=NcpmlZ;
    cpml_F_Z=NCELLZ-NcpmlZ;

    //exceitation plane
    inc_plane=NcpmlZ+inc_plane;
    //inc_plane=NCELLZ/2;
    NCELLZ+=1;
    NCELLY+=1;
    NCELLX+=1;
    cpml_x_lim=NCELLX-1;
    cpml_y_lim=NCELLY-1;
    cpml_z_lim=NCELLZ-1;
    //Lateral PML Boundaries
    Absorption = 0;
    Scattering = 0;

    if(Periodic_XY){
        NCELLX-=1;
        NCELLY-=1;
        cpml_x_lim=NCELLX;
        cpml_y_lim=NCELLY;
    }

    else if(Periodic_XZ){
        NCELLX-=1;
        NCELLZ-=1;
        cpml_x_lim=NCELLX;
        cpml_z_lim=NCELLZ;
    }

    else if(Periodic_YZ){
        NCELLY-=1;
        NCELLZ-=1;
        cpml_y_lim=NCELLY;
        cpml_z_lim=NCELLZ;
    }
    else{
      Absorption = 1;
      Scattering = 1;
    }

    if(PBC_CTW == 0)
    {
    NtfsfX=NtfsfY=NtfsfZ=10;
    if(Periodic_XY) NtfsfX = NtfsfY = 0;
    NtfsfX+=NcpmlX;
    NtfsfY+=NcpmlY;
    NtfsfZ+=NcpmlZ;
    printf("TFSF:%d\t%d\t%d\n",NtfsfX,NtfsfY,NtfsfZ);
    if(Absorption && Scattering){
      XSTARTAbs = NtfsfX + 5;
      YSTARTAbs = NtfsfY + 5;
      ZSTARTAbs = NtfsfZ + 5;
      XENDAbs = NCELLX-NtfsfX - 5;
      YENDAbs = NCELLY-NtfsfY - 5;
      ZENDAbs = NCELLZ-NtfsfZ - 5;

      printf("%d\t%d\n",XENDAbs ,XSTARTAbs);
      printf("%d\t%d\n",YENDAbs , YSTARTAbs);
      printf("%d\t%d\n",ZENDAbs , ZSTARTAbs);


      XSTARTSca = NtfsfX - 5;
      YSTARTSca = NtfsfY - 5;
      ZSTARTSca = NtfsfZ - 5;
      XENDSca = NCELLX-NtfsfX + 5;
      YENDSca = NCELLY-NtfsfY + 5;
      ZENDSca = NCELLZ-NtfsfZ + 5;

      XNEARAbs = XSTARTAbs;
      XFARAbs = XENDAbs-1;
      YNEARAbs = YSTARTAbs;
      YFARAbs = YENDAbs-1;
      ZNEARAbs = ZSTARTAbs;
      ZFARAbs = ZENDAbs-1;

      XNEARSca = XSTARTSca;
      XFARSca = XENDSca-1;
      YNEARSca = YSTARTSca;
      YFARSca = YENDSca-1;
      ZNEARSca = ZSTARTSca;
      ZFARSca = ZENDSca-1;

      printf("%d\t%d\n",XENDSca ,XSTARTSca);
      printf("%d\t%d\n",YENDSca , YSTARTSca);
      printf("%d\t%d\n",ZENDSca , ZSTARTSca);

    }
    if(Scattering){

    }
    }
    else if(PBC_CTW == 1){
      NtfsfX=NtfsfY=NtfsfZ=0;
    }
    //polynomial grading exponent (recommended 3 or 4)
    cpml_exp=3;
    //k_max
  //  max_stretch_factor_x=max_stretch_factor_y=max_stretch_factor_z=15;
    max_stretch_factor_x=max_stretch_factor_y=max_stretch_factor_z=5;

    //sigma_max (see eqn 7.66)
    // max_sigma_cpml_x=0.8*(cpml_exp+1)/(Z0*dx*pow(10,0.5));
    // max_sigma_cpml_y=0.8*(cpml_exp+1)/(Z0*dy*pow(10,0.5));
    // max_sigma_cpml_z=0.8*(cpml_exp+1)/(Z0*dz*pow(10,0.5));
    max_sigma_cpml_x=0.8*(cpml_exp+1)/(Z0*dx);
    max_sigma_cpml_y=0.8*(cpml_exp+1)/(Z0*dy);
    max_sigma_cpml_z=0.8*(cpml_exp+1)/(Z0*dz);


    //Suggested by CONVOLUTION PML CPML : ANEFFICIENT FDTD IMPLEMENTATION OFTHE CFS
    //� PML FOR ARBITRARY MEDIA by Roden et. Al
    max_alpha_x=max_alpha_y=max_alpha_z=0.05;

    //See eqn (7.79)
    exp_alpha_x=exp_alpha_y=exp_alpha_z=1;
     //PBC parameters
    period_x=(NCELLX-1)*dx;
    period_y=(NCELLY-1)*dy;

    //Fourier Analysis
    e_reflected=e_incident=0.0;


    first_medium_max = 1;
}

//Everything nessecary for CPML to work
void SETUP_CPML(void){

    int i;

    //set the ked_ and khd_ vectors to d_ (will further update later)
    for(i=0;i<NCELLX;i++){
        kedx[i]=dx;
        khdx[i]=dx;
    }
    for(i=0;i<NCELLY;i++){
        kedy[i]=dy;
        khdy[i]=dy;
    }
    for(i=0;i<NCELLZ;i++){
        kedz[i]=dz;
        khdz[i]=dz;
    }

   SETUP_CPML_X();
   printf("Setting Up CPMPLZ\n");
   SETUP_CPML_Z();
   printf("Done\n");
   SETUP_CPML_Y();

    return;
}

void SETUP_CPML_X(void){

    int i,j,k;
    real S,A,K;

    for(i=0;i<NcpmlX+1;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){

                psi_Ey_x_F[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
                psi_Ey_x_N[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
                psi_Ez_x_F[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
                psi_Ez_x_N[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;

        }
      }
    }

 for(i=0; i<NcpmlX;i++){
     for(j=0;j<NCELLY;j++){
         for(k=0;k<NCELLZ;k++){

                psi_Hy_x_F[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
                psi_Hy_x_N[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
                psi_Hz_x_F[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
                psi_Hz_x_N[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;


        }
      }
    }

    cudaMemcpy(psi_Ey_x_Fdev,psi_Ey_x_F,(NcpmlX+1)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ey_x_Ndev,psi_Ey_x_N,(NcpmlX+1)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ez_x_Fdev,psi_Ez_x_F,(NcpmlX+1)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ez_x_Ndev,psi_Ez_x_N,(NcpmlX+1)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);


    cudaMemcpy(psi_Hy_x_Fdev,psi_Hy_x_F,(NcpmlX)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hy_x_Ndev,psi_Hy_x_N,(NcpmlX)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hz_x_Fdev,psi_Hz_x_F,(NcpmlX)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hz_x_Ndev,psi_Hz_x_N,(NcpmlX)*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);

    //Y COMPONENT for E field
    for(i=0;i<NcpmlX+1;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow(((real)(NcpmlX-i)/((real)(NcpmlX))),cpml_exp)*max_sigma_cpml_x;
      // equation (7.79)
      A=max_alpha_x*pow((real)i/((real)(NcpmlX)),exp_alpha_x);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_x-1)*pow((real)(NcpmlX-i)/((real)(NcpmlX)),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      be_x_N[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)

    if(K*A+S==0.0){
        ce_x_N[i]=0.0;
      }
      else{
        ce_x_N[i]=S*(be_x_N[i]-1)/(K*(S+K*A));
      }
      kedx[i]=kedx[i]*K;
    }

    for(i=0;i<NcpmlX+1;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow(((real)(i)/((real)(NcpmlX))),cpml_exp)*(max_sigma_cpml_x);
      // equation (7.79)
      A=max_alpha_x*pow((real)(NcpmlX-i)/((real)(NcpmlX)),exp_alpha_x);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_x-1)*pow((real)(i)/((real)(NcpmlX)),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      be_x_F[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)
      if(K*A+S==0.0){
        ce_x_F[i]=0.0;
      }
      else{
        ce_x_F[i]=S*(be_x_F[i]-1)/(K*(S+K*A));
      }

      kedx[i+cpml_F_X]=kedx[i+cpml_F_X]*K;
      //printf("%f\t%f\t%f\n",be_z_F[i],ce_z_F[i],kedz[i+cpml_F_Z]);
    }

    //Y COMPONENT for H-field
    for(i=0;i<NcpmlX;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow((((real)NcpmlX-(real)i-0.5)/((real)NcpmlX)),cpml_exp)*max_sigma_cpml_x;
      // equation (7.79)
      A=max_alpha_x*pow(((real)i+0.5)/((real)NcpmlX),exp_alpha_x);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_x-1)*pow(((real)NcpmlX-(real)i-0.5)/((real)NcpmlX),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      bh_x_N[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)
      if(K*A+S==0.0){
        ch_x_N[i]=0;
      }
      else{
        ch_x_N[i]=S*(bh_x_N[i]-1)/(K*(S+K*A));
      }

      khdx[i]=khdx[i]*K;
    }
    /*bh_y_N[NcpmlY-1]=be_y_N[NcpmlY-1];
    ch_y_N[NcpmlY-1]=ce_y_N[NcpmlY-1];
    khdy[NcpmlY-1]=kedy[NcpmlY-1];*/

    for(i=0;i<NcpmlX;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow((((real)i+0.5)/((real)NcpmlX)),cpml_exp)*(max_sigma_cpml_x);
      // equation (7.79)
      A=max_alpha_x*pow(((real)NcpmlX-(real)i-0.5)/((real)NcpmlX),exp_alpha_x);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_x-1)*pow(((real)i+0.5)/((real)NcpmlX),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      bh_x_F[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)

    if(K*A+S==0.0){
        ch_x_F[i]=0;
      }
      else{
        ch_x_F[i]=S*(bh_x_F[i]-1)/(K*(S+K*A));
      }

      khdx[i+cpml_F_X]=khdx[i+cpml_F_X]*K;
    }
    /*bh_y_F[NcpmlY-1]=be_y_F[NcpmlY-1];
    ch_y_F[NcpmlY-1]=ce_y_F[NcpmlY-1];
    khdy[cpml_F_Y]=kedy[cpml_F_Y];*/
    cudaMemcpy(kedxdev,kedx,NCELLX*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(khdxdev,khdx,NCELLX*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ce_x_Ndev,ce_x_N,(NcpmlX+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ce_x_Fdev,ce_x_F,(NcpmlX+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(be_x_Ndev,be_x_N,(NcpmlX+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(be_x_Fdev,be_x_F,(NcpmlX+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ch_x_Ndev,ch_x_N,(NcpmlX)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ch_x_Fdev,ch_x_F,(NcpmlX)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(bh_x_Ndev,bh_x_N,(NcpmlX)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(bh_x_Fdev,bh_z_F,(NcpmlX)*sizeof(real),cudaMemcpyHostToDevice);
}

void SETUP_CPML_Y(void){

    int i,j,k;
    real S,A,K;

    for(i=0; i<NCELLX;i++){
        for(j=0;j<NcpmlY+1;j++){
            for(k=0;k<NCELLZ;k++){

                psi_Ex_y_F[ThreeDMap(i,j,k,NCELLZ,NcpmlY+1)]=0;
                psi_Ex_y_N[ThreeDMap(i,j,k,NCELLZ,NcpmlY+1)]=0;
                psi_Ez_y_F[ThreeDMap(i,j,k,NCELLZ,NcpmlY+1)]=0;
                psi_Ez_y_N[ThreeDMap(i,j,k,NCELLZ,NcpmlY+1)]=0;

        }
      }
    }

 for(i=0; i<NCELLX;i++){
     for(j=0;j<NcpmlY;j++){
         for(k=0;k<NCELLZ;k++){

                psi_Hx_y_F[ThreeDMap(i,j,k,NCELLZ,NcpmlY)]=0;
                psi_Hx_y_N[ThreeDMap(i,j,k,NCELLZ,NcpmlY)]=0;
                psi_Hz_y_F[ThreeDMap(i,j,k,NCELLZ,NcpmlY)]=0;
                psi_Hz_y_N[ThreeDMap(i,j,k,NCELLZ,NcpmlY)]=0;

        }
      }
    }
    cudaMemcpy(psi_Ex_y_Fdev,psi_Ex_y_F,(NcpmlY+1)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ex_y_Ndev,psi_Ex_y_N,(NcpmlY+1)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ez_y_Fdev,psi_Ez_y_F,(NcpmlY+1)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ez_y_Ndev,psi_Ez_y_N,(NcpmlY+1)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);


    cudaMemcpy(psi_Hx_y_Fdev,psi_Hx_y_F,(NcpmlY)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hx_y_Ndev,psi_Hx_y_N,(NcpmlY)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hz_y_Fdev,psi_Hz_y_F,(NcpmlY)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hz_y_Ndev,psi_Hz_y_N,(NcpmlY)*NCELLX*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    //Y COMPONENT for E field
    for(i=0;i<NcpmlY+1;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow(((real)(NcpmlY-i)/((real)(NcpmlY))),cpml_exp)*max_sigma_cpml_y;
      // equation (7.79)
      A=max_alpha_y*pow((real)i/((real)(NcpmlY)),exp_alpha_y);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_y-1)*pow((real)(NcpmlY-i)/((real)(NcpmlY)),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      be_y_N[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)

    if(K*A+S==0.0){
        ce_y_N[i]=0.0;
      }
      else{
        ce_y_N[i]=S*(be_y_N[i]-1)/(K*(S+K*A));
      }
      kedy[i]=kedy[i]*K;
    }

    for(i=0;i<NcpmlY+1;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow(((real)(i)/((real)(NcpmlY))),cpml_exp)*(max_sigma_cpml_y);
      // equation (7.79)
      A=max_alpha_y*pow((real)(NcpmlY-i)/((real)(NcpmlY)),exp_alpha_y);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_y-1)*pow((real)(i)/((real)(NcpmlY)),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      be_y_F[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)
      if(K*A+S==0.0){
        ce_y_F[i]=0.0;
      }
      else{
        ce_y_F[i]=S*(be_y_F[i]-1)/(K*(S+K*A));
      }

      kedy[i+cpml_F_Y]=kedy[i+cpml_F_Y]*K;
      //printf("%f\t%f\t%f\n",be_z_F[i],ce_z_F[i],kedz[i+cpml_F_Z]);
    }

    //Y COMPONENT for H-field
    for(i=0;i<NcpmlY;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow((((real)NcpmlY-(real)i-0.5)/((real)NcpmlY)),cpml_exp)*max_sigma_cpml_y;
      // equation (7.79)
      A=max_alpha_y*pow(((real)i+0.5)/((real)NcpmlY),exp_alpha_y);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_y-1)*pow(((real)NcpmlY-(real)i-0.5)/((real)NcpmlY),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      bh_y_N[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)
      if(K*A+S==0.0){
        ch_y_N[i]=0;
      }
      else{
        ch_y_N[i]=S*(bh_y_N[i]-1)/(K*(S+K*A));
      }

      khdy[i]=khdy[i]*K;
    }
    /*bh_y_N[NcpmlY-1]=be_y_N[NcpmlY-1];
    ch_y_N[NcpmlY-1]=ce_y_N[NcpmlY-1];
    khdy[NcpmlY-1]=kedy[NcpmlY-1];*/

    for(i=0;i<NcpmlY;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow((((real)i+0.5)/((real)NcpmlY)),cpml_exp)*(max_sigma_cpml_y);
      // equation (7.79)
      A=max_alpha_y*pow(((real)NcpmlY-(real)i-0.5)/((real)NcpmlY),exp_alpha_y);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_y-1)*pow(((real)i+0.5)/((real)NcpmlY),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      bh_y_F[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)

    if(K*A+S==0.0){
        ch_y_F[i]=0;
      }
      else{
        ch_y_F[i]=S*(bh_y_F[i]-1)/(K*(S+K*A));
      }

      khdy[i+cpml_F_Y]=khdy[i+cpml_F_Y]*K;
    }
    /*bh_y_F[NcpmlY-1]=be_y_F[NcpmlY-1];
    ch_y_F[NcpmlY-1]=ce_y_F[NcpmlY-1];
    khdy[cpml_F_Y]=kedy[cpml_F_Y];*/

    cudaMemcpy(kedydev,kedy,NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(khdydev,khdy,NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ce_y_Ndev,ce_y_N,(NcpmlY+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ce_y_Fdev,ce_y_F,(NcpmlY+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(be_y_Ndev,be_y_N,(NcpmlY+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(be_y_Fdev,be_y_F,(NcpmlY+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ch_y_Ndev,ch_y_N,(NcpmlY)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ch_y_Fdev,ch_y_F,(NcpmlY)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(bh_y_Ndev,bh_y_N,(NcpmlY)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(bh_y_Fdev,bh_y_F,(NcpmlY)*sizeof(real),cudaMemcpyHostToDevice);



}



void SETUP_CPML_Z(void){

    int i,j,k;
    real S,A,K;

    for(i=0; i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NcpmlZ+1;k++){

                psi_Ex_z_F[ThreeDMap(i,j,k,NcpmlZ+1,NCELLY)]=0.0;
                psi_Ex_z_N[ThreeDMap(i,j,k,NcpmlZ+1,NCELLY)]=0.0;
                psi_Ey_z_F[ThreeDMap(i,j,k,NcpmlZ+1,NCELLY)]=0.0;
                psi_Ey_z_N[ThreeDMap(i,j,k,NcpmlZ+1,NCELLY)]=0.0;


        }
      }
    }
// printf("here\n" );
 for(i=0; i<NCELLX;i++){
     for(j=0;j<NCELLY;j++){
         for(k=0;k<NcpmlZ;k++){

                psi_Hx_z_F[ThreeDMap(i,j,k,NcpmlZ,NCELLY)]=0.0;
                psi_Hx_z_N[ThreeDMap(i,j,k,NcpmlZ,NCELLY)]=0.0;
                psi_Hy_z_F[ThreeDMap(i,j,k,NcpmlZ,NCELLY)]=0.0;
                psi_Hy_z_N[ThreeDMap(i,j,k,NcpmlZ,NCELLY)]=0.0;


        }
      }
    }

    // printf("here\n" );
    cudaMemcpy(psi_Ex_z_Fdev,psi_Ex_z_F,(NcpmlZ+1)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ex_z_Ndev,psi_Ex_z_N,(NcpmlZ+1)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ey_z_Fdev,psi_Ey_z_F,(NcpmlZ+1)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Ey_z_Ndev,psi_Ey_z_N,(NcpmlZ+1)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);


    cudaMemcpy(psi_Hx_z_Fdev,psi_Hx_z_F,(NcpmlZ)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hx_z_Ndev,psi_Hx_z_N,(NcpmlZ)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hy_z_Fdev,psi_Hy_z_F,(NcpmlZ)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(psi_Hy_z_Ndev,psi_Hy_z_N,(NcpmlZ)*NCELLX*NCELLY*sizeof(real),cudaMemcpyHostToDevice);
    //Y COMPONENT for E field
    for(i=0;i<NcpmlZ+1;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow(((real)(NcpmlZ-i)/((real)(NcpmlZ))),cpml_exp)*max_sigma_cpml_z;
      // equation (7.79)
      A=max_alpha_z*pow((real)i/((real)(NcpmlZ)),exp_alpha_z);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_z-1)*pow((real)(NcpmlZ-i)/((real)(NcpmlZ)),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      be_z_N[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)

    if(K*A+S==0.0){
        ce_z_N[i]=0.0;
      }
      else{
        ce_z_N[i]=S*(be_z_N[i]-1)/(K*(S+K*A));
      }
      kedz[i]=kedz[i]*K;
    }
    //printf("here\n" );

    for(i=0;i<NcpmlZ+1;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow(((real)(i)/((real)(NcpmlZ))),cpml_exp)*(max_sigma_cpml_z);
      // equation (7.79)
      A=max_alpha_z*pow((real)(NcpmlZ-i)/((real)(NcpmlZ)),exp_alpha_z);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_z-1)*pow((real)(i)/((real)(NcpmlZ)),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      be_z_F[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)
      if(K*A+S==0.0){
        ce_z_F[i]=0.0;
      }
      else{
        ce_z_F[i]=S*(be_z_F[i]-1)/(K*(S+K*A));
      }

      kedz[i+cpml_F_Z]=kedz[i+cpml_F_Z]*K;
      //printf("%f\t%f\t%f\n",be_z_F[i],ce_z_F[i],kedz[i+cpml_F_Z]);
    }

    //Y COMPONENT for H-field
    for(i=0;i<NcpmlZ;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow((((real)NcpmlZ-(real)i-0.5)/((real)NcpmlZ)),cpml_exp)*max_sigma_cpml_z;
      // equation (7.79)
      A=max_alpha_z*pow(((real)i+0.5)/((real)NcpmlZ),exp_alpha_z);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_z-1)*pow(((real)NcpmlZ-(real)i-0.5)/((real)NcpmlZ),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      bh_z_N[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)
      if(K*A+S==0.0){
        ch_z_N[i]=0;
      }
      else{
        ch_z_N[i]=S*(bh_z_N[i]-1)/(K*(S+K*A));
      }

      khdz[i]=khdz[i]*K;
    }
    /*bh_y_N[NcpmlY-1]=be_y_N[NcpmlY-1];
    ch_y_N[NcpmlY-1]=ce_y_N[NcpmlY-1];
    khdy[NcpmlY-1]=kedy[NcpmlY-1];*/
    printf("here\n" );

    for(i=0;i<NcpmlZ;i++){
      //PML grading in the x direction (eqn 7.60a)
      S=pow((((real)i+0.5)/((real)NcpmlZ)),cpml_exp)*(max_sigma_cpml_z);
      // equation (7.79)
      A=max_alpha_z*pow(((real)NcpmlZ-(real)i-0.5)/((real)NcpmlZ),exp_alpha_z);
      // (eqn 7.60b)
      K=1+(max_stretch_factor_z-1)*pow(((real)i+0.5)/((real)NcpmlZ),cpml_exp);
      // (eqn 7.102) part of the c_w term. Also see (7.114a)
      bh_z_F[i]=exp(-(S/K+A)*(dt/EP0));
      //7.114(b)

    if(K*A+S==0.0){
        ch_z_F[i]=0;
      }
      else{
        ch_z_F[i]=S*(bh_z_F[i]-1)/(K*(S+K*A));
      }

      khdz[i+cpml_F_Z]=khdz[i+cpml_F_Z]*K;
    }
    /*bh_y_F[NcpmlY-1]=be_y_F[NcpmlY-1];
    ch_y_F[NcpmlY-1]=ce_y_F[NcpmlY-1];
    khdy[cpml_F_Y]=kedy[cpml_F_Y];*/
    printf("here\n" );
    cudaMemcpy(kedzdev,kedz,NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(khdzdev,khdz,NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ce_z_Ndev,ce_z_N,(NcpmlZ+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ce_z_Fdev,ce_z_F,(NcpmlZ+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(be_z_Ndev,be_z_N,(NcpmlZ+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(be_z_Fdev,be_z_F,(NcpmlZ+1)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ch_z_Ndev,ch_z_N,(NcpmlZ)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ch_z_Fdev,ch_z_F,(NcpmlZ)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(bh_z_Ndev,bh_z_N,(NcpmlZ)*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(bh_z_Fdev,bh_z_F,(NcpmlZ)*sizeof(real),cudaMemcpyHostToDevice);
}

void SETUP_Drude_CP(void){
  real FERMI_VELOCITY = 1.39e8 / 100.0; //for silver and gold (m/s)
  real D = 2.0;
  real C = 3.0*D/(D+2.0);
  real NL_COEFF = sqrt(C/D);

//  NL_COEFF = 1.0;
  real NONLOC;
    N_drude_poles=1;
    N_CP_poles=2;
    N_lorentz_poles = 2;
    int i,k,j,n;
    real something=0.0;
    Hydrodynamics = 0;
    WithMagField =0;
    WithConvection = 0;
    if(NONLOCAL == 3){
       Hydrodynamics = 1;
       WithMagField = 0;
       WithConvection = 0;
     }
     size_t extentD = NCELLX*NCELLY*NCELLZ*N_drude_poles*sizeof(real);
     size_t extentCP = NCELLX*NCELLY*NCELLZ*N_CP_poles*sizeof(real);
cudaError_t err;
     printf("HERE\n");
    //Allocate memory for parameter vectors
  err=  cudaMalloc(&C_1_cpdev,N_CP_poles*sizeof(real));
  err=  cudaMalloc(&C_2_cpdev,N_CP_poles*sizeof(real));
  err=  cudaMalloc(&C_3_cpdev,N_CP_poles*sizeof(real));
  err=  cudaMalloc(&C_4_cpdev,N_CP_poles*sizeof(real));
  err=  cudaMalloc(&C_5_cpdev,N_CP_poles*sizeof(real));
  if( cudaSuccess != err)
  {
      printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
      exit(-1);
  }
    C_1_cp=MALLOC1D(C_1_cp,N_CP_poles);
    C_2_cp=MALLOC1D(C_2_cp,N_CP_poles);
    C_3_cp=MALLOC1D(C_3_cp,N_CP_poles);
    C_4_cp=MALLOC1D(C_4_cp,N_CP_poles);
    C_5_cp=MALLOC1D(C_5_cp,N_CP_poles);
    C_cp=MALLOC1D(C_cp,N_CP_poles);

    err=  cudaMalloc(&d_1_ddev,N_drude_poles*sizeof(real));
    err=  cudaMalloc(&d_2_ddev,N_drude_poles*sizeof(real));
    err=  cudaMalloc(&d_3_ddev,N_drude_poles*sizeof(real));
    err=  cudaMalloc(&d_4_ddev,N_drude_poles*sizeof(real));
    err=  cudaMalloc(&d_5_ddev,N_drude_poles*sizeof(real));
    err=  cudaMalloc(&d_NLdev,N_drude_poles*sizeof(real));
    if( cudaSuccess != err)
    {
        printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
        exit(-1);
    }
    printf("HERE\n");

    d_1_d=MALLOC1D(d_1_d,N_drude_poles);
    d_2_d=MALLOC1D(d_2_d,N_drude_poles);
    d_3_d=MALLOC1D(d_3_d,N_drude_poles);
    d_4_d=MALLOC1D(d_4_d,N_drude_poles);
    d_5_d=MALLOC1D(d_5_d,N_drude_poles);
    d_d=MALLOC1D(d_d,N_drude_poles);
   d_NL = MALLOC1D(d_NL,N_drude_poles);

    psi_L = MALLOC1D(psi_L,N_lorentz_poles);
    psi_HD = MALLOC1D(psi_HD,N_drude_poles);
    alpha_HD1 = MALLOC1D(alpha_HD1,N_drude_poles);
    alpha_HD2 = MALLOC1D(alpha_HD2,N_drude_poles);
    alpha_L = MALLOC1D(alpha_L,N_lorentz_poles);
    eta_L = MALLOC1D(eta_L,N_lorentz_poles);
    eta_HD = MALLOC1D(eta_HD,N_drude_poles);

    w_D=MALLOC1D(w_D,N_drude_poles);
    gamma_d=MALLOC1D(gamma_d,N_drude_poles);

    A_cp=MALLOC1D(A_cp,N_CP_poles);
    OMEGA_cp=MALLOC1D(OMEGA_cp,N_CP_poles);
    phi_cp=MALLOC1D(phi_cp,N_CP_poles);
    GAMMA_cp=MALLOC1D(GAMMA_cp,N_CP_poles);

    d_eps_L = MALLOC1D(d_eps_L,N_lorentz_poles);
    delta_L = MALLOC1D(delta_L,N_lorentz_poles);
    omg_L = MALLOC1D(omg_L,N_lorentz_poles);
    printf("HERE\n");



    a_0_cp=MALLOC1D(a_0_cp,N_CP_poles);
    a_1_cp=MALLOC1D(a_1_cp,N_CP_poles);
    b_0_cp=MALLOC1D(b_0_cp,N_CP_poles);
    b_1_cp=MALLOC1D(b_1_cp,N_CP_poles);
    b_2_cp=MALLOC1D(b_2_cp,N_CP_poles);


    err=  cudaMalloc(&Px_cpdev,extentCP);
    err=  cudaMalloc(&Px_cp_ndev,extentCP);
    err=  cudaMalloc(&Px_cp_n_1dev,extentCP);
  err=    cudaMalloc(&Py_cpdev,extentCP);
    err=  cudaMalloc(&Py_cp_ndev,extentCP);
    err=  cudaMalloc(&Py_cp_n_1dev,extentCP);
    err=  cudaMalloc(&Pz_cpdev,extentCP);
    err=  cudaMalloc(&Pz_cp_ndev,extentCP);
    err=  cudaMalloc(&Pz_cp_n_1dev,extentCP);
    if( cudaSuccess != err)
    {
        printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
        exit(-1);
    }
    printf("HERE\n");

    Px_cp=MALLOC4D_Complex(Px_cp,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Px_cp_n=MALLOC4D_Complex(Px_cp_n,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Px_cp_n_1=MALLOC4D_Complex(Px_cp_n_1,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Py_cp=MALLOC4D_Complex(Py_cp,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Py_cp_n=MALLOC4D_Complex(Py_cp_n,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Py_cp_n_1=MALLOC4D_Complex(Py_cp_n_1,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Pz_cp=MALLOC4D_Complex(Pz_cp,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Pz_cp_n=MALLOC4D_Complex(Pz_cp_n,NCELLX,NCELLY,NCELLZ,N_CP_poles);
    Pz_cp_n_1=MALLOC4D_Complex(Pz_cp_n_1,NCELLX,NCELLY,NCELLZ,N_CP_poles);

    err= cudaMalloc(&Px_ddev,extentD);
    err= cudaMalloc(&Px_d_ndev,extentD);
    err= cudaMalloc(&Px_d_n_1dev,extentD);
    err= cudaMalloc(&Py_ddev,extentD);
    err= cudaMalloc(&Py_d_ndev,extentD);
    err= cudaMalloc(&Py_d_n_1dev,extentD);
    err= cudaMalloc(&Pz_ddev,extentD);
    err= cudaMalloc(&Pz_d_ndev,extentD);
    err= cudaMalloc(&Pz_d_n_1dev,extentD);
    err= cudaMalloc(&Px_d_n_2dev,extentD);
    err= cudaMalloc(&Py_d_n_2dev,extentD);
    err= cudaMalloc(&Pz_d_n_2dev,extentD);

    if( cudaSuccess != err)
    {
        printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
        exit(-1);
    }
    printf("HERE\n");


    Px_d=MALLOC4D_Complex(Px_d,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Px_d_n=MALLOC4D_Complex(Px_d_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Px_d_n_1=MALLOC4D_Complex(Px_d_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Py_d=MALLOC4D_Complex(Py_d,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Py_d_n=MALLOC4D_Complex(Py_d_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Py_d_n_1=MALLOC4D_Complex(Py_d_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Pz_d=MALLOC4D_Complex(Pz_d,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Pz_d_n=MALLOC4D_Complex(Pz_d_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Pz_d_n_1=MALLOC4D_Complex(Pz_d_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);

    Px_NL=MALLOC4D_Complex(Px_NL,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Px_NL_n=MALLOC4D_Complex(Px_NL_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Px_NL_n_1=MALLOC4D_Complex(Px_NL_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Py_NL=MALLOC4D_Complex(Py_NL,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Py_NL_n=MALLOC4D_Complex(Py_NL_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Py_NL_n_1=MALLOC4D_Complex(Py_NL_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Pz_NL=MALLOC4D_Complex(Pz_NL,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Pz_NL_n=MALLOC4D_Complex(Pz_NL_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Pz_NL_n_1=MALLOC4D_Complex(Pz_NL_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);

    Jx_NL=MALLOC4D_Complex(Jx_NL,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jx_NL_n=MALLOC4D_Complex(Jx_NL_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jx_NL_n_1=MALLOC4D_Complex(Jx_NL_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jy_NL=MALLOC4D_Complex(Jy_NL,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jy_NL_n=MALLOC4D_Complex(Jy_NL_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jy_NL_n_1=MALLOC4D_Complex(Jy_NL_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jz_NL=MALLOC4D_Complex(Jz_NL,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jz_NL_n=MALLOC4D_Complex(Jz_NL_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
    Jz_NL_n_1=MALLOC4D_Complex(Jz_NL_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);

    Jx_Lo=MALLOC4D_Complex(Jx_Lo,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jx_Lo_n=MALLOC4D_Complex(Jx_Lo_n,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jx_Lo_n_1=MALLOC4D_Complex(Jx_Lo_n_1,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jy_Lo=MALLOC4D_Complex(Jy_Lo,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jy_Lo_n=MALLOC4D_Complex(Jy_Lo_n,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jy_Lo_n_1=MALLOC4D_Complex(Jy_Lo_n_1,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jz_Lo=MALLOC4D_Complex(Jz_Lo,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jz_Lo_n=MALLOC4D_Complex(Jz_Lo_n,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);
    Jz_Lo_n_1=MALLOC4D_Complex(Jz_Lo_n_1,NCELLX,NCELLY,NCELLZ,N_lorentz_poles);

    printf("HERE\n");


    //SILVER
  if(material == 2){
    eps_inf=1.4447;
    w_D[0]=1.3280e16;
    gamma_d[0]=9.1269e13;

    A_cp[0]=-1.5951;
    phi_cp[0]=3.1288;
    OMEGA_cp[0]=8.2749e15;
    GAMMA_cp[0]=5.1770e15;

    A_cp[1]=0.25261;
    phi_cp[1]=-1.5066;
    OMEGA_cp[1]=6.1998e15;
    GAMMA_cp[1]=5.4126e14;

    if(NONLOCAL == 1) NONLOC = FERMI_VELOCITY*NL_COEFF;
    else NONLOC = 0;

    if(Hydrodynamics == 1){
      N_EQ = 5.541318196761572e28;
      printf("NEQ=%e\n",N_EQ);
       NONLOC = FERMI_VELOCITY*NL_COEFF;
    }

  }
    //Gold
  if(material == 3){
    eps_inf=1.1431;
    w_D[0]=1.3202e16;
    gamma_d[0]=1.0805e14;

    A_cp[0]=0.26698;
    phi_cp[0]=-1.2371;
    OMEGA_cp[0]=3.8711e15;
    GAMMA_cp[0]=4.4642e14;

    A_cp[1]=3.0834;
    phi_cp[1]=-1.0968;
    OMEGA_cp[1]=4.1684e15;
    GAMMA_cp[1]=2.3555e15;

    if(NONLOCAL == 1) NONLOC = FERMI_VELOCITY*NL_COEFF;
    else NONLOC = 0;
    if(Hydrodynamics == 1){
      N_EQ = 5.476415562682574e28;
       NONLOC = FERMI_VELOCITY*NL_COEFF;
    }
  }

  //Nonlocal SILVER
if(material == 4){
  // eps_inf=1.4447;
  // w_D[0]=1.3280e16;
  // gamma_d[0]=9.1269e13;
  //
  // A_cp[0]=-1.5951;
  // phi_cp[0]=3.1288;
  // OMEGA_cp[0]=8.2749e15;
  // GAMMA_cp[0]=5.1770e15;
  //
  // A_cp[1]=0.25261;
  // phi_cp[1]=-1.5066;
  // OMEGA_cp[1]=6.1998e15;
  // GAMMA_cp[1]=5.4126e14;
  //
  // NONLOC = FERMI_VELOCITY*NL_COEFF;
}
  //Nonlocal Gold
if(material == 5){


  // eps_inf=1.1431;
  // w_D[0]=1.3202e16;
  // gamma_d[0]=1.0805e14;
  //
  // A_cp[0]=0.26698;
  // phi_cp[0]=-1.2371;
  // OMEGA_cp[0]=3.8711e15;
  // GAMMA_cp[0]=4.4642e14;
  //
  // A_cp[1]=3.0834;
  // phi_cp[1]=-1.0968;
  // OMEGA_cp[1]=4.1684e15;
  // GAMMA_cp[1]=2.3555e15;
  //
  //
  //
  // NONLOC = FERMI_VELOCITY*NL_COEFF;
}
if(material == 6){

}
if(material == 7){
  eps_inf = 3.559;
  d_eps_L[0] = 2.912;
  d_eps_L[1] = 1.272;
  delta_L[0] = (2.0*PI/PLANKS)*1.541;
  delta_L[1] = (2.0*PI/PLANKS)*0.525;
  omg_L[0] = (2.0*PI/PLANKS)*4.693;
  omg_L[1] = (2.0*PI/PLANKS)*3.112;
  // delta_L[0] = 1.541;
  // delta_L[1] = 0.525;
  // omg_L[0] = 4.693;
  // omg_L[1] = 3.112;

  w_D[0] = (2.0*PI/PLANKS)*8.812;
  gamma_d[0] = (2.0*PI/PLANKS)*0.0752;
  // w_D[0] = 8.812;
  // gamma_d[0] = 0.0752;
  if(NONLOCAL==1) NONLOC = FERMI_VELOCITY*NL_COEFF;
  else NONLOC = 0;
}
printf("HERE\n");


  if(material >= 2 && material <=5){

    //set-up update coefficients
    for(i=0;i<N_CP_poles;i++){

        a_0_cp[i]=2.0*ep0*A_cp[i]*OMEGA_cp[i]*(OMEGA_cp[i]*cos(phi_cp[i])-GAMMA_cp[i]*sin(phi_cp[i]));
        a_1_cp[i]=-2.0*ep0*A_cp[i]*OMEGA_cp[i]*sin(phi_cp[i]);
        b_0_cp[i]=GAMMA_cp[i]*GAMMA_cp[i]+OMEGA_cp[i]*OMEGA_cp[i];
        b_1_cp[i]=2.0*GAMMA_cp[i];
        b_2_cp[i]=1.0;
        C_cp[i]=b_2_cp[i]/(dt*dt)+b_1_cp[i]/(2.0*dt)+b_0_cp[i]/4.0;
        C_1_cp[i]=(2*b_2_cp[i]/(dt*dt)-b_0_cp[i]/2.0)/C_cp[i];
        C_2_cp[i]=(b_1_cp[i]/(2.0*dt)-b_2_cp[i]/(dt*dt)-b_0_cp[i]/4.0)/C_cp[i];
        C_3_cp[i]=(a_0_cp[i]/4.0+a_1_cp[i]/(2.0*dt))/C_cp[i];
        C_4_cp[i]=a_0_cp[i]/(2.0*C_cp[i]);
        C_5_cp[i]=(a_0_cp[i]/4.0-a_1_cp[i]/(2.0*dt))/C_cp[i];

    }


    for(i=0;i<N_drude_poles;i++){
      if(Hydrodynamics==0){
        d_d[i]=1.0/(dt*dt)+gamma_d[i]/(2.0*dt);
        d_1_d[i]=2.0/(d_d[i]*dt*dt);
        d_2_d[i]=(gamma_d[i]/(2.0*dt)-1.0/(dt*dt))/d_d[i];
        d_3_d[i]=d_5_d[i]=ep0*w_D[i]*w_D[i]/(4.0*d_d[i]);
        d_4_d[i]=ep0*w_D[i]*w_D[i]/(2.0*d_d[i]);
        d_NL[i] = NONLOC*NONLOC/d_d[i];

        printf("%e,%e,%e,%e,%e,%e,%e,%e,%e\n",e0,me, C_E,d_d[i],d_1_d[i], d_2_d[i], d_3_d[i], d_4_d[i], d_5_d[i]);
        // printf("%e,%e,%e,%e,%e\n", cimag(d_1_d[i]), cimag(d_2_d[i]), cimag(d_3_d[i]), cimag(d_4_d[i]), cimag(d_5_d[i]));
        printf("%e\n",d_NL[i]);


        for(n=0;n<N_drude_poles;n++){
            something+=d_3_d[n];
        }
        for(n=0;n<N_CP_poles;n++){
            something+=C_3_cp[n];
        }
        C_E=(ep0*eps_inf+something);


        something=0.0;
        for(n=0;n<N_drude_poles;n++){
            something+=d_4_d[n];
        }
        for(n=0;n<N_CP_poles;n++){
            something+=C_4_cp[n];
        }
        C_E_1=-something+ep0*eps_inf;

        something=0.0;
        for(n=0;n<N_drude_poles;n++){
            something+=d_5_d[n];
        }
        for(n=0;n<N_CP_poles;n++){
            something+=C_5_cp[n];
        }
        C_E_2=something;

      }
      else{
        d_d[i] = 1.0/(dt*dt) + gamma_d[i]/(2.0*dt);
        d_1_d[i]=2.0/(d_d[i]*dt*dt);
        d_2_d[i]=(gamma_d[i]/(2.0*dt)-1.0/(dt*dt))/d_d[i];
        d_3_d[i]= -1.0*e0/(me)/d_d[i]/dt;
        d_4_d[i]= 1.0*e0/(me)/d_d[i]/dt;
        d_5_d[i] = -1.0/d_d[i];
        d_NL[i] = pow(N_EQ,1.0/3.0)*NONLOC*NONLOC/d_d[i];
        printf("%e,%e,%e,%e,%e,%e,%e,%e,%e\n",e0,me, C_E,d_d[i],d_1_d[i], d_2_d[i], d_3_d[i], d_4_d[i], d_5_d[i]);
        // printf("%e,%e,%e,%e,%e\n", cimag(d_1_d[i]), cimag(d_2_d[i]), cimag(d_3_d[i]), cimag(d_4_d[i]), cimag(d_5_d[i]));
        printf("%e\n",d_NL[i]);

        for(n=0;n<N_CP_poles;n++){
            something+=C_3_cp[n];
        }
        C_E=(ep0*eps_inf+something);


        something=0.0;

        for(n=0;n<N_CP_poles;n++){
            something+=C_4_cp[n];
        }
        C_E_1=-something+ep0*eps_inf;

        something=0.0;

        for(n=0;n<N_CP_poles;n++){
            something+=C_5_cp[n];
        }
        C_E_2=something;

      }

      printf("%d,%e,%e,%e,%e,%e,%e,%e,%e,%e\n",i,e0,me, C_E,d_d[i],d_1_d[i], d_2_d[i], d_3_d[i], d_4_d[i], d_5_d[i]);
      // printf("%e,%e,%e,%e,%e\n", cimag(d_1_d[i]), cimag(d_2_d[i]), cimag(d_3_d[i]), cimag(d_4_d[i]), cimag(d_5_d[i]));
      printf("%e\n",d_NL[i]);
    }

}



if(material >= 6){
    for(n=0;n<N_lorentz_poles;n++){
      alpha_L[n] = (2.0 - omg_L[n]*omg_L[n]*dt*dt)/(1.0 + delta_L[n]*dt);
      psi_L[n] = -1.0*(1.0 - delta_L[n]*dt)/(1.0 + delta_L[n]*dt);
      eta_L[n] = ep0*d_eps_L[n]*omg_L[n]*omg_L[n]*dt*dt/(1.0 + delta_L[n]*dt);
    }
    for(n=0;n<N_drude_poles;n++){
      alpha_HD1[n] = 4.0/(2.0 + gamma_d[n]*dt);
      alpha_HD2[n] = 2.0*dt*dt*NONLOC*NONLOC/(2.0 + gamma_d[n]*dt);
      psi_HD[n] = -1.0*(2.0 - gamma_d[n]*dt)/(2.0 + gamma_d[n]*dt);
      eta_HD[n] = 2.0*ep0*w_D[n]*w_D[n]*dt*dt/(2.0 + gamma_d[n]*dt);
    }
    C1_NL = eps_inf*ep0/dt;
    C2_NL = (1.0/(4.0*dt))*(eta_L[0]+eta_L[1]+eta_HD[0]);
    // printf("CNLs: %f\t%f\n",C1_NL,C2_NL);
    // printf("etas: %f\t%f\t%f\n",eta_L[0]/ep0,eta_L[1]/ep0,eta_HD[0]/ep0);
    printf("%e\n",alpha_HD2[0]);

}
printf("HERE\n");

    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                for(n=0;n<N_CP_poles;n++){
                    Px_cp[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Px_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Px_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Py_cp[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Py_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Py_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Pz_cp[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Pz_cp_n[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;
                    Pz_cp_n_1[FourDMap(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=0.0;

            }
          }
        }
    }
    printf("HERE1\n");

    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
           for(k=0;k<NCELLZ;k++){
             for(n=0;n<N_drude_poles;n++){
                Px_NL[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Px_NL_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Px_NL_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Py_NL[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Py_NL_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Py_NL_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Pz_NL[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Pz_NL_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Pz_NL_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;

                Px_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Px_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Px_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Py_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Py_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Py_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Pz_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Pz_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Pz_d_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;

                Jx_NL[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jx_NL_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jx_NL_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jy_NL[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jy_NL_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jy_NL_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jz_NL[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jz_NL_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
                Jz_NL_n_1[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=0.0;
            }
            }
          }
        }
        for(i=0;i<NCELLX;i++){
            for(j=0;j<NCELLY;j++){
               for(k=0;k<NCELLZ;k++){
                 for(n=0;n<N_lorentz_poles;n++){

                    Jx_Lo[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jx_Lo_n[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jx_Lo_n_1[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jy_Lo[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jy_Lo_n[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jy_Lo_n_1[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jz_Lo[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jz_Lo_n[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                    Jz_Lo_n_1[FourDMap(i,j,k,n,N_lorentz_poles,NCELLZ,NCELLY)]=0.0;
                }
                }
              }
            }


            printf("HERE1\n");


    FREE1D(w_D);
    FREE1D(gamma_d);
    FREE1D(A_cp);
    FREE1D(OMEGA_cp);
    FREE1D(phi_cp);
    FREE1D(GAMMA_cp);

    FREE1D(a_0_cp);
    FREE1D(a_1_cp);
    FREE1D(b_0_cp);
    FREE1D(b_1_cp);
    FREE1D(b_2_cp);
    printf("CE=%e\n",C_E);

cudaMemcpy(d_1_ddev,d_1_d,N_drude_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(d_2_ddev,d_2_d,N_drude_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(d_3_ddev,d_3_d,N_drude_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(d_4_ddev,d_4_d,N_drude_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(d_5_ddev,d_5_d,N_drude_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(d_NLdev,d_NL,N_drude_poles*sizeof(real),cudaMemcpyHostToDevice);

cudaMemcpy(C_1_cpdev,C_1_cp,N_CP_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(C_2_cpdev,C_2_cp,N_CP_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(C_3_cpdev,C_3_cp,N_CP_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(C_4_cpdev,C_4_cp,N_CP_poles*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(C_5_cpdev,C_5_cp,N_CP_poles*sizeof(real),cudaMemcpyHostToDevice);

    err=  cudaMemcpy(Px_cpdev,Px_cp,extentCP,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Px_cp_ndev,Px_cp_n,extentCP,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Px_cp_n_1dev,Px_cp_n_1,extentCP,cudaMemcpyHostToDevice);

    err=  cudaMemcpy(Py_cpdev,Py_cp,extentCP,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Py_cp_ndev,Py_cp_n,extentCP,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Py_cp_n_1dev,Py_cp_n_1,extentCP,cudaMemcpyHostToDevice);


    err=  cudaMemcpy(Pz_cpdev,Pz_cp,extentCP,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Pz_cp_ndev,Pz_cp_n,extentCP,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Pz_cp_n_1dev,Pz_cp_n_1,extentCP,cudaMemcpyHostToDevice);

    err=  cudaMemcpy(Px_ddev,Px_d,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Px_d_ndev,Px_d_n,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Px_d_n_1dev,Px_d_n_1,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Px_d_n_2dev,Px_NL,extentD,cudaMemcpyHostToDevice);


    err=  cudaMemcpy(Py_ddev,Py_d,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Py_d_ndev,Py_d_n,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Py_d_n_1dev,Py_d_n_1,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Py_d_n_2dev,Py_NL,extentD,cudaMemcpyHostToDevice);


    err=  cudaMemcpy(Pz_ddev,Pz_d,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Pz_d_ndev,Pz_d_n,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Pz_d_n_1dev,Pz_d_n_1,extentD,cudaMemcpyHostToDevice);
    err=  cudaMemcpy(Pz_d_n_2dev,Pz_NL,extentD,cudaMemcpyHostToDevice);

    // err=  cudaMemcpy(Py_cpdev,extentCP);
    // err=  cudaMemcpy(Py_cp_ndev,extentCP);
    // err=  cudaMemcpy(Py_cp_n_1dev,extentCP);
    // err=  cudaMemcpy(Pz_cpdev,extentCP);
    // err=  cudaMemcpy(Pz_cp_ndev,extentCP);
    // err=  cudaMemcpy(Pz_cp_n_1dev,extentCP);
    //
    // err= cudaMemcpy(Px_ddev,extentD);
    // err= cudaMemcpy(Px_d_ndev,extentD);
    // err= cudaMemcpy(Px_d_n_1dev,extentD);
    // err= cudaMemcpy(Py_ddev,extentD);
    // err= cudaMemcpy(Py_d_ndev,extentD);
    // err= cudaMemcpy(Py_d_n_1dev,extentD);
    // err= cudaMemcpy(Pz_ddev,extentD);
    // err= cudaMemcpy(Pz_d_ndev,extentD);
    // err= cudaMemcpy(Pz_d_n_1dev,extentD);
}



//1D memory allocation
real* MALLOC1D(real *grid, int SIZE){
  grid=(real *)malloc((SIZE)*sizeof(real));
  return grid;
}

real* MALLOC1D_double(real *grid, int SIZE){
  grid=(real *)malloc((SIZE)*sizeof(real));
  return grid;
}

comp* MALLOC1D_Complex(comp *grid, int SIZE){
  grid=(comp *)malloc((SIZE)*sizeof(comp));
  return grid;
}
real2* MALLOC1D_Real2(real2 *grid, int SIZE){
  grid=(real2 *)malloc((SIZE)*sizeof(real2));
  return grid;
}

double complex* MALLOC1D_Complex2(double complex *grid, int SIZE){
  grid=(double complex *)malloc((SIZE)*sizeof(double complex));
  return grid;
}

//2-D Memory Allocation
real* MALLOC2D(real *grid, int sizeX, int sizeZ){
    int i;
    grid=(real *)malloc((sizeX*sizeZ)*sizeof(real));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");
    }
    return grid;

}

real* MALLOC2D_double(real *grid, int sizeX, int sizeZ){
    int i;
    grid=(real *)malloc((sizeX*sizeZ)*sizeof(real));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }

    return grid;

}

comp* MALLOC2D_Complex(comp *grid, int sizeX, int sizeZ){
    int i;
    grid=(comp *)malloc((sizeX*sizeZ)*sizeof(comp));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }
    return grid;

}



// 3-D array memory allocation
real* MALLOC3D(real *grid, int sizeX, int sizeY, int sizeZ){
    int i,j;
    grid=(real *)malloc((sizeX*sizeY*sizeZ)*sizeof(real));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }

    return grid;

}

real* MALLOC3D_double(real *grid, int sizeX, int sizeY, int sizeZ){
    int i,j;
    grid=(real *)malloc((sizeX*sizeY*sizeZ)*sizeof(real));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }

    return grid;

}
real2* MALLOC3D_Real2(real2 *grid, int sizeX, int sizeY, int sizeZ){
    int i,j;
    grid=(real2 *)malloc((sizeX*sizeY*sizeZ)*sizeof(real2));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }

    return grid;

}

comp* MALLOC3D_Complex(comp *grid, int sizeX, int sizeY,int sizeZ){
    int i,j;
    grid=(comp *)malloc((sizeX*sizeY*sizeZ)*sizeof(comp));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }
    return grid;

}

double complex * MALLOC3D_Complex2(double complex *grid, int sizeX, int sizeY,int sizeZ){
    int i,j;
    grid=(double complex  *)malloc((sizeX*sizeY*sizeZ)*sizeof(double complex ));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }
    return grid;
}

int * MALLOC3D_int(int *grid,int sizeX, int sizeY, int sizeZ){
  int i,j;
  grid=(int *)malloc((sizeX*sizeY*sizeZ)*sizeof(int));
  if(!grid){
      printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

  }
  return grid;

}

//4-D Memory Allocation

real* MALLOC4D(real *grid, int sizeX, int sizeY, int sizeZ,int size4){
    int i,j,k;
    grid=(real *)malloc((sizeX*sizeY*sizeZ*size4)*sizeof(real));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }
      return grid;
}


real* MALLOC4D_double(real *grid, int sizeX, int sizeY, int sizeZ,int size4){
    int i,j,k;
    grid=(real *)malloc((sizeX*sizeY*sizeZ*size4)*sizeof(real));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }
      return grid;
}

comp* MALLOC4D_Complex(comp *grid, int sizeX, int sizeY,int sizeZ,int size4){
    int i,j,k;
    grid=(comp *)malloc((sizeX*sizeY*sizeZ*size4)*sizeof(comp));
    if(!grid){
        printf("\n\nERROR IN MEMORY ALLOCATION\n\n");

    }

    return grid;

}


//free arrays
void FREE1D(real *grid){
    free(grid);
}

void FREE1D_double(double *grid){
    free(grid);
}

void FREE1D_Complex(comp *grid){
    free(grid);
}

void FREE1D_Complex2(double complex *grid){
    free(grid);
}

void FREE2D(real *grid,int sizeX){
        free(grid);
    return;
}

void FREE2D_double(double *grid,int sizeX){

        free(grid);
    return;
}

void FREE2D_Complex(comp *grid,int sizeX){

    free(grid);
    return;
}

void FREE3D(real *grid,int sizeX,int sizeY){

    free(grid);
    return;
}


void FREE3D_double(double *grid,int sizeX,int sizeY){

    free(grid);
    return;
}


void FREE3D_Complex(comp *grid,int sizeX,int sizeY){
    int i,j;

    free(grid);
    return;
}

void FREE3D_Complex2(double complex *grid,int sizeX,int sizeY){

    free(grid);
    return;
}



void FREE4D(real *grid,int sizeX,int sizeY,int sizeZ){
    free(grid);
    return;
}

void FREE4D_double(double *grid,int sizeX,int sizeY,int sizeZ){

    free(grid);
    return;
}


void FREE4D_Complex(comp *grid,int sizeX,int sizeY,int sizeZ){

    free(grid);
    return;
}

//make the vector zero
real* ZERO_VECTORS2D(real *grid,int SizeX,int SizeY){
    int i,j;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){

                grid[TwoDMap(i,j,SizeY)]=0.0;

        }
    }
    return grid;
}

real* ZERO_VECTORS2D_double(real *grid,int SizeX,int SizeY){
    int i,j;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){

                grid[TwoDMap(i,j,SizeY)]=0.0;

        }
    }
    return grid;
}

real* ZERO_VECTORS3D(real *grid,int SizeX,int SizeY,int SizeZ){
    int i,j,k;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){
           for(k=0;k<SizeZ;k++) {
              grid[ThreeDMap(i,j,k,SizeZ,SizeY)]=0.0;
           }

        }
    }
    return grid;
}

real* ZERO_VECTORS3D_double(real *grid,int SizeX,int SizeY,int SizeZ){
    int i,j,k;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){
           for(k=0;k<SizeZ;k++) {
              grid[ThreeDMap(i,j,k,SizeZ,SizeY)]=0.0;
           }

        }
    }
    return grid;
}

real2* ZERO_VECTORS3D_Real2 (real2 *grid,int SizeX,int SizeY,int SizeZ){
    int i,j,k;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){
           for(k=0;k<SizeZ;k++) {
              grid[ThreeDMap(i,j,k,SizeZ,SizeY)]=0.0;
           }

        }
    }
    return grid;
}

comp * ZERO_VECTORS2D_Complex(comp *grid,int SizeX,int SizeY){
    int i,j;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){
                grid[TwoDMap(i,j,SizeY)]=0.0;

        }
    }
    return grid;
}

comp * ZERO_VECTORS3D_Complex(comp *grid,int SizeX,int SizeY,int SizeZ){
    int i,j,k;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){
                for(k=0;k<SizeZ;k++){
                    grid[ThreeDMap(i,j,k,SizeZ,SizeY)]=0.0;
                }
        }
    }
    return grid;
}

double complex * ZERO_VECTORS3D_Complex2(double complex *grid,int SizeX,int SizeY,int SizeZ){
    int i,j,k;
    for (i=0;i<SizeX;i++){
        for(j=0;j<SizeY;j++){
                for(k=0;k<SizeZ;k++){
                    grid[ThreeDMap(i,j,k,SizeZ,SizeY)]=0.0;
                }
        }
    }
    return grid;
}


real* ZERO_VECTORS1D(real *grid){
    int i;
    for(i=0;i<inc_Length;i++){
        grid[i]=0;
    }
    return grid;
}

real* ZERO_VECTORS1D_double(real *grid){
    int i;
    for(i=0;i<inc_Length;i++){
        grid[i]=0;
    }
    return grid;
}

comp* ZERO_VECTORS1D_Complex(comp *grid,int length){
    int i;
    for(i=0;i<length;i++){
        grid[i]=0;
    }
    return grid;
}

double complex* ZERO_VECTORS1D_Complex2(double complex *grid,int length){
    int i;
    for(i=0;i<length;i++){
        grid[i]=0;
    }
    return grid;
}



//define the electric conductivity
void DEF_SIGMA_E(void){
    int i,j,k;
    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                sigma_e[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
            }
        }
    }

    return;
}

//define the magnetic conductivity
void DEF_SIGMA_M(void){
    int i,j,k;
    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                sigma_m[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=0;
            }
        }
    }

    return;
}

//define the epsilon
void DEF_EPS(void){
    int i,j,k;
    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                    eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=EP0;

            }
        }
    }

    return;
}

//define the mu
void DEF_MU(void){
    int i,j,k;
    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                mu[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=MU0;
            }
        }
    }

    return;
}

void DIELECTRIC_SLAB(void){
    int i,j,k;
    real eps_r=10;
    //int location=inc_plane+ceil(NCELLZ/2)-2*NcpmlZ;
    int location=cpml_N_Z+200;
  //  int SLAB_SIZE=31;
    for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=location;k<NCELLZ;k++){
                   eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*eps_r;
            }
        }
    }
}

//define the update coefficients for E-field update equations on the H field components
real * DEF_UPDATE_COEFF_EonH(real *grid){
      int i,j,k;
      for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                grid[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(dt/(z0*eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)]))/((1+(dt*sigma_e[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(2*eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)])));
            }
        }
      }
    return grid;
}
//define the update coefficients for E-field update equations on the E field components
real * DEF_UPDATE_COEFF_EonE(real *grid){
      int i,j,k;
      for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                grid[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(1-(dt*sigma_e[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(2*eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)]))/(1+(dt*sigma_e[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(2*eps[ThreeDMap(i,j,k,NCELLZ,NCELLY)]));
            }
        }
      }
    return grid;
}

//define the update coefficients for H-field update equations on the E field components
real *DEF_UPDATE_COEFF_HonE(real *grid){
      int i,j,k;
      for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                grid[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(z0*dt/mu[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(1+(dt*sigma_m[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(2*mu[ThreeDMap(i,j,k,NCELLZ,NCELLY)]));
            }
        }
      }
    return grid;
}
//define the update coefficients for H-field update equations on the H field components
real *DEF_UPDATE_COEFF_HonH(real *grid){
      int i,j,k;
      for(i=0;i<NCELLX;i++){
        for(j=0;j<NCELLY;j++){
            for(k=0;k<NCELLZ;k++){
                grid[ThreeDMap(i,j,k,NCELLZ,NCELLY)]=(1-(dt*sigma_m[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(2*mu[ThreeDMap(i,j,k,NCELLZ,NCELLY)]))/(1+(dt*sigma_m[ThreeDMap(i,j,k,NCELLZ,NCELLY)])/(2*mu[ThreeDMap(i,j,k,NCELLZ,NCELLY)]));
            }
        }
      }
    return grid;
}

void SETUP_t_inc(void){
  int i;
  real n_inc = 1.0;
  real d_silver = (dispersive_slab - inc_plane)* dz;
  real arg, arg_sq;
  for(i=0;i<NUM_freq; i++){
    if(freq[i]>f_min){
      arg = (c0*k_x/(2*pi*freq[i]));
      arg_sq = arg*arg;
      arg = pow(1-arg_sq,0.5);
      t_inc[i] =  (d_silver*n_inc/(c0*arg));
      t_inc[i] =(double) ceil (t_inc[i] / dt + delay + width);
      // if(c0/freq[i] >= 675e-9){
      //     t_inc[i] += 1500;
      // }
      t_inc[i] = 0;
    }
    else{
      t_inc[i] = 1.0;
    }
  }
}
//
// void FT_source_calc(void){
//
//   int i,j,k,n,w;
//
//   int Tend_inc = Tend;
//   int NUM_freq_inc = Tend;
//
//   FILE *Input;
//   FILE *Input_w;
//
//   Input = fopen("Input_Pulse.txt", "w");
//   Input_w = fopen("Input_Pulse_w.txt", "w");
//
//   double max_omega = 2 * pi * freq[0];
//   double min_omega = 1;
//   double min_freq = min_omega/(2*pi);
//   double max_freq = max_omega/(2*pi);
//   double Tend_real = width_real*11;
//   double *omega2;
//   double *freq2;
//   double *time2;
//   double dw_2;
//
//   Ex_t_FT=MALLOC1D_Complex(Ex_t_FT,Tend_inc);
//   Hx_t_FT=MALLOC1D_Complex(Hx_t_FT,Tend_inc);
//   Ey_t_FT=MALLOC1D_Complex(Ey_t_FT,Tend_inc);
//   Hy_t_FT=MALLOC1D_Complex(Hy_t_FT,Tend_inc);
//   Ez_t_FT=MALLOC1D_Complex(Ez_t_FT,Tend_inc);
//   Hz_t_FT=MALLOC1D_Complex(Hz_t_FT,Tend_inc);
//
//   Ex_w_FT=MALLOC1D_Complex(Ex_w_FT,NUM_freq_inc);
//   Hx_w_FT=MALLOC1D_Complex(Hx_w_FT,NUM_freq_inc);
//   Ey_w_FT=MALLOC1D_Complex(Ey_w_FT,NUM_freq_inc);
//   Hy_w_FT=MALLOC1D_Complex(Hy_w_FT,NUM_freq_inc);
//   Ez_w_FT=MALLOC1D_Complex(Ez_w_FT,NUM_freq_inc);
//   Hz_w_FT=MALLOC1D_Complex(Hz_w_FT,NUM_freq_inc);
//
// omega2 = MALLOC1D_double(omega2, NUM_freq_inc);
// freq2 = MALLOC1D_double(freq2,NUM_freq_inc);
// time2 = MALLOC1D_double(time2,Tend_inc);
//
//   //printf("%e\n%e",freq[0], freq[NUM_freq_inc-1]);
//
//   dw_2 = (max_omega - min_omega)/(NUM_freq_inc - 1);
//
//
//   Hx_t_FT = ZERO_VECTORS1D_Complex(Hx_t_FT,Tend_inc);
//   Hy_t_FT = ZERO_VECTORS1D_Complex(Hy_t_FT,Tend_inc);
//   Hz_t_FT = ZERO_VECTORS1D_Complex(Hz_t_FT,Tend_inc);
//
//   Ex_t_FT = ZERO_VECTORS1D_Complex(Ex_t_FT,Tend_inc);
//   Ey_t_FT = ZERO_VECTORS1D_Complex(Ey_t_FT,Tend_inc);
//   Ez_t_FT = ZERO_VECTORS1D_Complex(Ez_t_FT,Tend_inc);
//
//   Hx_w_FT = ZERO_VECTORS1D_Complex(Hx_w_FT,NUM_freq_inc);
//   Hy_w_FT = ZERO_VECTORS1D_Complex(Hy_w_FT,NUM_freq_inc);
//   Hz_w_FT = ZERO_VECTORS1D_Complex(Hz_w_FT,NUM_freq_inc);
//
//   Ex_w_FT = ZERO_VECTORS1D_Complex(Ex_w_FT,NUM_freq_inc);
//   Ey_w_FT = ZERO_VECTORS1D_Complex(Ey_w_FT,NUM_freq_inc);
//   Ez_w_FT = ZERO_VECTORS1D_Complex(Ez_w_FT,NUM_freq_inc);
//
//   if(TEz){
//
//      for(n=0;n<Tend_inc;n++){
//         Ey_t_FT[n] = PULSE(n);
//         time2[n] = n*dt;
//       }
//
//       for(w=0;w<NUM_freq_inc;w++){
//         omega2[w] = min_omega + w*(max_omega - min_omega)/(NUM_freq_inc - 1);
//         freq2[w] = omega2[w]/(2*pi);
//       }
//
//
//       for(n=0;n<Tend_inc;n++){
//         for(w=0;w<NUM_freq_inc;w++){
//           Ey_w_FT[w] += Ey_t_FT[n] * cexp(-I*time2[n]*omega2[w])*dt;
//         }
//       }
//
//       for(w=0;w<NUM_freq_inc;w++){
//         Ey_w_FT[w] /= csqrt(2*pi);
//       }
//
//       for(w=0;w<NUM_freq_inc;w++){
//         Hx_w_FT[w] = Ey_w_FT[w]*csqrt((omega2[w]/c0)*(omega2[w]/c0) - (k_rho*k_rho))/(omega2[w]/c0);
//         fprintf(Input_w,"%e \t %e \t %e \n", omega2[w], cabs(Ey_w_FT[w]), cabs(Hx_w_FT[w]));
//       }
//       //printf("here\n");
//
//       for(n=0;n<Tend_inc;n++){
//           for(w=0;w<NUM_freq_inc;w++){
//           Hx_t_FT[n] += Hx_w_FT[w] * cexp(I*time2[n]*omega2[w])*dw_2;
//           Ey_t_FT[n] += Ey_w_FT[w] * cexp(I*time2[n]*omega2[w])*dw_2;
//
//         }
//         Ey_t_FT[n] /= csqrt(4*pi);
//         Hx_t_FT[n] /= csqrt(2*pi);
//         fprintf(Input,"%e \t %e \t %e \t %e\n", time2[n] , creal(Ey_t_FT[n]), creal(Hx_t_FT[n]), creal(PULSE(n)));
//       }
//       //printf("here\n");
//
//   }
//
//   fclose(Input);
//   fclose(Input_w);
//
//   free(omega2);
//   free(freq2);
//   free(time2);
// }
//


// void FT_source_calc(void){
//
//   int i,j,k,n,w;
//
//   // int Tend_inc = Tend;
//   // int NUM_freq_inc = 5000;
//   // if(Tend>250000) {
//   //   Tend_inc =  250000;
//   //   NUM_freq_inc = 5000;
//   // }
//   FILE *Input;
//   FILE *Input_w;
//
//   Input = fopen("Input_Pulse.txt", "w");
//   Input_w = fopen("Input_Pulse_w.txt", "w");
//
//   double Tend_real = Tend_inc*dt;
//   double *omega2;
//   double *freq2;
//   double *time2;
//   comp k_z,k_0;
//   double dw_2 = 2*pi/Tend_real;
//
//   Ex_t_FT=MALLOC1D_Complex(Ex_t_FT,Tend);
//   Hx_t_FT=MALLOC1D_Complex(Hx_t_FT,Tend);
//   Ey_t_FT=MALLOC1D_Complex(Ey_t_FT,Tend);
//   Hy_t_FT=MALLOC1D_Complex(Hy_t_FT,Tend);
//   Ez_t_FT=MALLOC1D_Complex(Ez_t_FT,Tend);
//   Hz_t_FT=MALLOC1D_Complex(Hz_t_FT,Tend);
//
//   Ex_w_FT=MALLOC1D_Complex(Ex_w_FT,NUM_freq_inc);
//   Hx_w_FT=MALLOC1D_Complex(Hx_w_FT,NUM_freq_inc);
//   Ey_w_FT=MALLOC1D_Complex(Ey_w_FT,NUM_freq_inc);
//   Hy_w_FT=MALLOC1D_Complex(Hy_w_FT,NUM_freq_inc);
//   Ez_w_FT=MALLOC1D_Complex(Ez_w_FT,NUM_freq_inc);
//   Hz_w_FT=MALLOC1D_Complex(Hz_w_FT,NUM_freq_inc);
//
// omega2 = MALLOC1D_double(omega2, NUM_freq_inc);
// freq2 = MALLOC1D_double(freq2,NUM_freq_inc);
// time2 = MALLOC1D_double(time2,Tend_inc);
//
//   //printf("%e\n%e",freq[0], freq[NUM_freq_inc-1]);
//
//   Hx_t_FT = ZERO_VECTORS1D_Complex(Hx_t_FT,Tend);
//   Hy_t_FT = ZERO_VECTORS1D_Complex(Hy_t_FT,Tend);
//   Hz_t_FT = ZERO_VECTORS1D_Complex(Hz_t_FT,Tend);
//
//   Ex_t_FT = ZERO_VECTORS1D_Complex(Ex_t_FT,Tend);
//   Ey_t_FT = ZERO_VECTORS1D_Complex(Ey_t_FT,Tend);
//   Ez_t_FT = ZERO_VECTORS1D_Complex(Ez_t_FT,Tend);
//
//   Hx_w_FT = ZERO_VECTORS1D_Complex(Hx_w_FT,NUM_freq_inc);
//   Hy_w_FT = ZERO_VECTORS1D_Complex(Hy_w_FT,NUM_freq_inc);
//   Hz_w_FT = ZERO_VECTORS1D_Complex(Hz_w_FT,NUM_freq_inc);
//
//   Ex_w_FT = ZERO_VECTORS1D_Complex(Ex_w_FT,NUM_freq_inc);
//   Ey_w_FT = ZERO_VECTORS1D_Complex(Ey_w_FT,NUM_freq_inc);
//   Ez_w_FT = ZERO_VECTORS1D_Complex(Ez_w_FT,NUM_freq_inc);
//
//   // FILE *SOURCE_IN;
//   // char filename[200];
//   // sprintf(filename,"Source.%d.txt",trials);
// 	// SOURCE_IN = fopen(filename,"r");
//   double Source_re1,Source_comp1,Source_re2,Source_comp2;
//
//   if(TEz){
//
//     printf("Define Function\n");
//      for(n=0;n<Tend_inc;n++){
//         //  fscanf(SOURCE_IN,"%lf %lf %lf %lf\n",&Source_re1, &Source_comp1, &Source_re2, &Source_comp2);
//         //  Ey_t_FT[n] = 1.0*(Source_re1 + I*Source_comp1);
//         //  //printf("%e\t%e\t%e\t%e\n",creal(Source_comp1),cimag(Source_comp1), creal(I*Source_comp1),cimag(I*Source_comp1));
//         //  Hx_t_FT[n] = -1.0*(Source_re2 + I*Source_comp2);
//         Ey_t_FT[n] = PULSE(n);
//       }
//
//       printf("DFT\n");
//
//       for(n=0;n<Tend_inc;n++){
//         for(w=0;w<NUM_freq_inc;w++){
//           Ey_w_FT[w] += Ey_t_FT[n] * cexp(-I*2*pi*n*w/Tend_inc);
//         }
//       }
//
//       printf("Transform\n");
//
//       for(w=0;w<NUM_freq_inc;w++){
//         Ey_w_FT[w] /= Tend_inc;
//         if(w>1){
//           k_z = csqrt((dw_2*w/c0)*(dw_2*w/c0) - (k_rho*k_rho));
//           k_0 = dw_2*w/c0;
//           Hx_w_FT[w] = (-1.0)*Ey_w_FT[w]*(k_z/k_0)*cexp(I*k_z*0.5*dz);
//           fprintf(Input_w,"%e \t %e \t %e \n", dw_2*w, cabs(Ey_w_FT[w]), cabs(Hx_w_FT[w]));
//         }
//       }
//       //printf("here\n");
//       printf("IDFT\n");
//
//       for(n=0;n<Tend_inc;n++){
//           for(w=1;w<NUM_freq_inc;w++){
//           Hx_t_FT[n] += Hx_w_FT[w] * cexp(I*2*pi*n*w/Tend_inc);
//         }
//
//         //fprintf(Input,"%e \t %e \t %e \t %e\n", time2[n] , creal(Ey_t_FT[n]), creal(Hx_t_FT[n]), creal(PULSE(n)));
//       }
//       printf("Output File\n");
//     for(n=0;n<Tend_inc;n++){
//       fprintf(Input,"%e \t %e \t %e\n", time2[n] , creal(Ey_t_FT[n]), creal(Hx_t_FT[n]));
//     }
//
//       i=0;
//       j=0;
//       //printf("here\n");
//       printf("Fourier Transform of Input Pulse\n");
//     for(n=0;n<Tend_inc;n++){
//       for(w=0;w<NUM_freq;w++){
//       //  for(i=0;i<NCELLX;i++){
//         //  for(j=0;j<NCELLY;j++){
//             E_Incident[w][i][j] += Ey_t_FT[n]*cexp(-I*2*pi*n*dt*freq[w]);//*cexp(-I*(i)*dx*k_x)*cexp(-I*(j+0.5)*dy*k_y);
//             H_Incident[w][i][j] += Hx_t_FT[n]*cexp(-I*2*pi*n*dt*freq[w]);//*cexp(-I*(i)*dx*k_x)*cexp(-I*(j+0.5)*dy*k_y);
//           //}
//         //}
//       }
//       printf("n=%d\n",n);
//
//     }
//
//   }
//
//   if(TMz){
//     printf("Define Function\n");
//     for(n=0;n<Tend_inc;n++){
//        Ex_t_FT[n] = PULSE(n);
//      }
//
//      printf("DFT\n");
//      for(n=0;n<Tend_inc;n++){
//        for(w=0;w<NUM_freq_inc;w++){
//          Ex_w_FT[w] += Ex_t_FT[n] * cexp(-I*2*pi*n*w/Tend_inc);
//        }
//      }
//
//
//      printf("Transform\n");
//      for(w=0;w<NUM_freq_inc;w++){
//        Ex_w_FT[w] /= Tend_inc;
//        if(w>1){
//          k_z = csqrt((dw_2*w/c0)*(dw_2*w/c0) - (k_rho*k_rho));
//          k_0 = dw_2*w/c0;
//          Hy_w_FT[w] = (-1.0)*Ex_w_FT[w]*(k_0/k_z)*cexp(I*k_z*0.5*dz);
//          fprintf(Input_w,"%e \t %e \t %e \n", dw_2*w, cabs(Ex_w_FT[w]), cabs(Hy_w_FT[w]));
//        }
//      }
//      //printf("here\n");
//      printf("IDFT\n");
//      for(n=0;n<Tend_inc;n++){
//          //Ey_t_FT[n] = 0.0;
//          for(w=1;w<NUM_freq_inc;w++){
//          //Hy_t_FT[n] += Hy_w_FT[w] * cexp(I*2*pi*n*w/Tend_inc);
//          Hy_t_FT[n] += Hy_w_FT[w] * cexp(I*2*pi*n*w/Tend_inc);
//
//        }
//
//        fprintf(Input,"%e \t %e \t %e \t %e\n", time2[n] , creal(Ex_t_FT[n]), creal(Hy_t_FT[n]), creal(PULSE(n)));
//      }
//      //printf("here\n");
//      printf("Fourier Transform of Input Pulse \n");
//      i=0;
//      j=0;
//    for(n=0;n<Tend_inc;n++){
//      for(w=0;w<NUM_freq;w++){
//        //for(i=0;i<NCELLX;i++){
//          //for(j=0;j<NCELLY;j++){
//            H_Incident[w][i][j] += Hy_t_FT[n]*cexp(-I*2*pi*n*dt*freq[w]);//*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//            E_Incident[w][i][j] += Ex_t_FT[n]*cexp(-I*2*pi*n*dt*freq[w]);//*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//          //}
//       // }
//      }
//    }
//   }
//   printf("HERE1\n");
//
//   fclose(Input);
//   fclose(Input_w);
//   //fclose(SOURCE_IN);
//   printf("HERE2\n");
//
//   free(omega2);
//   free(freq2);
//   free(time2);
//   printf("HERE3\n");
//
// }

void ALLOCATE_MEM(void){

     //allocate memory for 1-D incident fields

    e_inc=MALLOC1D_Complex(e_inc,inc_Length);
    h_inc=MALLOC1D_Complex(h_inc,inc_Length);
    ex_inc=MALLOC1D_Complex(ex_inc,inc_Length);
    ey_inc=MALLOC1D_Complex(ey_inc,inc_Length);
    ez_inc=MALLOC1D_Complex(ez_inc,inc_Length);
    hx_inc=MALLOC1D_Complex(hx_inc,inc_Length);
    hy_inc=MALLOC1D_Complex(hy_inc,inc_Length);
    hz_inc=MALLOC1D_Complex(hz_inc,inc_Length);

    //allocate memory for 3-D vectors
    ex=MALLOC3D_Complex(ex,NCELLX,NCELLY,NCELLZ);
    ey=MALLOC3D_Complex(ey,NCELLX,NCELLY,NCELLZ);
    ez=MALLOC3D_Complex(ez,NCELLX,NCELLY,NCELLZ);
//    Dx=MALLOC3D_Complex(Dx,NCELLX,NCELLY,NCELLZ);
//    Dy=MALLOC3D_Complex(Dy,NCELLX,NCELLY,NCELLZ);
//    Dz=MALLOC3D_Complex(Dz,NCELLX,NCELLY,NCELLZ);
    hx=MALLOC3D_Complex(hx,NCELLX,NCELLY,NCELLZ);
    hy=MALLOC3D_Complex(hy,NCELLX,NCELLY,NCELLZ);
    hz=MALLOC3D_Complex(hz,NCELLX,NCELLY,NCELLZ);


        size_t extent = NCELLX*NCELLY*NCELLZ * sizeof(real);
        size_t extent1 = NCELLX*NCELLY*NcpmlZ * sizeof(real);
        size_t extent2 = NCELLX*NcpmlY*NCELLZ * sizeof(real);
        size_t extent3 = NcpmlX*NCELLY*NCELLZ * sizeof(real);


        cudaMalloc(&hxdev,extent);
        cudaMalloc(&hydev,extent);
        cudaMalloc(&hzdev,extent);
        cudaMalloc(&exdev,extent);
        cudaMalloc(&ezdev,extent);
        cudaMalloc(&eydev,extent);

        cudaMalloc(&Chxhdev,extent);
        cudaMalloc(&Chxedev,extent);
        cudaMalloc(&Chyhdev,extent);
        cudaMalloc(&Chyedev,extent);
        cudaMalloc(&Chzhdev,extent);
        cudaMalloc(&Chzedev,extent);

        cudaMalloc(&Cexhdev,extent);
        cudaMalloc(&Cexedev,extent);
        cudaMalloc(&Ceyhdev,extent);
        cudaMalloc(&Ceyedev,extent);
        cudaMalloc(&Cezhdev,extent);
        cudaMalloc(&Cezedev,extent);

        cudaMalloc(&psi_Hx_z_Ndev,extent1);
        cudaMalloc(&psi_Hx_z_Fdev,extent1);
        cudaMalloc(&psi_Hx_y_Ndev,extent2);
        cudaMalloc(&psi_Hx_y_Fdev,extent2);

        cudaMalloc(&psi_Hy_z_Ndev,extent1);
        cudaMalloc(&psi_Hy_z_Fdev,extent1);
        cudaMalloc(&psi_Hy_x_Ndev,extent3);
        cudaMalloc(&psi_Hy_x_Fdev,extent3);

        cudaMalloc(&psi_Hz_y_Ndev,extent2);
        cudaMalloc(&psi_Hz_y_Fdev,extent2);
        cudaMalloc(&psi_Hz_x_Ndev,extent3);
        cudaMalloc(&psi_Hz_x_Fdev,extent3);

    //    size_t eextent = NCELLX*NCELLY*NCELLZ * sizeof(real);
        size_t eextent1 = NCELLX*NCELLY*(NcpmlZ+1) * sizeof(real);
        size_t eextent2 = NCELLX*(NcpmlY+1)*NCELLZ * sizeof(real);
        size_t eextent3 = (NcpmlX+1)*NCELLY*NCELLZ * sizeof(real);


        cudaMalloc(&psi_Ex_z_Ndev,eextent1);
        cudaMalloc(&psi_Ex_z_Fdev,eextent1);
        cudaMalloc(&psi_Ex_y_Ndev,eextent2);
        cudaMalloc(&psi_Ex_y_Fdev,eextent2);

        cudaMalloc(&psi_Ey_z_Ndev,eextent1);
        cudaMalloc(&psi_Ey_z_Fdev,eextent1);
        cudaMalloc(&psi_Ey_x_Ndev,eextent3);
        cudaMalloc(&psi_Ey_x_Fdev,eextent3);

        cudaMalloc(&psi_Ez_y_Ndev,eextent2);
        cudaMalloc(&psi_Ez_y_Fdev,eextent2);
        cudaMalloc(&psi_Ez_x_Ndev,eextent3);
        cudaMalloc(&psi_Ez_x_Fdev,eextent3);


        cudaMalloc(&khdydev,NCELLY*sizeof(real));
        cudaMalloc(&khdzdev,NCELLZ*sizeof(real));
        cudaMalloc(&khdxdev,NCELLX*sizeof(real));

        cudaMalloc(&kedydev,NCELLY*sizeof(real));
        cudaMalloc(&kedzdev,NCELLZ*sizeof(real));
        cudaMalloc(&kedxdev,NCELLX*sizeof(real));

        cudaMalloc(&bh_z_Ndev,NcpmlZ*sizeof(real));
        cudaMalloc(&bh_z_Fdev,NcpmlZ*sizeof(real));
        cudaMalloc(&ch_z_Ndev,NcpmlZ*sizeof(real));
        cudaMalloc(&ch_z_Fdev,NcpmlZ*sizeof(real));
        cudaMalloc(&bh_y_Ndev,NcpmlY*sizeof(real));
        cudaMalloc(&bh_y_Fdev,NcpmlY*sizeof(real));
        cudaMalloc(&ch_y_Ndev,NcpmlY*sizeof(real));
        cudaMalloc(&ch_y_Fdev,NcpmlY*sizeof(real));
        cudaMalloc(&bh_x_Ndev,NcpmlX*sizeof(real));
        cudaMalloc(&bh_x_Fdev,NcpmlX*sizeof(real));
        cudaMalloc(&ch_x_Ndev,NcpmlX*sizeof(real));
        cudaMalloc(&ch_x_Fdev,NcpmlX*sizeof(real));

        cudaMalloc(&be_z_Ndev,(NcpmlZ+1)*sizeof(real));
        cudaMalloc(&be_z_Fdev,(NcpmlZ+1)*sizeof(real));
        cudaMalloc(&ce_z_Ndev,(NcpmlZ+1)*sizeof(real));
        cudaMalloc(&ce_z_Fdev,(NcpmlZ+1)*sizeof(real));
        cudaMalloc(&be_y_Ndev,(NcpmlY+1)*sizeof(real));
        cudaMalloc(&be_y_Fdev,(NcpmlY+1)*sizeof(real));
        cudaMalloc(&ce_y_Ndev,(NcpmlY+1)*sizeof(real));
        cudaMalloc(&ce_y_Fdev,(NcpmlY+1)*sizeof(real));
        cudaMalloc(&be_x_Ndev,(NcpmlX+1)*sizeof(real));
        cudaMalloc(&be_x_Fdev,(NcpmlX+1)*sizeof(real));
        cudaMalloc(&ce_x_Ndev,(NcpmlX+1)*sizeof(real));
        cudaMalloc(&ce_x_Fdev,(NcpmlX+1)*sizeof(real));

        cudaMalloc(&ex_ndev,extent);
        cudaMalloc(&ez_ndev,extent);
        cudaMalloc(&ey_ndev,extent);

        cudaMalloc(&ex_n_1dev,extent);
        cudaMalloc(&ez_n_1dev,extent);
        cudaMalloc(&ey_n_1dev,extent);
        ex_n=MALLOC3D_Complex(ex_n,NCELLX,NCELLY,NCELLZ);
        ey_n=MALLOC3D_Complex(ey_n,NCELLX,NCELLY,NCELLZ);
        ez_n=MALLOC3D_Complex(ez_n,NCELLX,NCELLY,NCELLZ);

        ex_n_1=MALLOC3D_Complex(ex_n_1,NCELLX,NCELLY,NCELLZ);
        ey_n_1=MALLOC3D_Complex(ey_n_1,NCELLX,NCELLY,NCELLZ);
        ez_n_1=MALLOC3D_Complex(ez_n_1,NCELLX,NCELLY,NCELLZ);


    if(Hydrodynamics == 1){
      NDx=MALLOC3D_Complex(NDx,NCELLX,NCELLY,NCELLZ);
      NDy=MALLOC3D_Complex(NDy,NCELLX,NCELLY,NCELLZ);
      NDz=MALLOC3D_Complex(NDz,NCELLX,NCELLY,NCELLZ);
      cudaMalloc(&NDxdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      cudaMalloc(&NDydev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      cudaMalloc(&NDzdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      NDx_prev=MALLOC3D_Complex(NDx_prev,NCELLX,NCELLY,NCELLZ);
      NDy_prev=MALLOC3D_Complex(NDy_prev,NCELLX,NCELLY,NCELLZ);
      NDz_prev=MALLOC3D_Complex(NDz_prev,NCELLX,NCELLY,NCELLZ);
      cudaMalloc(&NDx_prevdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      cudaMalloc(&NDy_prevdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      cudaMalloc(&NDz_prevdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      hxPrev=MALLOC3D_Complex(hxPrev,NCELLX,NCELLY,NCELLZ);
      hyPrev=MALLOC3D_Complex(hyPrev,NCELLX,NCELLY,NCELLZ);
      hzPrev=MALLOC3D_Complex(hzPrev,NCELLX,NCELLY,NCELLZ);
      cudaMalloc(&hxPrevdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      cudaMalloc(&hyPrevdev,NCELLX*NCELLY*NCELLZ*sizeof(real));
      cudaMalloc(&hzPrevdev,NCELLX*NCELLY*NCELLZ*sizeof(real));


    }

    ExTransformNearZScaRe = MALLOC3D_Real2(ExTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    EyTransformNearZScaRe = MALLOC3D_Real2(EyTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    EzTransformNearZScaRe = MALLOC3D_Real2(EzTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    HxTransformNearZScaRe = MALLOC3D_Real2(HxTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    HyTransformNearZScaRe = MALLOC3D_Real2(HyTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    HzTransformNearZScaRe = MALLOC3D_Real2(HzTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);

    cudaMalloc(&ExTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
    cudaMalloc(&EyTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
    cudaMalloc(&HxTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
    cudaMalloc(&HyTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));

    ExTransformNearYScaRe = MALLOC3D_Real2(ExTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    EyTransformNearYScaRe = MALLOC3D_Real2(EyTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    EzTransformNearYScaRe = MALLOC3D_Real2(EzTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    HxTransformNearYScaRe = MALLOC3D_Real2(HxTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    HyTransformNearYScaRe = MALLOC3D_Real2(HyTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    HzTransformNearYScaRe = MALLOC3D_Real2(HzTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);

    cudaMalloc(&ExTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&EzTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HxTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HzTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));

    ExTransformNearXScaRe = MALLOC3D_Real2(ExTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    EyTransformNearXScaRe = MALLOC3D_Real2(EyTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    EzTransformNearXScaRe = MALLOC3D_Real2(EzTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    HxTransformNearXScaRe = MALLOC3D_Real2(HxTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    HyTransformNearXScaRe = MALLOC3D_Real2(HyTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    HzTransformNearXScaRe = MALLOC3D_Real2(HzTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

    cudaMalloc(&EyTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&EzTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HyTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HzTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));


        ExTransformNearZScaIm = MALLOC3D_Real2(ExTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        EyTransformNearZScaIm = MALLOC3D_Real2(EyTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        EzTransformNearZScaIm = MALLOC3D_Real2(EzTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        HxTransformNearZScaIm = MALLOC3D_Real2(HxTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        HyTransformNearZScaIm = MALLOC3D_Real2(HyTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        HzTransformNearZScaIm = MALLOC3D_Real2(HzTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);

        cudaMalloc(&ExTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
        cudaMalloc(&EyTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
        cudaMalloc(&HxTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
        cudaMalloc(&HyTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));

        ExTransformNearYScaIm = MALLOC3D_Real2(ExTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EyTransformNearYScaIm = MALLOC3D_Real2(EyTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EzTransformNearYScaIm = MALLOC3D_Real2(EzTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HxTransformNearYScaIm = MALLOC3D_Real2(HxTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HyTransformNearYScaIm = MALLOC3D_Real2(HyTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HzTransformNearYScaIm = MALLOC3D_Real2(HzTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

        cudaMalloc(&ExTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&EzTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HxTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HzTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));

        ExTransformNearXScaIm = MALLOC3D_Real2(ExTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EyTransformNearXScaIm = MALLOC3D_Real2(EyTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EzTransformNearXScaIm = MALLOC3D_Real2(EzTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HxTransformNearXScaIm = MALLOC3D_Real2(HxTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HyTransformNearXScaIm = MALLOC3D_Real2(HyTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HzTransformNearXScaIm = MALLOC3D_Real2(HzTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

        cudaMalloc(&EyTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&EzTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HyTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HzTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));


    ExTransformNearZAbsRe = MALLOC3D_Real2(ExTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EyTransformNearZAbsRe = MALLOC3D_Real2(EyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EzTransformNearZAbsRe = MALLOC3D_Real2(EzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HxTransformNearZAbsRe = MALLOC3D_Real2(HxTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HyTransformNearZAbsRe = MALLOC3D_Real2(HyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HzTransformNearZAbsRe = MALLOC3D_Real2(HzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);

    cudaMalloc(&ExTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&EyTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HxTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HyTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));

    ExTransformNearYAbsRe = MALLOC3D_Real2(ExTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformNearYAbsRe = MALLOC3D_Real2(EyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformNearYAbsRe = MALLOC3D_Real2(EzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformNearYAbsRe = MALLOC3D_Real2(HxTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformNearYAbsRe = MALLOC3D_Real2(HyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformNearYAbsRe = MALLOC3D_Real2(HzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&ExTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HxTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    ExTransformNearXAbsRe = MALLOC3D_Real2(ExTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformNearXAbsRe = MALLOC3D_Real2(EyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformNearXAbsRe = MALLOC3D_Real2(EzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformNearXAbsRe = MALLOC3D_Real2(HxTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformNearXAbsRe = MALLOC3D_Real2(HyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformNearXAbsRe = MALLOC3D_Real2(HzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&EyTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HyTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    ExTransformNearZAbsIm = MALLOC3D_Real2(ExTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EyTransformNearZAbsIm = MALLOC3D_Real2(EyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EzTransformNearZAbsIm = MALLOC3D_Real2(EzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HxTransformNearZAbsIm = MALLOC3D_Real2(HxTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HyTransformNearZAbsIm = MALLOC3D_Real2(HyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HzTransformNearZAbsIm = MALLOC3D_Real2(HzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);

    cudaMalloc(&ExTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&EyTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HxTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HyTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));

    ExTransformNearYAbsIm = MALLOC3D_Real2(ExTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformNearYAbsIm = MALLOC3D_Real2(EyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformNearYAbsIm = MALLOC3D_Real2(EzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformNearYAbsIm = MALLOC3D_Real2(HxTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformNearYAbsIm = MALLOC3D_Real2(HyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformNearYAbsIm = MALLOC3D_Real2(HzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&ExTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HxTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    ExTransformNearXAbsIm = MALLOC3D_Real2(ExTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformNearXAbsIm = MALLOC3D_Real2(EyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformNearXAbsIm = MALLOC3D_Real2(EzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformNearXAbsIm = MALLOC3D_Real2(HxTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformNearXAbsIm = MALLOC3D_Real2(HyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformNearXAbsIm = MALLOC3D_Real2(HzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&EyTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HyTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    // ExTransformFarZScaRe = MALLOC3D_Real2(ExTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // EyTransformFarZScaRe = MALLOC3D_Real2(EyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // EzTransformFarZScaRe = MALLOC3D_Real2(EzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // HxTransformFarZScaRe = MALLOC3D_Real2(HxTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // HyTransformFarZScaRe = MALLOC3D_Real2(HyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // HzTransformFarZScaRe = MALLOC3D_Real2(HzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    //
    // ExTransformFarYScaRe = MALLOC3D_Real2(ExTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EyTransformFarYScaRe = MALLOC3D_Real2(EyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EzTransformFarYScaRe = MALLOC3D_Real2(EzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HxTransformFarYScaRe = MALLOC3D_Real2(HxTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HyTransformFarYScaRe = MALLOC3D_Real2(HyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HzTransformFarYScaRe = MALLOC3D_Real2(HzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    //
    // ExTransformFarXScaRe = MALLOC3D_Real2(ExTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EyTransformFarXScaRe = MALLOC3D_Real2(EyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EzTransformFarXScaRe = MALLOC3D_Real2(EzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HxTransformFarXScaRe = MALLOC3D_Real2(HxTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HyTransformFarXScaRe = MALLOC3D_Real2(HyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HzTransformFarXScaRe = MALLOC3D_Real2(HzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    //
    //
    //
    // ExTransformFarZScaIm = MALLOC3D_Real2(ExTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // EyTransformFarZScaIm = MALLOC3D_Real2(EyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // EzTransformFarZScaIm = MALLOC3D_Real2(EzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // HxTransformFarZScaIm = MALLOC3D_Real2(HxTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // HyTransformFarZScaIm = MALLOC3D_Real2(HyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    // HzTransformFarZScaIm = MALLOC3D_Real2(HzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
    //
    // ExTransformFarYScaIm = MALLOC3D_Real2(ExTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EyTransformFarYScaIm = MALLOC3D_Real2(EyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EzTransformFarYScaIm = MALLOC3D_Real2(EzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HxTransformFarYScaIm = MALLOC3D_Real2(HxTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HyTransformFarYScaIm = MALLOC3D_Real2(HyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HzTransformFarYScaIm = MALLOC3D_Real2(HzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    //
    // ExTransformFarXScaIm = MALLOC3D_Real2(ExTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EyTransformFarXScaIm = MALLOC3D_Real2(EyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // EzTransformFarXScaIm = MALLOC3D_Real2(EzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HxTransformFarXScaIm = MALLOC3D_Real2(HxTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HyTransformFarXScaIm = MALLOC3D_Real2(HyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    // HzTransformFarXScaIm = MALLOC3D_Real2(HzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    //
    //
    // ExTransformFarZAbsRe = MALLOC3D_Real2(ExTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    // EyTransformFarZAbsRe = MALLOC3D_Real2(EyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    // EzTransformFarZAbsRe = MALLOC3D_Real2(EzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    // HxTransformFarZAbsRe = MALLOC3D_Real2(HxTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    // HyTransformFarZAbsRe = MALLOC3D_Real2(HyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    // HzTransformFarZAbsRe = MALLOC3D_Real2(HzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //
    // ExTransformFarYAbsRe = MALLOC3D_Real2(ExTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // EyTransformFarYAbsRe = MALLOC3D_Real2(EyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // EzTransformFarYAbsRe = MALLOC3D_Real2(EzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // HxTransformFarYAbsRe = MALLOC3D_Real2(HxTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // HyTransformFarYAbsRe = MALLOC3D_Real2(HyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // HzTransformFarYAbsRe = MALLOC3D_Real2(HzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //
    // ExTransformFarXAbsRe = MALLOC3D_Real2(ExTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // EyTransformFarXAbsRe = MALLOC3D_Real2(EyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // EzTransformFarXAbsRe = MALLOC3D_Real2(EzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // HxTransformFarXAbsRe = MALLOC3D_Real2(HxTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // HyTransformFarXAbsRe = MALLOC3D_Real2(HyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    // HzTransformFarXAbsRe = MALLOC3D_Real2(HzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //
    //
    //
    //     ExTransformFarZAbsIm = MALLOC3D_Real2(ExTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //     EyTransformFarZAbsIm = MALLOC3D_Real2(EyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //     EzTransformFarZAbsIm = MALLOC3D_Real2(EzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //     HxTransformFarZAbsIm = MALLOC3D_Real2(HxTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //     HyTransformFarZAbsIm = MALLOC3D_Real2(HyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //     HzTransformFarZAbsIm = MALLOC3D_Real2(HzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    //
    //     ExTransformFarYAbsIm = MALLOC3D_Real2(ExTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     EyTransformFarYAbsIm = MALLOC3D_Real2(EyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     EzTransformFarYAbsIm = MALLOC3D_Real2(EzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     HxTransformFarYAbsIm = MALLOC3D_Real2(HxTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     HyTransformFarYAbsIm = MALLOC3D_Real2(HyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     HzTransformFarYAbsIm = MALLOC3D_Real2(HzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //
    //     ExTransformFarXAbsIm = MALLOC3D_Real2(ExTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     EyTransformFarXAbsIm = MALLOC3D_Real2(EyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     EzTransformFarXAbsIm = MALLOC3D_Real2(EzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     HxTransformFarXAbsIm = MALLOC3D_Real2(HxTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     HyTransformFarXAbsIm = MALLOC3D_Real2(HyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //     HzTransformFarXAbsIm = MALLOC3D_Real2(HzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    //




    ExTransformFarZScaRe = MALLOC3D_Real2(ExTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    EyTransformFarZScaRe = MALLOC3D_Real2(EyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    EzTransformFarZScaRe = MALLOC3D_Real2(EzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    HxTransformFarZScaRe = MALLOC3D_Real2(HxTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    HyTransformFarZScaRe = MALLOC3D_Real2(HyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
    HzTransformFarZScaRe = MALLOC3D_Real2(HzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);

    cudaMalloc(&ExTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
    cudaMalloc(&EyTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
    cudaMalloc(&HxTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
    cudaMalloc(&HyTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));

    ExTransformFarYScaRe = MALLOC3D_Real2(ExTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    EyTransformFarYScaRe = MALLOC3D_Real2(EyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    EzTransformFarYScaRe = MALLOC3D_Real2(EzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    HxTransformFarYScaRe = MALLOC3D_Real2(HxTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    HyTransformFarYScaRe = MALLOC3D_Real2(HyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
    HzTransformFarYScaRe = MALLOC3D_Real2(HzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);

    cudaMalloc(&ExTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&EzTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HxTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HzTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));

    ExTransformFarXScaRe = MALLOC3D_Real2(ExTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    EyTransformFarXScaRe = MALLOC3D_Real2(EyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    EzTransformFarXScaRe = MALLOC3D_Real2(EzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    HxTransformFarXScaRe = MALLOC3D_Real2(HxTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    HyTransformFarXScaRe = MALLOC3D_Real2(HyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
    HzTransformFarXScaRe = MALLOC3D_Real2(HzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

    cudaMalloc(&EyTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&EzTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HyTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
    cudaMalloc(&HzTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));


        ExTransformFarZScaIm = MALLOC3D_Real2(ExTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        EyTransformFarZScaIm = MALLOC3D_Real2(EyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        EzTransformFarZScaIm = MALLOC3D_Real2(EzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        HxTransformFarZScaIm = MALLOC3D_Real2(HxTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        HyTransformFarZScaIm = MALLOC3D_Real2(HyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        HzTransformFarZScaIm = MALLOC3D_Real2(HzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);

        cudaMalloc(&ExTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
        cudaMalloc(&EyTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
        cudaMalloc(&HxTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));
        cudaMalloc(&HyTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1));

        ExTransformFarYScaIm = MALLOC3D_Real2(ExTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EyTransformFarYScaIm = MALLOC3D_Real2(EyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EzTransformFarYScaIm = MALLOC3D_Real2(EzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HxTransformFarYScaIm = MALLOC3D_Real2(HxTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HyTransformFarYScaIm = MALLOC3D_Real2(HyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HzTransformFarYScaIm = MALLOC3D_Real2(HzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

        cudaMalloc(&ExTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&EzTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HxTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HzTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1));

        ExTransformFarXScaIm = MALLOC3D_Real2(ExTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EyTransformFarXScaIm = MALLOC3D_Real2(EyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EzTransformFarXScaIm = MALLOC3D_Real2(EzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HxTransformFarXScaIm = MALLOC3D_Real2(HxTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HyTransformFarXScaIm = MALLOC3D_Real2(HyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HzTransformFarXScaIm = MALLOC3D_Real2(HzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

        cudaMalloc(&EyTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&EzTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HyTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));
        cudaMalloc(&HzTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1));


    ExTransformFarZAbsRe = MALLOC3D_Real2(ExTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EyTransformFarZAbsRe = MALLOC3D_Real2(EyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EzTransformFarZAbsRe = MALLOC3D_Real2(EzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HxTransformFarZAbsRe = MALLOC3D_Real2(HxTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HyTransformFarZAbsRe = MALLOC3D_Real2(HyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HzTransformFarZAbsRe = MALLOC3D_Real2(HzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);

    cudaMalloc(&ExTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&EyTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HxTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HyTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));

    ExTransformFarYAbsRe = MALLOC3D_Real2(ExTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformFarYAbsRe = MALLOC3D_Real2(EyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformFarYAbsRe = MALLOC3D_Real2(EzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformFarYAbsRe = MALLOC3D_Real2(HxTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformFarYAbsRe = MALLOC3D_Real2(HyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformFarYAbsRe = MALLOC3D_Real2(HzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&ExTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HxTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    ExTransformFarXAbsRe = MALLOC3D_Real2(ExTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformFarXAbsRe = MALLOC3D_Real2(EyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformFarXAbsRe = MALLOC3D_Real2(EzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformFarXAbsRe = MALLOC3D_Real2(HxTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformFarXAbsRe = MALLOC3D_Real2(HyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformFarXAbsRe = MALLOC3D_Real2(HzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&EyTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HyTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    ExTransformFarZAbsIm = MALLOC3D_Real2(ExTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EyTransformFarZAbsIm = MALLOC3D_Real2(EyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    EzTransformFarZAbsIm = MALLOC3D_Real2(EzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HxTransformFarZAbsIm = MALLOC3D_Real2(HxTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HyTransformFarZAbsIm = MALLOC3D_Real2(HyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
    HzTransformFarZAbsIm = MALLOC3D_Real2(HzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);

    cudaMalloc(&ExTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&EyTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HxTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));
    cudaMalloc(&HyTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1));

    ExTransformFarYAbsIm = MALLOC3D_Real2(ExTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformFarYAbsIm = MALLOC3D_Real2(EyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformFarYAbsIm = MALLOC3D_Real2(EzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformFarYAbsIm = MALLOC3D_Real2(HxTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformFarYAbsIm = MALLOC3D_Real2(HyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformFarYAbsIm = MALLOC3D_Real2(HzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&ExTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HxTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));

    ExTransformFarXAbsIm = MALLOC3D_Real2(ExTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EyTransformFarXAbsIm = MALLOC3D_Real2(EyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    EzTransformFarXAbsIm = MALLOC3D_Real2(EzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HxTransformFarXAbsIm = MALLOC3D_Real2(HxTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HyTransformFarXAbsIm = MALLOC3D_Real2(HyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
    HzTransformFarXAbsIm = MALLOC3D_Real2(HzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

    cudaMalloc(&EyTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&EzTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HyTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));
    cudaMalloc(&HzTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1));




        ExTransformNearZScaRe = ZERO_VECTORS3D_Real2(ExTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
        EyTransformNearZScaRe = ZERO_VECTORS3D_Real2(EyTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
        EzTransformNearZScaRe = ZERO_VECTORS3D_Real2(EzTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
        HxTransformNearZScaRe = ZERO_VECTORS3D_Real2(HxTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
        HyTransformNearZScaRe = ZERO_VECTORS3D_Real2(HyTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
        HzTransformNearZScaRe = ZERO_VECTORS3D_Real2(HzTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);


        cudaMemcpy(ExTransformNearZScaRedev,ExTransformNearZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(EyTransformNearZScaRedev,EyTransformNearZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HxTransformNearZScaRedev,HxTransformNearZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HyTransformNearZScaRedev,HyTransformNearZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);

        ExTransformNearYScaRe = ZERO_VECTORS3D_Real2(ExTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
        EyTransformNearYScaRe = ZERO_VECTORS3D_Real2(EyTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
        EzTransformNearYScaRe = ZERO_VECTORS3D_Real2(EzTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
        HxTransformNearYScaRe = ZERO_VECTORS3D_Real2(HxTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
        HyTransformNearYScaRe = ZERO_VECTORS3D_Real2(HyTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
        HzTransformNearYScaRe = ZERO_VECTORS3D_Real2(HzTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);

        cudaMemcpy(ExTransformNearYScaRedev,ExTransformNearYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(EzTransformNearYScaRedev,EzTransformNearYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HxTransformNearYScaRedev,HxTransformNearYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HzTransformNearYScaRedev,HzTransformNearYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

        ExTransformNearXScaRe = ZERO_VECTORS3D_Real2(ExTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EyTransformNearXScaRe = ZERO_VECTORS3D_Real2(EyTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        EzTransformNearXScaRe = ZERO_VECTORS3D_Real2(EzTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HxTransformNearXScaRe = ZERO_VECTORS3D_Real2(HxTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HyTransformNearXScaRe = ZERO_VECTORS3D_Real2(HyTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        HzTransformNearXScaRe = ZERO_VECTORS3D_Real2(HzTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

        cudaMemcpy(EyTransformNearXScaRedev,EyTransformNearXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(EzTransformNearXScaRedev,EzTransformNearXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HyTransformNearXScaRedev,HyTransformNearXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HzTransformNearXScaRedev,HzTransformNearXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

            ExTransformNearZScaIm = ZERO_VECTORS3D_Real2(ExTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
            EyTransformNearZScaIm = ZERO_VECTORS3D_Real2(EyTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
            EzTransformNearZScaIm = ZERO_VECTORS3D_Real2(EzTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
            HxTransformNearZScaIm = ZERO_VECTORS3D_Real2(HxTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
            HyTransformNearZScaIm = ZERO_VECTORS3D_Real2(HyTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
            HzTransformNearZScaIm = ZERO_VECTORS3D_Real2(HzTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);

            cudaMemcpy(ExTransformNearZScaImdev,ExTransformNearZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(EyTransformNearZScaImdev,EyTransformNearZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HxTransformNearZScaImdev,HxTransformNearZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HyTransformNearZScaImdev,HyTransformNearZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);

            ExTransformNearYScaIm = ZERO_VECTORS3D_Real2(ExTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            EyTransformNearYScaIm = ZERO_VECTORS3D_Real2(EyTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            EzTransformNearYScaIm = ZERO_VECTORS3D_Real2(EzTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            HxTransformNearYScaIm = ZERO_VECTORS3D_Real2(HxTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            HyTransformNearYScaIm = ZERO_VECTORS3D_Real2(HyTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            HzTransformNearYScaIm = ZERO_VECTORS3D_Real2(HzTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

            cudaMemcpy(ExTransformNearYScaImdev,ExTransformNearYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(EzTransformNearYScaImdev,EzTransformNearYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HxTransformNearYScaImdev,HxTransformNearYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HzTransformNearYScaImdev,HzTransformNearYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

            ExTransformNearXScaIm = ZERO_VECTORS3D_Real2(ExTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            EyTransformNearXScaIm = ZERO_VECTORS3D_Real2(EyTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            EzTransformNearXScaIm = ZERO_VECTORS3D_Real2(EzTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            HxTransformNearXScaIm = ZERO_VECTORS3D_Real2(HxTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            HyTransformNearXScaIm = ZERO_VECTORS3D_Real2(HyTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
            HzTransformNearXScaIm = ZERO_VECTORS3D_Real2(HzTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

            cudaMemcpy(EyTransformNearXScaImdev,EyTransformNearXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(EzTransformNearXScaImdev,EzTransformNearXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HyTransformNearXScaImdev,HyTransformNearXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HzTransformNearXScaImdev,HzTransformNearXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

        // ExTransformNearZAbsRe = ZERO_VECTORS3D_Real2(ExTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // EyTransformNearZAbsRe = ZERO_VECTORS3D_Real2(EyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // EzTransformNearZAbsRe = ZERO_VECTORS3D_Real2(EzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HxTransformNearZAbsRe = ZERO_VECTORS3D_Real2(HxTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HyTransformNearZAbsRe = ZERO_VECTORS3D_Real2(HyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HzTransformNearZAbsRe = ZERO_VECTORS3D_Real2(HzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //
        // ExTransformNearYAbsRe = ZERO_VECTORS3D_Real2(ExTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EyTransformNearYAbsRe = ZERO_VECTORS3D_Real2(EyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EzTransformNearYAbsRe = ZERO_VECTORS3D_Real2(EzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HxTransformNearYAbsRe = ZERO_VECTORS3D_Real2(HxTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HyTransformNearYAbsRe = ZERO_VECTORS3D_Real2(HyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HzTransformNearYAbsRe = ZERO_VECTORS3D_Real2(HzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //
        // ExTransformNearXAbsRe = ZERO_VECTORS3D_Real2(ExTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EyTransformNearXAbsRe = ZERO_VECTORS3D_Real2(EyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EzTransformNearXAbsRe = ZERO_VECTORS3D_Real2(EzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HxTransformNearXAbsRe = ZERO_VECTORS3D_Real2(HxTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HyTransformNearXAbsRe = ZERO_VECTORS3D_Real2(HyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HzTransformNearXAbsRe = ZERO_VECTORS3D_Real2(HzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //
        //
        // ExTransformNearZAbsIm = ZERO_VECTORS3D_Real2(ExTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // EyTransformNearZAbsIm = ZERO_VECTORS3D_Real2(EyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // EzTransformNearZAbsIm = ZERO_VECTORS3D_Real2(EzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HxTransformNearZAbsIm = ZERO_VECTORS3D_Real2(HxTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HyTransformNearZAbsIm = ZERO_VECTORS3D_Real2(HyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HzTransformNearZAbsIm = ZERO_VECTORS3D_Real2(HzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //
        // ExTransformNearYAbsIm = ZERO_VECTORS3D_Real2(ExTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EyTransformNearYAbsIm = ZERO_VECTORS3D_Real2(EyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EzTransformNearYAbsIm = ZERO_VECTORS3D_Real2(EzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HxTransformNearYAbsIm = ZERO_VECTORS3D_Real2(HxTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HyTransformNearYAbsIm = ZERO_VECTORS3D_Real2(HyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HzTransformNearYAbsIm = ZERO_VECTORS3D_Real2(HzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //
        // ExTransformNearXAbsIm = ZERO_VECTORS3D_Real2(ExTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EyTransformNearXAbsIm = ZERO_VECTORS3D_Real2(EyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EzTransformNearXAbsIm = ZERO_VECTORS3D_Real2(EzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HxTransformNearXAbsIm = ZERO_VECTORS3D_Real2(HxTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HyTransformNearXAbsIm = ZERO_VECTORS3D_Real2(HyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HzTransformNearXAbsIm = ZERO_VECTORS3D_Real2(HzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);













        ExTransformNearZAbsRe = ZERO_VECTORS3D_Real2(ExTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
        EyTransformNearZAbsRe = ZERO_VECTORS3D_Real2(EyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
        EzTransformNearZAbsRe = ZERO_VECTORS3D_Real2(EzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
        HxTransformNearZAbsRe = ZERO_VECTORS3D_Real2(HxTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
        HyTransformNearZAbsRe = ZERO_VECTORS3D_Real2(HyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
        HzTransformNearZAbsRe = ZERO_VECTORS3D_Real2(HzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);


        cudaMemcpy(ExTransformNearZAbsRedev,ExTransformNearZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(EyTransformNearZAbsRedev,EyTransformNearZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HxTransformNearZAbsRedev,HxTransformNearZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HyTransformNearZAbsRedev,HyTransformNearZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);

        ExTransformNearYAbsRe = ZERO_VECTORS3D_Real2(ExTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
        EyTransformNearYAbsRe = ZERO_VECTORS3D_Real2(EyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
        EzTransformNearYAbsRe = ZERO_VECTORS3D_Real2(EzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
        HxTransformNearYAbsRe = ZERO_VECTORS3D_Real2(HxTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
        HyTransformNearYAbsRe = ZERO_VECTORS3D_Real2(HyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
        HzTransformNearYAbsRe = ZERO_VECTORS3D_Real2(HzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);

        cudaMemcpy(ExTransformNearYAbsRedev,ExTransformNearYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(EzTransformNearYAbsRedev,EzTransformNearYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HxTransformNearYAbsRedev,HxTransformNearYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HzTransformNearYAbsRedev,HzTransformNearYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);

        ExTransformNearXAbsRe = ZERO_VECTORS3D_Real2(ExTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        EyTransformNearXAbsRe = ZERO_VECTORS3D_Real2(EyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        EzTransformNearXAbsRe = ZERO_VECTORS3D_Real2(EzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        HxTransformNearXAbsRe = ZERO_VECTORS3D_Real2(HxTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        HyTransformNearXAbsRe = ZERO_VECTORS3D_Real2(HyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        HzTransformNearXAbsRe = ZERO_VECTORS3D_Real2(HzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

        cudaMemcpy(EyTransformNearXAbsRedev,EyTransformNearXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(EzTransformNearXAbsRedev,EzTransformNearXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HyTransformNearXAbsRedev,HyTransformNearXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
        cudaMemcpy(HzTransformNearXAbsRedev,HzTransformNearXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);

            ExTransformNearZAbsIm = ZERO_VECTORS3D_Real2(ExTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
            EyTransformNearZAbsIm = ZERO_VECTORS3D_Real2(EyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
            EzTransformNearZAbsIm = ZERO_VECTORS3D_Real2(EzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
            HxTransformNearZAbsIm = ZERO_VECTORS3D_Real2(HxTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
            HyTransformNearZAbsIm = ZERO_VECTORS3D_Real2(HyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
            HzTransformNearZAbsIm = ZERO_VECTORS3D_Real2(HzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);

            cudaMemcpy(ExTransformNearZAbsImdev,ExTransformNearZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(EyTransformNearZAbsImdev,EyTransformNearZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HxTransformNearZAbsImdev,HxTransformNearZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HyTransformNearZAbsImdev,HyTransformNearZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);

            ExTransformNearYAbsIm = ZERO_VECTORS3D_Real2(ExTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            EyTransformNearYAbsIm = ZERO_VECTORS3D_Real2(EyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            EzTransformNearYAbsIm = ZERO_VECTORS3D_Real2(EzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            HxTransformNearYAbsIm = ZERO_VECTORS3D_Real2(HxTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            HyTransformNearYAbsIm = ZERO_VECTORS3D_Real2(HyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            HzTransformNearYAbsIm = ZERO_VECTORS3D_Real2(HzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

            cudaMemcpy(ExTransformNearYAbsImdev,ExTransformNearYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(EzTransformNearYAbsImdev,EzTransformNearYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HxTransformNearYAbsImdev,HxTransformNearYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HzTransformNearYAbsImdev,HzTransformNearYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);

            ExTransformNearXAbsIm = ZERO_VECTORS3D_Real2(ExTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            EyTransformNearXAbsIm = ZERO_VECTORS3D_Real2(EyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            EzTransformNearXAbsIm = ZERO_VECTORS3D_Real2(EzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            HxTransformNearXAbsIm = ZERO_VECTORS3D_Real2(HxTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            HyTransformNearXAbsIm = ZERO_VECTORS3D_Real2(HyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
            HzTransformNearXAbsIm = ZERO_VECTORS3D_Real2(HzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

            cudaMemcpy(EyTransformNearXAbsImdev,EyTransformNearXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(EzTransformNearXAbsImdev,EzTransformNearXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HyTransformNearXAbsImdev,HyTransformNearXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
            cudaMemcpy(HzTransformNearXAbsImdev,HzTransformNearXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);










        // ExTransformFarZScaRe = ZERO_VECTORS3D_Real2(ExTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // EyTransformFarZScaRe = ZERO_VECTORS3D_Real2(EyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // EzTransformFarZScaRe = ZERO_VECTORS3D_Real2(EzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // HxTransformFarZScaRe = ZERO_VECTORS3D_Real2(HxTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // HyTransformFarZScaRe = ZERO_VECTORS3D_Real2(HyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // HzTransformFarZScaRe = ZERO_VECTORS3D_Real2(HzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        //
        // ExTransformFarYScaRe = ZERO_VECTORS3D_Real2(ExTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EyTransformFarYScaRe = ZERO_VECTORS3D_Real2(EyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EzTransformFarYScaRe = ZERO_VECTORS3D_Real2(EzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HxTransformFarYScaRe = ZERO_VECTORS3D_Real2(HxTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HyTransformFarYScaRe = ZERO_VECTORS3D_Real2(HyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HzTransformFarYScaRe = ZERO_VECTORS3D_Real2(HzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        //
        // ExTransformFarXScaRe = ZERO_VECTORS3D_Real2(ExTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EyTransformFarXScaRe = ZERO_VECTORS3D_Real2(EyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EzTransformFarXScaRe = ZERO_VECTORS3D_Real2(EzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HxTransformFarXScaRe = ZERO_VECTORS3D_Real2(HxTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HyTransformFarXScaRe = ZERO_VECTORS3D_Real2(HyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HzTransformFarXScaRe = ZERO_VECTORS3D_Real2(HzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        //
        //
        //
        // ExTransformFarZScaIm = ZERO_VECTORS3D_Real2(ExTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // EyTransformFarZScaIm = ZERO_VECTORS3D_Real2(EyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // EzTransformFarZScaIm = ZERO_VECTORS3D_Real2(EzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // HxTransformFarZScaIm = ZERO_VECTORS3D_Real2(HxTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // HyTransformFarZScaIm = ZERO_VECTORS3D_Real2(HyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        // HzTransformFarZScaIm = ZERO_VECTORS3D_Real2(HzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
        //
        // ExTransformFarYScaIm = ZERO_VECTORS3D_Real2(ExTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EyTransformFarYScaIm = ZERO_VECTORS3D_Real2(EyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EzTransformFarYScaIm = ZERO_VECTORS3D_Real2(EzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HxTransformFarYScaIm = ZERO_VECTORS3D_Real2(HxTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HyTransformFarYScaIm = ZERO_VECTORS3D_Real2(HyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HzTransformFarYScaIm = ZERO_VECTORS3D_Real2(HzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        //
        // ExTransformFarXScaIm = ZERO_VECTORS3D_Real2(ExTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EyTransformFarXScaIm = ZERO_VECTORS3D_Real2(EyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // EzTransformFarXScaIm = ZERO_VECTORS3D_Real2(EzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HxTransformFarXScaIm = ZERO_VECTORS3D_Real2(HxTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HyTransformFarXScaIm = ZERO_VECTORS3D_Real2(HyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        // HzTransformFarXScaIm = ZERO_VECTORS3D_Real2(HzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
        //
        //
        // ExTransformFarZAbsRe = ZERO_VECTORS3D_Real2(ExTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // EyTransformFarZAbsRe = ZERO_VECTORS3D_Real2(EyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // EzTransformFarZAbsRe = ZERO_VECTORS3D_Real2(EzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HxTransformFarZAbsRe = ZERO_VECTORS3D_Real2(HxTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HyTransformFarZAbsRe = ZERO_VECTORS3D_Real2(HyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        // HzTransformFarZAbsRe = ZERO_VECTORS3D_Real2(HzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //
        // ExTransformFarYAbsRe = ZERO_VECTORS3D_Real2(ExTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EyTransformFarYAbsRe = ZERO_VECTORS3D_Real2(EyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EzTransformFarYAbsRe = ZERO_VECTORS3D_Real2(EzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HxTransformFarYAbsRe = ZERO_VECTORS3D_Real2(HxTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HyTransformFarYAbsRe = ZERO_VECTORS3D_Real2(HyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HzTransformFarYAbsRe = ZERO_VECTORS3D_Real2(HzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //
        // ExTransformFarXAbsRe = ZERO_VECTORS3D_Real2(ExTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EyTransformFarXAbsRe = ZERO_VECTORS3D_Real2(EyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // EzTransformFarXAbsRe = ZERO_VECTORS3D_Real2(EzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HxTransformFarXAbsRe = ZERO_VECTORS3D_Real2(HxTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HyTransformFarXAbsRe = ZERO_VECTORS3D_Real2(HyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        // HzTransformFarXAbsRe = ZERO_VECTORS3D_Real2(HzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //
        //
        //
        //     ExTransformFarZAbsIm = ZERO_VECTORS3D_Real2(ExTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //     EyTransformFarZAbsIm = ZERO_VECTORS3D_Real2(EyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //     EzTransformFarZAbsIm = ZERO_VECTORS3D_Real2(EzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //     HxTransformFarZAbsIm = ZERO_VECTORS3D_Real2(HxTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //     HyTransformFarZAbsIm = ZERO_VECTORS3D_Real2(HyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //     HzTransformFarZAbsIm = ZERO_VECTORS3D_Real2(HzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
        //
        //     ExTransformFarYAbsIm = ZERO_VECTORS3D_Real2(ExTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     EyTransformFarYAbsIm = ZERO_VECTORS3D_Real2(EyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     EzTransformFarYAbsIm = ZERO_VECTORS3D_Real2(EzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     HxTransformFarYAbsIm = ZERO_VECTORS3D_Real2(HxTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     HyTransformFarYAbsIm = ZERO_VECTORS3D_Real2(HyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     HzTransformFarYAbsIm = ZERO_VECTORS3D_Real2(HzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //
        //     ExTransformFarXAbsIm = ZERO_VECTORS3D_Real2(ExTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     EyTransformFarXAbsIm = ZERO_VECTORS3D_Real2(EyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     EzTransformFarXAbsIm = ZERO_VECTORS3D_Real2(EzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     HxTransformFarXAbsIm = ZERO_VECTORS3D_Real2(HxTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     HyTransformFarXAbsIm = ZERO_VECTORS3D_Real2(HyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
        //     HzTransformFarXAbsIm = ZERO_VECTORS3D_Real2(HzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);



                ExTransformFarZScaRe = ZERO_VECTORS3D_Real2(ExTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
                EyTransformFarZScaRe = ZERO_VECTORS3D_Real2(EyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
                EzTransformFarZScaRe = ZERO_VECTORS3D_Real2(EzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
                HxTransformFarZScaRe = ZERO_VECTORS3D_Real2(HxTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
                HyTransformFarZScaRe = ZERO_VECTORS3D_Real2(HyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);
                HzTransformFarZScaRe = ZERO_VECTORS3D_Real2(HzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca +1,YENDSca - YSTARTSca+1);


                cudaMemcpy(ExTransformFarZScaRedev,ExTransformFarZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(EyTransformFarZScaRedev,EyTransformFarZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(HxTransformFarZScaRedev,HxTransformFarZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(HyTransformFarZScaRedev,HyTransformFarZScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);

                ExTransformFarYScaRe = ZERO_VECTORS3D_Real2(ExTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
                EyTransformFarYScaRe = ZERO_VECTORS3D_Real2(EyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
                EzTransformFarYScaRe = ZERO_VECTORS3D_Real2(EzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
                HxTransformFarYScaRe = ZERO_VECTORS3D_Real2(HxTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
                HyTransformFarYScaRe = ZERO_VECTORS3D_Real2(HyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);
                HzTransformFarYScaRe = ZERO_VECTORS3D_Real2(HzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca +1,ZENDSca - ZSTARTSca+1);

                cudaMemcpy(ExTransformFarYScaRedev,ExTransformFarYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(EzTransformFarYScaRedev,EzTransformFarYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(HxTransformFarYScaRedev,HxTransformFarYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(HzTransformFarYScaRedev,HzTransformFarYScaRe,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

                ExTransformFarXScaRe = ZERO_VECTORS3D_Real2(ExTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                EyTransformFarXScaRe = ZERO_VECTORS3D_Real2(EyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                EzTransformFarXScaRe = ZERO_VECTORS3D_Real2(EzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                HxTransformFarXScaRe = ZERO_VECTORS3D_Real2(HxTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                HyTransformFarXScaRe = ZERO_VECTORS3D_Real2(HyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                HzTransformFarXScaRe = ZERO_VECTORS3D_Real2(HzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

                cudaMemcpy(EyTransformFarXScaRedev,EyTransformFarXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(EzTransformFarXScaRedev,EzTransformFarXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(HyTransformFarXScaRedev,HyTransformFarXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                cudaMemcpy(HzTransformFarXScaRedev,HzTransformFarXScaRe,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

                    ExTransformFarZScaIm = ZERO_VECTORS3D_Real2(ExTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
                    EyTransformFarZScaIm = ZERO_VECTORS3D_Real2(EyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
                    EzTransformFarZScaIm = ZERO_VECTORS3D_Real2(EzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
                    HxTransformFarZScaIm = ZERO_VECTORS3D_Real2(HxTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
                    HyTransformFarZScaIm = ZERO_VECTORS3D_Real2(HyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);
                    HzTransformFarZScaIm = ZERO_VECTORS3D_Real2(HzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,YENDSca - YSTARTSca + 1);

                    cudaMemcpy(ExTransformFarZScaImdev,ExTransformFarZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(EyTransformFarZScaImdev,EyTransformFarZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(HxTransformFarZScaImdev,HxTransformFarZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(HyTransformFarZScaImdev,HyTransformFarZScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyHostToDevice);

                    ExTransformFarYScaIm = ZERO_VECTORS3D_Real2(ExTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    EyTransformFarYScaIm = ZERO_VECTORS3D_Real2(EyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    EzTransformFarYScaIm = ZERO_VECTORS3D_Real2(EzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    HxTransformFarYScaIm = ZERO_VECTORS3D_Real2(HxTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    HyTransformFarYScaIm = ZERO_VECTORS3D_Real2(HyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    HzTransformFarYScaIm = ZERO_VECTORS3D_Real2(HzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

                    cudaMemcpy(ExTransformFarYScaImdev,ExTransformFarYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(EzTransformFarYScaImdev,EzTransformFarYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(HxTransformFarYScaImdev,HxTransformFarYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(HzTransformFarYScaImdev,HzTransformFarYScaIm,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);

                    ExTransformFarXScaIm = ZERO_VECTORS3D_Real2(ExTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    EyTransformFarXScaIm = ZERO_VECTORS3D_Real2(EyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    EzTransformFarXScaIm = ZERO_VECTORS3D_Real2(EzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    HxTransformFarXScaIm = ZERO_VECTORS3D_Real2(HxTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    HyTransformFarXScaIm = ZERO_VECTORS3D_Real2(HyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);
                    HzTransformFarXScaIm = ZERO_VECTORS3D_Real2(HzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca + 1 ,ZENDSca - ZSTARTSca + 1);

                    cudaMemcpy(EyTransformFarXScaImdev,EyTransformFarXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(EzTransformFarXScaImdev,EzTransformFarXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(HyTransformFarXScaImdev,HyTransformFarXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);
                    cudaMemcpy(HzTransformFarXScaImdev,HzTransformFarXScaIm,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyHostToDevice);



                            ExTransformFarZAbsRe = ZERO_VECTORS3D_Real2(ExTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
                            EyTransformFarZAbsRe = ZERO_VECTORS3D_Real2(EyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
                            EzTransformFarZAbsRe = ZERO_VECTORS3D_Real2(EzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
                            HxTransformFarZAbsRe = ZERO_VECTORS3D_Real2(HxTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
                            HyTransformFarZAbsRe = ZERO_VECTORS3D_Real2(HyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);
                            HzTransformFarZAbsRe = ZERO_VECTORS3D_Real2(HzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,YENDAbs - YSTARTAbs+1);


                            cudaMemcpy(ExTransformFarZAbsRedev,ExTransformFarZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(EyTransformFarZAbsRedev,EyTransformFarZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(HxTransformFarZAbsRedev,HxTransformFarZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(HyTransformFarZAbsRedev,HyTransformFarZAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);

                            ExTransformFarYAbsRe = ZERO_VECTORS3D_Real2(ExTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
                            EyTransformFarYAbsRe = ZERO_VECTORS3D_Real2(EyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
                            EzTransformFarYAbsRe = ZERO_VECTORS3D_Real2(EzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
                            HxTransformFarYAbsRe = ZERO_VECTORS3D_Real2(HxTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
                            HyTransformFarYAbsRe = ZERO_VECTORS3D_Real2(HyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);
                            HzTransformFarYAbsRe = ZERO_VECTORS3D_Real2(HzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs +1,ZENDAbs - ZSTARTAbs+1);

                            cudaMemcpy(ExTransformFarYAbsRedev,ExTransformFarYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(EzTransformFarYAbsRedev,EzTransformFarYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(HxTransformFarYAbsRedev,HxTransformFarYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(HzTransformFarYAbsRedev,HzTransformFarYAbsRe,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);

                            ExTransformFarXAbsRe = ZERO_VECTORS3D_Real2(ExTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                            EyTransformFarXAbsRe = ZERO_VECTORS3D_Real2(EyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                            EzTransformFarXAbsRe = ZERO_VECTORS3D_Real2(EzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                            HxTransformFarXAbsRe = ZERO_VECTORS3D_Real2(HxTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                            HyTransformFarXAbsRe = ZERO_VECTORS3D_Real2(HyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                            HzTransformFarXAbsRe = ZERO_VECTORS3D_Real2(HzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

                            cudaMemcpy(EyTransformFarXAbsRedev,EyTransformFarXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(EzTransformFarXAbsRedev,EzTransformFarXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(HyTransformFarXAbsRedev,HyTransformFarXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                            cudaMemcpy(HzTransformFarXAbsRedev,HzTransformFarXAbsRe,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);

                                ExTransformFarZAbsIm = ZERO_VECTORS3D_Real2(ExTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
                                EyTransformFarZAbsIm = ZERO_VECTORS3D_Real2(EyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
                                EzTransformFarZAbsIm = ZERO_VECTORS3D_Real2(EzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
                                HxTransformFarZAbsIm = ZERO_VECTORS3D_Real2(HxTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
                                HyTransformFarZAbsIm = ZERO_VECTORS3D_Real2(HyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);
                                HzTransformFarZAbsIm = ZERO_VECTORS3D_Real2(HzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,YENDAbs - YSTARTAbs + 1);

                                cudaMemcpy(ExTransformFarZAbsImdev,ExTransformFarZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(EyTransformFarZAbsImdev,EyTransformFarZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(HxTransformFarZAbsImdev,HxTransformFarZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(HyTransformFarZAbsImdev,HyTransformFarZAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyHostToDevice);

                                ExTransformFarYAbsIm = ZERO_VECTORS3D_Real2(ExTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                EyTransformFarYAbsIm = ZERO_VECTORS3D_Real2(EyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                EzTransformFarYAbsIm = ZERO_VECTORS3D_Real2(EzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                HxTransformFarYAbsIm = ZERO_VECTORS3D_Real2(HxTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                HyTransformFarYAbsIm = ZERO_VECTORS3D_Real2(HyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                HzTransformFarYAbsIm = ZERO_VECTORS3D_Real2(HzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

                                cudaMemcpy(ExTransformFarYAbsImdev,ExTransformFarYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(EzTransformFarYAbsImdev,EzTransformFarYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(HxTransformFarYAbsImdev,HxTransformFarYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(HzTransformFarYAbsImdev,HzTransformFarYAbsIm,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);

                                ExTransformFarXAbsIm = ZERO_VECTORS3D_Real2(ExTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                EyTransformFarXAbsIm = ZERO_VECTORS3D_Real2(EyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                EzTransformFarXAbsIm = ZERO_VECTORS3D_Real2(EzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                HxTransformFarXAbsIm = ZERO_VECTORS3D_Real2(HxTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                HyTransformFarXAbsIm = ZERO_VECTORS3D_Real2(HyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);
                                HzTransformFarXAbsIm = ZERO_VECTORS3D_Real2(HzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs + 1 ,ZENDAbs - ZSTARTAbs + 1);

                                cudaMemcpy(EyTransformFarXAbsImdev,EyTransformFarXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(EzTransformFarXAbsImdev,EzTransformFarXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(HyTransformFarXAbsImdev,HyTransformFarXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);
                                cudaMemcpy(HzTransformFarXAbsImdev,HzTransformFarXAbsIm,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyHostToDevice);



    Ex_Reflected = MALLOC3D_Complex2(Ex_Reflected,NUM_freq,NCELLX,NCELLY);
    Hx_Reflected = MALLOC3D_Complex2(Hx_Reflected,NUM_freq,NCELLX,NCELLY);
    Ey_Reflected = MALLOC3D_Complex2(Ey_Reflected,NUM_freq,NCELLX,NCELLY);
    Hy_Reflected = MALLOC3D_Complex2(Hy_Reflected,NUM_freq,NCELLX,NCELLY);

    Ex_Transmitted  = MALLOC3D_Complex2(Ex_Transmitted ,NUM_freq,NCELLX,NCELLY);
    Hx_Transmitted  = MALLOC3D_Complex2(Hx_Transmitted ,NUM_freq,NCELLX,NCELLY);
    Ey_Transmitted  = MALLOC3D_Complex2(Ey_Transmitted ,NUM_freq,NCELLX,NCELLY);
    Hy_Transmitted  = MALLOC3D_Complex2(Hy_Transmitted ,NUM_freq,NCELLX,NCELLY);

    E_Incident = MALLOC3D_Complex2(E_Incident,NUM_freq,NCELLX,NCELLY);
    H_Incident = MALLOC3D_Complex2(H_Incident,NUM_freq,NCELLX,NCELLY);



    Cexe=MALLOC3D(Cexe,NCELLX,NCELLY,NCELLZ);
    Cexh=MALLOC3D(Cexh,NCELLX,NCELLY,NCELLZ);
    Ceye=MALLOC3D(Ceye,NCELLX,NCELLY,NCELLZ);
    Ceyh=MALLOC3D(Ceyh,NCELLX,NCELLY,NCELLZ);
    Ceze=MALLOC3D(Ceze,NCELLX,NCELLY,NCELLZ);
    Cezh=MALLOC3D(Cezh,NCELLX,NCELLY,NCELLZ);

    Chxe=MALLOC3D(Chxe,NCELLX,NCELLY,NCELLZ);
    Chxh=MALLOC3D(Chxh,NCELLX,NCELLY,NCELLZ);
    Chye=MALLOC3D(Chye,NCELLX,NCELLY,NCELLZ);
    Chyh=MALLOC3D(Chyh,NCELLX,NCELLY,NCELLZ);
    Chze=MALLOC3D(Chze,NCELLX,NCELLY,NCELLZ);
    Chzh=MALLOC3D(Chzh,NCELLX,NCELLY,NCELLZ);

    eps=MALLOC3D(eps,NCELLX,NCELLY,NCELLZ);
    mu=MALLOC3D(mu,NCELLX,NCELLY,NCELLZ);
    sigma_e=MALLOC3D(sigma_e,NCELLX,NCELLY,NCELLZ);
    sigma_m=MALLOC3D(sigma_m,NCELLX,NCELLY,NCELLZ);
    //
    psi_Ex_y_N=MALLOC3D_Complex(psi_Ex_y_N,NCELLX,NcpmlY+1,NCELLZ);
    psi_Ez_y_N=MALLOC3D_Complex(psi_Ez_y_N,NCELLX,NcpmlY+1,NCELLZ);
    psi_Ex_y_F=MALLOC3D_Complex(psi_Ex_y_F,NCELLX,NcpmlY+1,NCELLZ);
    psi_Ez_y_F=MALLOC3D_Complex(psi_Ez_y_F,NCELLX,NcpmlY+1,NCELLZ);

    psi_Ex_z_N=MALLOC3D_Complex(psi_Ex_z_N,NCELLX,NCELLY,NcpmlZ+1);
    psi_Ey_z_N=MALLOC3D_Complex(psi_Ey_z_N,NCELLX,NCELLY,NcpmlZ+1);
    psi_Ex_z_F=MALLOC3D_Complex(psi_Ex_z_F,NCELLX,NCELLY,NcpmlZ+1);
    psi_Ey_z_F=MALLOC3D_Complex(psi_Ey_z_F,NCELLX,NCELLY,NcpmlZ+1);

    psi_Ey_x_N=MALLOC3D_Complex(psi_Ey_x_N,NcpmlX+1,NCELLY,NCELLZ);
    psi_Ez_x_N=MALLOC3D_Complex(psi_Ez_x_N,NcpmlX+1,NCELLY,NCELLZ);
    psi_Ey_x_F=MALLOC3D_Complex(psi_Ey_x_F,NcpmlX+1,NCELLY,NCELLZ);
    psi_Ez_x_F=MALLOC3D_Complex(psi_Ez_x_F,NcpmlX+1,NCELLY,NCELLZ);

    psi_Hx_y_F=MALLOC3D_Complex(psi_Hx_y_F,NCELLX,NcpmlY,NCELLZ);
    psi_Hz_y_F=MALLOC3D_Complex(psi_Hz_y_F,NCELLX,NcpmlY,NCELLZ);
    psi_Hx_y_N=MALLOC3D_Complex(psi_Hx_y_N,NCELLX,NcpmlY,NCELLZ);
    psi_Hz_y_N=MALLOC3D_Complex(psi_Hz_y_N,NCELLX,NcpmlY,NCELLZ);

    psi_Hx_z_F=MALLOC3D_Complex(psi_Hx_z_F,NCELLX,NCELLY,NcpmlZ);
    psi_Hy_z_F=MALLOC3D_Complex(psi_Hy_z_F,NCELLX,NCELLY,NcpmlZ);
    psi_Hx_z_N=MALLOC3D_Complex(psi_Hx_z_N,NCELLX,NCELLY,NcpmlZ);
    psi_Hy_z_N=MALLOC3D_Complex(psi_Hy_z_N,NCELLX,NCELLY,NcpmlZ);

    psi_Hz_x_F=MALLOC3D_Complex(psi_Hz_x_F,NcpmlX,NCELLY,NCELLZ);
    psi_Hy_x_F=MALLOC3D_Complex(psi_Hy_x_F,NcpmlX,NCELLY,NCELLZ);
    psi_Hy_x_N=MALLOC3D_Complex(psi_Hy_x_N,NcpmlX,NCELLY,NCELLZ);
    psi_Hz_x_N=MALLOC3D_Complex(psi_Hz_x_N,NcpmlX,NCELLY,NCELLZ);

    kedx=MALLOC1D(kedx,NCELLX);
    kedy=MALLOC1D(kedy,NCELLY);
    kedz=MALLOC1D(kedz,NCELLZ);
    khdx=MALLOC1D(khdx,NCELLX);
    khdy=MALLOC1D(khdy,NCELLY);
    khdz=MALLOC1D(khdz,NCELLZ);

    be_x_N=MALLOC1D(be_x_N,NcpmlX+1);
    be_y_N=MALLOC1D(be_y_N,NcpmlY+1);
    be_z_N=MALLOC1D(be_z_N,NcpmlZ+1);
    bh_x_N=MALLOC1D(bh_x_N,NcpmlX);
    bh_y_N=MALLOC1D(bh_y_N,NcpmlY);
    bh_z_N=MALLOC1D(bh_z_N,NcpmlZ);

    ce_x_N=MALLOC1D(ce_x_N,NcpmlX+1);
    ce_y_N=MALLOC1D(ce_y_N,NcpmlY+1);
    ce_z_N=MALLOC1D(ce_z_N,NcpmlZ+1);
    ch_x_N=MALLOC1D(ch_x_N,NcpmlX);
    ch_y_N=MALLOC1D(ch_y_N,NcpmlY);
    ch_z_N=MALLOC1D(ch_z_N,NcpmlZ);

    be_x_F=MALLOC1D(be_x_F,NcpmlX+1);
    be_y_F=MALLOC1D(be_y_F,NcpmlY+1);
    be_z_F=MALLOC1D(be_z_F,NcpmlZ+1);
    bh_x_F=MALLOC1D(bh_x_F,NcpmlX);
    bh_y_F=MALLOC1D(bh_y_F,NcpmlY);
    bh_z_F=MALLOC1D(bh_z_F,NcpmlZ);

    ce_x_F=MALLOC1D(ce_x_F,NcpmlX+1);
    ce_y_F=MALLOC1D(ce_y_F,NcpmlY+1);
    ce_z_F=MALLOC1D(ce_z_F,NcpmlZ+1);
    ch_x_F=MALLOC1D(ch_x_F,NcpmlX);
    ch_y_F=MALLOC1D(ch_y_F,NcpmlY);
    ch_z_F=MALLOC1D(ch_z_F,NcpmlZ);

    sigma_e_x=MALLOC1D(sigma_e_x,NcpmlX);
    sigma_e_y=MALLOC1D(sigma_e_y,NcpmlY);
    sigma_e_z=MALLOC1D(sigma_e_z,NcpmlZ);
    sigma_h_x=MALLOC1D(sigma_h_x,NcpmlX);
    sigma_h_y=MALLOC1D(sigma_h_y,NcpmlY);
    sigma_h_z=MALLOC1D(sigma_h_z,NcpmlZ);

    t_inc = MALLOC1D_double(t_inc, NUM_freq);

    E_incident=MALLOC1D_Complex2(E_incident,NUM_freq);
    E_reflected=MALLOC1D_Complex2(E_reflected,NUM_freq);
    E_transmitted=MALLOC1D_Complex2(E_transmitted,NUM_freq);


return ;
}

void FREE_MEM(void){

    FREE1D_Complex2(E_incident);
    FREE1D_Complex2(E_reflected);

    FREE1D_Complex(e_inc);
    FREE1D_Complex(h_inc);

    FREE1D_Complex(ex_inc);
    FREE1D_Complex(ey_inc);
    FREE1D_Complex(ez_inc);
    FREE1D_Complex(hx_inc);
    FREE1D_Complex(hy_inc);
    FREE1D_Complex(hz_inc);

    FREE3D_Complex(ez,NCELLX,NCELLY);
    FREE3D_Complex(ey,NCELLX,NCELLY);
    FREE3D_Complex(ex,NCELLX,NCELLY);
    FREE3D_Complex(hx,NCELLX,NCELLY);
    FREE3D_Complex(hy,NCELLX,NCELLY);
    FREE3D_Complex(hz,NCELLX,NCELLY);
    FREE3D_Complex(hxPrev,NCELLX,NCELLY);
    FREE3D_Complex(hyPrev,NCELLX,NCELLY);
    FREE3D_Complex(hzPrev,NCELLX,NCELLY);
    // FREE3D_Complex(ExTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformNearZScaRe,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformNearYScaRe,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EyTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EzTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HxTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HyTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HzTransformNearXScaRe,NUM_freq,YENDSca - YSTARTSca);
    //
    // FREE3D_Complex(ExTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformFarZScaRe,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformFarYScaRe,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HxTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HyTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HzTransformFarXScaRe,NUM_freq,YENDSca - YSTARTSca);
    //
    //
    //
    //
    // FREE3D_Complex(ExTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformNearZScaIm,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformNearYScaIm,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EyTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EzTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HxTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HyTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HzTransformNearXScaIm,NUM_freq,YENDSca - YSTARTSca);
    //
    // FREE3D_Complex(ExTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformFarZScaIm,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(EzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HxTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HyTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca);
    // FREE3D_Complex(HzTransformFarYScaIm,NUM_freq,XENDSca - XSTARTSca);
    //
    // FREE3D_Complex(ExTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(EzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HxTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HyTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca);
    // FREE3D_Complex(HzTransformFarXScaIm,NUM_freq,YENDSca - YSTARTSca);
    //
    //
    //
    //
    //
    //
    // FREE3D_Complex(ExTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformNearZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformNearYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HxTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HyTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HzTransformNearXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    //
    // FREE3D_Complex(ExTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformFarZAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformFarYAbsRe,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HxTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HyTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HzTransformFarXAbsRe,NUM_freq,YENDAbs - YSTARTAbs);
    //
    // FREE3D_Complex(ExTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformNearZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformNearYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HxTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HyTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HzTransformNearXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    //
    // FREE3D_Complex(ExTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformFarZAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(EzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HxTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HyTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    // FREE3D_Complex(HzTransformFarYAbsIm,NUM_freq,XENDAbs - XSTARTAbs);
    //
    // FREE3D_Complex(ExTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(EzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HxTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HyTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);
    // FREE3D_Complex(HzTransformFarXAbsIm,NUM_freq,YENDAbs - YSTARTAbs);

    FREE4D_Complex(Pz_cp,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Py_cp,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Px_cp,NCELLX,NCELLY,NCELLZ);

    FREE4D_Complex(Pz_cp_n,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Py_cp_n,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Px_cp_n,NCELLX,NCELLY,NCELLZ);

    FREE4D_Complex(Pz_cp_n_1,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Py_cp_n_1,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Px_cp_n_1,NCELLX,NCELLY,NCELLZ);

    FREE4D_Complex(Px_d,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Px_d_n,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Px_d_n_1,NCELLX,NCELLY,NCELLZ);

    FREE4D_Complex(Py_d,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Py_d_n,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Py_d_n_1,NCELLX,NCELLY,NCELLZ);

    FREE4D_Complex(Pz_d,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Pz_d_n,NCELLX,NCELLY,NCELLZ);
    FREE4D_Complex(Pz_d_n_1,NCELLX,NCELLY,NCELLZ);


    FREE3D(Cexe,NCELLX,NCELLY);
    FREE3D(Ceye,NCELLX,NCELLY);
    FREE3D(Ceze,NCELLX,NCELLY);
    FREE3D(Chxe,NCELLX,NCELLY);
    FREE3D(Chye,NCELLX,NCELLY);
    FREE3D(Chze,NCELLX,NCELLY);
    FREE3D(Cexh,NCELLX,NCELLY);
    FREE3D(Ceyh,NCELLX,NCELLY);
    FREE3D(Cezh,NCELLX,NCELLY);
    FREE3D(Chxh,NCELLX,NCELLY);
    FREE3D(Chyh,NCELLX,NCELLY);
    FREE3D(Chzh,NCELLX,NCELLY);

    FREE3D(sigma_e,NCELLX,NCELLY);
    FREE3D(sigma_m,NCELLX,NCELLY);
    FREE3D(eps,NCELLX,NCELLY);
    FREE3D(mu,NCELLX,NCELLY);

    FREE3D_Complex(psi_Ex_y_N,NCELLX,NcpmlY+1);
    FREE3D_Complex(psi_Ex_z_N,NCELLX,NCELLY);
    FREE3D_Complex(psi_Ey_x_N,NcpmlX+1,NCELLY);
    FREE3D_Complex(psi_Ey_z_N,NCELLX,NCELLY);
    FREE3D_Complex(psi_Ez_y_N,NCELLX,NcpmlY+1);
    FREE3D_Complex(psi_Ez_x_N,NcpmlX+1,NCELLY);
    FREE3D_Complex(psi_Hx_z_N,NCELLX,NCELLY);
    FREE3D_Complex(psi_Hx_y_N,NCELLX,NcpmlY);
    FREE3D_Complex(psi_Hy_x_N,NcpmlX,NCELLY);
    FREE3D_Complex(psi_Hy_z_N,NCELLX,NCELLY);
    FREE3D_Complex(psi_Hz_x_N,NcpmlX,NCELLY);
    FREE3D_Complex(psi_Hz_y_N,NCELLX,NcpmlY);

    FREE3D_Complex(psi_Ex_y_F,NCELLX,NcpmlY+1);
    FREE3D_Complex(psi_Ex_z_F,NCELLX,NCELLY);
    FREE3D_Complex(psi_Ey_x_F,NcpmlX+1,NCELLY);
    FREE3D_Complex(psi_Ey_z_F,NCELLX,NCELLY);
    FREE3D_Complex(psi_Ez_y_F,NCELLX,NcpmlY+1);
    FREE3D_Complex(psi_Ez_x_F,NcpmlX+1,NCELLY);
    FREE3D_Complex(psi_Hx_z_F,NCELLX,NCELLY);
    FREE3D_Complex(psi_Hx_y_F,NCELLX,NcpmlY);
    FREE3D_Complex(psi_Hy_x_F,NcpmlX,NCELLY);
    FREE3D_Complex(psi_Hy_z_F,NCELLX,NCELLY);
    FREE3D_Complex(psi_Hz_x_F,NcpmlX,NCELLY);
    FREE3D_Complex(psi_Hz_y_F,NCELLX,NcpmlY);

    FREE3D_Complex2(E_Incident,NUM_freq,NCELLX);
    FREE3D_Complex2(Ex_Reflected,NUM_freq,NCELLX);
    FREE3D_Complex2(H_Incident,NUM_freq,NCELLX);
    FREE3D_Complex2(Hx_Reflected,NUM_freq,NCELLX);
    FREE3D_Complex2(Ey_Reflected,NUM_freq,NCELLX);
    FREE3D_Complex2(Hy_Reflected,NUM_freq,NCELLX);
    FREE3D_Complex2(Hx_Transmitted,NUM_freq,NCELLX);
    FREE3D_Complex2(Ex_Transmitted,NUM_freq,NCELLX);
    FREE3D_Complex2(Ey_Transmitted,NUM_freq,NCELLX);
    FREE3D_Complex2(Hy_Transmitted,NUM_freq,NCELLX);




    FREE1D(be_x_N);
    FREE1D(be_y_N);
    FREE1D(be_z_N);
    FREE1D(bh_x_N);
    FREE1D(bh_y_N);
    FREE1D(bh_z_N);

    FREE1D(ce_x_N);
    FREE1D(ce_y_N);
    FREE1D(ce_z_N);
    FREE1D(ch_x_N);
    FREE1D(ch_y_N);
    FREE1D(ch_z_N);

    FREE1D(be_x_F);
    FREE1D(be_y_F);
    FREE1D(be_z_F);
    FREE1D(bh_x_F);
    FREE1D(bh_y_F);
    FREE1D(bh_z_F);

    FREE1D(ce_x_F);
    FREE1D(ce_y_F);
    FREE1D(ce_z_F);
    FREE1D(ch_x_F);
    FREE1D(ch_y_F);
    FREE1D(ch_z_F);


    FREE1D(kedx);
    FREE1D(kedy);
    FREE1D(kedz);
    FREE1D(khdx);
    FREE1D(khdy);
    FREE1D(khdz);

    FREE1D(sigma_e_x);
    FREE1D(sigma_e_y);
    FREE1D(sigma_e_z);
    FREE1D(sigma_h_x);
    FREE1D(sigma_h_y);
    FREE1D(sigma_h_z);

  //  FREE1D_double(t_inc);

    return;
}
