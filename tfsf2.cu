#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "extern_var.h"

void SETUP_TFSF(void){

    //length of the incident field vectors
    //inc_Length=4*ceil(sqrt(NCELLX*NCELLX+NCELLY*NCELLY+NCELLZ*NCELLZ));
    //polarization and incident angles in degrees
    inc_phi=0;
    inc_theta=0;
    polar_psi=0;
    //convert to radians
    inc_phi*=pi/180.0;
    inc_theta*=pi/180.0;
    polar_psi*=pi/180.0 ;

    d_1D=dx;//*pow(pow(sin(inc_theta),4)*(pow(sin(inc_phi),4)+pow(cos(inc_phi),4))+pow(cos(inc_theta),4),0.5);
   // printf("%e\n",d_1D);
    /*printf("\n\n\n ex: %f\n",(cos(polar_psi)*sin(inc_phi)-sin(polar_psi)*cos(inc_theta)*cos(inc_phi)));
    printf("ey: %f\n",(-cos(polar_psi)*cos(inc_phi)-sin(polar_psi)*cos(inc_theta)*sin(inc_phi)));
    printf("ez: %f\n",(sin(polar_psi)*sin(inc_theta)));
    printf("hx: %f\n",(sin(polar_psi)*sin(inc_phi)+cos(polar_psi)*cos(inc_theta)*cos(inc_phi)));
    printf("hy: %f\n",(-sin(polar_psi)*cos(inc_phi)+cos(polar_psi)*cos(inc_theta)*sin(inc_phi)));
    printf("hz: %f\n\n\n",(-cos(polar_psi)*sin(inc_theta)));*/

    //where 1-D array meets 3-D array
    m0=floor(inc_Length/2);

    i_0=NtfsfX;
    j_0=NtfsfY;
    k_0=NtfsfZ;

    //incident pulse parameter;

    real BW=BandWidth;
    width_real = 6/(BW*2*pi);
    width=ceil((6/(BW*2*pi))/dt);
    delay=6*width;

    printf("delay = %f \t width = %f\n", delay, width);
  // f_0=k_rho*c0/(2*pi)+BW/1.3;
  // f_0;
  //f_0=k_rho*c0/(2*pi)+BW/2;

  printf("f_0 = %e\n",f_0);
	//f_0 = 6283200.000000*c0/(2*pi)+BW/2;
   // f_0=1e10;
    f_min=k_rho*c0/(2*pi);


    //first order Mur ABC buffers
    e1=e2=0;

}

// void TFSF_CORRECT(void){
//
//    UPDATE_h_inc();
//    // printf("HERE1");
//
//   if(PBC_CTW == 0){
//     if(Periodic_XY==1) {
//
//       CORRECT_Z();
//     }
//     else{
//       CORRECT_X();
//
//       CORRECT_Y();
//
//       CORRECT_Z();
//
//     //  printf("Here");
//
//     }
//
//   }
//
// else{
//   // SOURCE_IN_B();
//   // SOURCE_IN_E();
//   SOURCE_IN();
// }
//     UPDATE_e_inc();
//     //SOURCE_IN();
//
// }
//update 1-D incident fields
//
// void UPDATE_e_inc(void){
//     real pulse;
//     real factor = c0*dt;
//     int i,m;
//
//     //1st Order Mur ABC buffers
//     e1=e_inc[1];
//
//     e2=e_inc[inc_Length-2];
//
//     //#pragma omp parallel for private(i)// schedule(guided)
//     for(i=1;i<inc_Length-1;i++){
//         e_inc[i]=e_inc[i]-(dt/(z0*ep0*d_1D))*(h_inc[i]-h_inc[i-1]);
//     }
//   //  printf("%f\t%f\n",(dt/(z0*ep0*d_1D))*(h_inc[i]-h_inc[i-1]),((c0*dt-d_1D)/(c0*dt+d_1D)));
//     //1st order Mur ABC
//
//     e_inc[0]=e1+((factor-d_1D)/(factor+d_1D))*(e_inc[1]-e_inc[0]);
//
//     e_inc[inc_Length-1]=e2+((factor-d_1D)/(factor+d_1D))*(e_inc[inc_Length-2]-e_inc[inc_Length-1]);
//
//     //introduce source
//     #ifdef DOUBLECOMPLEX
//     pulse=exp(-pow((real)(t-delay)/(real)width,2)/2.0)*cexp(I*2*pi*f_0*(t)*dt);
//     #endif
//     #ifdef DOUBLEPRECISION
//     pulse=exp(-pow((real)(t-delay)/(real)width,2)/2.0)*sin(2*pi*f_0*(t)*dt);
//     #endif
//     for(m=0;m<NUM_freq;m++){
//       E_incident[m] += pulse*cexp(I*2*pi*t*dt*freq[m]);
//     }
//     e_inc[m0-50]+=pulse;
//
//   //  printf("%f\n",pulse);
// }
//
// void UPDATE_h_inc(void){
//     int i;
//
// //    //#pragma omp parallel for private(i)
//     for(i=0;i<inc_Length-1;i++){
//         h_inc[i]=h_inc[i]-(z0*dt/(mu0*d_1D))*(e_inc[i+1]-e_inc[i]);
//     }
// }
//
// //correcting
// void CORRECT_Y(void){
//  int i,j,k;
//  real d,d_prime,d_2_prime,e_inc_d,h_inc_d,e_x_inc,e_z_inc,h_z_inc,h_x_inc;
// ////#pragma omp parallel for collapse(2)
//  for(i=NtfsfX;i<=NCELLX-NtfsfX;i++){
//     for(k=NtfsfZ;k<=NCELLZ-NtfsfZ;k++){
//         j=NtfsfY;
//         d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_prime=d-(int)d;
//         e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//         e_z_inc=e_inc_d*(sin(polar_psi)*sin(inc_theta));
//
//         if(k != NCELLZ-NtfsfZ){
//            hx[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]*e_z_inc/dy;
//         }
//
//         j=NtfsfY;
//         d=(i-i_0+0.5)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_prime=d-(int)d;
//         e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//         e_x_inc=e_inc_d*(cos(polar_psi)*sin(inc_phi)-sin(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//          if(i != NCELLX-NtfsfX ){
//            hz[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]-=Chze[ThreeDMap(i,j-1,k,NCELLZ,NCELLY)]*e_x_inc/dy;
//         }
//
//         j=NtfsfY;
//         d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_2_prime=d+0.5;
//         d_prime=d_2_prime-(int)(d_2_prime);
//         h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//         h_z_inc=h_inc_d*(-cos(polar_psi)*sin(inc_theta));
//
//         if(i != NCELLX-NtfsfX){
//             ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dy;
//         }
//
//
//         j=NtfsfY;
//         d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j-0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_2_prime=d+0.5;
//         d_prime=d_2_prime-(int)(d_2_prime);
//         h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//         h_x_inc=h_inc_d*(sin(polar_psi)*sin(inc_phi)+cos(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//         if(k != NCELLZ-NtfsfZ){
//             ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dy;
//         }
//
//
//         j=NCELLY-NtfsfY;
//         d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_prime=d-(int)d;
//         e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//         e_z_inc=e_inc_d*(sin(polar_psi)*sin(inc_theta));
//
//         if(k != NCELLZ-NtfsfZ){
//            hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Chxe[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_z_inc/dy;
//         }
//
//
//         j=NCELLY-NtfsfY;
//         d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_prime=d-(int)d;
//         e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//         e_x_inc=e_inc_d*(cos(polar_psi)*sin(inc_phi)-sin(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//         if(i != NCELLX-NtfsfX ){
//            hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_x_inc/dy;
//         }
//
//
//         j=NCELLY-NtfsfY;
//         d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_2_prime=d+0.5;
//         d_prime=d_2_prime-(int)(d_2_prime);
//         h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//         h_z_inc=h_inc_d*(-cos(polar_psi)*sin(inc_theta));
//
//         if(i != NCELLX-NtfsfX){
//             ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dy;
//         }
//
//         j=NCELLY-NtfsfY;
//         d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//         d=(dx/d_1D)*d;
//         d_2_prime=d+0.5;
//         d_prime=d_2_prime-(int)(d_2_prime);
//         h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//         h_x_inc=h_inc_d*(sin(polar_psi)*sin(inc_phi)+cos(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//         if(k != NCELLZ-NtfsfZ){
//             ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dy;
//         }
//
//     }
//  }
// }
//
// void CORRECT_X(void){
//  int i,j,k;
//  real d,d_prime,d_2_prime,e_inc_d,h_inc_d,e_y_inc,e_z_inc,h_z_inc,h_y_inc;
//  ////#pragma omp parallel for collapse(2)
//  for(j=NtfsfY;j<=NCELLY-NtfsfY;j++){
//         for(k=NtfsfZ;k<=NCELLZ-NtfsfZ; k++){
//
//             i=NtfsfX;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_z_inc=e_inc_d*(sin(polar_psi)*sin(inc_theta));
//
//             if(k != NCELLZ-NtfsfZ){
//                 hy[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]*e_z_inc/dx;
//             }
//
//
//             i=NtfsfX;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_y_inc=e_inc_d*(-cos(polar_psi)*cos(inc_phi)-sin(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(j != NCELLY-NtfsfY ){
//                 hz[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMap(i-1,j,k,NCELLZ,NCELLY)]*e_y_inc/dx;
//             }
//
//
//             i=NtfsfX;
//             d=(i-0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_z_inc=h_inc_d*(-cos(polar_psi)*sin(inc_theta));
//
//             if(j != NCELLY-NtfsfY){
//                 ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dx;
//             }
//
//
//             i=NtfsfX;
//             d=(i-0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_y_inc=h_inc_d*(-sin(polar_psi)*cos(inc_phi)+cos(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(k != NCELLZ-NtfsfZ){
//                 ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dx;
//             }
//
//
//
//             i=NCELLX-NtfsfX;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_z_inc=e_inc_d*(sin(polar_psi)*sin(inc_theta));
//
//             if(k != NCELLZ-NtfsfZ){
//                 hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Chye[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_z_inc/dx;
//             }
//
//
//             i=NCELLX-NtfsfX;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_y_inc=e_inc_d*(-cos(polar_psi)*cos(inc_phi)-sin(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(j != NCELLY-NtfsfY ){
//                 hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Chze[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_y_inc/dx;
//             }
//
//
//             i=NCELLX-NtfsfX;
//             d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_z_inc=h_inc_d*(-cos(polar_psi)*sin(inc_theta));
//
//             if(j != NCELLY-NtfsfY){
//                 ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dx;
//             }
//
//
//             i=NCELLX-NtfsfX;
//             d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_y_inc=h_inc_d*(-sin(polar_psi)*cos(inc_phi)+cos(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(k != NCELLZ-NtfsfZ){
//                 ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dx;
//             }
//
//     }
//  }
// }
//
// void CORRECT_Z(void){
//  int i,j,k;
//  real d,d_prime,d_2_prime,e_inc_d,h_inc_d,e_y_inc,e_x_inc,h_x_inc,h_y_inc;
// // //#pragma omp parallel for collapse(2)
//  for(i=NtfsfX;i<=NCELLX-NtfsfX;i++){
//     for(j=NtfsfY;j<=NCELLY-NtfsfY;j++){
// //    if(t==2001)printf("%d\t%d\n",i,j);
//             k=NtfsfZ;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//           //  if(t==2001) printf("%f,%d",d,(int)d);
//
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//             e_y_inc=e_inc_d*(-cos(polar_psi)*cos(inc_phi)-sin(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(j != NCELLY-NtfsfY){
//                 hx[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]-=Chxe[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]*e_y_inc/dz;
//             }
//
//
//             k=NtfsfZ;
//             d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_x_inc=e_inc_d*(cos(polar_psi)*sin(inc_phi)-sin(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//             if(i != NCELLX-NtfsfX ){
//                 hy[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]+=Chye[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]*e_x_inc/dz;
//             }
//
//
//
//             k=NtfsfZ;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_x_inc=h_inc_d*(sin(polar_psi)*sin(inc_phi)+cos(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//             if(j != NCELLY-NtfsfY){
//                 ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dz;
//             }
//
//
//             k=NtfsfZ;
//             d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k-0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_y_inc=h_inc_d*(-sin(polar_psi)*cos(inc_phi)+cos(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(i != NCELLX-NtfsfX ){
//                 ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dz;
//             }
//
// if(Periodic_XY == 0){
//             k=NCELLZ-NtfsfZ;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_y_inc=e_inc_d*(-cos(polar_psi)*cos(inc_phi)-sin(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//             if(j != NCELLY-NtfsfY){
//                 hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_y_inc/dz;
//             }
//             if(i != NCELLX-NtfsfX ){
//                 hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_x_inc/dz;
//             }
//
//
//             k=NCELLZ-NtfsfZ;
//             d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_prime=d-(int)d;
//             e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
//
//             e_x_inc=e_inc_d*(cos(polar_psi)*sin(inc_phi)-sin(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//             if(i != NCELLX-NtfsfX ){
//                 hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*e_x_inc/dz;
//             }
//
//
//             k=NCELLZ-NtfsfZ;
//             d=(i-i_0)*sin(inc_theta)*cos(inc_phi)+(j+0.5-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_x_inc=h_inc_d*(sin(polar_psi)*sin(inc_phi)+cos(polar_psi)*cos(inc_theta)*cos(inc_phi));
//
//             if(j != NCELLY-NtfsfY){
//                 ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dz;
//             }
//
//
//             k=NCELLZ-NtfsfZ;
//             d=(i+0.5-i_0)*sin(inc_theta)*cos(inc_phi)+(j-j_0)*sin(inc_theta)*sin(inc_phi)+(k+0.5-k_0)*cos(inc_theta);
//             d=(dx/d_1D)*d;
//             d_2_prime=d+0.5;
//             d_prime=d_2_prime-(int)(d_2_prime);
//             h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];
//
//             h_y_inc=h_inc_d*(-sin(polar_psi)*cos(inc_phi)+cos(polar_psi)*cos(inc_theta)*sin(inc_phi));
//
//              if(i != NCELLX-NtfsfX ){
//                 ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dz;
//             }
// }
//     }
//  }
// }
//
// void SOURCE_IN(void){
//   int i,j,k;
//   k=inc_plane;
//
//   ////#pragma omp parallel for collapse(2) private(i,j) schedule(guided)
//   for(i=0;i<NCELLX-1;i++){
//        for(j=0;j<NCELLY-1;j++){
//            if(TEz){
//              //hy[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)] += (Chye[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]/dz)*PULSE(t)*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//              hx[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)] -= (Chxe[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]/dz)*PULSE(t)*cexp(-I*(i)*dx*k_x)*cexp(-I*(j+0.5)*dy*k_y);
//
//            }
//            if(TMz){
//              ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)] -= (Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]/dz)*PULSE(t)*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//            }
//
//        }
//
//    }
// }
//
// void SOURCE_IN_B(void){
//    int i,j,k;
//    k=inc_plane;
//
//    ////#pragma omp parallel for collapse(2) private(i,j) schedule(guided)
//    for(i=0;i<NCELLX;i++){
//         for(j=0;j<NCELLY;j++){
//             if(TEz){
//               hy[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)] += (Chye[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]/dz)*Ex_t_FT[t]*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//               hx[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)] -= (Chxe[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]/dz)*Ey_t_FT[t]*cexp(-I*(i)*dx*k_x)*cexp(-I*(j+0.5)*dy*k_y);
//
//             }
//             if(TMz){
//               hy[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)] -= (Chye[ThreeDMap(i,j,k-1,NCELLZ,NCELLY)]/dz)*Ex_t_FT[t]*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//             }
//
//         }
//
//     }
// }
//
// void SOURCE_IN_E(void){
//   int i,j,k;
//   k=inc_plane;
//
//
//   ////#pragma omp parallel for collapse(2) private(i,j) schedule(guided)
//   for(i=0;i<NCELLX;i++){
//        for(j=0;j<NCELLY;j++){
//            if(TEz){
//               ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)] -= (Ceyh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]/dz)*0.5*(Hx_t_FT[t]+Hx_t_FT[t+1])*cexp(-I*(i)*dx*k_x)*cexp(-I*(j+0.5)*dy*k_y);
//            }
//            if(TMz){
//               ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)] += (Ceyh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]/dz)*0.5*(Hx_t_FT[t+1]+Hx_t_FT[t])*cexp(-I*(i)*dx*k_x)*cexp(-I*(j+0.5)*dy*k_y);
//               ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)] += (Cexh[ThreeDMap(i,j,k,NCELLZ,NCELLY)]/dz)*0.5*(Hy_t_FT[t+1]+Hy_t_FT[t])*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//            }
//
//        }
//
//    }
// }
//
// comp PULSE(int t){
//   comp pulse;
//   //pulse=exp(-pow((real)(t-delay)/(real)width,2)/2.0)*cexp(I*2*pi*f_0*(t)*dt);
//   if(t < 3*(delay+2*width)) pulse = exp(-pow((real)(t-delay)/(real)width,2)/2.0)*cos(2*pi*f_0*(t)*dt);
//   else pulse = 0;
//   // if(t>3*(delay+2*width)){
//   //   return 0;
//   // }
// return (pulse);
// }
