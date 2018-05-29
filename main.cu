#include <stdlib.h>
#include <stdio.h>
#include "global_var.h"
#include <math.h>
#include <time.h>
#include <cuda.h>
// #include <cuda_runtime.h>
// #include <device_launch_parameters.h>
//#include<conio.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <cuda_runtime_api.h>
// #include <cuda.h>
#include "cuda_profiler_api.h"
#include <cudaProfiler.h>
__device__  int ThreeDMapD(int i,int j,int k,int SizeZ,int SizeY){
 int num = k + SizeZ*j +SizeY*SizeZ*i;
 return num;
}


__device__  int FourDMapD(int i,int j,int k,int n,int SizeN,int SizeZ,int SizeY){
 int num = n + SizeN*( k + SizeZ*j +SizeY*SizeZ*i);
 return num;
}

__device__  int TwoDMapD(int i,int j,int size){
 int num = j + i*size;
 return num;
}




    __global__ void ScattAbs(real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,int NUM_freq,int t,real dt,real* freq,real pi,int XSTARTAbs,int XENDAbs,int YSTARTAbs,int YENDAbs,int ZSTARTAbs,int ZENDAbs,int XSTARTSca,int XENDSca,int YSTARTSca,int YENDSca,int ZSTARTSca,int ZENDSca,int XNEARAbs,int XFARAbs,int YNEARAbs,int YFARAbs,int ZNEARAbs,int ZFARAbs,
    int XNEARSca,int XFARSca,int YNEARSca,int YFARSca,int ZNEARSca,int ZFARSca,real2 *ExTransformNearZAbsRe,real2 *ExTransformNearZAbsIm,real2 *EyTransformNearZAbsRe,real2 *EyTransformNearZAbsIm,real2 *HxTransformNearZAbsRe,real2 *HxTransformNearZAbsIm,real2 *HyTransformNearZAbsRe,real2 *HyTransformNearZAbsIm,
    real2 *ExTransformFarZAbsRe,real2 *ExTransformFarZAbsIm,real2 *EyTransformFarZAbsRe,real2 *EyTransformFarZAbsIm,real2 *HxTransformFarZAbsRe,real2 *HxTransformFarZAbsIm,real2 *HyTransformFarZAbsRe,real2 *HyTransformFarZAbsIm,
    real2 *ExTransformNearYAbsRe,real2 *ExTransformNearYAbsIm,real2 *EzTransformNearYAbsRe,real2 *EzTransformNearYAbsIm,real2 *HxTransformNearYAbsRe,real2 *HxTransformNearYAbsIm,real2 *HzTransformNearYAbsRe,real2 *HzTransformNearYAbsIm,
    real2 *ExTransformFarYAbsRe,real2 *ExTransformFarYAbsIm,real2 *EzTransformFarYAbsRe,real2 *EzTransformFarYAbsIm,real2 *HxTransformFarYAbsRe,real2 *HxTransformFarYAbsIm,real2 *HzTransformFarYAbsRe,real2 *HzTransformFarYAbsIm,
    real2 *EyTransformNearXAbsRe,real2 *EyTransformNearXAbsIm,real2 *EzTransformNearXAbsRe,real2 *EzTransformNearXAbsIm,real2 *HyTransformNearXAbsRe,real2 *HyTransformNearXAbsIm,real2 *HzTransformNearXAbsRe,real2 *HzTransformNearXAbsIm,
    real2 *EyTransformFarXAbsRe,real2 *EyTransformFarXAbsIm,real2 *EzTransformFarXAbsRe,real2 *EzTransformFarXAbsIm,real2 *HyTransformFarXAbsRe,real2 *HyTransformFarXAbsIm,real2 *HzTransformFarXAbsRe,real2 *HzTransformFarXAbsIm,
    real2 *ExTransformNearZScaRe,real2 *ExTransformNearZScaIm,real2 *EyTransformNearZScaRe,real2 *EyTransformNearZScaIm,real2 *HxTransformNearZScaRe,real2 *HxTransformNearZScaIm,real2 *HyTransformNearZScaRe,real2 *HyTransformNearZScaIm,
    real2 *ExTransformFarZScaRe,real2 *ExTransformFarZScaIm,real2 *EyTransformFarZScaRe,real2 *EyTransformFarZScaIm,real2 *HxTransformFarZScaRe,real2 *HxTransformFarZScaIm,real2 *HyTransformFarZScaRe,real2 *HyTransformFarZScaIm,
    real2 *ExTransformNearYScaRe,real2 *ExTransformNearYScaIm,real2 *EzTransformNearYScaRe,real2 *EzTransformNearYScaIm,real2 *HxTransformNearYScaRe,real2 *HxTransformNearYScaIm,real2 *HzTransformNearYScaRe,real2 *HzTransformNearYScaIm,
    real2 *ExTransformFarYScaRe,real2 *ExTransformFarYScaIm,real2 *EzTransformFarYScaRe,real2 *EzTransformFarYScaIm,real2 *HxTransformFarYScaRe,real2 *HxTransformFarYScaIm,real2 *HzTransformFarYScaRe,real2 *HzTransformFarYScaIm,
    real2 *EyTransformNearXScaRe,real2 *EyTransformNearXScaIm,real2 *EzTransformNearXScaRe,real2 *EzTransformNearXScaIm,real2 *HyTransformNearXScaRe,real2 *HyTransformNearXScaIm,real2 *HzTransformNearXScaRe,real2 *HzTransformNearXScaIm,
    real2 *EyTransformFarXScaRe,real2 *EyTransformFarXScaIm,real2 *EzTransformFarXScaRe,real2 *EzTransformFarXScaIm,real2 *HyTransformFarXScaRe,real2 *HyTransformFarXScaIm,real2 *HzTransformFarXScaRe,real2 *HzTransformFarXScaIm,int NCELLX,int NCELLY,int NCELLZ){


    int freq_count,i,j,k,II,JJ,KK;
    real TransVecERe;
    real TransVecHRe;
    real TransVecEIm;
    real TransVecHIm;

    int idx = blockDim.x * blockIdx.x + threadIdx.x;
  //  freq_count = blockIdx.y;

    i = idx / (NCELLZ*NCELLY);
    j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
    k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
    for(freq_count=0;freq_count<NUM_freq;freq_count++){
    //if(freq_count < NUM_freq){
    #ifdef FlOATPRECISION
    TransVecERe = cosf(2.0*pi*(t+1.0)*dt*freq[freq_count]);
    TransVecEIm = sinf(2.0*pi*(t+1.0)*dt*freq[freq_count]);
    TransVecHRe = cosf(2.0*pi*(t+0.5)*dt*freq[freq_count]);
    TransVecHIm = sinf(2.0*pi*(t+0.5)*dt*freq[freq_count]);
    #endif

    #ifdef DOUBLEPRECISION
    TransVecERe = cos((real2)2.0*pi*(t+1.0)*dt*freq[freq_count]);
    TransVecEIm = sin((real2)2.0*pi*(t+1.0)*dt*freq[freq_count]);
    TransVecHRe = cos((real2)2.0*pi*(t+0.5)*dt*freq[freq_count]);
    TransVecHIm = sin((real2)2.0*pi*(t+0.5)*dt*freq[freq_count]);
    #endif
    //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
      // for(i=XSTARTAbs;i<XENDAbs;i++){
      //   for(j=YSTARTAbs;j<YENDAbs;j++){
          if(i>=XSTARTAbs && i<XENDAbs && j>=YSTARTAbs && j<YENDAbs && k == ZNEARAbs){
        //  printf("%d\t%d\t%d\t%d\n",i,j,XENDAbs,YENDAbs);
          //k = ZNEARAbs;
          II = i - XSTARTAbs;
          JJ = j - YSTARTAbs;

          ExTransformNearZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2) ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformNearZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
        //  EzTransformNearZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformNearZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformNearZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
      //    EzTransformNearZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformNearZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformNearZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
        //  HzTransformNearZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformNearZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformNearZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HzTransformNearZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }
          if(i>=XSTARTAbs && i<XENDAbs && j>=YSTARTAbs && j<YENDAbs && k == ZFARAbs){
          // k = ZFARAbs;
          II = i - XSTARTAbs;
          JJ = j - YSTARTAbs;

          ExTransformFarZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformFarZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
      //    EzTransformFarZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformFarZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] +=(real2) ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformFarZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
    //      EzTransformFarZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformFarZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformFarZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
      //    HzTransformFarZAbsRe[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformFarZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformFarZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HzTransformFarZAbsIm[ThreeDMapD(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
    }

      //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
      // for(i=XSTARTAbs;i<XENDAbs;i++){
      //   for(k=ZSTARTAbs;k<ZENDAbs;k++){
          if(i>=XSTARTAbs && i<XENDAbs && k>=ZSTARTAbs && k<ZENDAbs && j == YNEARAbs){

          II = i - XSTARTAbs;
          KK = k - ZSTARTAbs;

          // j=YNEARAbs;

          ExTransformNearYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
      //    EyTransformNearYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformNearYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformNearYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
    //      EyTransformNearYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformNearYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] +=(real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformNearYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] +=(real2) hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
      //    HyTransformNearYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformNearYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] +=(real2) hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformNearYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HyTransformNearYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformNearYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }

        if(i>=XSTARTAbs && i<XENDAbs && k>=ZSTARTAbs && k<ZENDAbs && j == YFARAbs){
          II = i - XSTARTAbs;
          KK = k - ZSTARTAbs;
          // j=YFARAbs;

          ExTransformFarYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
    //      EyTransformFarYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformFarYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformFarYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
    //      EyTransformFarYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformFarYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformFarYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
    //      HyTransformFarYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformFarYAbsRe[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformFarYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HyTransformFarYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformFarYAbsIm[ThreeDMapD(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }
    //  }

      //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
      if(j>=YSTARTAbs && j<YENDAbs && k>=ZSTARTAbs && k<ZENDAbs && i == XNEARAbs){

      // for(j=YSTARTAbs;j<YENDAbs;j++){
      //   for(k=ZSTARTAbs;k<ZENDAbs;k++){
          KK = k - ZSTARTAbs;
          JJ = j - YSTARTAbs;

          // i=XNEARAbs;

      //    ExTransformNearXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformNearXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] +=(real2) ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformNearXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

    //        ExTransformNearXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformNearXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformNearXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

    //      HxTransformNearXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformNearXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] +=(real2) hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformNearXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

      //    HxTransformNearXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformNearXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformNearXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }

          if(j>=YSTARTAbs && j<YENDAbs && k>=ZSTARTAbs && k<ZENDAbs && i == XFARAbs){
            KK = k - ZSTARTAbs;
            JJ = j - YSTARTAbs;
          // i=XFARAbs;

      //    ExTransformFarXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformFarXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformFarXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

    //      ExTransformFarXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformFarXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformFarXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

      //    HxTransformFarXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformFarXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] +=(real2) hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformFarXAbsRe[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

        //  HxTransformFarXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformFarXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformFarXAbsIm[ThreeDMapD(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }
    //  }

    if(i>=XSTARTSca && i<XENDSca && j>=YSTARTSca && j<YENDSca && k == ZNEARSca){
      // for(i=XSTARTSca;i<XENDSca;i++){
      //   for(j=YSTARTSca;j<YENDSca;j++){
      //     k = ZNEARSca;
          II = i - XSTARTSca;
          JJ = j - YSTARTSca;

          ExTransformNearZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformNearZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
      //    EzTransformNearZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformNearZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformNearZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
      //    EzTransformNearZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformNearZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformNearZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
      //    HzTransformNearZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformNearZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformNearZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HzTransformNearZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
    }
    if(i>=XSTARTSca && i<XENDSca && j>=YSTARTSca && j<YENDSca && k == ZFARSca){
      II = i - XSTARTSca;
      JJ = j - YSTARTSca;
          // k = ZFARSca;

          ExTransformFarZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformFarZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
      //    EzTransformFarZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformFarZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformFarZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
      //    EzTransformFarZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformFarZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformFarZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
      //    HzTransformFarZScaRe[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformFarZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformFarZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HzTransformFarZScaIm[ThreeDMapD(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }
    //  }

      //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
      if(i>=XSTARTSca && i<XENDSca && k>=ZSTARTSca && k<ZENDSca && j == YNEARSca){
      // for(i=XSTARTSca;i<XENDSca;i++){
      //   for(k=ZSTARTSca;k<ZENDSca;k++){

          // j=YNEARSca;
          II = i - XSTARTSca;
          KK = k - ZSTARTSca;
          ExTransformNearYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
        //  EyTransformNearYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformNearYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformNearYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
      //    EyTransformNearYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformNearYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformNearYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
        //  HyTransformNearYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformNearYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformNearYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HyTransformNearYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformNearYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
    }
    if(i>=XSTARTSca && i<XENDSca && k>=ZSTARTSca && k<ZENDSca && j == YFARSca){

          // j=YFARSca;
          II = i - XSTARTSca;
          KK = k - ZSTARTSca;
          ExTransformFarYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
      //    EyTransformFarYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformFarYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

          ExTransformFarYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
      //    EyTransformFarYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformFarYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

          HxTransformFarYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
      //    HyTransformFarYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformFarYScaRe[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

          HxTransformFarYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //    HyTransformFarYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformFarYScaIm[ThreeDMapD(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] +=(real2) hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
      //  }
      }

      //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
      if(j>=YSTARTSca && j<YENDSca && k>=ZSTARTSca && k<ZENDSca && i == XNEARSca){

      // for(j=YSTARTSca;j<YENDSca;j++){
      //   for(k=ZSTARTSca;k<ZENDSca;k++){
          KK = k - ZSTARTSca;
          JJ = j - YSTARTSca;
          // i=XNEARSca;

      //    ExTransformNearXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformNearXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] +=(real2) ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformNearXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] +=(real2) ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

    //        ExTransformNearXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformNearXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] +=(real2) ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformNearXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

    //        HxTransformNearXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformNearXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformNearXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

    //        HxTransformNearXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformNearXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformNearXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }
          if(j>=YSTARTSca && j<YENDSca && k>=ZSTARTSca && k<ZENDSca && i == XFARSca){

          // i=XFARSca;
          KK = k - ZSTARTSca;
          JJ = j - YSTARTSca;
    //      ExTransformFarXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EyTransformFarXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;
          EzTransformFarXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecERe;

      //    ExTransformFarXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EyTransformFarXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;
          EzTransformFarXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecEIm;

      //    HxTransformFarXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HyTransformFarXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;
          HzTransformFarXScaRe[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHRe;

      //    HxTransformFarXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HyTransformFarXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
          HzTransformFarXScaIm[ThreeDMapD(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += (real2)hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*TransVecHIm;
        }
    //  }
    }
    return;
    }

//correcting
__global__ void CORRECT_Y(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real *e_inc,real *h_inc,real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,real inc_theta,
  real inc_phi,real polar_psi,real polar_theta,real dx,real dy,real dz,real dt,int i_0,int j_0,int k_0,real d_1D,int m0,real *Cexe,real *Ceye,real *Ceze,real *Cexh,real *Ceyh,real *Cezh,real *Chxe,real *Chye,real *Chze,real *Chxh,real *Chyh,real *Chzh,int Periodic_XY){

 int i,j,k;

 int idx = blockDim.x * blockIdx.x + threadIdx.x;

 i = idx / (NCELLZ*NCELLY);
 j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
 k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
 real d,d_prime,d_2_prime,e_inc_d,h_inc_d,e_x_inc,e_z_inc,h_z_inc,h_x_inc;

////#pragma omp parallel for collapse(2)
 // for(i=NtfsfX;i<=NCELLX-NtfsfX;i++){
 //    for(k=NtfsfZ;k<=NCELLZ-NtfsfZ;k++){
      if(i>=NtfsfX && i<=NCELLX-NtfsfX && k>=NtfsfZ && k<=NCELLZ-NtfsfZ){
        // j=NtfsfY;
        if(j==NtfsfY){
        d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_prime=d-(int)d;
        e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

        e_z_inc=e_inc_d*(sin((float)polar_psi)*sin((float)inc_theta));

        if(k != NCELLZ-NtfsfZ){
           hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]*e_z_inc/dy;
        }

        // j=NtfsfY;
        d=(i-i_0+0.5)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_prime=d-(int)d;
        e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

        e_x_inc=e_inc_d*(cos((float)polar_psi)*sin((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

         if(i != NCELLX-NtfsfX ){
           hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]-=Chze[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]*e_x_inc/dy;
        }

        // j=NtfsfY;
        d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_2_prime=d+0.5;
        d_prime=d_2_prime-(int)(d_2_prime);
        h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

        h_z_inc=h_inc_d*(-cos((float)polar_psi)*sin((float)inc_theta));

        if(i != NCELLX-NtfsfX){
            ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dy;
        }


        // j=NtfsfY;
        d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_2_prime=d+0.5;
        d_prime=d_2_prime-(int)(d_2_prime);
        h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

        h_x_inc=h_inc_d*(sin((float)polar_psi)*sin((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

        if(k != NCELLZ-NtfsfZ){
            ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dy;
        }

      }
      if(j == NCELLY-NtfsfY){
        // j=NCELLY-NtfsfY;
        d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_prime=d-(int)d;
        e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

        e_z_inc=e_inc_d*(sin((float)polar_psi)*sin((float)inc_theta));

        if(k != NCELLZ-NtfsfZ){
           hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_z_inc/dy;
        }


        // j=NCELLY-NtfsfY;
        d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_prime=d-(int)d;
        e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

        e_x_inc=e_inc_d*(cos((float)polar_psi)*sin((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

        if(i != NCELLX-NtfsfX ){
           hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_x_inc/dy;
        }


        // j=NCELLY-NtfsfY;
        d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_2_prime=d+0.5;
        d_prime=d_2_prime-(int)(d_2_prime);
        h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

        h_z_inc=h_inc_d*(-cos((float)polar_psi)*sin((float)inc_theta));

        if(i != NCELLX-NtfsfX){
            ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dy;
        }

        // j=NCELLY-NtfsfY;
        d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
        d=(dx/d_1D)*d;
        d_2_prime=d+0.5;
        d_prime=d_2_prime-(int)(d_2_prime);
        h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

        h_x_inc=h_inc_d*(sin((float)polar_psi)*sin((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

        if(k != NCELLZ-NtfsfZ){
            ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dy;
        }
}
    // }
 }
}

__global__ void CORRECT_X(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real *e_inc,real *h_inc,real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,real inc_theta,
  real inc_phi,real polar_psi,real polar_theta,real dx,real dy,real dz,real dt,int i_0,int j_0,int k_0,real d_1D,int m0,real *Cexe,real *Ceye,real *Ceze,real *Cexh,real *Ceyh,real *Cezh,real *Chxe,real *Chye,real *Chze,real *Chxh,real *Chyh,real *Chzh,int Periodic_XY){
 int i,j,k;
 ////#pragma omp parallel for collapse(2)
 int idx = blockDim.x * blockIdx.x + threadIdx.x;

 i = idx / (NCELLZ*NCELLY);
 j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
 k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
 real d,d_prime,d_2_prime,e_inc_d,h_inc_d,e_y_inc,e_z_inc,h_z_inc,h_y_inc;

 // for(j=NtfsfY;j<=NCELLY-NtfsfY;j++){
 //        for(k=NtfsfZ;k<=NCELLZ-NtfsfZ; k++){
          if(j>=NtfsfY && j<=NCELLY-NtfsfY && k>=NtfsfZ && k<=NCELLZ-NtfsfZ){
            if(i==NtfsfX){
            // i=NtfsfX;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_z_inc=e_inc_d*(sin((float)polar_psi)*sin((float)inc_theta));

            if(k != NCELLZ-NtfsfZ){
                hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]*e_z_inc/dx;
            }


            // i=NtfsfX;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_y_inc=e_inc_d*(-cos((float)polar_psi)*cos((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(j != NCELLY-NtfsfY ){
                hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]*e_y_inc/dx;
            }


            // i=NtfsfX;
            d=(i-0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_z_inc=h_inc_d*(-cos((float)polar_psi)*sin((float)inc_theta));

            if(j != NCELLY-NtfsfY){
                ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dx;
            }


            // i=NtfsfX;
            d=(i-0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_y_inc=h_inc_d*(-sin((float)polar_psi)*cos((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(k != NCELLZ-NtfsfZ){
                ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dx;
            }


}
          if(i==NCELLX-NtfsfX){
            // i=NCELLX-NtfsfX;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_z_inc=e_inc_d*(sin((float)polar_psi)*sin((float)inc_theta));

            if(k != NCELLZ-NtfsfZ){
                hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_z_inc/dx;
            }


            // i=NCELLX-NtfsfX;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_y_inc=e_inc_d*(-cos((float)polar_psi)*cos((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(j != NCELLY-NtfsfY ){
                hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_y_inc/dx;
            }


            // i=NCELLX-NtfsfX;
            d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_z_inc=h_inc_d*(-cos((float)polar_psi)*sin((float)inc_theta));

            if(j != NCELLY-NtfsfY){
                ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_z_inc/dx;
            }


            // i=NCELLX-NtfsfX;
            d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_y_inc=h_inc_d*(-sin((float)polar_psi)*cos((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(k != NCELLZ-NtfsfZ){
                ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dx;
            }
}
    // }
 }
}

__global__ void CORRECT_Z(int NtfsfX,int NtfsfY,int NtfsfZ,int NCELLX,int NCELLY,int NCELLZ,real *e_inc,real *h_inc,real *ex,real *ey,real *ez,real *hx,real *hy,real *hz,real inc_theta,
  real inc_phi,real polar_psi,real polar_theta,real dx,real dy,real dz,real dt,int i_0,int j_0,int k_0,real d_1D,int m0,real *Cexe,real *Ceye,real *Ceze,real *Cexh,real *Ceyh,real *Cezh,real *Chxe,real *Chye,real *Chze,real *Chxh,real *Chyh,real *Chzh,int Periodic_XY){

 int i,j,k;
// //#pragma omp parallel for collapse(2)
 // for(i=NtfsfX;i<=NCELLX-NtfsfX;i++){
 //    for(j=NtfsfY;j<=NCELLY-NtfsfY;j++){

      int idx = blockDim.x * blockIdx.x + threadIdx.x;

      i = idx / (NCELLZ*NCELLY);
      j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
      k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
      real d,d_prime,d_2_prime,e_inc_d,h_inc_d,e_y_inc,e_x_inc,h_x_inc,h_y_inc;

       if(j>=NtfsfY && j<=NCELLY-NtfsfY && i>=NtfsfX && i<=NCELLX-NtfsfX){

//    if(t==2001)printf("%d\t%d\n",i,j);
          if(k==NtfsfZ){
            // k=NtfsfZ;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
          //  if(t==2001) printf("%f,%d",d,(int)d);

            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];
            e_y_inc=e_inc_d*(-cos((float)polar_psi)*cos((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(j != NCELLY-NtfsfY){
                hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]-=Chxe[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]*e_y_inc/dz;
            }


            // k=NtfsfZ;
            d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_x_inc=e_inc_d*(cos((float)polar_psi)*sin((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

            if(i != NCELLX-NtfsfX ){
                hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]+=Chye[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]*e_x_inc/dz;
            }



            // k=NtfsfZ;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_x_inc=h_inc_d*(sin((float)polar_psi)*sin((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

            if(j != NCELLY-NtfsfY){
                ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dz;
            }


            // k=NtfsfZ;
            d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_y_inc=h_inc_d*(-sin((float)polar_psi)*cos((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(i != NCELLX-NtfsfX ){
                ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dz;
            }
          }

if(Periodic_XY == 0){
            if(k==NCELLZ-NtfsfZ){
            // k=NCELLZ-NtfsfZ;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_y_inc=e_inc_d*(-cos((float)polar_psi)*cos((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

            if(j != NCELLY-NtfsfY){
                hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_y_inc/dz;
            }
            // if(i != NCELLX-NtfsfX ){
            //     hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_x_inc/dz;
            // }


            // k=NCELLZ-NtfsfZ;
            d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_prime=d-(int)d;
            e_inc_d=(1-d_prime)*e_inc[m0+(int)d]+d_prime*e_inc[m0+(int)d+1];

            e_x_inc=e_inc_d*(cos((float)polar_psi)*sin((float)inc_phi)-sin((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

            if(i != NCELLX-NtfsfX ){
                hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*e_x_inc/dz;
            }


            // k=NCELLZ-NtfsfZ;
            d=(i-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j+0.5-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_x_inc=h_inc_d*(sin((float)polar_psi)*sin((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*cos((float)inc_phi));

            if(j != NCELLY-NtfsfY){
                ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_x_inc/dz;
            }


            // k=NCELLZ-NtfsfZ;
            d=(i+0.5-i_0)*sin((float)inc_theta)*cos((float)inc_phi)+(j-j_0)*sin((float)inc_theta)*sin((float)inc_phi)+(k+0.5-k_0)*cos((float)inc_theta);
            d=(dx/d_1D)*d;
            d_2_prime=d+0.5;
            d_prime=d_2_prime-(int)(d_2_prime);
            h_inc_d=(1-d_prime)*h_inc[m0-1+(int)d_2_prime]+d_prime*h_inc[m0+(int)d_2_prime];

            h_y_inc=h_inc_d*(-sin((float)polar_psi)*cos((float)inc_phi)+cos((float)polar_psi)*cos((float)inc_theta)*sin((float)inc_phi));

             if(i != NCELLX-NtfsfX ){
                ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*h_y_inc/dz;
            }
}
    }
 }
}



  __global__ void UPDATE_hx(real *hx,real *hxPrev,real *ez,real *ey,real *Chxh,real *Chxe,real *psi_Hx_z_N,real *psi_Hx_z_F,real *psi_Hx_y_N,real *psi_Hx_y_F,real *khdy,real
    *khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlY,int Hydrodynamics){
//void UPDATE_hx(void){
// cudaProfilerStart();

      int i,j,k,j2,k2;
      comp Curl_E;
      // __shared__ double Ey[110][110][110];
      // __shared__ double Ez[110][110][110];

      int idx = blockDim.x * blockIdx.x + threadIdx.x;

      i = idx / (NCELLZ*NCELLY);
      j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
      k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
      // Ey[i][j][k] = ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
      // Ez[i][j][k] = ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];


      if(Periodic_XY){
  	////#pragma omp parallel for collapse(3) private(i,j,k,Curl_E,j2,k2) // schedule(static)
  	// for(k=0;k<NCELLZ;k++){
  	// 	for(i=0;i<NCELLX;i++){
  	//         for(j=0;j<NCELLY;j++){

              if(i<NCELLX && j<NCELLY && k<NCELLZ){
  	           //     for(k=0;k<NCELLZ-1;k++){
  									 //if(i==1) printf("%d %d %d\n",i,j,k);
  		                	if(j==NCELLY-1){
  												#ifdef DOUBLECOMPLEX
  		                				Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]*cexp(-I*k_y*period_y)-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
  												#endif
  												#ifdef DOUBLEPRECISION
  		                				Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
  												#endif
  		                   	    hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chxh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;

  		                	}
  		                	 else{
  														Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
  		                   	    hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chxh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  		                	}
  										 //Z-CPML
  										 if(k<cpml_N_Z && i<cpml_x_lim && j<cpml_y_lim){
  											 	//Near Z-PML
  											 		psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  													hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
  										 }
  										 if(k>=cpml_F_Z && j<cpml_y_lim && i<cpml_x_lim){
  											 //Far Z-PML
  											 		k2 = k - cpml_F_Z;
  													psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  													hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
  										 }
                     }
  	    //             }
  	    //     }
  	    // }
      }

      else{
  			//  //#pragma omp target device(0) MapD(Chxe[:NCELLX-1][:NCELLY-1][:NCELLZ-1],Chxh[:NCELLX-1][:NCELLY-1][:NCELLZ-1],ez[:NCELLX-1][:NCELLY-1][:NCELLZ-1],ey[:NCELLX-1][:NCELLY-1][:NCELLZ-1],khdy[:NCELLY-1],khdz[:NCELLZ-1],bh_z_N[:NcpmlZ-1],bh_z_F[:NcpmlZ-1],ch_z_N[:NcpmlZ-1],ch_z_F[:NcpmlZ-1],bh_y_N[:NcpmlY-1],bh_y_F[:NcpmlY-1],ch_y_N[:NcpmlY-1],ch_y_F[:NcpmlY-1]) 		MapD(tofrom:hx[:NCELLX-1][:NCELLY-1][:NCELLZ-1],psi_Hx_z_N[:NCELLX-1][:NCELLY-1][:cpml_N_Z-1],psi_Hx_z_F[:NCELLX-1][:NCELLY-1][:cpml_N_Z-1],psi_Hx_y_N[:NCELLX-1][:cpml_N_Y-1][:NCELLZ-1],psi_Hx_y_F[:NCELLX-1][:cpml_N_Y-1][:NCELLZ-1])
  			//  {
  			// //#pragma omp parallel for collapse(3) private(i,j,k,Curl_E,j2,k2) // schedule(static)
  	    // for(i=0;i<NCELLX;i++){
  	    //     for(j=0;j<NCELLY-1;j++){
  	    //             for(k=0;k<NCELLZ-1;k++){
                      if(i<NCELLX && j<(NCELLY-1) && k<(NCELLZ-1)){
                        if(Hydrodynamics >= 1) hxPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  	                    Curl_E=(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k]-(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j];
  	                    hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chxh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  											//Z-CPML
  											if(k<cpml_N_Z && i<cpml_x_lim && j<cpml_y_lim){
  												 //Near Z-PML
  													 psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  													 hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
  											}
  											if(k>=cpml_F_Z && j<cpml_y_lim && i<cpml_x_lim){
  												//Far Z-PML
  													 k2 = k - cpml_F_Z;
  													 psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ey[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  													 hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
  											}
  										   //Y- PML
  											if(j<cpml_N_Y && i<cpml_x_lim && j<cpml_y_lim){
  													 psi_Hx_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]=bh_y_N[j]*psi_Hx_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]+ch_y_N[j]*(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
  													 hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)];
  											}
  											if(j>=cpml_F_Y && i<cpml_x_lim && k<cpml_z_lim){
  													j2 = j - cpml_F_Y;
  													psi_Hx_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]=bh_y_F[j2]*psi_Hx_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]+ch_y_F[j2]*(ez[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
  													hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chxe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hx_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)];
  											}
                      }
  	  //               }
  	  //       }
  	  //  // }
  		// }
    }
    // cudaProfilerStop();

    return;
  }

  __global__   void UPDATE_hy(real *hy,real *hyPrev,real *ez,real *ex,real *Chyh,real *Chye,real *psi_Hy_z_N,real *psi_Hy_z_F,real *psi_Hy_x_N,real *psi_Hy_x_F,real *khdx,real
    *khdz,real *bh_z_N,real *bh_z_F,real *ch_z_N,real *ch_z_F,real *bh_x_N,real *bh_x_F,real *ch_x_N,real *ch_x_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_Z,int cpml_F_Z,int cpml_N_X,int cpml_F_X,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlZ,int NcpmlX,int Hydrodynamics){

      int i,j,k,n,i2,k2;
      comp Curl_E;
      int idx = blockDim.x * blockIdx.x + threadIdx.x;


      i = idx / (NCELLZ*NCELLY);
      j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
      k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
      if(Periodic_XY){
  		////#pragma omp parallel for collapse(3) private(Curl_E,i,i2,j,k,k2) // schedule(static)
  		// for(k=0;k<NCELLZ;k++){
  		// for(i=0;i<NCELLX;i++){
  	  //       for(j=0;j<NCELLY;j++){
              if(i<NCELLX && j<NCELLY && k<NCELLZ){
  	             //   for(k=0;k<NCELLZ-1;k++){
  	                	if(i==NCELLX-1){
  											#ifdef DOUBLECOMPLEX
  	                		 Curl_E=(ez[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]*cexp(-I*k_x*period_x)-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
  											 #endif
  											 #ifdef DOUBLEPRECISION
  											 Curl_E=(ez[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
  											 #endif
  	                    	 hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  	                	}
  	                	else{
  	                		 Curl_E=(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
  	                    	 hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  	                }
  										//Z-PML
  										if(k<cpml_N_Z && j<cpml_y_lim && k<cpml_z_lim){
  											psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  											hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
  										}
  										if(k>=cpml_F_Z && j<cpml_y_lim && k<cpml_z_lim){
  											k2 = k - cpml_F_Z;
  											psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  											hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
  										}
                    }
  	    //             }
  	    //     }
  	    // }
      }

      else{
  			////#pragma omp parallel for collapse(3) private(Curl_E,i,i2,j,k,k2) // schedule(static)
  	    // for(i=0;i<NCELLX-1;i++){
  	    //     for(j=0;j<NCELLY;j++){
  	    //             for(k=0;k<NCELLZ-1;k++){
                      if(i<(NCELLX-1) && j<NCELLY && k<(NCELLZ-1)){
                        if(Hydrodynamics >= 1)hyPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  	                    Curl_E=(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i]-(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdz[k];
  	                    hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;

  											//Z-PML
  											if(k<cpml_N_Z && j<cpml_y_lim && k<cpml_z_lim){
  												psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]=bh_z_N[k]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)]+ch_z_N[k]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_N[ThreeDMapD(i,j,k,NcpmlZ,NCELLY)];
  											}
  											if(k>=cpml_F_Z && j<cpml_y_lim && k<cpml_z_lim){
  												k2 = k - cpml_F_Z;
  												psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]=bh_z_F[k2]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)]+ch_z_F[k2]*(ex[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dz;
  												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_z_F[ThreeDMapD(i,j,k2,NcpmlZ,NCELLY)];
  											}
  											//X-PML
  											if(i<cpml_N_X && j<cpml_y_lim && k<cpml_z_lim){
  												psi_Hy_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=bh_x_N[i]*psi_Hy_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ch_x_N[i]*(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
  												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  											}
  											if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){
  												i2 = i - cpml_F_X;
  												psi_Hy_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=bh_x_F[i2]*psi_Hy_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ch_x_F[i2]*(ez[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
  												hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hy_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
  											}
                      }
  	    //             }
  	    //     }
  	    // }
      }
    return;
  }

  //void UPDATE_hz(void){
  __global__   void UPDATE_hz(real *hz,real *hzPrev,real *ey,real *ex,real *Chzh,real *Chze,real *psi_Hz_x_N,real *psi_Hz_x_F,real *psi_Hz_y_N,real *psi_Hz_y_F,real *khdx,real
    *khdy,real *bh_x_N,real *bh_x_F,real *ch_x_N,real *ch_x_F,real *bh_y_N,real *bh_y_F,real *ch_y_N,real *ch_y_F,int NCELLX,int NCELLY,int NCELLZ,int Periodic_XY,real dx,real dy,real dz,real dt,int cpml_N_X,int cpml_F_X,int cpml_N_Y,int cpml_F_Y,int cpml_z_lim,int cpml_y_lim,int cpml_x_lim,int NcpmlY,int NcpmlX,int Hydrodynamics){

      int i,j,k,j2,i2;
      comp Curl_E;
      int idx = blockDim.x * blockIdx.x + threadIdx.x;

      i = idx / (NCELLZ*NCELLY);
      j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
      k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
      if(Periodic_XY){
    //  //#pragma omp parallel for collapse(3) private(Curl_E,i,j,k,j2,i2) // schedule(static)
  		// for(k=0;k<NCELLZ;k++){
  	  //   for(i=0;i<NCELLX;i++){
  	  //       for(j=0;j<NCELLY;j++){
          if(i<NCELLX && j<NCELLY && k<NCELLZ){
  	          //      for(k=0;k<NCELLZ;k++){
  	                	if(i==NCELLX-1 || j== NCELLY-1){
  	                		if(i==NCELLX-1 && j==NCELLY-1){
  	                			//printf("%d,%d,%d\n",i,j,k);
  												#ifdef DOUBLECOMPLEX
  	                			 Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]*cexp(-I*k_y*period_y)-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]*cexp(-I*k_x*period_x)-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  												 #endif
  												 #ifdef DOUBLEPRECISION
  													Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  													#endif
  	                   			 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  	                		}
  	                		else if(i==NCELLX-1){
  												#ifdef DOUBLECOMPLEX
  	                			 Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]*cexp(-I*k_x*period_x)-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  												 #endif
  												 #ifdef DOUBLEPRECISION
  												 Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(0,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  												 #endif
  	                    		 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  	                		}
  	                		else{
  												#ifdef DOUBLECOMPLEX
  	                			 Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]*cexp(-I*k_y*period_y)-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  												 #endif
  												 #ifdef DOUBLEPRECISION
  												 Curl_E=(ex[ThreeDMapD(i,0,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  												 #endif
  	                    		 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  	                		}
  	                	}
  	                	else{
  							 Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  	                    	 hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  	               	}
                  }
  	    //             }
  	    //     }
  	    // }
      }

      else{
  		//	//#pragma omp parallel for collapse(3) private(Curl_E,i,j,k,j2,i2) // schedule(static)
  	    // for(i=0;i<NCELLX-1;i++){
  	    //     for(j=0;j<NCELLY-1;j++){
  	    //             for(k=0;k<NCELLZ;k++){
                      if(i<(NCELLX-1) && j<(NCELLY-1) && k<NCELLZ){
                        if(Hydrodynamics>=1)  hzPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] =  hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  	                    Curl_E=(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdy[j]-(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/khdx[i];
  	                    hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Chzh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_E;
  											//X-PML
  											if(i<cpml_N_X && j<cpml_y_lim && k<cpml_z_lim){
  												psi_Hz_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=bh_x_N[i]*psi_Hz_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ch_x_N[i]*(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
  												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
  											}
  											if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){
  												i2 = i - cpml_F_X;
  												psi_Hz_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=bh_x_F[i2]*psi_Hz_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ch_x_F[i2]*(ey[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dx;
  												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
  											}
  											//Y-PML
  											if(j<cpml_N_Y && i<cpml_x_lim && k<cpml_z_lim){
  												psi_Hz_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]=bh_y_N[j]*psi_Hz_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)]+ch_y_N[j]*(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
  												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY)];
  											}
  											if(j>=cpml_F_Y && i<cpml_x_lim && k<cpml_z_lim){
  												j2 = j - cpml_F_Y;
  												psi_Hz_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]=bh_y_F[j2]*psi_Hz_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)]+ch_y_F[j2]*(ex[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)])/dy;
  												hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Chze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Hz_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY)];
  											}
                      }
  	    //             }
  	    //     }
  	    // }
      }
    return;
  }



    __global__ void UPDATE_ex(real *ex,real *ex_n,real *ex_n_1,real *hy,real *hz,real *Cexe,real *Cexh,real *kedy,real *kedz,int *mat_matrix,int *mat_matrixX,int first_medium_max,real *psi_Ex_z_N,
    real *psi_Ex_z_F,real *psi_Ex_y_N,real *psi_Ex_y_F,real *Px_cp,real *Px_cp_n,real *Px_cp_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Px_d_n_2,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Py_d_n_2,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *Pz_d_n_2,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
    real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,
    real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_Y,int cpml_F_Y,int cpml_N_Z,int cpml_F_Z,int NcpmlY,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY,
  real *NDx,real *NDy,real *NDz,real *NDx_prev,real *NDy_prev,real *NDz_prev,real e0,real N_EQ){
        int i,j,k,n,k2,j2;
        comp Curl_H, Div_Grad=0.0,J_T,dummy_var;
    		comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,NxHold;
    			    comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
    					real INV_DX = 1.0/dx;
    					real INV_DY = 1.0/dy;
    					real INV_DZ = 1.0/dz;

              int idx = blockDim.x * blockIdx.x + threadIdx.x;

              i = idx / (NCELLZ*NCELLY);
              j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
              k = idx - i*NCELLZ*NCELLY - j*NCELLZ;

        if(Periodic_XY){
    			////#pragma omp parallel for collapse(3) private(Curl_H,i,j,j2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
    			// for(k=1;k<NCELLZ-1;k++){
    			// 	for(i=0;i<NCELLX;i++){
    		  //       for(j=0;j<NCELLY;j++){
                  if(i<NCELLX && j<NCELLY && k>0 && k<(NCELLZ-1)){
    	                //for(k=1;k<NCELLZ-1;k++){
    	                	if(j==0){
    											#ifdef DOUBLECOMPLEX
    												Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)]*cexp(I*k_y*period_y))/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
    											#endif
    											#ifdef DOUBLEPRECISION
    												Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)])/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
    											#endif

    	                	}
    	                	else{
    	                		Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];
    	                	}

    										if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){


    											//  Div_Grad = Calc_DIV_GRADx(i,j,k);
                          Div_Grad = 0.0;

                          C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

                          for(n=0;n<N_drude_poles;n++){
                              C_P_1+=(d_1_d[n]-1)*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                              C_P_3+=(d_2_d[n])*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                             C_P_NL += d_NL[n]*Div_Grad;
                          }
                          for(n=0;n<N_CP_poles;n++){
                              C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                              C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                          }
                          ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                          ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                          ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4 - C_P_NL);


    												//printf("%e\n",Div_Grad);
    												//Z-CPML
    												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    													//Near-Z-PML
    													psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
    												}
    												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    													k2 = k - cpml_F_Z ;
    													psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=(1/C_E)*dt*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    												}

    											 	for(n=0;n<N_CP_poles;n++){
    											 							Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    											 	}

    												for(n=0;n<N_drude_poles;n++){
    													Px_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;

    												}
                            for(n=0;n<N_drude_poles;n++){
                              Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d_n[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                              Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d[FourDMap(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];


                            }
                            for(n=0;n<N_CP_poles;n++){
                              Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                              Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];

                            }
    										}

    										else{
    												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Cexe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
    												//Z-CPML
    												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    													//Near-Z-PML
    													psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												}
    												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    													k2 = k - cpml_F_Z ;
    													psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
    												}
    										}
                      }



    	  //       }
    	  //   }
        // }
    	}
        //No PBCs
        else{
    			////#pragma omp target teams distribute parallel for collapse(3) schedule(static,1) private(Curl_H,i,j,j2,k,k2,n,dummy_var,J_T,Div_Grad)
    		//	//#pragma omp parallel for collapse(3) private(Curl_H,i,j,j2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
    	    // for(i=0;i<NCELLX-1;i++){
    	    //     for(j=1;j<NCELLY-1;j++){
    	    //             for(k=1;k<NCELLZ-1;k++){
                        if(i<(NCELLX-1) && j>0 && j<(NCELLY-1) && k>0 && k<(NCELLZ-1)){
    	                    Curl_H=(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j]-(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k];

    											if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){


    													if(Hydrodynamics == 0)
    													{
    														//Div_Grad = Calc_DIV_GRADx(i,j,k);
                                Div_Grad = 0.0;
    														// CP_D_ex(i,j,k,Curl_H,Div_Grad);

                                 C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

                                 for(n=0;n<N_drude_poles;n++){
                                     C_P_1+=(d_1_d[n]-1)*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                     C_P_3+=(d_2_d[n])*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                    C_P_NL += d_NL[n]*Div_Grad;
                                 }
                                 for(n=0;n<N_CP_poles;n++){
                                     C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                     C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                 }
                                 ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                 ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                 ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4 - C_P_NL);




    														for(n=0;n<N_CP_poles;n++){
    																				Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    														}

    														for(n=0;n<N_drude_poles;n++){
    																				Px_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;

    														}

                                for(n=0;n<N_drude_poles;n++){
                                  Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                  Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                }
                                for(n=0;n<N_CP_poles;n++){
                                  Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                  Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                }
    													}
    													else{



    														Vx1 = Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)];
    														Vx2 = Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)];
    														Nx1 = NDx_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)];
    														Nx2 = NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)];

    														Vy1 = 0.5*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Vy2 = 0.5*(Py_d_n[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Ny1 = 0.5*(NDy_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    														Ny2 = 0.5*(NDy_prev[ThreeDMapD(i+1,j-1,k,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);

    														Vz1 = 0.5*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Vz2 = 0.5*(Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Nz1 = 0.5*(NDz_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    														Nz2 = 0.5*(NDz_prev[ThreeDMapD(i+1,j,k-1,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);

    														NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] - 2.0*dt*INV_DX*(0.5*(Nx1*Vx1-Nx2*Vx2) + (Ny1*Vy1-Ny2*Vy2) + (Nz1*Vz1-Nz2*Vz2) + (0.5*(Vx1-Vx2) + (Vy1-Vy2) + (Vz1-Vz2))*N_EQ);



    												    C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;


    												    for(n=0;n<N_CP_poles;n++){
    												        C_P_2+=(C_1_cp[n]-1)*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
    												        C_P_4+=(C_2_cp[n])*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
    												    }
    												    ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												    ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												    ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0 + dt*Px_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*e0*(N_EQ + NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]) + C_E_1*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] -C_E_2*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_2-C_P_4);

    														for(n=0;n<N_CP_poles;n++){
    																				Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    														}
                                for(n=0;n<N_drude_poles;n++){
                                  Px_d_n_2[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                  Px_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                  Px_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Px_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                }
                                for(n=0;n<N_CP_poles;n++){
                                  Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                  Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                }

                                NxHold = NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                NDx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NxHold;

    													}

    												//	printf("%e\n",Div_Grad);
    													//Z-CPML
    													// if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    													// 	//Near-Z-PML
    													// 	psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													// 	ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    													// }
    													// if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    													// 	k2 = k - cpml_F_Z ;
    													// 	psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													// 	ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=(1/C_E)*dt*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    													// }



    											}

    	                    else{
    	                        ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Cexe[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
    	                    }
    											//Z-CPML
    											if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    												//Near-Z-PML
    												psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
    											}
    											if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){ //Far Z PML
    												k2 = k - cpml_F_Z;
    												psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    												// if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
    												// 	ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=(1/C_E)*dt*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    												// }
    												// else{
    													ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
    												//}
    											}
    											//Y PML
    											if(j<cpml_N_Y+1 && i<cpml_x_lim && k<cpml_z_lim){ //Near Y PML
    												psi_Ex_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]=be_y_N[j]*psi_Ex_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]+ce_y_N[j]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
    												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)];
    											}
    											if(j>=cpml_F_Y && i<cpml_x_lim && k<cpml_z_lim){ //Far Y PML
    												j2 = j - cpml_F_Y;
    												psi_Ex_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]=be_y_F[j2]*psi_Ex_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]+ce_y_F[j2]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
    												ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cexh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ex_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)];
    											}
    									// 		if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
    									// 				for(n=0;n<N_CP_poles;n++){
    									// 										Px_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Px_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Px_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    									// 				}
    									// 				for(n=0;n<N_drude_poles;n++){
    									// 										Px_d[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=d_1_d[n]*Px_d_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+d_2_d[n]*Px_d_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+d_3_d[n]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ex_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    									// 				}
    	                // }
                    }
    	//         }
    	//     }
    	// }
    }

      return;
    }






    __global__ void UPDATE_ey(real *ey,real *ey_n,real *ey_n_1,real *hx,real *hz,real *Ceye,real *Ceyh,real *kedx,real *kedz,int *mat_matrix,int *mat_matrixY,int first_medium_max,real *psi_Ey_z_N,
    real *psi_Ey_z_F,real *psi_Ey_x_N,real *psi_Ey_x_F,real *Py_cp,real *Py_cp_n,real *Py_cp_n_1, real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Px_d_n_2,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Py_d_n_2,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *Pz_d_n_2,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
    real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_z_N,real *ce_z_F,real *be_z_N,real *be_z_F,real *ce_x_N,real *ce_x_F,real *be_x_N,real *be_x_F,
    real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_X,int cpml_F_X,int cpml_N_Z,int cpml_F_Z,int NcpmlX,int NcpmlZ,real C_E_1,real C_E_2,int Periodic_XY,
    real *NDx,real *NDy,real *NDz,real *NDx_prev,real *NDy_prev,real *NDz_prev,real e0,real N_EQ){
        int i,j,k,i2,k2,n;
        comp Curl_H,Div_Grad=0.0,dummy_var,J_T;
    		comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
    		real INV_DX = 1.0/dx;
    		real INV_DY = 1.0/dy;
    		real INV_DZ = 1.0/dz;
    		comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,NyHold;
        int idx = blockDim.x * blockIdx.x + threadIdx.x;

        i = idx / (NCELLZ*NCELLY);
        j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
        k = idx - i*NCELLZ*NCELLY - j*NCELLZ;

        if(Periodic_XY){
        	////#pragma omp parallel for collapse(3) private(Curl_H,i,j,i2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
    			// for(k=0;k<NCELLZ-1;k++){
    			//   for(i=0;i<NCELLX;i++){
    		  //       for(j=0;j<NCELLY;j++){
                  if(i<NCELLX && j<NCELLY && k>0 && k<(NCELLZ-1)){
    	               // for(k=1;k<NCELLZ-1;k++){
    	                	if(i==0){
    											#ifdef DOUBLECOMPLEX
    	                		Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]*cexp(I*period_x*k_x))/kedx[i];
    											#endif
    											#ifdef DOUBLEPRECISION
    											Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)])/kedx[i];
    											#endif
    											  //printf("%d,%d,%d \t %f\t%f\t%f\n",i,j,k,creal(Curl_H),creal(hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]), creal(hz[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]));

    	                	}

    	                	else{
    	                		Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i];
    	                	}
    										if(mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] <6){

    												//Div_Grad = Calc_DIV_GRADy(i,j,k);
                            Div_Grad = 0.0;
    											// printf("%e\n",d_NL[0]*Div_Grad);
    												// CP_D_ey(i,j,k,Curl_H,Div_Grad);
                            C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
                            //printf("here");
                            for(n=0;n<N_drude_poles;n++){
                                C_P_1+=(d_1_d[n]-1.0)*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                C_P_3+=(d_2_d[n])*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                C_P_NL += d_NL[n]*Div_Grad;
                            }
                            for(n=0;n<N_CP_poles;n++){
                                C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                            }
                            ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                            ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                            ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);


    											//	printf("%e\n",ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    												//Z-CPML
    												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    													//Here we are in the near Z-PML
    													psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
    												}
    												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    													//Here we are in the far Z-PML
    														k2 = k - cpml_F_Z ;
    														psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    															ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    												}

    														for(n=0;n<N_CP_poles;n++){
    																				Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    														}

    														//if(Hydrodynamics == 0){
    														for(n=0;n<N_drude_poles;n++){
    															Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;

    														}

                                for(n=0;n<N_drude_poles;n++){
                                  Py_d_n_2[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                  Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                  Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];


                                }
                                for(n=0;n<N_CP_poles;n++){

                                  Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                  Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];

                                }
    													//}

    										}

    										else{
    												ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
    												//Z-CPML
    												if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    													//Here we are in the near Z-PML
    													psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
    												}
    												if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    													//Here we are in the far Z-PML
    														k2 = k - cpml_F_Z ;
    														psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    														if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
    															ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    														}
    														else{
    															ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
    														}
    												}
    										}

                      }

    	      //       }
    	      //   }
       	    // }
        }

        else{
    		//	//#pragma omp parallel for collapse(3) private(Curl_H,i,j,i2,k,k2,n,dummy_var,J_T,Div_Grad) // schedule(static)
    	    // for(i=1;i<NCELLX-1;i++){
    	    //     for(j=0;j<NCELLY-1;j++){
    	    //             for(k=1;k<NCELLZ-1;k++){
                        if(i>0 && i<(NCELLX-1) && j<(NCELLY-1) && k>0 && k<(NCELLZ-1)){
    	                    Curl_H=(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/kedz[k]-(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i];
    											if(mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] <6){

    												if(Hydrodynamics == 0){
    													//Div_Grad = Calc_DIV_GRADy(i,j,k);
                              Div_Grad= 0.0;
    												// printf("%e\n",d_NL[0]*Div_Grad);
    													// CP_D_ey(i,j,k,Curl_H,Div_Grad);
                              C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;
                          		//printf("here");
                              for(n=0;n<N_drude_poles;n++){
                                  C_P_1+=(d_1_d[n]-1.0)*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                  C_P_3+=(d_2_d[n])*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                          				C_P_NL += d_NL[n]*Div_Grad;
                              }
                              for(n=0;n<N_CP_poles;n++){
                                  C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                  C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                              }
                              ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                              ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                              ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);

    												//	printf("%e\n",ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    													//Z-CPML
    													// if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    													// 	//Here we are in the near Z-PML
    													// 	psi_Ey_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													// 	ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    													// }
    													// if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    													// 	//Here we are in the far Z-PML
    													// 		k2 = k - cpml_F_Z ;
    													// 		psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    													// 			ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    													// }

    															for(n=0;n<N_CP_poles;n++){
    																					Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    															}

    															for(n=0;n<N_drude_poles;n++){
    																Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;

    															}
                                  for(n=0;n<N_drude_poles;n++){

                                    Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                    Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];


                                  }
                                  for(n=0;n<N_CP_poles;n++){
                                    // Py_cp_n_2[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];

                                    Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                    Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];

                                  }
    													}
    													else{

    														Vy1 = Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)];
    														Vy2 = Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)];
    														Ny1 = NDy_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)];
    														Ny2 = NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)];

    														Vx1 = 0.5*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Vx2 = 0.5*(Px_d_n[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Nx1 = 0.5*(NDx_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    														Nx2 = 0.5*(NDx_prev[ThreeDMapD(i-1,j+1,k,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);

    														Vz1 = 0.5*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Vz2 = 0.5*(Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
    														Nz1 = 0.5*(NDz_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    														Nz2 = 0.5*(NDz_prev[ThreeDMapD(i,j+1,k-1,NCELLZ,NCELLY)] + NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);

    														 NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] - 2.0*dt*INV_DX*((Nx1*Vx1-Nx2*Vx2) + 0.5*(Ny1*Vy1-Ny2*Vy2) + (Nz1*Vz1-Nz2*Vz2) + ((Vx1-Vx2) + 0.5*(Vy1-Vy2) + (Vz1-Vz2))*N_EQ);


    											    C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

    											    for(n=0;n<N_CP_poles;n++){
    											        C_P_2+=(C_1_cp[n]-1.0)*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
    											        C_P_4+=(C_2_cp[n])*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
    											    }
    											    ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    											    ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    											    ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0 + dt*Py_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*e0*(N_EQ + NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]) + C_E_1*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_2-C_P_4);
    													for(n=0;n<N_CP_poles;n++){
    																			Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ey_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    													}

                              for(n=0;n<N_CP_poles;n++){
                                Py_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                Py_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Py_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];

                              }
                              for(n=0;n<N_drude_poles;n++){
                                Py_d_n_2[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                Py_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                Py_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Py_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];


                              }


                                                              NyHold = NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                                              NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                                              NDy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NyHold;

    											}
    										}


    										else{
    												ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceye[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;

    										}

    										if(k<cpml_N_Z+1 && i<cpml_x_lim && j<cpml_y_lim){
    											//Here we are in the near Z-PML
    											psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]=be_z_N[k]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)]+ce_z_N[k]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    											ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_N[ThreeDMapD(i,j,k,NcpmlZ+1,NCELLY)];
    										}
    										if(k>=cpml_F_Z && i<cpml_x_lim && j<cpml_y_lim){
    												//Here we are in the far Z-PML
    												k2 = k - cpml_F_Z;
    												psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]=be_z_F[k2]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]+ce_z_F[k2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])/dz;
    												// if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
    												// 	ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=(1/C_E)*dt*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)]/z0;
    												// }
    												// else{
    													ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_z_F[ThreeDMapD(i,j,k2,NcpmlZ+1,NCELLY)];
    												//}
    										}
    										//X-CPML
    										if(i<cpml_N_X+1 && j<cpml_y_lim && k<cpml_z_lim){
    											//Here we are in the near-X-PML
    											psi_Ey_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_x_N[i]*psi_Ey_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_x_N[i]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
    											ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    										}
    										if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){
    											//Here we are in the far-X-PML
    											i2 = i - cpml_F_X;
    											psi_Ey_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=be_x_F[i2]*psi_Ey_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ce_x_F[i2]*(hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
    											ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Ceyh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ey_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
    										}


                      }

    								}
    	  //       }
        //
    	  //   }
        // }

      return;
    }


    //void UPDATE_ez(void){
      __global__ void UPDATE_ez(real *ez,real *ez_n,real *ez_n_1,real *hx,real *hy,real *Ceze,real *Cezh,real *kedx,real *kedy,int *mat_matrix,int *mat_matrixZ,int first_medium_max,real *psi_Ez_y_N,
      real *psi_Ez_y_F,real *psi_Ez_x_N,real *psi_Ez_x_F,real *Pz_cp,real *Pz_cp_n,real *Pz_cp_n_1,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Px_d_n_2,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Py_d_n_2,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *Pz_d_n_2,real *C_1_cp,real *C_2_cp,real *C_3_cp,real *C_4_cp,real *C_5_cp,real *d_1_d,
      real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,real C_E,real z0,int N_CP_poles,int N_drude_poles,real *ce_y_N,real *ce_y_F,real *be_y_N,real *be_y_F,real *ce_x_N,real *ce_x_F,real *be_x_N,real *be_x_F,
      real dx,real dy,real dz,real dt,int NCELLX,int NCELLY,int NCELLZ,int Hydrodynamics,int cpml_x_lim,int cpml_y_lim,int cpml_z_lim,int cpml_N_X,int cpml_F_X,int cpml_N_Y,int cpml_F_Y,int NcpmlX,int NcpmlY,real C_E_1,real C_E_2,int Periodic_XY,
    real *NDx,real *NDy,real *NDz,real *NDx_prev,real *NDy_prev,real *NDz_prev,real e0,real N_EQ){

        int i,j,k;//,i2,j2,n;
        // comp Curl_H,Div_Grad=0.0,dummy_var,J_T;
    		// comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
    		// comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2;

    		real INV_DX = 1.0/dx;
    		real INV_DY = 1.0/dy;
    		real INV_DZ = 1.0/dz;
        int idx = blockDim.x * blockIdx.x + threadIdx.x;
        comp Curl_H,Div_Grad=0.0,dummy_var,J_T;
        comp C_P_1,C_P_2,C_P_3,C_P_4,C_P_NL;
        comp Vx1,Vx2,Vy1,Vy2,Vz1,Vz2,Nx1,Nx2,Ny1,Ny2,Nz1,Nz2,NzHold;
        int i2,j2,k2,n;
        i = idx / (NCELLZ*NCELLY);
        j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
        k = idx - i*NCELLZ*NCELLY - j*NCELLZ;

        if(Periodic_XY){
    	//	//#pragma omp parallel for collapse(3) private(Curl_H,i,j,k,i2,j2,n,dummy_var,J_T,Div_Grad) // schedule(static)
    		// for(k=0;k<NCELLZ-1;k++){
    		// for(i=0;i<NCELLX;i++){
    		//         for(j=0;j<NCELLY;j++){
                  if(i<NCELLX,j<NCELLY,k<NCELLZ){
    		           //     for(k=0;k<NCELLZ-1;k++){
    											//printf("Thread %d, ready to work\n",omp_get_thread_num());

    		                	if(i==0 || j==0){
    		                		if(i==0 && j==0){
    													#ifdef DOUBLECOMPLEX
    				        						Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]*cexp(I*k_x*period_x))/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)]*cexp(I*k_y*period_y))/kedy[j];
    														#endif
    														#ifdef DOUBLEPRECISION
    														Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)])/kedy[j];
    														#endif


    			               	    }
    			               	    else if(i==0){
    													#ifdef DOUBLECOMPLEX
    				                 	Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)]*cexp(I*k_x*period_x))/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
    													#endif
    													#ifdef DOUBLEPRECISION
    													Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(NCELLX-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
    													#endif


    				                }
    				                else{
    													#ifdef DOUBLECOMPLEX
    				                	Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)]*cexp(I*k_y*period_y))/kedy[j];
    													#endif
    													#ifdef DOUBLEPRECISION
    													Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,NCELLY-1,k,NCELLZ,NCELLY)])/kedy[j];
    													#endif


    				                }
    		        			}

    		        			else{
    		        				Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];

    		        			}
    									if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]<6){

    											// Div_Grad  = Calc_DIV_GRADz(i,j,k);
                          Div_Grad = 0.0;

    											// CP_D_ez(i,j,k,Curl_H,Div_Grad);
                          C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

                          for(n=0;n<N_drude_poles;n++){
                              C_P_1+=(d_1_d[n]-1)*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                              C_P_3+=(d_2_d[n])*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                              C_P_NL += d_NL[n]*Div_Grad;
                          }
                          for(n=0;n<N_CP_poles;n++){
                              C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                              C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                          }
                          ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                          ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                          ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);


    												for(n=0;n<N_CP_poles;n++){
    																		Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												}

    												for(n=0;n<N_drude_poles;n++){
    													Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
    												}

                            for(n=0;n<N_drude_poles;n++){

                              Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                              Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];

                            }
                            for(n=0;n<N_CP_poles;n++){

                              Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                              Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                            }
    									}

    									else{
    											ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
    									}

                    }
    		      //     }
    		      //   }
    					// }
        }

        else{
    	//		//#pragma omp parallel for collapse(3) private(Curl_H,i,j,k,i2,j2,n,dummy_var,J_T,Div_Grad) // schedule(static)
    	    // for(i=1;i<NCELLX-1;i++){
    	    //     for(j=1;j<NCELLY-1;j++){
    	    //             for(k=0;k<NCELLZ-1;k++){
                        if(i>0 && i<(NCELLX-1) && j>0 && j<(NCELLY-1) && k<(NCELLZ-1)){
    	                    Curl_H=(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/kedx[i]-(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/kedy[j];
    											if(mat_matrixZ[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max && mat_matrixZ[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]<6){

    												if(Hydrodynamics == 0){
    			 									 //Div_Grad  = Calc_DIV_GRADz(i,j,k);
                             Div_Grad = 0.0;
    			 									// CP_D_ez(i,j,k,Curl_H,Div_Grad);
                            C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

                            // Vz1 = d_1_d[0];
                            // Vy2 = d_2_d[0];
                            // Vy1 = d_3_d[0];
                            // Vx1 = d_4_d[0];
                            // Vx2 = d_5_d[0];


                            for(n=0;n<N_drude_poles;n++){
                                C_P_1+=(d_1_d[n]-1)*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                C_P_3+=(d_2_d[n])*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                C_P_NL += d_NL[n]*Div_Grad;
                            }
                            for(n=0;n<N_CP_poles;n++){
                                C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                            }
                            ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                            ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                            ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0+C_E_1*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_1-C_P_2-C_P_3-C_P_4-C_P_NL);


    			 										 for(n=0;n<N_CP_poles;n++){
    			 																 Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    			 										 }
    			 										 for(n=0;n<N_drude_poles;n++){
    			 											 Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_NL[n]*Div_Grad;
    			 										 }

                               for(n=0;n<N_drude_poles;n++){

                                 Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                                 Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];

                               }
                               for(n=0;n<N_CP_poles;n++){

                                 Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                                 Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                               }

    											 }
    											 else{

    												 C_P_1=C_P_2=C_P_3=C_P_4=C_P_NL=0.0;

    												 Vz1 = Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)];
    												 Vz2 = Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)];
    												 Nz1 = NDz_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)];
    												 Nz2 = NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)];

    												 Vx1 = 0.5*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    												 Vx2 = 0.5*(Px_d_n[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    												 Nx1 = 0.5*(NDx_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    												 Nx2 = 0.5*(NDx_prev[ThreeDMapD(i-1,j,k+1,NCELLZ,NCELLY)] + NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);

    												 Vy1 = 0.5*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    												 Vy2 = 0.5*(Py_d_n[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
    												 Ny1 = 0.5*(NDy_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
    												 Ny2 = 0.5*(NDy_prev[ThreeDMapD(i,j-1,k+1,NCELLZ,NCELLY)] + NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);

    												  NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] - 2.0*dt*INV_DX*((Nx1*Vx1-Nx2*Vx2) + (Ny1*Vy1-Ny2*Vy2) + 0.5*(Vz1-Vz2) + ((Vx1-Vx2) + (Vy1-Vy2) + 0.5*(Vz1-Vz2))*N_EQ);


    												 for(n=0;n<N_CP_poles;n++){
    														 C_P_2+=(C_1_cp[n]-1)*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
    														 C_P_4+=(C_2_cp[n])*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
    												 }
    												 ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												 ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												 ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=(1/C_E)*(dt*Curl_H/z0 + dt*Pz_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*e0*(NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + N_EQ)+C_E_1*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_E_2*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-C_P_2-C_P_4);

    												 for(n=0;n<N_CP_poles;n++){
    																		 Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    												 }

                             for(n=0;n<N_drude_poles;n++){
                               Pz_d_n_2[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                               Pz_d_n_1[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];
                               Pz_d_n[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)] = Pz_d[FourDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)];

                             }
                             for(n=0;n<N_CP_poles;n++){

                               Pz_cp_n_1[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                               Pz_cp_n[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)] = Pz_cp[FourDMapD(i,j,k,n,N_CP_poles,NCELLZ,NCELLY)];
                             }

                                                             NzHold = NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                                             NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                                                             NDz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] = NzHold;

    			 							 }
                       }


    	                    else{
    	                        ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=Ceze[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*Curl_H;
    	                    }
    											//Y CPML
    											if(j<cpml_N_Y && i<cpml_x_lim && k<cpml_z_lim){ //Near Y PML
    												psi_Ez_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]=be_y_N[j]*psi_Ez_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)]+ce_y_N[j]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
    												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_y_N[ThreeDMapD(i,j,k,NCELLZ,NcpmlY+1)];
    											}
    											if(j>=cpml_F_Y && i<cpml_x_lim &&  k<cpml_z_lim){
    												j2 = j - cpml_F_Y;
    												psi_Ez_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]=be_y_F[j2]*psi_Ez_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)]+ce_y_F[j2]*(hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])/dy;
    												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_y_F[ThreeDMapD(i,j2,k,NCELLZ,NcpmlY+1)];
    											}
    											//X PML
    											if(i<cpml_N_X+1 && j<cpml_y_lim && k<cpml_z_lim){//Near X-PML
    												psi_Ez_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]=be_x_N[i]*psi_Ez_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ce_x_N[i]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
    												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_x_N[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    											}
    											if(i>=cpml_F_X && j<cpml_y_lim && k<cpml_z_lim){//far X-PML
    												i2 = i - cpml_F_X;
    												psi_Ez_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]=be_x_F[i2]*psi_Ez_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)]+ce_x_F[i2]*(hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])/dx;
    												ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+=Cezh[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]*psi_Ez_x_F[ThreeDMapD(i2,j,k,NCELLZ,NCELLY)];
    											}
    											// if(mat_matrix[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] > first_medium_max){
    											// 		 for(n=0;n<N_CP_poles;n++){
    											// 				Pz_cp[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=C_1_cp[n]*Pz_cp_n[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+C_2_cp[n]*Pz_cp_n_1[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+C_3_cp[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_4_cp[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+C_5_cp[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    											// 				}
    											// 		for(n=0;n<N_drude_poles;n++){
    											// 				Pz_d[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]=d_1_d[n]*Pz_d_n[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_2_d[n]*Pz_d_n_1[ThreeDMapD(i,j,k,n,N_drude_poles,NCELLZ,NCELLY)]+d_3_d[n]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_4_d[n]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+d_5_d[n]*ez_n_1[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
    											// 				}
    											// }

    	                }
          //           }
    	    //     }
          //
    	    // }
        }

      return;
    }



    __global__ void UpdateHydroPx(real *ex,real *ex_n,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Px_d_n_2,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Py_d_n_2,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *Pz_d_n_2,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
			real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me,real z0){
      int i,j,k;

      real Vx1,Vx2,Vx3,Vy1,Vy2,Vy3,Vz1,Vz2,Vz3,Hx1,Hz1,Hy1,Hx2,Hz2,Hy2,Ex1,Ey1,Ez1,VdotGrad,VdotGrad2,VdotGrad3,DivV,VcrossH,VcrossH2,Pressure,ND1,ND2,ND3,Grad_Div,Grad_Div2;
      real INV_DX,INV_DY,INV_DZ;
      INV_DX = 1.0/dx;
      INV_DY = 1.0/dy;
      INV_DZ = 1.0/dz;
      Grad_Div = 0.0;
      Grad_Div2 =0.0;

      ////#pragma omp parallel for collapse(3)  // schedule(static)
      // for(i=0;i<NCELLX-1;i++){
      //     for(j=1;j<NCELLY-1;j++){
      //             for(k=1;k<NCELLZ-1;k++){

      int idx = blockDim.x * blockIdx.x + threadIdx.x;

      i = idx / (NCELLZ*NCELLY);
      j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
      k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
                    if(i<NCELLX-1 && j>0 && j<NCELLY-1 && k>0 && k<NCELLZ-1){
                    if(mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]> first_medium && mat_matrixX[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){
                      ND1 = N_EQ + NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];
                      //
                      // Grad_Div = Calc_DIV_GRADx(i,j,k);
                      // Grad_Div2 = Calc_DIV_GRADx2(i,j,k);


                      Grad_Div2 = INV_DX*INV_DX*(Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])
                                 + INV_DX*INV_DY*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i+1,j-1,k,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])
                                 + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i+1,j,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i+1,j,k-1,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);

                                 if(i==0 && j==0){
                                     Grad_Div = INV_DX*INV_DX*(Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DY*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i+1,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                                 }

                                 else if(i==NCELLX-1 && j==0){
                                     Grad_Div = INV_DX*INV_DX*(Px_d_n[FourDMapD(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DY*(Py_d_n[FourDMapD(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(0,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(0,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                                 }
                                 else if(i==NCELLX-1){
                                     Grad_Div = INV_DX*INV_DX*(Px_d_n[FourDMapD(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DY*(Py_d_n[FourDMapD(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(0,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(0,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(0,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                                 }
                                 else if(i==0){
                                     Grad_Div = INV_DX*INV_DX*(Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DY*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                                 }
                                 else if(j==0){
                                     Grad_Div = INV_DX*INV_DX*(Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DY*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i+1,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                                 }
                                 else{
                                     Grad_Div = INV_DX*INV_DX*(Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DY*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                                + INV_DX*INV_DZ*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                                 }
                      // Grad_Div = 0.0;
                      // Grad_Div2 = 0.0;

                      // Vx1 = d_1_d[0];
                      // Vx2 = d_2_d[0];
                      // Vy1 = d_3_d[0];
                      // Vy2 = d_4_d[0];
                      // Vz1 = d_5_d[0];
                      //




                      Vx1 = 0.5 * (Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vx2 = 0.5 * (Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vx2 = 0.5 * (Px_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);


                      Vy1 = 0.25 * (Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vy2 = 0.25 * (Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vy3 = 0.25 * (Py_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i+1,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);

                      //Vy1 = 0.5 * (Vy1 + Vy2);
                      Vz1 = 0.25 * (Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vz2 = 0.25 * (Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vz3 = 0.25 * (Pz_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i+1,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);

                      //Vz1 = 0.5 * (Vz1 + Vz2);
                      Hy2 = 0.5 * (hyPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hyPrev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);
                      Hz2 = 0.5 * (hzPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hzPrev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);
                      Hy1 = 0.5 * (hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hy[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);
                      Hz1 = 0.5 * (hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hz[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);

                      if(WithConvection==1) {
                        VdotGrad = 0.5*(Vx1*(Px_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                        VdotGrad2 = 0.5*(Vx2*(Px_d_n_1[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Px_d_n_1[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Px_d_n_1[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                        VdotGrad3= 0.5*(Vx3*(Px_d_n_2[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_2[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy3*(Px_d_n_2[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_2[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz3*(Px_d_n_2[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_2[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);

                        // VdotGrad = 1.0*(Vx1*(Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                        // VdotGrad2 = 1.0*(Vx2*(Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Px_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);

                        // VdotGrad = 0.5*(3.0*VdotGrad - 4.0*VdotGrad2 + VdotGrad3)/dt;
                        VdotGrad = (VdotGrad-VdotGrad2)/dt;
                      }
                      else VdotGrad = 0.0;

                      if(WithMagField==1){
                        VcrossH = Vy1*Hz1 - Vz1*Hy1;
                        VcrossH2 = Vy2*Hz2 - Vz2*Hy2;
                        VcrossH = (VcrossH - VcrossH2)/dt;
                      }
                      else VcrossH = 0.0;

                         Px_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] = d_1_d[0]*Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_2_d[0]*Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] +d_3_d[0]*ex[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]
                         + d_4_d[0]*ex_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_5_d[0]*(VdotGrad + (mu0*e0/me)*VcrossH/z0) + d_NL[0]*(Grad_Div + Grad_Div2/N_EQ)/powf(ND1,1.0/3.0);


                    }

             }
              //   }
              // }
    }

    __global__ void UpdateHydroPy(real *ey,real *ey_n,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Px_d_n_2,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Py_d_n_2,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *Pz_d_n_2,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
			real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me,real z0){
      int i,j,k;
      real Vx1,Vx2,Vx3,Vy1,Vy2,Vy3,Vz1,Vz2,Vz3,Hx1,Hz1,Hy1,Hx2,Hz2,Hy2,Ex1,Ey1,Ez1,VdotGrad,VdotGrad2,VdotGrad3,DivV,VcrossH,VcrossH2,Pressure,ND1,ND2,ND3,Grad_Div,Grad_Div2;
      real INV_DX,INV_DY,INV_DZ;
      INV_DX = 1.0/dx;
      INV_DY = 1.0/dy;
      INV_DZ = 1.0/dz;
      Grad_Div = 0.0;
      Grad_Div2 =0.0;
      ////#pragma omp parallel for collapse(3) // schedule(static)
      // for(i=1;i<NCELLX-1;i++){
      //     for(j=0;j<NCELLY-1;j++){
      //             for(k=1;k<NCELLZ-1;k++){

  int idx = blockDim.x * blockIdx.x + threadIdx.x;

  i = idx / (NCELLZ*NCELLY);
  j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
  k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
                if(j<NCELLY-1 && i>0 && i<NCELLX-1 && k>0 && k<NCELLZ-1){
                    if(mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]> first_medium && mat_matrixY[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){

                      // Grad_Div = Calc_DIV_GRADy(i,j,k);
                      // Grad_Div2 = Calc_DIV_GRADy2(i,j,k);


                      Grad_Div2 = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])
                                 + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i-1,j+1,k,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)])
                                 + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j+1,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j+1,k-1,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);

                      if(i==0 && j==0){
                          Grad_Div = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(NCELLX-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                        }

                        else if(i==0 && j==NCELLY-1){
                          Grad_Div = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(NCELLX-1,0,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,0,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                          }
                          else if(j==NCELLY-1){
                            Grad_Div = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,0,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,0,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,0,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                          }
                        else if(j==0){
                          Grad_Div = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                        }
                        else if(i==0){
                          Grad_Div = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(NCELLX-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                        }
                        else{
                          Grad_Div = INV_DY*INV_DY*(Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DX*INV_DY*(Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                     + INV_DY*INV_DZ*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]+Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                        }

                      //
                      // Grad_Div = 0.0;
                      // Grad_Div2 = 0.0;

                      ND1 = N_EQ + NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];

                      Vy1 = 0.5 * (Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vy2 = 0.5 * (Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vy3 = 0.5 * (Py_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);


                      Vx1 = 0.25 * (Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vx2 = 0.25 * (Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vx3 = 0.25 * (Px_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i-1,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)]);

                     //	Vx1 = 0.5 * (Vx1 + Vx2);
                      Vz1 = 0.25 * (Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vz2 = 0.25 * (Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);
                      Vz3 = 0.25 * (Pz_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i,j+1,k-1,0,N_drude_poles,NCELLZ,NCELLY)]);

                      //Vz1 = 0.5 * (Vz1 + Vz2);
                      Hx1 = 0.5 * (hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hx[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);
                      Hz1 = 0.5 * (hz[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hz[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);
                      Hx2 = 0.5 * (hxPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hxPrev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)]);
                      Hz2 = 0.5 * (hzPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hzPrev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);


                       if(WithMagField==1){
                         VcrossH = Vz1*Hx1 - Vx1*Hz1;
                         VcrossH2 = Vz2*Hx2 - Vx2*Hz2;
                         VcrossH = (VcrossH - VcrossH2)/dt;
                       }
                       else VcrossH = 0.0;
                       if(WithConvection==1) {
                         VdotGrad = 0.5*(Vx1*(Py_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Py_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                         VdotGrad2 = 0.5*(Vx2*(Py_d_n_1[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Py_d_n_1[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Py_d_n_1[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                         VdotGrad3 = 0.5*(Vx3*(Py_d_n_2[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_2[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy3*(Py_d_n_2[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_2[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz3*(Py_d_n_2[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_2[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);

                         // VdotGrad = 1.0*(Vx1*(Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                         // VdotGrad2 = 1.0*(Vx2*(Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Py_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);

                         // VdotGrad = (VdotGrad - VdotGrad2)/dt;
                         // VdotGrad = 0.5*(3.0*VdotGrad - 4.0*VdotGrad2 + VdotGrad3)/dt;
                         VdotGrad = (VdotGrad-VdotGrad2)/dt;


                       }
                       else VdotGrad = 0.0;

                     Py_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] = d_1_d[0]*Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_2_d[0]*Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]
                     + d_3_d[0]*ey[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_4_d[0]*ey_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ d_5_d[0]*(VdotGrad + (mu0*e0/me)*VcrossH/z0) + d_NL[0]*(Grad_Div + Grad_Div2/N_EQ)/powf(ND1,1.0/3.0);

                    }


                //   }
                // }
              }
    }

    __global__ void UpdateHydroPz(real *ez,real *ez_n,real *Px_d,real *Px_d_n,real *Px_d_n_1,real *Px_d_n_2,real *Py_d,real *Py_d_n,real *Py_d_n_1,real *Py_d_n_2,real *Pz_d,real *Pz_d_n,real *Pz_d_n_1,real *Pz_d_n_2,real *NDx,real *NDx_prev,real *NDy,real *NDy_prev,real *NDz,real *NDz_prev,
      real *hx,real *hy,real *hz,real *hxPrev,real *hyPrev,real *hzPrev,int WithConvection,int WithMagField,real N_EQ,int *mat_matrixX,int *mat_matrixY,int *mat_matrixZ,real dt,real dx,real dy,real dz,int NCELLX,int NCELLY,int NCELLZ,int first_medium,real *d_1_d,real *d_2_d,real *d_3_d,real *d_4_d,real *d_5_d,real *d_NL,int N_drude_poles,real mu0,real e0, real me,real z0){
      int i,j,k;
      real Vx1,Vx2,Vx3,Vy1,Vy2,Vy3,Vz1,Vz2,Vz3,Hx1,Hz1,Hy1,Hx2,Hz2,Hy2,Ex1,Ey1,Ez1,VdotGrad,VdotGrad2,VdotGrad3,DivV,VcrossH,VcrossH2,Pressure,ND1,ND2,ND3,Grad_Div,Grad_Div2;
      real INV_DX,INV_DY,INV_DZ;
      INV_DX = 1.0/dx;
      INV_DY = 1.0/dy;
      INV_DZ = 1.0/dz;
      Grad_Div = 0.0;
      Grad_Div2 =0.0;
      ////#pragma omp parallel for collapse(3) // schedule(static)
      // for(i=1;i<NCELLX-1;i++){
      //     for(j=1;j<NCELLY-1;j++){
      //             for(k=0;k<NCELLZ-1;k++){

                      int idx = blockDim.x * blockIdx.x + threadIdx.x;

                      i = idx / (NCELLZ*NCELLY);
                      j = (idx - i*NCELLZ*NCELLY) / NCELLZ;
                      k = idx - i*NCELLZ*NCELLY - j*NCELLZ;
                    if(k<NCELLZ-1 && i>0 && i<NCELLX-1 && j>0 && j<NCELLY-1){

                    if(mat_matrixZ[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]> first_medium && mat_matrixZ[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] < 6){


                        Grad_Div2 = INV_DZ*INV_DZ*(Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)]*NDz_prev[ThreeDMapD(i,j,k-1,NCELLZ,NCELLY)])
                                   + INV_DZ*INV_DY*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j-1,k+1,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDy_prev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)])
                                   + INV_DX*INV_DZ*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i,j,k+1,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i-1,j,k+1,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]*NDx_prev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)]);

                          if(i==0 && j==0){
                            Grad_Div = INV_DZ*INV_DZ*(Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DZ*INV_DY*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,NCELLY-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DX*INV_DZ*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(NCELLX-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                          }
                          else if(i==0){
                            Grad_Div = INV_DZ*INV_DZ*(Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DZ*INV_DY*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DX*INV_DZ*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(NCELLX-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(NCELLX-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                          }
                          else if(j==0){
                            Grad_Div = INV_DZ*INV_DZ*(Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DZ*INV_DY*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,NCELLY-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,NCELLY-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DX*INV_DZ*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                          }
                          else{
                            Grad_Div = INV_DZ*INV_DZ*(Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - 2.0*Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DZ*INV_DY*(Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Py_d_n[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])
                                       + INV_DX*INV_DZ*(Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]-Px_d_n[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)]+Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);

                          }


                      // Grad_Div = 0.0;
                      // Grad_Div2 =0.0;
                      ND1 = N_EQ + NDz_prev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)];

                     Vz1 = 0.5 * (Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                     Vz2 = 0.5 * (Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                     Vz3 = 0.5 * (Pz_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Pz_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);


                     Vx1 = 0.25 * (Px_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                     Vx2 = 0.25 * (Px_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                     Vx3 = 0.25 * (Px_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i-1,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Px_d_n_2[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)]);

                     //Vx1 = 0.5 * (Vx1 + Vx2);
                     Vy1 = 0.25 * (Py_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                     Vy2 = 0.25 * (Py_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);
                     Vy3 = 0.25 * (Py_d_n_2[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i,j-1,k+1,0,N_drude_poles,NCELLZ,NCELLY)] + Py_d_n_2[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)]);

                     //	Vy1 = 0.5 * (Vy1 + Vy2);
                     Hy1 = 0.5 * (hy[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)] + hy[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
                     Hx1 = 0.5 * (hx[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hx[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);
                     Hy2 = 0.5 * (hyPrev[ThreeDMapD(i-1,j,k,NCELLZ,NCELLY)] + hyPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]);
                     Hx2 = 0.5 * (hxPrev[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + hxPrev[ThreeDMapD(i,j-1,k,NCELLZ,NCELLY)]);

                     if(WithMagField==1){
                      VcrossH = Vx1*Hy1 - Vy1*Hx1;
                      VcrossH2 = Vx2*Hy2 - Vy2*Hx2;
                      VcrossH = (VcrossH - VcrossH2)/dt;
                     }
                     else VcrossH = 0.0;

                     if(WithConvection==1){
                      VdotGrad = 0.5*(Vx1*(Pz_d_n[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Pz_d_n[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Pz_d_n[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                      VdotGrad2 = 0.5*(Vx2*(Pz_d_n_1[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Pz_d_n_1[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Pz_d_n_1[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                      VdotGrad3 = 0.5*(Vx3*(Pz_d_n_2[FourDMapD(i+1,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_2[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy3*(Pz_d_n_2[FourDMapD(i,j+1,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_2[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz3*(Pz_d_n_2[FourDMapD(i,j,k+1,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_2[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);

                      // VdotGrad = 1.0*(Vx1*(Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy1*(Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz1*(Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);
                      // VdotGrad2 = 1.0*(Vx2*(Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMapD(i-1,j,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DX + Vy2*(Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMapD(i,j-1,k,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DY  + Vz2*(Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] - Pz_d_n_1[FourDMapD(i,j,k-1,0,N_drude_poles,NCELLZ,NCELLY)])*INV_DZ);

                      // VdotGrad = (VdotGrad - VdotGrad2)/dt;
                      // VdotGrad = 0.5*(3.0*VdotGrad - 4.0*VdotGrad2 + VdotGrad3)/dt;
                      VdotGrad = (VdotGrad-VdotGrad2)/dt;


                     }
                     else VdotGrad = 0.0;

                     Pz_d[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] = d_1_d[0]*Pz_d_n[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)] + d_2_d[0]*Pz_d_n_1[FourDMapD(i,j,k,0,N_drude_poles,NCELLZ,NCELLY)]
                     + d_3_d[0]*ez[ThreeDMapD(i,j,k,NCELLZ,NCELLY)] + d_4_d[0]*ez_n[ThreeDMapD(i,j,k,NCELLZ,NCELLY)]+ d_5_d[0]*(VdotGrad + (mu0*e0/me)*VcrossH/z0) + d_NL[0]*(Grad_Div + Grad_Div2/N_EQ)/powf(ND1,1.0/3.0);
                  //
                  //   }
                  //
                  //
                  // }
                }
              }
    }



  __global__  void UPDATE_e_inc(real* e_inc,real* h_inc,int inc_Length,real d_1D,real c0,real dt,real z0,real ep0,int t,real delay,real width,real pi,real f_0,int m0){
        real pulse;
        real factor = c0*dt;
        real e1,e2;
        int m;

        //1st Order Mur ABC buffers


        int i = blockDim.x * blockIdx.x + threadIdx.x;
        if(i<inc_Length){
        e1=e_inc[1];

        e2=e_inc[inc_Length-2];
      }
        //#pragma omp parallel for private(i)// schedule(guided)
    //    for(i=1;i<inc_Length-1;i++){
    __syncthreads();
    if(i>0 && i<inc_Length-1){
            e_inc[i]=e_inc[i]-(dt/(z0*ep0*d_1D))*(h_inc[i]-h_inc[i-1]);
      }
    //    }
      //  printf("%f\t%f\n",(dt/(z0*ep0*d_1D))*(h_inc[i]-h_inc[i-1]),((c0*dt-d_1D)/(c0*dt+d_1D)));
        //1st order Mur ABC
      //  __syncthreads();
      __syncthreads();
      //
      // if(i==inc_Length-1){
      //   e_inc[0]=e1+((factor-d_1D)/(factor+d_1D))*(e_inc[1]-e_inc[0]);
      //
      //   e_inc[inc_Length-1]=e2+((factor-d_1D)/(factor+d_1D))*(e_inc[inc_Length-2]-e_inc[inc_Length-1]);
      //
      //
      //   }
      //   if(i==m0-50){
      //     //introduce source
      //     #ifdef DOUBLECOMPLEX
      //     pulse=exp(-pow((real)(t-delay)/(real)width,2)/2.0)*cexp(I*2*pi*f_0*(t)*dt);
      //     #endif
      //     #ifndef DOUBLECOMPLEX
      //     pulse=exp(-powf((real)(t-delay)/(real)width,2)/2.0)*sin(2*pi*f_0*(t)*dt);
      //     #endif
      //   //  pulse =0.0;
      //     e_inc[m0-50]+=pulse;
      //   }
      //  printf("%f\n",pulse);
    }

  __global__  void UPDATE_h_inc(real* e_inc,real* h_inc,int inc_Length,real d_1D,real c0,real dt,real z0,real ep0,int t,real delay,real width,real pi,real f_0,real mu0){
        // int i;

    //    //#pragma omp parallel for private(i)
      //  for(i=0;i<inc_Length-1;i++){
      int i = blockDim.x * blockIdx.x + threadIdx.x;
      if(i<inc_Length-1){
            h_inc[i]=h_inc[i]-(z0*dt/(mu0*d_1D))*(e_inc[i+1]-e_inc[i]);
        }
    }

comp *Incident_spec;

int main(void)
{

  real  StaticBuild = 1000.0;
  real Transmit;
  int CALC_REFL = 0;
  // if(argc<=1){
  //   printf("Please Input max and min trials \n");
  //   exit(1);
  // }
  //
  //  int trials, max_trials, min_trials, num_trials;
   clock_t begin, End;
   real time_spent;
   int i;
// printf("Max threads:%d\n", omp_get_max_threads() );
//   //  printf("argc=%d\n", argc);
//   //  printf("min_trials=%d \t max_trials=%d\n",(int)atoi(argv[1]),(int)atoi(argv[1]));
//   //  min_trials = atoi(argv[1]);
//   //  max_trials = atoi(argv[2]);
//
//   //omp_set_num_threads(100);
//    #pragma omp parallel
//    {
//         if(omp_get_thread_num() == 0){
//            printf("Number of parallel cores: %d\n", omp_get_num_threads());
//         }
// //         printf("Thread %d, ready to work\n",omp_get_thread_num());
// //
// // if(omp_in_parallel()) printf("In Parallel Construct\n");
// 	//else printf("NOT\n");
//     }
//
//     printf("NUMBER OF DEVICES = %d\n",omp_get_num_devices());
    begin=clock();

    READ_DATA_FILE();

    if(StaticField == 1){
      printf("Static Field\n");
    }

    #ifdef DOUBLEPRECISION
    printf("Double Precision, Real Number Simulation\n");
    PBC_CTW = 0;
    if(num_trials != 1) {
      printf("Error, Only One Trial Allowed with Double Precision\n");
    exit(-1) ;
  }
    #endif
    #ifdef DOUBLECOMPLEX
    printf("Complex Double Precision, Complex Number Simulation\n");
    PBC_CTW = 1;
    if(num_trials == 1) printf("Double Check number of trials. Only 1 selected.\n");
    #endif


    Periodic_XY=0;
    Periodic_XZ=0;
    Periodic_YZ=0;

    if(Periodic_XY == 1){
      CALC_REFL = 1;
    }

    if(TE_TM == 0){
      TEz=1;
      TMz=0;
      printf("TE Polarization");
    }
    else if(TE_TM == 1){
      TEz=0;
      TMz=1;
      printf("TM Polarization");

    }
    else{
      printf("Invalid Polarization!!\n");
      exit(-1);
    }



    real polar_theta=0;
    // num_trials=15;
    polar_psi=polar_theta*(3.14159265359)/180;
    //Infinite dispersive slab (1-yes 0-no)
    //inf_disp_slab=1;
    printf("Setup Constants\n");
    SETUP_CONST();
    printf("Setup Material Matrix\n");
    MATERIAL_MATRIX();
    printf("Entering Trials Iterations\n");

    for(trials=min_trials;trials<=max_trials;trials++){

    k_rho=trials*3.1416e7/num_trials;
    k_x=k_rho*cos(polar_psi);
    k_y=k_rho*sin(polar_psi);
    printf("kx=%f, ky=%f\n", k_x,k_y);

    freq=MALLOC1D_double(freq,NUM_freq);
    real *freqdev;
    cudaMalloc(&freqdev,NUM_freq*sizeof(real));
real lam_max,lam_min;
lam_max = 700e-9;
lam_min = 200e-9;

  if(WL_or_freq == 1){ //Wavelength plot
    for(i=0;i<NUM_freq;i++){
      freq[i] = c0/(lam_min + i*(lam_max- lam_min)/(NUM_freq-1.0));
      }
      // if(NONLOCAL >= 3){
      //   // freq[0] = f_0;
      //   freq[1] = f_0*2.0;
      //   if(NUM_freq >= 3) freq[2] = f_0*3.0;
      // }
  }
  else{ //frequency plot
    for(i=0;i<NUM_freq;i++){
    freq[i] = c0/700e-9 + i*(c0/200e-9 - c0/700e-9)/(NUM_freq - 1.0);
  }
  }
cudaMemcpy(freqdev,freq,sizeof(real)*NUM_freq,cudaMemcpyHostToDevice);

    SETUP_Drude_CP();

    SETUP_SNAPSHOT();

    ALLOCATE_MEM();

    SETUP_TFSF();
    printf("fmin = %e, WL_max = %e\n",f_min, c0/f_min);


    //make the field vectors all zeros
    ex=ZERO_VECTORS3D_Complex(ex,NCELLX,NCELLY,NCELLZ);
    ey=ZERO_VECTORS3D_Complex(ey,NCELLX,NCELLY,NCELLZ);
    ez=ZERO_VECTORS3D_Complex(ez,NCELLX,NCELLY,NCELLZ);
    hx=ZERO_VECTORS3D_Complex(hx,NCELLX,NCELLY,NCELLZ);
    hy=ZERO_VECTORS3D_Complex(hy,NCELLX,NCELLY,NCELLZ);
    hz=ZERO_VECTORS3D_Complex(hz,NCELLX,NCELLY,NCELLZ);
cudaMemcpy(exdev,ex,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(eydev,ey,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(ezdev,ez,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(hxdev,hx,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(hydev,hy,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
cudaMemcpy(hzdev,hz,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);


printf("E and H fields zeroed\n");


    if(Hydrodynamics == 1){
    hxPrev=ZERO_VECTORS3D_Complex(hxPrev,NCELLX,NCELLY,NCELLZ);
    hyPrev=ZERO_VECTORS3D_Complex(hyPrev,NCELLX,NCELLY,NCELLZ);
    hzPrev=ZERO_VECTORS3D_Complex(hzPrev,NCELLX,NCELLY,NCELLZ);
    cudaMemcpy(hxPrevdev,hxPrev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(hyPrevdev,hyPrev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(hzPrevdev,hzPrev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    NDx=ZERO_VECTORS3D_Complex(NDx,NCELLX,NCELLY,NCELLZ);
    NDy=ZERO_VECTORS3D_Complex(NDy,NCELLX,NCELLY,NCELLZ);
    NDz=ZERO_VECTORS3D_Complex(NDz,NCELLX,NCELLY,NCELLZ);
    cudaMemcpy(NDxdev,NDx,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(NDydev,NDy,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(NDzdev,NDz,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    NDx_prev=ZERO_VECTORS3D_Complex(NDx_prev,NCELLX,NCELLY,NCELLZ);
    NDy_prev=ZERO_VECTORS3D_Complex(NDy_prev,NCELLX,NCELLY,NCELLZ);
    NDz_prev=ZERO_VECTORS3D_Complex(NDz_prev,NCELLX,NCELLY,NCELLZ);
    cudaMemcpy(NDx_prevdev,NDx_prev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(NDy_prevdev,NDy_prev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(NDz_prevdev,NDz_prev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
  }

    ex_n=ZERO_VECTORS3D_Complex(ex_n,NCELLX,NCELLY,NCELLZ);
    ey_n=ZERO_VECTORS3D_Complex(ey_n,NCELLX,NCELLY,NCELLZ);
    ez_n=ZERO_VECTORS3D_Complex(ez_n,NCELLX,NCELLY,NCELLZ);
    cudaMemcpy(ex_ndev,ex_n,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ey_ndev,ey_n,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ez_ndev,ez_n,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    ex_n_1=ZERO_VECTORS3D_Complex(ex_n_1,NCELLX,NCELLY,NCELLZ);
    ey_n_1=ZERO_VECTORS3D_Complex(ey_n_1,NCELLX,NCELLY,NCELLZ);
    ez_n_1=ZERO_VECTORS3D_Complex(ez_n_1,NCELLX,NCELLY,NCELLZ);
    cudaMemcpy(ex_n_1dev,ex_n_1,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ey_n_1dev,ey_n_1,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(ez_n_1dev,ez_n_1,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    printf("E and H fields zeroed\n");

    E_incident=ZERO_VECTORS1D_Complex2(E_incident,NUM_freq);
    E_reflected=ZERO_VECTORS1D_Complex2(E_reflected,NUM_freq);
    E_transmitted=ZERO_VECTORS1D_Complex2(E_transmitted,NUM_freq);


    Incident_spec=MALLOC1D_Complex(Incident_spec,NUM_freq);
    Incident_spec=ZERO_VECTORS1D_Complex(Incident_spec,NUM_freq);

    Ex_Reflected = ZERO_VECTORS3D_Complex2(Ex_Reflected, NUM_freq, NCELLX,NCELLY);
    Hx_Reflected = ZERO_VECTORS3D_Complex2(Hx_Reflected, NUM_freq, NCELLX,NCELLY);
    Ey_Reflected = ZERO_VECTORS3D_Complex2(Ey_Reflected, NUM_freq, NCELLX,NCELLY);
    Hy_Reflected = ZERO_VECTORS3D_Complex2(Hy_Reflected, NUM_freq, NCELLX,NCELLY);

    Ex_Transmitted = ZERO_VECTORS3D_Complex2(Ex_Transmitted, NUM_freq, NCELLX,NCELLY);
    Hx_Transmitted = ZERO_VECTORS3D_Complex2(Hx_Transmitted, NUM_freq, NCELLX,NCELLY);
    Ey_Transmitted = ZERO_VECTORS3D_Complex2(Ey_Transmitted, NUM_freq, NCELLX,NCELLY);
    Hy_Transmitted = ZERO_VECTORS3D_Complex2(Hy_Transmitted, NUM_freq, NCELLX,NCELLY);

    E_Incident = ZERO_VECTORS3D_Complex2(E_Incident, NUM_freq, NCELLX, NCELLY);
    H_Incident = ZERO_VECTORS3D_Complex2(H_Incident, NUM_freq, NCELLX, NCELLY);

    printf("Spectral Vectors Zeroed\n");
     e_inc=ZERO_VECTORS1D_Complex(e_inc,inc_Length);
     h_inc=ZERO_VECTORS1D_Complex(h_inc,inc_Length);
    //
    // ex_inc=ZERO_VECTORS1D_Complex(ex_inc,inc_Length);
    // ey_inc=ZERO_VECTORS1D_Complex(ey_inc,inc_Length);
    // ez_inc=ZERO_VECTORS1D_Complex(ez_inc,inc_Length);
    // hx_inc=ZERO_VECTORS1D_Complex(hx_inc,inc_Length);
    // hy_inc=ZERO_VECTORS1D_Complex(hy_inc,inc_Length);
    // hz_inc=ZERO_VECTORS1D_Complex(hz_inc,inc_Length);

    //Define permitivity, permeability, and conductivities (electric and magnetic)
    DEF_EPS();
    DEF_MU();
    DEF_SIGMA_E();
    DEF_SIGMA_M();
  //  DIELECTRIC_SLAB();

    //Define update coefficients
    printf("Setting Up update coefficients for E\n");
    Cexe=DEF_UPDATE_COEFF_EonE(Cexe);
    Cexh=DEF_UPDATE_COEFF_EonH(Cexh);
    Ceye=DEF_UPDATE_COEFF_EonE(Ceye);
    Ceyh=DEF_UPDATE_COEFF_EonH(Ceyh);
    Ceze=DEF_UPDATE_COEFF_EonE(Ceze);
    Cezh=DEF_UPDATE_COEFF_EonH(Cezh);
    cudaMemcpy(Cexedev,Cexe,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Cexhdev,Cexh,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Ceyedev,Ceye,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Ceyhdev,Ceyh,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Cezedev,Ceze,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Cezhdev,Cezh,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    printf("Setting Up update coefficients for H\n");

    Chxe=DEF_UPDATE_COEFF_HonE(Chxe);
//    printf("here\n");
    Chxh=DEF_UPDATE_COEFF_HonH(Chxh);
  //  printf("here1\n");
    Chye=DEF_UPDATE_COEFF_HonE(Chye);
  //  printf("here2\n");
    Chyh=DEF_UPDATE_COEFF_HonH(Chyh);
  //  printf("here3\n");
    Chze=DEF_UPDATE_COEFF_HonE(Chze);
  //  printf("here4\n");
    Chzh=DEF_UPDATE_COEFF_HonH(Chzh);

    cudaMemcpy(Chxedev,Chxe,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Chxhdev,Chxh,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Chyedev,Chye,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Chyhdev,Chyh,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Chzedev,Chze,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
    cudaMemcpy(Chzhdev,Chzh,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);



    printf("Setting Up CPML\n");
    SETUP_CPML();
  //  SETUP_t_inc();
  //  FT_source_calc();
trigger = 0;

// real FOURIERVEC_RE[NUM_freq*Tend];
// real FOURIERVEC_IE[NUM_freq*Tend];
// real FOURIERVEC_RH[NUM_freq*Tend];
// real FOURIERVEC_IH[NUM_freq*Tend];
// // real *FOURIERVEC_REdev,*FOURIERVEC_RHdev,*FOURIERVEC_IHdev,*FOURIERVEC_IEdev;
// // cudaMalloc(FOURIERVEC_REdev,sizeof()
// int mm;
//
// for(mm=0;mm<NUM_freq;mm++){
// for(t=0;t<Tend;t++){
//   FOURIERVEC_RE[TwoDMap(mm,t,Tend)] = cos(2.0*pi*(t+1)*dt*freq[mm]);
//   FOURIERVEC_IE[TwoDMap(mm,t,Tend)] = sin(2.0*pi*(t+1)*dt*freq[mm]);
//   FOURIERVEC_RH[TwoDMap(mm,t,Tend)] = cos(2.0*pi*(t+0.5)*dt*freq[mm]);
//   FOURIERVEC_IH[TwoDMap(mm,t,Tend)] = sin(2.0*pi*(t+0.5)*dt*freq[mm]);
// }
// }
comp pulse;
// if(StaticField == 0){
FILE *SourceFile;
SourceFile = fopen("SourceFile.txt","w");
 int m;
for(t=0;t<Tend;t++){
  #ifdef DOUBLECOMPLEX
  pulse=exp(-pow((real)(t-delay)/(real)width,2)/2.0)*cexp(I*2*pi*f_0*(t)*dt);
  #endif
  #ifndef DOUBLECOMPLEX
  pulse=MAX_AMP*exp(-powf((real)(t-delay)/(real)width,2)/2.0)*sin(2*pi*f_0*(t)*dt);
  #endif
fprintf(SourceFile,"%e\t%e\n",t*dt,pulse);
for(m=0;m<NUM_freq;m++){
  E_incident[m] += pulse*cexp(I*2*pi*t*dt*freq[m]);
}
}
fclose(SourceFile);
// }

cudaMalloc(&e_incdev,inc_Length*sizeof(comp));
cudaMalloc(&h_incdev,inc_Length*sizeof(comp));
cudaMemcpy(e_incdev,e_inc,inc_Length*sizeof(comp),cudaMemcpyHostToDevice);
cudaMemcpy(h_incdev,h_inc,inc_Length*sizeof(comp),cudaMemcpyHostToDevice);

//   printf("Malloc H sources\n");
//   Hx_source = MALLOC3D_Complex(Hx_source,NCELLX,NCELLY,Tend_inc);
//   Hy_source = MALLOC3D_Complex(Hy_source,NCELLX,NCELLY,Tend_inc);
// //  Hz_source = MALLOC3D_Complex(Hz_source,NCELLX,NCELLY,Tend_inc);
//   printf("Malloc E sources\n");
//   Ex_source = MALLOC3D_Complex(Ex_source,NCELLX,NCELLY,Tend_inc);
//   Ey_source = MALLOC3D_Complex(Ey_source,NCELLX,NCELLY,Tend_inc);
// //  Ez_source = MALLOC3D_Complex(Ez_source,NCELLX,NCELLY,Tend_inc);
//
//   printf("Zero H sources\n");
//   Hx_source = ZERO_VECTORS3D_Complex(Hx_source,NCELLX,NCELLY,Tend_inc);
//   Hy_source = ZERO_VECTORS3D_Complex(Hy_source,NCELLX,NCELLY,Tend_inc);
// //  Hz_source = ZERO_VECTORS3D_Complex(Hz_source,NCELLX,NCELLY,Tend_inc);
//
//   printf("Zero E Sources\n");
//   Ex_source = ZERO_VECTORS3D_Complex(Ex_source,NCELLX,NCELLY,Tend_inc);
//   Ey_source = ZERO_VECTORS3D_Complex(Ey_source,NCELLX,NCELLY,Tend_inc);
// //  Ez_source = ZERO_VECTORS3D_Complex(Ez_source,NCELLX,NCELLY,Tend_inc);
  //printf("Setting Up Source");


//if(PBC_CTW == 1) SOURCE_SETUP();
printf("Source Setup Finished \n");
Source=fopen("source.txt","w");
trigger = 1;
//    int m;
//    for(m=0;m<NcpmlZ+1;m++){
//
//        printf("%e\t%e\n",kedz[m],kedz[NCELLZ-m-1]);
//    }

FILE *Test;
Test = fopen("Test.txt","w");
Test_offset;


    char filename[100];
    FILE* Reflected;
    //Reflected=fopen(filename,"w");
    real Reflectivity;
    int i;
    real Amplitude;
real *EzInf,*JzInf,*EzInfdev,*JzInfdev;
cudaMalloc(&EzInfdev,sizeof(real));
cudaMalloc(&JzInfdev,sizeof(real));
EzInf = (real *)malloc(sizeof(real));
JzInf = (real *)malloc(sizeof(real));


    centerx = floor(NCELLX/2);
    centery = floor(NCELLY/2);
    int Number;
  	int threadsPerBlock = 350;
    Number = NCELLX * NCELLY *NCELLZ;
  	int blocksPerGrid = Number/threadsPerBlock + 1;
    dim3 blocksPerGrid2(blocksPerGrid,NUM_freq);
    int threadsPerBlock1D = threadsPerBlock;
    int blocksPerGrid1D = inc_Length/threadsPerBlock1D + 1;
    //Time Stepping
    real factor = c0*dt;

     cudaError_t err = cudaSuccess;
      cudaProfilerStart();
      err = cudaGetLastError();
      if( err != cudaSuccess)
      {
          printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
          exit(-1);
      }

    FILE *OUTFile;
    OUTFile = fopen("OUT.txt","w");
    for(t=0;t<=Tend;t++)
    {
        printf("%d\t%d\t%e\t%e\n",t,trials,ey[ThreeDMap((int)NCELLX/2,(int)NCELLY/2,(int)NCELLZ/2,NCELLZ,NCELLY)],pulse);

        if(isnan(ex[ThreeDMap((int)NCELLX/2,(int)NCELLY/2,(int)NCELLZ/2 - int(NCELLZ/4),NCELLZ,NCELLY)]) ) {
          fprintf(OUTFile, "NAN Detected\n");
          break;
        }
        fprintf(OUTFile,"%e\t%e\t%e\t%e\n",t*dt,ex[ThreeDMap((int)NCELLX/2,(int)NCELLY/2,(int)NCELLZ/2 ,NCELLZ,NCELLY)],ey[ThreeDMap((int)NCELLX/2,(int)NCELLY/2,(int)NCELLZ/2,NCELLZ,NCELLY)],ez[ThreeDMap((int)NCELLX/2,(int)NCELLY/2,(int)NCELLZ/2,NCELLZ,NCELLY)]);
        //UPDATE_B();
// printf("here\n");
      UPDATE_hx <<<blocksPerGrid, threadsPerBlock>>> (hxdev,hxPrevdev,ezdev,eydev,Chxhdev,Chxedev,psi_Hx_z_Ndev,psi_Hx_z_Fdev,psi_Hx_y_Ndev,psi_Hx_y_Fdev,khdydev,khdzdev,bh_z_Ndev,bh_z_Fdev,ch_z_Ndev,ch_z_Fdev,bh_y_Ndev,bh_y_Fdev,ch_y_Ndev,ch_y_Fdev,NCELLX,NCELLY,NCELLZ,Periodic_XY,dx,dy,dz,dt,cpml_N_Z,cpml_F_Z,cpml_N_Y,cpml_F_Y,cpml_z_lim,cpml_y_lim,cpml_x_lim,NcpmlZ,NcpmlY,Hydrodynamics);
        err = cudaGetLastError();
      if( cudaSuccess != err)
      {
          printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
          exit(-1);
      }
      // printf("here\n");

      // err = cudaDeviceSynchronize();
      // if( cudaSuccess != err)
      // {
      //     printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
      //     exit(-1);
      // }
      //

      UPDATE_hy <<<blocksPerGrid, threadsPerBlock>>> (hydev,hyPrevdev,ezdev,exdev,Chyhdev,Chyedev,psi_Hy_z_Ndev,psi_Hy_z_Fdev,psi_Hy_x_Ndev,psi_Hy_x_Fdev,khdxdev,khdzdev,bh_z_Ndev,bh_z_Fdev,ch_z_Ndev,ch_z_Fdev,bh_x_Ndev,bh_x_Fdev,ch_x_Ndev,ch_x_Fdev,NCELLX,NCELLY,NCELLZ,Periodic_XY,dx,dy,dz,dt,cpml_N_Z,cpml_F_Z,cpml_N_X,cpml_F_X,cpml_z_lim,cpml_y_lim,cpml_x_lim,NcpmlZ,NcpmlX,Hydrodynamics);
        err = cudaGetLastError();

      if( cudaSuccess != err)
      {
          printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
          exit(-1);
      }

      // err = cudaDeviceSynchronize();
      //
      // if( cudaSuccess != err)
      // {
      //     printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
      //     exit(-1);
      // }

       UPDATE_hz <<<blocksPerGrid, threadsPerBlock>>> (hzdev,hzPrevdev,eydev,exdev,Chzhdev,Chzedev,psi_Hz_x_Ndev,psi_Hz_x_Fdev,psi_Hz_y_Ndev,psi_Hz_y_Fdev,khdxdev,khdydev,bh_x_Ndev,bh_x_Fdev,ch_x_Ndev,ch_x_Fdev,bh_y_Ndev,bh_y_Fdev,ch_y_Ndev,ch_y_Fdev,NCELLX,NCELLY,NCELLZ,Periodic_XY,dx,dy,dz,dt,cpml_N_X,cpml_F_X,cpml_N_Y,cpml_F_Y,cpml_z_lim,cpml_y_lim,cpml_x_lim,NcpmlY,NcpmlX,Hydrodynamics);
        err = cudaGetLastError();
      if( cudaSuccess != err)
      {
          printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
          exit(-1);
      }

    // err = cudaDeviceSynchronize();
        if( cudaSuccess != err)
        {
            printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
            exit(-1);
        }

      //   cudaMemcpy(ex,exdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
      //   cudaMemcpy(ey,eydev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
      //   cudaMemcpy(ez,ezdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
      // cudaMemcpy(hx,hxdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
      // cudaMemcpy(hy,hydev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
      // cudaMemcpy(hz,hzdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
     // TFSF_CORRECT();
     UPDATE_h_inc<<<blocksPerGrid1D, threadsPerBlock1D>>>(e_incdev,h_incdev,inc_Length,d_1D,c0,dt,z0,ep0,t,delay,width,pi,f_0,mu0);
     err = cudaGetLastError();
   if( cudaSuccess != err)
   {
       printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
       exit(-1);
   }
     CORRECT_X<<<blocksPerGrid, threadsPerBlock>>>(NtfsfX,NtfsfY,NtfsfZ,NCELLX,NCELLY,NCELLZ,e_incdev,h_incdev,exdev,eydev,ezdev,hxdev,hydev,hzdev,inc_theta,inc_phi,polar_psi,polar_theta,dx,dy,dz,dt,i_0,j_0,k_0,d_1D,m0,Cexedev,Ceyedev,Cezedev,Cexhdev,Ceyhdev,Cezhdev,Chxedev,Chyedev,Chzedev,Chxhdev,Chyhdev,Chzhdev,Periodic_XY);
     err = cudaGetLastError();
   if( cudaSuccess != err)
   {
       printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
       exit(-1);
   }
     CORRECT_Y<<<blocksPerGrid, threadsPerBlock>>>(NtfsfX,NtfsfY,NtfsfZ,NCELLX,NCELLY,NCELLZ,e_incdev,h_incdev,exdev,eydev,ezdev,hxdev,hydev,hzdev,inc_theta,inc_phi,polar_psi,polar_theta,dx,dy,dz,dt,i_0,j_0,k_0,d_1D,m0,Cexedev,Ceyedev,Cezedev,Cexhdev,Ceyhdev,Cezhdev,Chxedev,Chyedev,Chzedev,Chxhdev,Chyhdev,Chzhdev,Periodic_XY);
     err = cudaGetLastError();
   if( cudaSuccess != err)
   {
       printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
       exit(-1);
   }

     CORRECT_Z<<<blocksPerGrid, threadsPerBlock>>>(NtfsfX,NtfsfY,NtfsfZ,NCELLX,NCELLY,NCELLZ,e_incdev,h_incdev,exdev,eydev,ezdev,hxdev,hydev,hzdev,inc_theta,inc_phi,polar_psi,polar_theta,dx,dy,dz,dt,i_0,j_0,k_0,d_1D,m0,Cexedev,Ceyedev,Cezedev,Cexhdev,Ceyhdev,Cezhdev,Chxedev,Chyedev,Chzedev,Chxhdev,Chyhdev,Chzhdev,Periodic_XY);
     err = cudaGetLastError();
   if( cudaSuccess != err)
   {
       printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
       exit(-1);
   }


     UPDATE_e_inc<<<blocksPerGrid1D, threadsPerBlock1D>>>(e_incdev,h_incdev,inc_Length,d_1D,c0,dt,z0,ep0,t,delay,width,pi,f_0,m0);
     //
     err = cudaGetLastError();
   if( cudaSuccess != err)
   {
       printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
       exit(-1);
   }

   cudaMemcpy(e_inc,e_incdev,inc_Length*sizeof(comp),cudaMemcpyDeviceToHost);

     e_inc[0]=e1+((factor-d_1D)/(factor+d_1D))*(e_inc[1]-e_inc[0]);

     e_inc[inc_Length-1]=e2+((factor-d_1D)/(factor+d_1D))*(e_inc[inc_Length-2]-e_inc[inc_Length-1]);

     e1=e_inc[1];

     e2=e_inc[inc_Length-2];


     if(StaticField == 0){
       //introduce source
       #ifdef DOUBLECOMPLEX
       pulse=exp(-pow((real)(t-delay)/(real)width,2)/2.0)*cexp(I*2*pi*f_0*(t)*dt);
       #endif
       #ifndef DOUBLECOMPLEX
       pulse=MAX_AMP*exp(-powf((real)(t-delay)/(real)width,2)/2.0)*sin(2*pi*f_0*(t)*dt);
       #endif
     //  pulse =0.0;
       e_inc[m0-50]+=pulse;
     }
     else{
          if(t<StaticBuild) pulse = (real) MAX_AMP * t/StaticBuild;
          else pulse = MAX_AMP * 1.0;
          e_inc[m0-50]+=pulse;

     }
       cudaMemcpy(e_incdev,e_inc,inc_Length*sizeof(comp),cudaMemcpyHostToDevice);

   // cudaDeviceSynchronize();

     // cudaMemcpy(exdev,ex,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
     // cudaMemcpy(eydev,ey,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
     // cudaMemcpy(ezdev,ez,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
     // cudaMemcpy(hxdev,hx,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
     // cudaMemcpy(hydev,hy,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
     // cudaMemcpy(hzdev,hz,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyHostToDevice);
     if(Hydrodynamics == 1){
     UpdateHydroPx<<<blocksPerGrid, threadsPerBlock>>>(exdev,ex_ndev,Px_ddev,Px_d_ndev,Px_d_n_1dev,Px_d_n_2dev,Py_ddev,Py_d_ndev,Py_d_n_1dev,Py_d_n_2dev,Pz_ddev,Pz_d_ndev,Pz_d_n_1dev,Pz_d_n_2dev,NDxdev,NDx_prevdev,NDydev,NDy_prevdev,NDzdev,NDz_prevdev,
       hxdev,hydev,hzdev,hxPrevdev,hyPrevdev,hzPrevdev,WithConvection,WithMagField,N_EQ,mat_matrixXdev,mat_matrixYdev,mat_matrixZdev,dt,dx,dy,dz,NCELLX,NCELLY,NCELLZ,first_medium,d_1_ddev,d_2_ddev,d_3_ddev,d_4_ddev,d_5_ddev,d_NLdev,N_drude_poles, mu0, e0,  me,z0);
       UpdateHydroPy<<<blocksPerGrid, threadsPerBlock>>>(eydev,ey_ndev,Px_ddev,Px_d_ndev,Px_d_n_1dev,Px_d_n_2dev,Py_ddev,Py_d_ndev,Py_d_n_1dev,Py_d_n_2dev,Pz_ddev,Pz_d_ndev,Pz_d_n_1dev,Pz_d_n_2dev,NDxdev,NDx_prevdev,NDydev,NDy_prevdev,NDzdev,NDz_prevdev,
         hxdev,hydev,hzdev,hxPrevdev,hyPrevdev,hzPrevdev,WithConvection,WithMagField,N_EQ,mat_matrixXdev,mat_matrixYdev,mat_matrixZdev,dt,dx,dy,dz,NCELLX,NCELLY,NCELLZ,first_medium,d_1_ddev,d_2_ddev,d_3_ddev,d_4_ddev,d_5_ddev,d_NLdev,N_drude_poles, mu0, e0,  me,z0);

         UpdateHydroPz<<<blocksPerGrid, threadsPerBlock>>>(ezdev,ez_ndev,Px_ddev,Px_d_ndev,Px_d_n_1dev,Px_d_n_2dev,Py_ddev,Py_d_ndev,Py_d_n_1dev,Py_d_n_2dev,Pz_ddev,Pz_d_ndev,Pz_d_n_1dev,Pz_d_n_2dev,NDxdev,NDx_prevdev,NDydev,NDy_prevdev,NDzdev,NDz_prevdev,
           hxdev,hydev,hzdev,hxPrevdev,hyPrevdev,hzPrevdev,WithConvection,WithMagField,N_EQ,mat_matrixXdev,mat_matrixYdev,mat_matrixZdev,dt,dx,dy,dz,NCELLX,NCELLY,NCELLZ,first_medium,d_1_ddev,d_2_ddev,d_3_ddev,d_4_ddev,d_5_ddev,d_NLdev,N_drude_poles, mu0, e0,  me,z0);
}


        UPDATE_ex <<<blocksPerGrid, threadsPerBlock>>> (exdev,ex_ndev,ex_n_1dev,hydev,hzdev,Cexedev,Cexhdev,kedydev,kedzdev,mat_matrixdev,mat_matrixXdev,first_medium_max,psi_Ex_z_Ndev,psi_Ex_z_Fdev,psi_Ex_y_Ndev,psi_Ex_y_Fdev,Px_cpdev,Px_cp_ndev,Px_cp_n_1dev,Px_ddev,Px_d_ndev,Px_d_n_1dev,Px_d_n_2dev,Py_ddev,Py_d_ndev,Py_d_n_1dev,Py_d_n_2dev,Pz_ddev,Pz_d_ndev,Pz_d_n_1dev,Pz_d_n_2dev,
       C_1_cpdev,C_2_cpdev,C_3_cpdev,C_4_cpdev,C_5_cpdev,d_1_ddev,d_2_ddev,d_3_ddev,d_4_ddev,d_5_ddev,d_NLdev,C_E,z0,N_CP_poles,N_drude_poles,ce_z_Ndev,ce_z_Fdev,be_z_Ndev,be_z_Fdev,ce_y_Ndev,ce_y_Fdev,be_y_Ndev,be_y_Fdev,dx,dy,dz,dt,NCELLX,NCELLY,NCELLZ,
       Hydrodynamics,cpml_x_lim,cpml_y_lim,cpml_z_lim,cpml_N_Y,cpml_F_Y,cpml_N_Z,cpml_F_Z,NcpmlY,NcpmlZ, C_E_1,C_E_2,Periodic_XY,NDxdev,NDydev,NDzdev,NDx_prevdev,NDy_prevdev,NDz_prevdev, e0, N_EQ);
       err = cudaGetLastError();
     if( cudaSuccess != err)
     {
         printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
         exit(-1);
     }
        //UPDATE_ey();
        UPDATE_ey <<<blocksPerGrid, threadsPerBlock>>> (eydev,ey_ndev,ey_n_1dev,hxdev,hzdev,Ceyedev,Ceyhdev,kedxdev,kedzdev,mat_matrixdev,mat_matrixYdev,first_medium_max,psi_Ey_z_Ndev,psi_Ey_z_Fdev,psi_Ey_x_Ndev,psi_Ey_x_Fdev,Py_cpdev,Py_cp_ndev,Py_cp_n_1dev,Px_ddev,Px_d_ndev,Px_d_n_1dev,Px_d_n_2dev,Py_ddev,Py_d_ndev,Py_d_n_1dev,Py_d_n_2dev,Pz_ddev,Pz_d_ndev,Pz_d_n_1dev,Pz_d_n_2dev,
       C_1_cpdev,C_2_cpdev,C_3_cpdev,C_4_cpdev,C_5_cpdev,d_1_ddev,d_2_ddev,d_3_ddev,d_4_ddev,d_5_ddev,d_NLdev,C_E,z0,N_CP_poles,N_drude_poles,ce_z_Ndev,ce_z_Fdev,be_z_Ndev,be_z_Fdev,ce_y_Ndev,ce_x_Fdev,be_x_Ndev,be_x_Fdev,dx,dy,dz,dt,NCELLX,NCELLY,NCELLZ,
       Hydrodynamics,cpml_x_lim,cpml_y_lim,cpml_z_lim,cpml_N_X,cpml_F_X,cpml_N_Z,cpml_F_Z,NcpmlX,NcpmlZ, C_E_1,C_E_2,Periodic_XY,NDxdev,NDydev,NDzdev,NDx_prevdev,NDy_prevdev,NDz_prevdev, e0, N_EQ);
       err = cudaGetLastError();
     if( cudaSuccess != err)
     {
         printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
         exit(-1);
     }

        // UPDATE_ez();
        UPDATE_ez <<<blocksPerGrid, threadsPerBlock>>> (ezdev,ez_ndev,ez_n_1dev,hxdev,hydev,Cezedev,Cezhdev,kedxdev,kedydev,mat_matrixdev,mat_matrixZdev,first_medium_max,psi_Ez_y_Ndev,psi_Ez_y_Fdev,psi_Ez_x_Ndev,psi_Ez_x_Fdev,Pz_cpdev,Pz_cp_ndev,Pz_cp_n_1dev,Px_ddev,Px_d_ndev,Px_d_n_1dev,Px_d_n_2dev,Py_ddev,Py_d_ndev,Py_d_n_1dev,Py_d_n_2dev,Pz_ddev,Pz_d_ndev,Pz_d_n_1dev,Pz_d_n_2dev,
       C_1_cpdev,C_2_cpdev,C_3_cpdev,C_4_cpdev,C_5_cpdev,d_1_ddev,d_2_ddev,d_3_ddev,d_4_ddev,d_5_ddev,d_NLdev,C_E,z0,N_CP_poles,N_drude_poles,ce_y_Ndev,ce_y_Fdev,be_y_Ndev,be_y_Fdev,ce_x_Ndev,ce_x_Fdev,be_x_Ndev,be_x_Fdev,dx,dy,dz,dt,NCELLX,NCELLY,NCELLZ,
       Hydrodynamics,cpml_x_lim,cpml_y_lim,cpml_z_lim,cpml_N_X,cpml_F_X,cpml_N_Y,cpml_F_Y,NcpmlX,NcpmlY, C_E_1,C_E_2,Periodic_XY,NDxdev,NDydev,NDzdev,NDx_prevdev,NDy_prevdev,NDz_prevdev, e0, N_EQ);

       err = cudaGetLastError();
     if( cudaSuccess != err)
     {
         printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
         exit(-1);
     }

//      cudaMemcpy(EzInf,EzInfdev,sizeof(real),cudaMemcpyDeviceToHost);
//      cudaMemcpy(JzInf,JzInfdev,sizeof(real),cudaMemcpyDeviceToHost);
// printf("%e\t%e\n",JzInf[0],EzInf[0]);

// cudaDeviceSynchronize();


  //  Fourier_Transform();
 if(StaticField == 0 && t % 4 == 0){
  ScattAbs <<<blocksPerGrid, threadsPerBlock>>>  (exdev,eydev,ezdev,hxdev,hydev,hzdev,NUM_freq,t,dt,freqdev,pi,XSTARTAbs,XENDAbs,YSTARTAbs,YENDAbs,ZSTARTAbs,ZENDAbs,XSTARTSca,XENDSca,YSTARTSca,YENDSca,ZSTARTSca,ZENDSca,XNEARAbs,XFARAbs,YNEARAbs,YFARAbs,ZNEARAbs,ZFARAbs,
   XNEARSca,XFARSca,YNEARSca,YFARSca,ZNEARSca,ZFARSca,ExTransformNearZAbsRedev,ExTransformNearZAbsImdev,EyTransformNearZAbsRedev,EyTransformNearZAbsImdev,HxTransformNearZAbsRedev,HxTransformNearZAbsImdev,HyTransformNearZAbsRedev,HyTransformNearZAbsImdev,
   ExTransformFarZAbsRedev,ExTransformFarZAbsImdev,EyTransformFarZAbsRedev,EyTransformFarZAbsImdev,HxTransformFarZAbsRedev,HxTransformFarZAbsImdev,HyTransformFarZAbsRedev,HyTransformFarZAbsImdev,
   ExTransformNearYAbsRedev,ExTransformNearYAbsImdev,EzTransformNearYAbsRedev,EzTransformNearYAbsImdev,HxTransformNearYAbsRedev,HxTransformNearYAbsImdev,HzTransformNearYAbsRedev,HzTransformNearYAbsImdev,
   ExTransformFarYAbsRedev,ExTransformFarYAbsImdev,EzTransformFarYAbsRedev,EzTransformFarYAbsImdev,HxTransformFarYAbsRedev,HxTransformFarYAbsImdev,HzTransformFarYAbsRedev,HzTransformFarYAbsImdev,
   EyTransformNearXAbsRedev,EyTransformNearXAbsImdev,EzTransformNearXAbsRedev,EzTransformNearXAbsImdev,HyTransformNearXAbsRedev,HyTransformNearXAbsImdev,HzTransformNearXAbsRedev,HzTransformNearXAbsImdev,
   EyTransformFarXAbsRedev,EyTransformFarXAbsImdev,EzTransformFarXAbsRedev,EzTransformFarXAbsImdev,HyTransformFarXAbsRedev,HyTransformFarXAbsImdev,HzTransformFarXAbsRedev,HzTransformFarXAbsImdev,
   ExTransformNearZScaRedev,ExTransformNearZScaImdev,EyTransformNearZScaRedev,EyTransformNearZScaImdev,HxTransformNearZScaRedev,HxTransformNearZScaImdev,HyTransformNearZScaRedev,HyTransformNearZScaImdev,
   ExTransformFarZScaRedev,ExTransformFarZScaImdev,EyTransformFarZScaRedev,EyTransformFarZScaImdev,HxTransformFarZScaRedev,HxTransformFarZScaImdev,HyTransformFarZScaRedev,HyTransformFarZScaImdev,
   ExTransformNearYScaRedev,ExTransformNearYScaImdev,EzTransformNearYScaRedev,EzTransformNearYScaImdev,HxTransformNearYScaRedev,HxTransformNearYScaImdev,HzTransformNearYScaRedev,HzTransformNearYScaImdev,
   ExTransformFarYScaRedev,ExTransformFarYScaImdev,EzTransformFarYScaRedev,EzTransformFarYScaImdev,HxTransformFarYScaRedev,HxTransformFarYScaImdev,HzTransformFarYScaRedev,HzTransformFarYScaImdev,
   EyTransformNearXScaRedev,EyTransformNearXScaImdev,EzTransformNearXScaRedev,EzTransformNearXScaImdev,HyTransformNearXScaRedev,HyTransformNearXScaImdev,HzTransformNearXScaRedev,HzTransformNearXScaImdev,
   EyTransformFarXScaRedev,EyTransformFarXScaImdev,EzTransformFarXScaRedev,EzTransformFarXScaImdev,HyTransformFarXScaRedev,HyTransformFarXScaImdev,HzTransformFarXScaRedev,HzTransformFarXScaImdev,NCELLX,NCELLY,NCELLZ);
}
// cudaDeviceSynchronize();
        //
        // #ifdef DOUBLECOMPLEX
        // fprintf(Test, "%e\t%e\t%e\t%e\t%e\t%e\n",creal(ex[centerx][centery][inc_plane - Test_offset]),creal(ey[centerx][centery][inc_plane - Test_offset]),creal(ez[centerx][centery][inc_plane - Test_offset]),creal(hx[centerx][centery][inc_plane - Test_offset]),creal(hy[centerx][centery][inc_plane - Test_offset]),creal(hz[centerx][centery][inc_plane - Test_offset] ));
        // #endif
        //
        // #ifdef DOUBLEPRECISION
        // fprintf(Test, "%e\t%e\t%e\t%e\t%e\t%e\n",ex[centerx][centery][25],ey[centerx][centery][25],ez[centerx][centery][25],hx[centerx][centery][25],hy[centerx][centery][25],hz[centerx][centery][25]);
        // #endif
        // cudaMemcpy(ex,exdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
        cudaMemcpy(ey,eydev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
        // cudaMemcpy(ez,ezdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);

      if(t%t_skip==0){
          //snapshot_count++;
//          //  SNAPSHOT_2D();
//            //aux field:
      if(Snap_in == 1) SNAPSHOT_1D();
      if(Snap_in == 2){
        cudaDeviceSynchronize();

          // cudaMemcpy(ex,NDxdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
          // cudaMemcpy(ey,NDydev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
          // cudaMemcpy(ez,NDzdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);

          cudaMemcpy(ex,Px_ddev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
          cudaMemcpy(ey,Py_ddev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
          cudaMemcpy(ez,Pz_ddev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);

          cudaMemcpy(NDx,NDxdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
          cudaMemcpy(NDy,NDydev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
          cudaMemcpy(NDz,NDzdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);

           SNAPSHOT_2D();
           SNAPSHOT_2D_N();

           cudaDeviceSynchronize();

      }
      if(Snap_in == 3){
        cudaDeviceSynchronize();
        cudaMemcpy(ex,exdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
        cudaMemcpy(ey,eydev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);
        cudaMemcpy(ez,ezdev,NCELLX*NCELLY*NCELLZ*sizeof(real),cudaMemcpyDeviceToHost);

        // cudaMemcpy(e_inc,e_incdev,inc_Length*sizeof(comp),cudaMemcpyDeviceToHost);
          err = cudaGetLastError();
        if( cudaSuccess != err)
        {
            printf( "Cuda error: %s.\n",cudaGetErrorString( err) );
            exit(-1);
        }
cudaDeviceSynchronize();
        SNAPSHOT_1D();
        SNAPSHOT_2D();
      }
    }


  if(CALC_REFL==1 && t>0 && t%(Tend/10) == 0 && StaticField == 0){
	Reflectance_XZ();
	 if(TMz) sprintf(filename,"ReflectivityTM.%d.%d.txt",trials,t);
   if(TEz) sprintf(filename,"ReflectivityTE.%d.%d.txt",trials,t);

    Reflected=fopen(filename,"w");
    double theta;
    printf("%f",f_min);
    for(i=1;i<NUM_freq;i++){

      if(PBC_CTW == 1){
      if(freq[i]>f_min){

       Amplitude = creal(E_incident[i]);
        if(TEz) Reflectivity = -1*creal(E_reflected[i]);
        if(TMz) Reflectivity = 1*creal(E_reflected[i]);
        theta = asin(c0*k_x/(2*pi*freq[i]));
        theta = 180*theta/pi;

              if(WL_or_freq == 1)  fprintf(Reflected,"%e\t%e\t%e\t%e\t%e\t%e\n",c0/freq[i],(Reflectivity/Amplitude),Amplitude,Reflectivity,(Reflectivity/(NCELLX*dx*NCELLZ*dz))/(cabs(Incident_spec[i])*cabs(Incident_spec[i])),cabs(Incident_spec[i])*cabs(Incident_spec[i]));
              else fprintf(Reflected,"%e\t%e\t%e\t%e\t%e\t%e\n",freq[i],(Reflectivity/Amplitude),Amplitude,Reflectivity,t_inc[i],theta);

        }
        else{
            Amplitude=0.0;
            Reflectivity=1.0;
            theta = 0.0;
            if(WL_or_freq == 1) fprintf(Reflected,"%e\t%f\t%f\t%f\t%e\t%e\n",c0/freq[i],Reflectivity,Amplitude,Reflectivity,t_inc[i],theta);
            else fprintf(Reflected,"%e\t%f\t%f\t%f\t%e\t%e\n",freq[i],Reflectivity,Amplitude,Reflectivity,t_inc[i],theta);
        }
    }

  else{
    Amplitude = sqrt(creal(E_incident[i])*creal(E_incident[i]) + cimag(E_incident[i])*cimag(E_incident[i]));
    Amplitude = Amplitude*Amplitude;
    if(TEz) Reflectivity = 1*creal(E_reflected[i]);
    else if(TMz) Reflectivity = 1*creal(E_reflected[i]);
    else Reflectivity = creal(E_reflected[i]);
    Transmit = creal(E_transmitted[i]);
    if(WL_or_freq == 1)  fprintf(Reflected,"%e\t%e\t%e\t%e\t%e\t%e\n",c0/freq[i],(Reflectivity/Amplitude),(Transmit/Amplitude),Amplitude,Reflectivity,Transmit);
    else fprintf(Reflected,"%e\t%e\t%e\t%e\t%e\t%e\n",freq[i],(Reflectivity/Amplitude),(Transmit/Amplitude),Amplitude,Reflectivity,Transmit);
  }
}
    fclose(Reflected);

    }


  if((Scattering == 1 || Absorption == 1) && t%(Tend/10) == 0 && StaticField == 0){
    cudaDeviceSynchronize();

    cudaMemcpy(ExTransformNearZScaRe,ExTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(EyTransformNearZScaRe,EyTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HxTransformNearZScaRe,HxTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HyTransformNearZScaRe,HyTransformNearZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);

    cudaMemcpy(ExTransformNearYScaRe,ExTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(EzTransformNearYScaRe,EzTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HxTransformNearYScaRe,HxTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HzTransformNearYScaRe,HzTransformNearYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);



  cudaMemcpy(EyTransformNearXScaRe,EyTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EzTransformNearXScaRe,EzTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HyTransformNearXScaRe,HyTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HzTransformNearXScaRe,HzTransformNearXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

  cudaMemcpy(ExTransformNearZScaIm,ExTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EyTransformNearZScaIm,EyTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HxTransformNearZScaIm,HxTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HyTransformNearZScaIm,HyTransformNearZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);

  cudaMemcpy(ExTransformNearYScaIm,ExTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EzTransformNearYScaIm,EzTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HxTransformNearYScaIm,HxTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HzTransformNearYScaIm,HzTransformNearYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

  cudaMemcpy(EyTransformNearXScaIm,EyTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EzTransformNearXScaIm,EzTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HyTransformNearXScaIm,HyTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HzTransformNearXScaIm,HzTransformNearXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

  cudaMemcpy(ExTransformNearZAbsRe,ExTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EyTransformNearZAbsRe,EyTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HxTransformNearZAbsRe,HxTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HyTransformNearZAbsRe,HyTransformNearZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);

  cudaMemcpy(ExTransformNearYAbsRe,ExTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EzTransformNearYAbsRe,EzTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HxTransformNearYAbsRe,HxTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HzTransformNearYAbsRe,HzTransformNearYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);

  cudaMemcpy(EyTransformNearXAbsRe,EyTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(EzTransformNearXAbsRe,EzTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HyTransformNearXAbsRe,HyTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
  cudaMemcpy(HzTransformNearXAbsRe,HzTransformNearXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);


    cudaMemcpy(ExTransformNearZAbsIm,ExTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(EyTransformNearZAbsIm,EyTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HxTransformNearZAbsIm,HxTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HyTransformNearZAbsIm,HyTransformNearZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);

    cudaMemcpy(ExTransformNearYAbsIm,ExTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(EzTransformNearYAbsIm,EzTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HxTransformNearYAbsIm,HxTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HzTransformNearYAbsIm,HzTransformNearYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);

    cudaMemcpy(EyTransformNearXAbsIm,EyTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(EzTransformNearXAbsIm,EzTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HyTransformNearXAbsIm,HyTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
    cudaMemcpy(HzTransformNearXAbsIm,HzTransformNearXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarZScaRe,ExTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EyTransformFarZScaRe,EyTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarZScaRe,HxTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarZScaRe,HyTransformFarZScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarYScaRe,ExTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarYScaRe,EzTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarYScaRe,HxTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarYScaRe,HzTransformFarYScaRedev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(EyTransformFarXScaRe,EyTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarXScaRe,EzTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarXScaRe,HyTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarXScaRe,HzTransformFarXScaRedev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarZScaIm,ExTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EyTransformFarZScaIm,EyTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarZScaIm,HxTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarZScaIm,HyTransformFarZScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(YENDSca-YSTARTSca+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarYScaIm,ExTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarYScaIm,EzTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarYScaIm,HxTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarYScaIm,HzTransformFarYScaImdev,sizeof(real2)*NUM_freq*(XENDSca-XSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(EyTransformFarXScaIm,EyTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarXScaIm,EzTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarXScaIm,HyTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarXScaIm,HzTransformFarXScaImdev,sizeof(real2)*NUM_freq*(YENDSca-YSTARTSca+1)*(ZENDSca-ZSTARTSca+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarZAbsRe,ExTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EyTransformFarZAbsRe,EyTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarZAbsRe,HxTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarZAbsRe,HyTransformFarZAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarYAbsRe,ExTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarYAbsRe,EzTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarYAbsRe,HxTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarYAbsRe,HzTransformFarYAbsRedev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);


      cudaMemcpy(EyTransformFarXAbsRe,EyTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarXAbsRe,EzTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarXAbsRe,HyTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarXAbsRe,HzTransformFarXAbsRedev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);


      cudaMemcpy(ExTransformFarZAbsIm,ExTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EyTransformFarZAbsIm,EyTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarZAbsIm,HxTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarZAbsIm,HyTransformFarZAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(YENDAbs-YSTARTAbs+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(ExTransformFarYAbsIm,ExTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarYAbsIm,EzTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HxTransformFarYAbsIm,HxTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarYAbsIm,HzTransformFarYAbsImdev,sizeof(real2)*NUM_freq*(XENDAbs-XSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);

      cudaMemcpy(EyTransformFarXAbsIm,EyTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(EzTransformFarXAbsIm,EzTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HyTransformFarXAbsIm,HyTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaMemcpy(HzTransformFarXAbsIm,HzTransformFarXAbsImdev,sizeof(real2)*NUM_freq*(YENDAbs-YSTARTAbs+1)*(ZENDAbs-ZSTARTAbs+1),cudaMemcpyDeviceToHost);
      cudaDeviceSynchronize();

    CalculateAbsScatt();
  }


}
cudaProfilerStop();
fclose(OUTFile);
fclose(Test);
    fclose(Source);
    printf("\n TOTAL SNAPSHOTS=%d",snapshot_count);




//    FREE_MEM();

    }
    End=clock();
    time_spent=(real)(End-begin)/CLOCKS_PER_SEC;
    printf("Time Spent= %e",time_spent);
    return 0;
}


void Reflectance_XZ()
{
    real  P_r,P_t,P_r_inc, Area, dsurf1, dsurf2,k_rho,k_z;
    real u,v,u2,v2;
    double complex  zz,zz_trans;
    int i, k, m;
    double complex cEx, cEy, ccHx, ccHy;
    double complex cEx_trans, cEy_trans, ccHx_trans, ccHy_trans;

    double complex cEx_inc, cEy_inc, ccHx_inc, cc_Hz_inc;

    FILE *Spectrum;
    Spectrum = fopen("Spectrum_ref.txt","w");

    for (m = 0; m < NUM_freq; m++) {
    		Area = 0.0;
    		P_r = 0.0;
        P_r_inc = 0.0;



      //  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

    		for (i = 0; i < NCELLX; i++) {
    			dsurf1 = dx;
    			for (k = 0; k < NCELLY; k++) {
    			   dsurf2 = dy;

      				//XZ-NEAR = Reflected
      				if(PBC_CTW == 1 || PBC_CTW == 0) Area += dsurf1*dsurf2;
            //  if(PBC_CTW == 0) Area += 1;


      				// u = ReTFEx_XZ[m][i][0][k];
      				// v = ImTFEx_XZ[m][i][0][k];
              //
      				// cEx = u + I_UNIT*v;
              if(PBC_CTW == 1){
                u = creal(Ey_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)]-E_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                v = cimag(Ey_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)]-E_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
              }
              else{
                u = creal(Ey_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                v = cimag(Ey_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                u2 = creal(Ex_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                v2 = cimag(Ex_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
              }


      				cEy = u + I*v;
              cEx = u2 + I*v2;

              u = creal(E_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
              v = cimag(E_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);

              cEy_inc = u + I*v;

              if(!PBC_CTW){
                u = creal(Ey_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                v = cimag(Ey_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                u2 = creal(Ex_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
                v2 = cimag(Ex_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
              }

              cEy_trans = u + I*v;
              cEx_trans = u2 + I*v2;

              //
      				// u = ReTFHx_XZ[m][i][0][k];// + ReTFHx_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
      				// v = ImTFHx_XZ[m][i][0][k];// + ImTFHx_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
              //
      				// ccHx = u - I_UNIT*v;
              if(PBC_CTW == 1){
                u = creal(Hx_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)] - H_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);// + ReTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
        				v = cimag(Hx_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)] - H_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);// + ImTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
              }
              else{
                u = creal(Hx_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ReTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
                v = cimag(Hx_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ImTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
                u2 = creal(Hy_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ReTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
                v2 = cimag(Hy_Reflected[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ImTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
              }


      				ccHx = u - I*v;
              ccHy = u2 - I*v2;
              u = creal(H_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);
              v = cimag(H_Incident[ThreeDMap(m,i,k,NCELLY,NCELLX)]);

              ccHx_inc = u - I*v;

              if(!PBC_CTW){
                u = creal(Hx_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ReTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
                v = cimag(Hx_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ImTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
                u2 = creal(Hy_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ReTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
                v2 = cimag(Hy_Transmitted[ThreeDMap(m,i,k,NCELLY,NCELLX)] );// + ImTFHz_xz_n[ThreeDMap(m,i,k,NCELLY,NCELLX)];
              }
              ccHx_trans = u - I*v;
              ccHy_trans = u2 - I*v2;

      				zz = cEy*ccHx - cEx*ccHy;
              zz_trans = cEy_trans*ccHx_trans - cEx_trans*ccHy_trans;


      				P_r += creal(zz)*dsurf1*dsurf2;
              P_t += creal(zz_trans)*dsurf1*dsurf2;



              zz = cEy_inc*ccHx_inc;

              P_r_inc += creal(zz)*dsurf1*dsurf2;

    				}
    		}

        //printf("%e\n", P_r/Area);
       E_reflected[m] = P_r/Area;
       if(!PBC_CTW) E_transmitted[m] = P_t/Area;
       if(PBC_CTW) E_incident[m] = P_r_inc/Area;

    }
    //Area_loc = Area;
  fclose(Spectrum);
	return;
}

// void SOURCE_SETUP(void){
//
//     int *** mat_matrix_safe;
//     printf("Creating Second MAT Matrix\n");
//
//     mat_matrix_safe = MALLOC3D_int(mat_matrix_safe,NCELLX,NCELLY,NCELLZ);
//
//     double complex ex_tot=0.0;
//     double complex ex_source=0.0;
//     double complex ey_tot=0.0;
//     double complex ey_source=0.0;
//     double complex ez_tot=0.0;
//     double complex ez_source=0.0;
//
//     double complex hx_tot=0.0;
//     double complex hx_source=0.0;
//     double complex hy_tot=0.0;
//     double complex hy_source=0.0;
//     double complex hz_tot=0.0;
//     double complex hz_source=0.0;
//
//
//     int n,w,i,j,k=inc_plane;
//
//     printf("Switching Material Matrices\n");
//
//     for(i=0;i<NCELLX;i++){
//       for(j=0;j<NCELLY;j++){
//         for(k=0;k<NCELLZ;k++){
//
//           mat_matrix_safe[i][j][k] = mat_matrix[i][j][k];
//
//           mat_matrix[i][j][k] = first_medium;
//         }
//       }
//     }
//
//     FILE *Test2;
//
//    Test2 = fopen("Test2.txt","w");
//
//    k = inc_plane;
//
//     for(t=0;t<=Tend_inc;t++)
//     {
//       printf("Source Setup\t%d\n",t);
//
//       UPDATE_B();
//     //  printf("B updated\n");
//       SOURCE_IN();
//
//
//       UPDATE_E();
//     //  printf("E updated\n");
//
//             for(i=0;i<NCELLX;i++){
//               for(j=0;j<NCELLY;j++){
//             //   //  printf("%d\n",Spect_loc);
//             //     Ex_source[i][j][t] = ex[i][j][inc_plane-Spect_loc];
//             //     Ey_source[i][j][t] = ey[i][j][inc_plane-Spect_loc];
//             // //    Ez_source[i][j][t] = ez[i][j][inc_plane-Spect_loc];
//             //
//             //     Hx_source[i][j][t] = hx[i][j][inc_plane-Spect_loc];
//             //     Hy_source[i][j][t] = hy[i][j][inc_plane-Spect_loc];
//             // //    Hz_source[i][j][t] = hz[i][j][inc_plane-Spect_loc];
//
//                 ex_tot = (ex[i][j][inc_plane-Spect_loc])*cexp(I*(i+0.5)*dx*k_x)*cexp(I*j*dy*k_y);
//                 ey_tot = (ey[i][j][inc_plane-Spect_loc])*cexp(I*i*dx*k_x)*cexp(I*(j+0.5)*dy*k_y);
//
//                 hx_tot = (hx[i][j][inc_plane-Spect_loc])*cexp(I*(i)*dx*k_x)*cexp(I*(j+0.5)*dy*k_y);
//                 hy_tot = (hy[i][j][inc_plane-Spect_loc])*cexp(I*(i+0.5)*dx*k_x)*cexp(I*(j)*dy*k_y);
//
//
//                 for(w=0;w<NUM_freq;w++){
//                   if(TEz){
//                     H_Incident[w][i][j] += hx_tot*cexp(-I*2*pi*t*dt*freq[w]);//*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//                     E_Incident[w][i][j] += ey_tot*cexp(-I*2*pi*t*dt*freq[w]);//*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//                   }
//                   if(TMz){
//                     H_Incident[w][i][j] += hy_tot*cexp(-I*2*pi*t*dt*freq[w]);//*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//                     E_Incident[w][i][j] += ex_tot*cexp(-I*2*pi*t*dt*freq[w]);//*cexp(-I*(i+0.5)*dx*k_x)*cexp(-I*(j)*dy*k_y);
//                   }
//                 // if(i == 0 && j ==0) Incident_spec[w] += PULSE(t)*cexp(-I*2*pi*t*dt*freq[w]);
//                 }
//
//
//               }
//             }
//       fprintf(Test2, "%e\t%e\t%e\t%e\t%e\t%e\n",creal(ex[centerx][centery][inc_plane - Test_offset]),creal(ey[centerx][centery][inc_plane - Test_offset]),creal(ez[centerx][centery][inc_plane - Test_offset]),creal(hx[centerx][centery][inc_plane - Test_offset]),creal(hy[centerx][centery][inc_plane - Test_offset]),creal(hz[centerx][centery][inc_plane - Test_offset]));
//
//       if(t%t_skip==0 && t>NUM_freq_inc ){
//             snapshot_count++;
// //          //  SNAPSHOT_2D();
// //            //aux field:
//       if(Snap_in == 1) SNAPSHOT_1D();
// //
// //
// //
//     }
//     }
//
//     FREE3D_Complex(ez,NCELLX,NCELLY);
//     FREE3D_Complex(ey,NCELLX,NCELLY);
//     FREE3D_Complex(ex,NCELLX,NCELLY);
//     FREE3D_Complex(hx,NCELLX,NCELLY);
//     FREE3D_Complex(hy,NCELLX,NCELLY);
//     FREE3D_Complex(hz,NCELLX,NCELLY);
//
//     FREE3D_Complex(psi_Ex_y_N,NCELLX,NcpmlY+1);
//     FREE3D_Complex(psi_Ex_z_N,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Ey_x_N,NcpmlX+1,NCELLY);
//     FREE3D_Complex(psi_Ey_z_N,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Ez_y_N,NCELLX,NcpmlY+1);
//     FREE3D_Complex(psi_Ez_x_N,NcpmlX+1,NCELLY);
//     FREE3D_Complex(psi_Hx_z_N,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Hx_y_N,NCELLX,NcpmlY);
//     FREE3D_Complex(psi_Hy_x_N,NcpmlX,NCELLY);
//     FREE3D_Complex(psi_Hy_z_N,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Hz_x_N,NcpmlX,NCELLY);
//     FREE3D_Complex(psi_Hz_y_N,NCELLX,NcpmlY);
//
//     FREE3D_Complex(psi_Ex_y_F,NCELLX,NcpmlY+1);
//     FREE3D_Complex(psi_Ex_z_F,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Ey_x_F,NcpmlX+1,NCELLY);
//     FREE3D_Complex(psi_Ey_z_F,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Ez_y_F,NCELLX,NcpmlY+1);
//     FREE3D_Complex(psi_Ez_x_F,NcpmlX+1,NCELLY);
//     FREE3D_Complex(psi_Hx_z_F,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Hx_y_F,NCELLX,NcpmlY);
//     FREE3D_Complex(psi_Hy_x_F,NcpmlX,NCELLY);
//     FREE3D_Complex(psi_Hy_z_F,NCELLX,NCELLY);
//     FREE3D_Complex(psi_Hz_x_F,NcpmlX,NCELLY);
//     FREE3D_Complex(psi_Hz_y_F,NCELLX,NcpmlY);
//
//     FREE4D_Complex(Pz_cp,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Py_cp,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Px_cp,NCELLX,NCELLY,NCELLZ);
//
//     FREE4D_Complex(Pz_cp_n,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Py_cp_n,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Px_cp_n,NCELLX,NCELLY,NCELLZ);
//
//     FREE4D_Complex(Pz_cp_n_1,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Py_cp_n_1,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Px_cp_n_1,NCELLX,NCELLY,NCELLZ);
//
//     FREE4D_Complex(Px_d,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Px_d_n,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Px_d_n_1,NCELLX,NCELLY,NCELLZ);
//
//     FREE4D_Complex(Py_d,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Py_d_n,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Py_d_n_1,NCELLX,NCELLY,NCELLZ);
//
//     FREE4D_Complex(Pz_d,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Pz_d_n,NCELLX,NCELLY,NCELLZ);
//     FREE4D_Complex(Pz_d_n_1,NCELLX,NCELLY,NCELLZ);
//
//     psi_Ex_y_N=MALLOC3D_Complex(psi_Ex_y_N,NCELLX,NcpmlY+1,NCELLZ);
//     psi_Ez_y_N=MALLOC3D_Complex(psi_Ez_y_N,NCELLX,NcpmlY+1,NCELLZ);
//     psi_Ex_y_F=MALLOC3D_Complex(psi_Ex_y_F,NCELLX,NcpmlY+1,NCELLZ);
//     psi_Ez_y_F=MALLOC3D_Complex(psi_Ez_y_F,NCELLX,NcpmlY+1,NCELLZ);
//
//     psi_Ex_z_N=MALLOC3D_Complex(psi_Ex_z_N,NCELLX,NCELLY,NcpmlZ+1);
//     psi_Ey_z_N=MALLOC3D_Complex(psi_Ey_z_N,NCELLX,NCELLY,NcpmlZ+1);
//     psi_Ex_z_F=MALLOC3D_Complex(psi_Ex_z_F,NCELLX,NCELLY,NcpmlZ+1);
//     psi_Ey_z_F=MALLOC3D_Complex(psi_Ey_z_F,NCELLX,NCELLY,NcpmlZ+1);
//
//     psi_Ey_x_N=MALLOC3D_Complex(psi_Ey_x_N,NcpmlX+1,NCELLY,NCELLZ);
//     psi_Ez_x_N=MALLOC3D_Complex(psi_Ez_x_N,NcpmlX+1,NCELLY,NCELLZ);
//     psi_Ey_x_F=MALLOC3D_Complex(psi_Ey_x_F,NcpmlX+1,NCELLY,NCELLZ);
//     psi_Ez_x_F=MALLOC3D_Complex(psi_Ez_x_F,NcpmlX+1,NCELLY,NCELLZ);
//
//     psi_Hx_y_F=MALLOC3D_Complex(psi_Hx_y_F,NCELLX,NcpmlY,NCELLZ);
//     psi_Hz_y_F=MALLOC3D_Complex(psi_Hz_y_F,NCELLX,NcpmlY,NCELLZ);
//     psi_Hx_y_N=MALLOC3D_Complex(psi_Hx_y_N,NCELLX,NcpmlY,NCELLZ);
//     psi_Hz_y_N=MALLOC3D_Complex(psi_Hz_y_N,NCELLX,NcpmlY,NCELLZ);
//
//     psi_Hx_z_F=MALLOC3D_Complex(psi_Hx_z_F,NCELLX,NCELLY,NcpmlZ);
//     psi_Hy_z_F=MALLOC3D_Complex(psi_Hy_z_F,NCELLX,NCELLY,NcpmlZ);
//     psi_Hx_z_N=MALLOC3D_Complex(psi_Hx_z_N,NCELLX,NCELLY,NcpmlZ);
//     psi_Hy_z_N=MALLOC3D_Complex(psi_Hy_z_N,NCELLX,NCELLY,NcpmlZ);
//
//     psi_Hz_x_F=MALLOC3D_Complex(psi_Hz_x_F,NcpmlX,NCELLY,NCELLZ);
//     psi_Hy_x_F=MALLOC3D_Complex(psi_Hy_x_F,NcpmlX,NCELLY,NCELLZ);
//     psi_Hy_x_N=MALLOC3D_Complex(psi_Hy_x_N,NcpmlX,NCELLY,NCELLZ);
//     psi_Hz_x_N=MALLOC3D_Complex(psi_Hz_x_N,NcpmlX,NCELLY,NCELLZ);
//
//
//     ex=MALLOC3D_Complex(ex,NCELLX,NCELLY,NCELLZ);
//     ey=MALLOC3D_Complex(ey,NCELLX,NCELLY,NCELLZ);
//     ez=MALLOC3D_Complex(ez,NCELLX,NCELLY,NCELLZ);
// //    Dx=MALLOC3D_Complex(Dx,NCELLX,NCELLY,NCELLZ);
// //    Dy=MALLOC3D_Complex(Dy,NCELLX,NCELLY,NCELLZ);
// //    Dz=MALLOC3D_Complex(Dz,NCELLX,NCELLY,NCELLZ);
//     hx=MALLOC3D_Complex(hx,NCELLX,NCELLY,NCELLZ);
//     hy=MALLOC3D_Complex(hy,NCELLX,NCELLY,NCELLZ);
//     hz=MALLOC3D_Complex(hz,NCELLX,NCELLY,NCELLZ);
//
//     Px_cp=MALLOC4D_Complex(Px_cp,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Px_cp_n=MALLOC4D_Complex(Px_cp_n,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Px_cp_n_1=MALLOC4D_Complex(Px_cp_n_1,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Py_cp=MALLOC4D_Complex(Py_cp,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Py_cp_n=MALLOC4D_Complex(Py_cp_n,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Py_cp_n_1=MALLOC4D_Complex(Py_cp_n_1,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Pz_cp=MALLOC4D_Complex(Pz_cp,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Pz_cp_n=MALLOC4D_Complex(Pz_cp_n,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//     Pz_cp_n_1=MALLOC4D_Complex(Pz_cp_n_1,NCELLX,NCELLY,NCELLZ,N_CP_poles);
//
//     Px_d=MALLOC4D_Complex(Px_d,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Px_d_n=MALLOC4D_Complex(Px_d_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Px_d_n_1=MALLOC4D_Complex(Px_d_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Py_d=MALLOC4D_Complex(Py_d,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Py_d_n=MALLOC4D_Complex(Py_d_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Py_d_n_1=MALLOC4D_Complex(Py_d_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Pz_d=MALLOC4D_Complex(Pz_d,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Pz_d_n=MALLOC4D_Complex(Pz_d_n,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//     Pz_d_n_1=MALLOC4D_Complex(Pz_d_n_1,NCELLX,NCELLY,NCELLZ,N_drude_poles);
//
//
//     ex=ZERO_VECTORS3D_Complex(ex,NCELLX,NCELLY,NCELLZ);
//     ey=ZERO_VECTORS3D_Complex(ey,NCELLX,NCELLY,NCELLZ);
//     ez=ZERO_VECTORS3D_Complex(ez,NCELLX,NCELLY,NCELLZ);
//     hx=ZERO_VECTORS3D_Complex(hx,NCELLX,NCELLY,NCELLZ);
//     hy=ZERO_VECTORS3D_Complex(hy,NCELLX,NCELLY,NCELLZ);
//     hz=ZERO_VECTORS3D_Complex(hz,NCELLX,NCELLY,NCELLZ);
//
//  //    for(i=0;i<NCELLX;i++){
//  //      for(j=0;j<NCELLY;j++){
//  //        for(k=0;k<NCELLZ;k++){
//  //          ex[i][j][k] = 0.0;
//  //          ey[i][j][k] = 0.0;
//  //          ez[i][j][k] = 0.0;
//  //          hx[i][j][k] = 0.0;
//  //          hy[i][j][k] = 0.0;
//  //          hz[i][j][k] = 0.0;
//  //
//  //        }
//  //      }
//  //    }
//  //
//  //    for(i=0; i<NCELLX;i++){
//  //        for(j=0;j<NCELLY;j++){
//  //            for(k=0;k<NcpmlZ+1;k++){
//  //
//  //                psi_Ex_z_F[i][j][k]=0.0;
//  //                psi_Ex_z_N[i][j][k]=0.0;
//  //                psi_Ey_z_F[i][j][k]=0.0;
//  //                psi_Ey_z_N[i][j][k]=0.0;
//  //
//  //
//  //        }
//  //      }
//  //    }
//  //
//  // for(i=0; i<NCELLX;i++){
//  //     for(j=0;j<NCELLY;j++){
//  //         for(k=0;k<NcpmlZ;k++){
//  //
//  //                psi_Hx_z_F[i][j][k]=0.0;
//  //                psi_Hx_z_N[i][j][k]=0.0;
//  //                psi_Hy_z_F[i][j][k]=0.0;
//  //                psi_Hy_z_N[i][j][k]=0.0;
//  //
//  //
//  //        }
//  //      }
//  //    }
//
//  for(i=0; i<NCELLX;i++){
//      for(j=0;j<NCELLY;j++){
//          for(k=0;k<NcpmlZ+1;k++){
//
//              psi_Ex_z_F[i][j][k]=0.0;
//              psi_Ex_z_N[i][j][k]=0.0;
//              psi_Ey_z_F[i][j][k]=0.0;
//              psi_Ey_z_N[i][j][k]=0.0;
//
//
//      }
//    }
//  }
//  printf("Here\n");
//
// for(i=0; i<NCELLX;i++){
//   for(j=0;j<NCELLY;j++){
//       for(k=0;k<NcpmlZ;k++){
//
//              psi_Hx_z_F[i][j][k]=0.0;
//              psi_Hx_z_N[i][j][k]=0.0;
//              psi_Hy_z_F[i][j][k]=0.0;
//              psi_Hy_z_N[i][j][k]=0.0;
//
//
//      }
//    }
//  }
//  printf("Here\n");
//
//  for(i=0;i<NCELLX;i++){
//      for(j=0;j<NCELLY;j++){
//          for(k=0;k<NCELLZ;k++){
//              for(n=0;n<N_CP_poles;n++){
//                  Px_cp[i][j][k][n]=0.0;
//                  Px_cp_n[i][j][k][n]=0.0;
//                  Px_cp_n_1[i][j][k][n]=0.0;
//                  Py_cp[i][j][k][n]=0.0;
//                  Py_cp_n[i][j][k][n]=0.0;
//                  Py_cp_n_1[i][j][k][n]=0.0;
//                  Pz_cp[i][j][k][n]=0.0;
//                  Pz_cp_n[i][j][k][n]=0.0;
//                  Pz_cp_n_1[i][j][k][n]=0.0;
//
//                  printf("%d,%d,%d,%d\n",i,j,k,n);
//
//          }
//        }
//      }
//  }
//  printf("Here\n");
//
//  for(i=0;i<NCELLX;i++){
//      for(j=0;j<NCELLY;j++){
//         for(k=0;k<NCELLZ;k++){
//           for(n=0;n<N_drude_poles;n++){
//              Px_d[i][j][k][n]=0.0;
//              Px_d_n[i][j][k][n]=0.0;
//              Px_d_n_1[i][j][k][n]=0.0;
//              Py_d[i][j][k][n]=0.0;
//              Py_d_n[i][j][k][n]=0.0;
//              Py_d_n_1[i][j][k][n]=0.0;
//              Pz_d[i][j][k][n]=0.0;
//              Pz_d_n[i][j][k][n]=0.0;
//              Pz_d_n_1[i][j][k][n]=0.0;
//          }
//        }
//      }
//  }
//  printf("Here\n");
//
//     for(i=0;i<NCELLX;i++){
//       for(j=0;j<NCELLY;j++){
//         for(k=0;k<NCELLZ;k++){
//
//           mat_matrix[i][j][k] = mat_matrix_safe[i][j][k];
//
//         }
//       }
//     }
//     printf("Here\n");
//
//     FILE *Spectrum;
//     Spectrum = fopen("Spectrum.txt","w");
//
//     i=1;
//     j=1;
//
//     for(w=0;w<NUM_freq;w++){
//       fprintf(Spectrum,"%e\t%e\t%e\n",freq[w],cabs(E_Incident[w][1][1]),cabs(H_Incident[w][1][1]));
//     }
//
//     fclose(Spectrum);
//
//     fclose(Test2);
//
// }



__host__ __device__  int ThreeDMap(int i,int j,int k,int SizeZ,int SizeY){
  int num = k + SizeZ*j +SizeY*SizeZ*i;
  return num;
}


__host__ __device__  int FourDMap(int i,int j,int k,int n,int SizeN,int SizeZ,int SizeY){
  int num = n + SizeN*( k + SizeZ*j +SizeY*SizeZ*i);
  return num;
}

__host__ __device__  int TwoDMap(int i,int j,int size){
  int num = j + i*size;
  return num;
}

//
// int ThreeDMap(int i,int j,int k,int SizeZ,int SizeY){
//   int num = k + SizeZ*j +SizeY*SizeZ*i;
//   return num;
// }
//
//
// int FourDMap(int i,int j,int k,int n,int SizeN,int SizeZ,int SizeY){
//   int num = n + SizeN*( k + SizeZ*j +SizeY*SizeZ*i);
//   return num;
// }
//
// int TwoDMap(int i,int j,int size){
//   int num = j + i*size;
//   return num;
// }
