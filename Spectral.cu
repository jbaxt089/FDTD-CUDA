#include <stdio.h>
#include <stdlib.h>
#include "extern_var.h"
#include <complex.h>
 #include <math.h>


 void Fourier_Transform(void){
    int freq_count,i,j,II,JJ,KK;
    int k=inc_plane - Spect_loc;
    int k2 = NCELLZ - NcpmlZ -5;
    if(PBC_CTW == 0) {
	k = NtfsfZ - 5;
	k2 = NCELLZ-NtfsfZ-5;

}

#ifdef FlOATPRECISION
  float complex TransExpE[NUM_freq];
  float complex TransExpH[NUM_freq];
#else
  double complex TransExpE[NUM_freq];
  double complex TransExpH[NUM_freq];
#endif


comp TransVecERe[NUM_freq];
comp TransVecHRe[NUM_freq];
comp TransVecEIm[NUM_freq];
comp TransVecHIm[NUM_freq];

  if(Periodic_XY){

    for(freq_count=0;freq_count<NUM_freq;freq_count++){
      TransExpE[freq_count] = cexp(-I*2*pi*t*dt*freq[freq_count]);
      //#pragma omp parallel for private(i,j) collapse(2)
            for(j=0;j<NCELLY;j++){
                for(i=0;i<NCELLX;i++){
                    if(TEz){

                      #ifdef DOUBLECOMPLEX
                      comp  ex_tot = (ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*(i+0.5)*dx*k_x)*cexp(I*j*dy*k_y);
                      comp  ey_tot = (ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*i*dx*k_x)*cexp(I*(j+0.5)*dy*k_y);

                      comp  hx_tot = (hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*(i)*dx*k_x)*cexp(I*(j+0.5)*dy*k_y);
                      comp  hy_tot = (hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*(i+0.5)*dx*k_x)*cexp(I*(j)*dy*k_y);
                        #endif

                        #ifndef DOUBLECOMPLEX
                        comp ex_tot = ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
                        comp ey_tot = ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
                        comp hx_tot = hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
                        comp hy_tot = hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)];

                        comp ex_tot2 = ex[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        comp ey_tot2 = ey[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        comp hx_tot2 = hx[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        comp hy_tot2 = hy[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        #endif



                         Ey_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ey_tot*TransExpE[freq_count] ;
                         Ex_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ex_tot*TransExpE[freq_count] ;
                         //Ez_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ez_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);

                         Hx_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot*TransExpE[freq_count];
                         Hy_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot*TransExpE[freq_count] ;


                         Ey_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ey_tot2*TransExpE[freq_count] ;
                         Ex_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ex_tot2*TransExpE[freq_count] ;
                         //Ez_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ez_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);

                         Hx_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot2*TransExpE[freq_count] ;
                         Hy_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot2*TransExpE[freq_count] ;

                      }

                    else if(TMz){

                      #ifdef DOUBLECOMPLEX
                      comp  ex_tot = (ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*(i+0.5)*dx*k_x)*cexp(I*j*dy*k_y);
                      comp  ey_tot = (ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*i*dx*k_x)*cexp(I*(j+0.5)*dy*k_y);

                      comp  hx_tot = (hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*(i)*dx*k_x)*cexp(I*(j+0.5)*dy*k_y);
                      comp  hy_tot = (hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)])*cexp(I*(i+0.5)*dx*k_x)*cexp(I*(j)*dy*k_y);
                        #endif

                        #ifndef DOUBLECOMPLEX
                        comp ex_tot = ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
                        comp ey_tot = ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
                        comp hx_tot = hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)];
                        comp hy_tot = hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)];

                        comp ex_tot2 = ex[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        comp ey_tot2 = ey[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        comp hx_tot2 = hx[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        comp hy_tot2 = hy[ThreeDMap(i,j,k2,NCELLZ,NCELLY)];
                        #endif


                         Ey_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ey_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);
                         Ex_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ex_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);
                         //Ez_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ez_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);

                         Hx_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);
                         Hy_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);


                         Ey_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ey_tot2*cexp(-I*2*pi*t*dt*freq[freq_count]);
                         Ex_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ex_tot2*cexp(-I*2*pi*t*dt*freq[freq_count]);
                         //Ez_Reflected[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += ez_tot*cexp(-I*2*pi*t*dt*freq[freq_count]);

                         Hx_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot2*cexp(-I*2*pi*t*dt*freq[freq_count]);
                         Hy_Transmitted[ThreeDMap(freq_count,i,j,NCELLY,NCELLX)] += hx_tot2*cexp(-I*2*pi*t*dt*freq[freq_count]);


                    }
                }
            }

    }
}
else{
    if(Absorption == 1){







      for(freq_count=0;freq_count<NUM_freq;freq_count++){
        TransVecERe[freq_count] = cos(2.0*pi*(t+1.0)*dt*freq[freq_count]);
        TransVecEIm[freq_count] = sin(2.0*pi*(t+1.0)*dt*freq[freq_count]);
        TransVecHRe[freq_count] = cos(2.0*pi*(t+0.5)*dt*freq[freq_count]);
        TransVecHIm[freq_count] = sin(2.0*pi*(t+0.5)*dt*freq[freq_count]);

        //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
            for(i=XSTARTAbs;i<XENDAbs;i++){
              for(j=YSTARTAbs;j<YENDAbs;j++){
              //  printf("%d\t%d\t%d\t%d\n",i,j,XENDAbs,YENDAbs);
                k = ZNEARAbs;
                II = i - XSTARTAbs;
                JJ = j - YSTARTAbs;

                ExTransformNearZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformNearZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
              //  EzTransformNearZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformNearZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformNearZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
            //    EzTransformNearZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformNearZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformNearZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
              //  HzTransformNearZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformNearZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformNearZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HzTransformNearZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];

                k = ZFARAbs;

                ExTransformFarZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformFarZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
            //    EzTransformFarZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformFarZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformFarZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
          //      EzTransformFarZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformFarZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformFarZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
            //    HzTransformFarZAbsRe[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformFarZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformFarZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HzTransformFarZAbsIm[ThreeDMap(freq_count,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
              }
            }

            //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
            for(i=XSTARTAbs;i<XENDAbs;i++){
              for(k=ZSTARTAbs;k<ZENDAbs;k++){
                II = i - XSTARTAbs;
                KK = k - ZSTARTAbs;

                j=YNEARAbs;

                ExTransformNearYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
            //    EyTransformNearYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformNearYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformNearYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
          //      EyTransformNearYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformNearYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformNearYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
            //    HyTransformNearYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformNearYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformNearYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HyTransformNearYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformNearYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];

                j=YFARAbs;

                ExTransformFarYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
          //      EyTransformFarYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformFarYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformFarYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
          //      EyTransformFarYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformFarYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformFarYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
          //      HyTransformFarYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformFarYAbsRe[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformFarYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HyTransformFarYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformFarYAbsIm[ThreeDMap(freq_count,II,KK,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
              }
            }

            //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
            for(j=YSTARTAbs;j<YENDAbs;j++){
              for(k=ZSTARTAbs;k<ZENDAbs;k++){
                KK = k - ZSTARTAbs;
                JJ = j - YSTARTAbs;

                i=XNEARAbs;

            //    ExTransformNearXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformNearXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformNearXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

        //        ExTransformNearXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformNearXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformNearXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

          //      HxTransformNearXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformNearXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformNearXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

            //    HxTransformNearXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformNearXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformNearXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];

                i=XFARAbs;

            //    ExTransformFarXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformFarXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformFarXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

          //      ExTransformFarXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformFarXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformFarXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

            //    HxTransformFarXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformFarXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformFarXAbsRe[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

              //  HxTransformFarXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformFarXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformFarXAbsIm[ThreeDMap(freq_count,JJ,KK,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
              }
            }
    }
  }

    if(Scattering == 1){
      for(freq_count=0;freq_count<NUM_freq;freq_count++){
        //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
            for(i=XSTARTSca;i<XENDSca;i++){
              for(j=YSTARTSca;j<YENDSca;j++){
                k = ZNEARSca;
                II = i - XSTARTSca;
                JJ = j - YSTARTSca;

                ExTransformNearZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformNearZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
            //    EzTransformNearZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformNearZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformNearZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
            //    EzTransformNearZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformNearZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformNearZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
            //    HzTransformNearZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformNearZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformNearZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HzTransformNearZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];

                k = ZFARSca;

                ExTransformFarZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformFarZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
            //    EzTransformFarZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformFarZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformFarZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
            //    EzTransformFarZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformFarZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformFarZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
            //    HzTransformFarZScaRe[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformFarZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformFarZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HzTransformFarZScaIm[ThreeDMap(freq_count,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
              }
            }

            //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
            for(i=XSTARTSca;i<XENDSca;i++){
              for(k=ZSTARTSca;k<ZENDSca;k++){

                j=YNEARSca;
                II = i - XSTARTSca;
                KK = k - ZSTARTSca;
                ExTransformNearYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
              //  EyTransformNearYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformNearYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformNearYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
            //    EyTransformNearYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformNearYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformNearYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
              //  HyTransformNearYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformNearYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformNearYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HyTransformNearYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformNearYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];

                j=YFARSca;

                ExTransformFarYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
            //    EyTransformFarYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformFarYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

                ExTransformFarYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
            //    EyTransformFarYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformFarYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

                HxTransformFarYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
            //    HyTransformFarYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformFarYScaRe[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

                HxTransformFarYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
            //    HyTransformFarYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformFarYScaIm[ThreeDMap(freq_count,II,KK,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
              }
            }

            //#pragma omp parallel for private(i,j,k,II,JJ,KK) collapse(2)
            for(j=YSTARTSca;j<YENDSca;j++){
              for(k=ZSTARTSca;k<ZENDSca;k++){
                KK = k - ZSTARTSca;
                JJ = j - YSTARTSca;
                i=XNEARSca;

            //    ExTransformNearXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformNearXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformNearXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

        //        ExTransformNearXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformNearXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformNearXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

        //        HxTransformNearXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformNearXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformNearXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

        //        HxTransformNearXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformNearXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformNearXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];

                i=XFARSca;

          //      ExTransformFarXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EyTransformFarXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];
                EzTransformFarXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecERe[freq_count];

            //    ExTransformFarXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ex[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EyTransformFarXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ey[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];
                EzTransformFarXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += ez[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecEIm[freq_count];

            //    HxTransformFarXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HyTransformFarXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];
                HzTransformFarXScaRe[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHRe[freq_count];

            //    HxTransformFarXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hx[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HyTransformFarXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hy[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
                HzTransformFarXScaIm[ThreeDMap(freq_count,JJ,KK,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)] += hz[ThreeDMap(i,j,k,NCELLZ,NCELLY)]*TransVecHIm[freq_count];
              }
            }
    }
}
}
    return;
}


void CalculateAbsScatt(void){
  real  P_r,P_t,P_r_inc, Area, dsurf1, dsurf2,k_rho,k_z;
  real u,v,u2,v2;
//  #ifdef FlOATPRECISION
  float complex  zz,zz_trans;
  int i, k,j,II,JJ,KK, m;
  float complex cEx, cEy, ccHx, ccHy;
  float complex cEx_trans, cEy_trans, ccHx_trans, ccHy_trans;
  float complex cEx_inc, cEy_inc, ccHx_inc, cc_Hz_inc;
  int freq_count;
  real2 *QNEAXsca,*QFARXsca,*QNEAYsca,*QFARYsca,*QNEAZsca,*QFARZsca;
  real2 *QNEAXabs,*QFARXabs,*QNEAYabs,*QFARYabs,*QNEAZabs,*QFARZabs;
  float Amplitude;
  // #else
  // double complex  zz,zz_trans;
  // int i, k,j,II,JJ,KK, m;
  // double complex cEx, cEy, ccHx, ccHy;
  // double complex cEx_trans, cEy_trans, ccHx_trans, ccHy_trans;
  // double complex cEx_inc, cEy_inc, ccHx_inc, cc_Hz_inc;
  // int freq_count;
  // double *QNEAXsca,*QFARXsca,*QNEAYsca,*QFARYsca,*QNEAZsca,*QFARZsca;
  // double *QNEAXabs,*QFARXabs,*QNEAYabs,*QFARYabs,*QNEAZabs,*QFARZabs;
  // double Amplitude;
  // #endif
  QNEAXsca = MALLOC1D_Real2(QNEAXsca,NUM_freq);
  QNEAYsca = MALLOC1D_Real2(QNEAYsca,NUM_freq);
  QNEAZsca = MALLOC1D_Real2(QNEAZsca,NUM_freq);
  QFARXsca = MALLOC1D_Real2(QFARXsca,NUM_freq);
  QFARYsca = MALLOC1D_Real2(QFARYsca,NUM_freq);
  QFARZsca = MALLOC1D_Real2(QFARZsca,NUM_freq);
  QNEAXabs = MALLOC1D_Real2(QNEAXabs,NUM_freq);
  QNEAYabs = MALLOC1D_Real2(QNEAYabs,NUM_freq);
  QNEAZabs = MALLOC1D_Real2(QNEAZabs,NUM_freq);
  QFARXabs = MALLOC1D_Real2(QFARXabs,NUM_freq);
  QFARYabs = MALLOC1D_Real2(QFARYabs,NUM_freq);
  QFARZabs = MALLOC1D_Real2(QFARZabs,NUM_freq);




  //
  // if(Absorption == 1){
  //   for(freq_count=0;freq_count<NUM_freq;freq_count++){
  //     m = freq_count;
  //     // TransVecERe[freq_count] = cos(2.0*pi*(t+1.0)*dt*freq[freq_count]);
  //     // TransVecEIm[freq_count] = sin(2.0*pi*(t+1.0)*dt*freq[freq_count]);
  //     // TransVecHRe[freq_count] = cos(2.0*pi*(t+0.5)*dt*freq[freq_count]);
  //     // TransVecHIm[freq_count] = sin(2.0*pi*(t+0.5)*dt*freq[freq_count]);
  //
  //     //#pragma omp parallel for private(i,j) collapse(2)
  //         for(i=XSTARTAbs;i<XENDAbs;i++){
  //           for(j=YSTARTAbs;j<YENDAbs;j++){
  //           //  printf("%d\t%d\t%d\t%d\n",i,j,XENDAbs,YENDAbs);
  //             II = i - XSTARTAbs;
  //             JJ = j - YSTARTAbs;
  //
  //             ExTransformNearZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             EyTransformNearZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             EzTransformNearZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //
  //             ExTransformNearZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             EyTransformNearZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             EzTransformNearZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //
  //             HxTransformNearZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             HyTransformNearZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             HzTransformNearZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //
  //             HxTransformNearZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             HyTransformNearZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             HzTransformNearZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //
  //
  //             ExTransformFarZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             EyTransformFarZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             EzTransformFarZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //
  //             ExTransformFarZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             EyTransformFarZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             EzTransformFarZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //
  //             HxTransformFarZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             HyTransformFarZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //             HzTransformFarZAbsRe[freq_count][II][JJ] /= (Amplitude);
  //
  //             HxTransformFarZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             HyTransformFarZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //             HzTransformFarZAbsIm[freq_count][II][JJ] /= (Amplitude);
  //           }
  //         }
  //
  //         //#pragma omp parallel for private(i,k) collapse(2)
  //         for(i=XSTARTAbs;i<XENDAbs;i++){
  //           for(k=ZSTARTAbs;k<ZENDAbs;k++){
  //             II = i - XSTARTAbs;
  //             KK = k - ZSTARTAbs;
  //
  //
  //             ExTransformNearYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             EyTransformNearYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             EzTransformNearYAbsRe[freq_count][II][KK] /= (Amplitude);
  //
  //             ExTransformNearYAbsIm[freq_count][II][KK]/= (Amplitude);
  //             EyTransformNearYAbsIm[freq_count][II][KK]/= (Amplitude);
  //             EzTransformNearYAbsIm[freq_count][II][KK]/= (Amplitude);
  //
  //             HxTransformNearYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             HyTransformNearYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             HzTransformNearYAbsRe[freq_count][II][KK] /= (Amplitude);
  //
  //             HxTransformNearYAbsIm[freq_count][II][KK] /= (Amplitude);
  //             HyTransformNearYAbsIm[freq_count][II][KK] /= (Amplitude);
  //             HzTransformNearYAbsIm[freq_count][II][KK]/= (Amplitude);
  //
  //             ExTransformFarYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             EyTransformFarYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             EzTransformFarYAbsRe[freq_count][II][KK] /= (Amplitude);
  //
  //             ExTransformFarYAbsIm[freq_count][II][KK] /= (Amplitude);
  //             EyTransformFarYAbsIm[freq_count][II][KK] /= (Amplitude);
  //             EzTransformFarYAbsIm[freq_count][II][KK] /= (Amplitude);
  //
  //             HxTransformFarYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             HyTransformFarYAbsRe[freq_count][II][KK] /= (Amplitude);
  //             HzTransformFarYAbsRe[freq_count][II][KK] /= (Amplitude);
  //
  //             HxTransformFarYAbsIm[freq_count][II][KK] /= (Amplitude);
  //             HyTransformFarYAbsIm[freq_count][II][KK] /= (Amplitude);
  //             HzTransformFarYAbsIm[freq_count][II][KK] /= (Amplitude);
  //           }
  //         }
  //
  //         //#pragma omp parallel for private(j,k) collapse(2)
  //         for(j=YSTARTAbs;j<YENDAbs;j++){
  //           for(k=ZSTARTAbs;k<ZENDAbs;k++){
  //             KK = k - ZSTARTAbs;
  //             JJ = j - YSTARTAbs;
  //
  //
  //             ExTransformNearXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             EyTransformNearXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             EzTransformNearXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //             ExTransformNearXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //             EyTransformNearXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //             EzTransformNearXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //
  //             HxTransformNearXAbsRe[freq_count][JJ][KK]/= (Amplitude);
  //             HyTransformNearXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             HzTransformNearXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //             HxTransformNearXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //             HyTransformNearXAbsIm[freq_count][JJ][KK]/= (Amplitude);
  //             HzTransformNearXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //
  //
  //             ExTransformFarXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             EyTransformFarXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             EzTransformFarXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //             ExTransformFarXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //             EyTransformFarXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //             EzTransformFarXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //
  //             HxTransformFarXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             HyTransformFarXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //             HzTransformFarXAbsRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //             HxTransformFarXAbsIm[freq_count][JJ][KK]/= (Amplitude);
  //             HyTransformFarXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //             HzTransformFarXAbsIm[freq_count][JJ][KK] /= (Amplitude);
  //           }
  //         }
  // }
  // }
  //
  //
  //   if(Scattering == 1){
  //     for(freq_count=0;freq_count<NUM_freq;freq_count++){
  //       // TransVecERe[freq_count] = cos(2.0*pi*(t+1.0)*dt*freq[freq_count]);
  //       // TransVecEIm[freq_count] = sin(2.0*pi*(t+1.0)*dt*freq[freq_count]);
  //       // TransVecHRe[freq_count] = cos(2.0*pi*(t+0.5)*dt*freq[freq_count]);
  //       // TransVecHIm[freq_count] = sin(2.0*pi*(t+0.5)*dt*freq[freq_count]);
  //
  //       //#pragma omp parallel for private(i,j) collapse(2)
  //           for(i=XSTARTSca;i<XENDSca;i++){
  //             for(j=YSTARTSca;j<YENDSca;j++){
  //             //  printf("%d\t%d\t%d\t%d\n",i,j,XENDSca,YENDSca);
  //               II = i - XSTARTSca;
  //               JJ = j - YSTARTSca;
  //
  //               ExTransformNearZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               EyTransformNearZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               EzTransformNearZScaRe[freq_count][II][JJ] /= (Amplitude);
  //
  //               ExTransformNearZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               EyTransformNearZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               EzTransformNearZScaIm[freq_count][II][JJ] /= (Amplitude);
  //
  //               HxTransformNearZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               HyTransformNearZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               HzTransformNearZScaRe[freq_count][II][JJ] /= (Amplitude);
  //
  //               HxTransformNearZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               HyTransformNearZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               HzTransformNearZScaIm[freq_count][II][JJ] /= (Amplitude);
  //
  //
  //               ExTransformFarZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               EyTransformFarZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               EzTransformFarZScaRe[freq_count][II][JJ] /= (Amplitude);
  //
  //               ExTransformFarZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               EyTransformFarZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               EzTransformFarZScaIm[freq_count][II][JJ] /= (Amplitude);
  //
  //               HxTransformFarZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               HyTransformFarZScaRe[freq_count][II][JJ] /= (Amplitude);
  //               HzTransformFarZScaRe[freq_count][II][JJ] /= (Amplitude);
  //
  //               HxTransformFarZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               HyTransformFarZScaIm[freq_count][II][JJ] /= (Amplitude);
  //               HzTransformFarZScaIm[freq_count][II][JJ] /= (Amplitude);
  //             }
  //           }
  //
  //           //#pragma omp parallel for private(i,k) collapse(2)
  //           for(i=XSTARTSca;i<XENDSca;i++){
  //             for(k=ZSTARTSca;k<ZENDSca;k++){
  //               II = i - XSTARTSca;
  //               KK = k - ZSTARTSca;
  //
  //
  //               ExTransformNearYScaRe[freq_count][II][KK] /= (Amplitude);
  //               EyTransformNearYScaRe[freq_count][II][KK] /= (Amplitude);
  //               EzTransformNearYScaRe[freq_count][II][KK] /= (Amplitude);
  //
  //               ExTransformNearYScaIm[freq_count][II][KK]/= (Amplitude);
  //               EyTransformNearYScaIm[freq_count][II][KK]/= (Amplitude);
  //               EzTransformNearYScaIm[freq_count][II][KK]/= (Amplitude);
  //
  //               HxTransformNearYScaRe[freq_count][II][KK] /= (Amplitude);
  //               HyTransformNearYScaRe[freq_count][II][KK] /= (Amplitude);
  //               HzTransformNearYScaRe[freq_count][II][KK] /= (Amplitude);
  //
  //               HxTransformNearYScaIm[freq_count][II][KK] /= (Amplitude);
  //               HyTransformNearYScaIm[freq_count][II][KK] /= (Amplitude);
  //               HzTransformNearYScaIm[freq_count][II][KK]/= (Amplitude);
  //
  //               ExTransformFarYScaRe[freq_count][II][KK] /= (Amplitude);
  //               EyTransformFarYScaRe[freq_count][II][KK] /= (Amplitude);
  //               EzTransformFarYScaRe[freq_count][II][KK] /= (Amplitude);
  //
  //               ExTransformFarYScaIm[freq_count][II][KK] /= (Amplitude);
  //               EyTransformFarYScaIm[freq_count][II][KK] /= (Amplitude);
  //               EzTransformFarYScaIm[freq_count][II][KK] /= (Amplitude);
  //
  //               HxTransformFarYScaRe[freq_count][II][KK] /= (Amplitude);
  //               HyTransformFarYScaRe[freq_count][II][KK] /= (Amplitude);
  //               HzTransformFarYScaRe[freq_count][II][KK] /= (Amplitude);
  //
  //               HxTransformFarYScaIm[freq_count][II][KK] /= (Amplitude);
  //               HyTransformFarYScaIm[freq_count][II][KK] /= (Amplitude);
  //               HzTransformFarYScaIm[freq_count][II][KK] /= (Amplitude);
  //             }
  //           }
  //
  //           //#pragma omp parallel for private(j,k) collapse(2)
  //           for(j=YSTARTSca;j<YENDSca;j++){
  //             for(k=ZSTARTSca;k<ZENDSca;k++){
  //               KK = k - ZSTARTSca;
  //               JJ = j - YSTARTSca;
  //
  //
  //               ExTransformNearXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               EyTransformNearXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               EzTransformNearXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //               ExTransformNearXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //               EyTransformNearXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //               EzTransformNearXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //
  //               HxTransformNearXScaRe[freq_count][JJ][KK]/= (Amplitude);
  //               HyTransformNearXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               HzTransformNearXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //               HxTransformNearXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //               HyTransformNearXScaIm[freq_count][JJ][KK]/= (Amplitude);
  //               HzTransformNearXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //
  //
  //               ExTransformFarXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               EyTransformFarXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               EzTransformFarXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //               ExTransformFarXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //               EyTransformFarXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //               EzTransformFarXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //
  //               HxTransformFarXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               HyTransformFarXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //               HzTransformFarXScaRe[freq_count][JJ][KK] /= (Amplitude);
  //
  //               HxTransformFarXScaIm[freq_count][JJ][KK]/= (Amplitude);
  //               HyTransformFarXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //               HzTransformFarXScaIm[freq_count][JJ][KK] /= (Amplitude);
  //             }
  //           }
  //   }
  //   }
  //



  for (m = 0; m < NUM_freq; m++) {

      Area = 0.0;
      P_r = 0.0;
      P_r_inc = 0.0;
      P_t = 0.0;
      Amplitude = (float)cabs(E_incident[m]);
      // Amplitude = 1.0;
      // if(t==0) Amplitude = 1.0;

    //  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

    for(i=XSTARTSca;i<XENDSca;i++){
      dsurf1 = dx;
      for(j=YSTARTSca;j<YENDSca;j++){
        dsurf2 = dy;
        II = i - XSTARTSca;
        JJ = j - YSTARTSca;


              u = EyTransformNearZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);
              v = EyTransformNearZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);
              u2 = ExTransformNearZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);
              v2 = ExTransformNearZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);

            cEy = u + I*v;
            cEx = u2 + I*v2;

            u = EyTransformFarZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);
            v = EyTransformFarZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);
            u2 = ExTransformFarZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);
            v2 = ExTransformFarZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)] / (Amplitude);


            cEy_trans = u + I*v;
            cEx_trans = u2 + I*v2;

              u = HxTransformNearZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
              v = HxTransformNearZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
              u2 = HyTransformNearZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
              v2 = HyTransformNearZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);



            ccHx = u - I*v;
            ccHy = u2 - I*v2;


              u = HxTransformFarZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
              v = HxTransformFarZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
              u2 = HyTransformFarZScaRe[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
              v2 = HyTransformFarZScaIm[ThreeDMap(m,II,JJ,YENDSca-YSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);// + ImTFHz_xz_n[m][i][K];

            ccHx_trans = u - I*v;
            ccHy_trans = u2 - I*v2;

            zz = cEx*ccHy - cEy*ccHx;
            zz_trans = cEx_trans*ccHy_trans - cEy_trans*ccHx_trans;


            P_r += creal(zz)*dsurf1*dsurf2;
            P_t += creal(zz_trans)*dsurf1*dsurf2;


          }
      }
    QNEAZsca[m] = P_r;
    QFARZsca[m] = P_t;

    Area = 0.0;
    P_r = 0.0;
    P_r_inc = 0.0;
    P_t = 0.0;



  //  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

  for(i=XSTARTSca;i<XENDSca;i++){
    dsurf1 = dx;
    for(j=ZSTARTSca;j<ZENDSca;j++){
      dsurf2 = dy;
      II = i - XSTARTSca;
      JJ = j - ZSTARTSca;


            u = ExTransformNearYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            v = ExTransformNearYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            u2 = EzTransformNearYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            v2 = EzTransformNearYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);

          cEy = u + I*v;
          cEx = u2 + I*v2;

          u = ExTransformFarYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
          v = ExTransformFarYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
          u2 = EzTransformFarYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
          v2 = EzTransformFarYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);


          cEy_trans = u + I*v;
          cEx_trans = u2 + I*v2;

            u = HzTransformNearYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            v = HzTransformNearYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            u2 = HxTransformNearYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            v2 = HxTransformNearYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);



          ccHx = u - I*v;
          ccHy = u2 - I*v2;


            u = HzTransformFarYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            v = HzTransformFarYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            u2 = HxTransformFarYScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);
            v2 = HxTransformFarYScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,XENDSca-XSTARTSca+1)]/(Amplitude);// + ImTFHz_xz_n[m][i][K];

          ccHx_trans = u - I*v;
          ccHy_trans = u2 - I*v2;

          zz = cEx*ccHy - cEy*ccHx;
          zz_trans = cEx_trans*ccHy_trans - cEy_trans*ccHx_trans;

          P_r += creal(zz)*dsurf1*dsurf2;
          P_t += creal(zz_trans)*dsurf1*dsurf2;


        }
    }
  QNEAYsca[m] = P_r;
  QFARYsca[m] = P_t;

  Area = 0.0;
  P_r = 0.0;
  P_r_inc = 0.0;
  P_t = 0.0;



//  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

for(i=YSTARTSca;i<YENDSca;i++){
  dsurf1 = dx;
  for(j=ZSTARTSca;j<ZENDSca;j++){
    dsurf2 = dy;
    II = i - YSTARTSca;
    JJ = j - ZSTARTSca;


          u = EzTransformNearXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          v = EzTransformNearXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          u2 = EyTransformNearXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          v2 = EyTransformNearXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);

        cEy = u + I*v;
        cEx = u2 + I*v2;

        u = EzTransformFarXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
        v = EzTransformFarXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
        u2 = EyTransformFarXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
        v2 = EyTransformFarXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);


        cEy_trans = u + I*v;
        cEx_trans = u2 + I*v2;

          u = HyTransformNearXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          v = HyTransformNearXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          u2 = HzTransformNearXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          v2 = HzTransformNearXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);



        ccHx = u - I*v;
        ccHy = u2 - I*v2;


          v = HyTransformFarXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          u = HyTransformFarXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          u2 = HzTransformFarXScaRe[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);
          v2 = HzTransformFarXScaIm[ThreeDMap(m,II,JJ,ZENDSca-ZSTARTSca + 1,YENDSca-YSTARTSca+1)]/(Amplitude);// + ImTFHz_xz_n[m][i][K];

        ccHx_trans = u - I*v;
        ccHy_trans = u2 - I*v2;

        zz = cEx*ccHy - cEy*ccHx;
        zz_trans = cEx_trans*ccHy_trans - cEy_trans*ccHx_trans;

        P_r += creal(zz)*dsurf1*dsurf2;
        P_t += creal(zz_trans)*dsurf1*dsurf2;


      }
  }
QNEAXsca[m] = P_r;
QFARXsca[m] = P_t;













Area = 0.0;
P_r = 0.0;
P_r_inc = 0.0;
P_t = 0.0;



//  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

for(i=XSTARTAbs;i<XENDAbs;i++){
dsurf1 = dx;
for(j=YSTARTAbs;j<YENDAbs;j++){
  dsurf2 = dy;
  II = i - XSTARTAbs;
  JJ = j - YSTARTAbs;


        u = EyTransformNearZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        v = EyTransformNearZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        u2 = ExTransformNearZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        v2 = ExTransformNearZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);

      cEy = u + I*v;
      cEx = u2 + I*v2;

      u = EyTransformFarZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v = EyTransformFarZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      u2 = ExTransformFarZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v2 = ExTransformFarZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);


      cEy_trans = u + I*v;
      cEx_trans = u2 + I*v2;

        v = HxTransformNearZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        u = HxTransformNearZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        u2 = HyTransformNearZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        v2 = HyTransformNearZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);



      ccHx = u - I*v;
      ccHy = u2 - I*v2;


        u = HxTransformFarZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        v = HxTransformFarZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        u2 = HyTransformFarZAbsRe[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
        v2 = HyTransformFarZAbsIm[ThreeDMap(m,II,JJ,YENDAbs-YSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);// + ImTFHz_xz_n[m][i][K];

      ccHx_trans = u - I*v;
      ccHy_trans = u2 - I*v2;

      zz = cEx*ccHy - cEy*ccHx;
      zz_trans = cEx_trans*ccHy_trans - cEy_trans*ccHx_trans;


      P_r += creal(zz)*dsurf1*dsurf2;
      P_t += creal(zz_trans)*dsurf1*dsurf2;


    }
}
QNEAZabs[m] = P_r;
QFARZabs[m] = P_t;

Area = 0.0;
P_r = 0.0;
P_r_inc = 0.0;
P_t = 0.0;



//  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

for(i=XSTARTAbs;i<XENDAbs;i++){
dsurf1 = dx;
for(j=ZSTARTAbs;j<ZENDAbs;j++){
dsurf2 = dy;
II = i - XSTARTAbs;
JJ = j - ZSTARTAbs;


      u = ExTransformNearYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v = ExTransformNearYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      u2 = EzTransformNearYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v2 = EzTransformNearYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);

    cEy = u + I*v;
    cEx = u2 + I*v2;

    v = ExTransformFarYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
    u = ExTransformFarYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
    v2 = EzTransformFarYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
    u2 = EzTransformFarYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);


    cEy_trans = u + I*v;
    cEx_trans = u2 + I*v2;

      u = HzTransformNearYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v = HzTransformNearYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      u2 = HxTransformNearYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v2 = HxTransformNearYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);



    ccHx = u - I*v;
    ccHy = u2 - I*v2;


      u = HzTransformFarYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v = HzTransformFarYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      u2 = HxTransformFarYAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);
      v2 = HxTransformFarYAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,XENDAbs-XSTARTAbs+1)]/(Amplitude);// + ImTFHz_xz_n[m][i][K];

    ccHx_trans = u - I*v;
    ccHy_trans = u2 - I*v2;

    zz = cEx*ccHy - cEy*ccHx;
    zz_trans = cEx_trans*ccHy_trans - cEy_trans*ccHx_trans;

    P_r += creal(zz)*dsurf1*dsurf2;
    P_t += creal(zz_trans)*dsurf1*dsurf2;


  }
}
QNEAYabs[m] = P_r;
QFARYabs[m] = P_t;

Area = 0.0;
P_r = 0.0;
P_r_inc = 0.0;
P_t = 0.0;



//  fprintf(Spectrum,"%e\t%e\t%e\t%e\t%e\n",freq[m],creal(E_Reflected[m][3][3]),cimag(E_Reflected[m][3][3]),creal(H_Reflected[m][3][3]),cimag(H_Reflected[m][3][3]));

for(i=YSTARTAbs;i<YENDAbs;i++){
dsurf1 = dx;
for(j=ZSTARTAbs;j<ZENDAbs;j++){
dsurf2 = dy;
II = i - YSTARTAbs;
JJ = j - ZSTARTAbs;


    u = EzTransformNearXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    v = EzTransformNearXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    u2 = EyTransformNearXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    v2 = EyTransformNearXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);

  cEy = u + I*v;
  cEx = u2 + I*v2;

  v = EzTransformFarXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
  u = EzTransformFarXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
  u2 = EyTransformFarXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
  v2 = EyTransformFarXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);


  cEy_trans = u + I*v;
  cEx_trans = u2 + I*v2;

    u = HyTransformNearXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    v = HyTransformNearXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    u2 = HzTransformNearXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    v2 = HzTransformNearXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);





  ccHx = u - I*v;
  ccHy = u2 - I*v2;


    u = HyTransformFarXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    v = HyTransformFarXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    u2 = HzTransformFarXAbsRe[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);
    v2 = HzTransformFarXAbsIm[ThreeDMap(m,II,JJ,ZENDAbs-ZSTARTAbs + 1,YENDAbs-YSTARTAbs+1)]/(Amplitude);// + ImTFHz_xz_n[m][i][K];

  ccHx_trans = u - I*v;
  ccHy_trans = u2 - I*v2;

  zz = cEx*ccHy - cEy*ccHx;
  zz_trans = cEx_trans*ccHy_trans - cEy_trans*ccHx_trans;

  P_r += creal(zz)*dsurf1*dsurf2;
  P_t += creal(zz_trans)*dsurf1*dsurf2;


}
}
QNEAXabs[m] = P_r;
QFARXabs[m] = P_t;






  }










FILE *SCATABS;

char filename[100];
FILE *Snap;

static char name[10]={'S','c','a','t','A','b','s'};

sprintf(filename,"%s.%d.txt",name,t);
SCATABS = fopen(filename,"w");
double ScaTot,AbsTot,ExtTot;
for(m=0;m<NUM_freq;m++){

  ScaTot = -QNEAXsca[m] + QFARXsca[m] - QNEAYsca[m] + QFARYsca[m] - QNEAZsca[m] + QFARZsca[m];
  AbsTot = -QNEAXabs[m] + QFARXabs[m] - QNEAYabs[m] + QFARYabs[m] - QNEAZabs[m] + QFARZabs[m];
  if(nano_sphere == 1){
    // printf("%e\n",pi*nano_sphere_radius*nano_sphere_radius);
  ScaTot = ScaTot/(pi*nano_sphere_radius*nano_sphere_radius);
  AbsTot = -AbsTot/(pi*nano_sphere_radius*nano_sphere_radius);
  ExtTot = AbsTot + ScaTot;
}

  fprintf(SCATABS,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",c0/freq[m],ScaTot,AbsTot,ExtTot,QNEAXsca[m],QFARXsca[m] ,QNEAYsca[m] ,QFARYsca[m] ,QNEAZsca[m],QFARZsca[m],QNEAXabs[m],QFARXabs[m],QNEAYabs[m],QFARYabs[m],QNEAZabs[m],QFARZabs[m],cabs(E_incident[m]),creal(E_incident[m]),cimag(E_incident[m]));

}
fclose(SCATABS);


}
