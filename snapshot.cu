#include <stdio.h>
#include <stdlib.h>
#include "extern_var.h"
#include <math.h>

//pre-name of snapshot files
static char snapshot_name[10]={'s','i','m'};
static char aux_snapshot_name[10]={'a','u','x','_','s','i','m'};

void SETUP_SNAPSHOT(void){
    //how many spacial steps do you want to skip in each direction?
    // x_skip = 1;
    // y_skip = 1;
    // z_skip = 1;
    // if(NCELLX > 30) x_skip=2;
    // if(NCELLY > 30) y_skip=1;
    // if(NCELLZ >100) z_skip=1;
    //how many time steps in between snapshots?
    //t_skip=10;
    //counter
    snapshot_count=0;

}

void SNAPSHOT_2D(){
    int i,j,k;
    char filename[100];
    real ETOT;
    //Is it the appropriate time step?
    FILE *Snap;

    //fprintf(Snap,"%d\t%d\n",NCELLX/x_skip,NCELLY/y_skip);
    int XZ = 0;
    int XY = 0;
    int YZ = 1;
    if(XZ == 1){
      sprintf(filename,"videofiles/XZ_%s.%d.txt",snapshot_name,snapshot_count);
      Snap=fopen(filename,"w");
    for(i=0;i<NCELLX;i+=x_skip){
       for(k=0;k<NCELLZ;k+=z_skip){
           // for(k=0;k<NCELLZ;k+=z_skip){

                    #ifdef DOUBLECOMPLEX
                    fprintf(Snap,"%f,",creal(ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]));
                    #endif
                    #ifndef DOUBLECOMPLEX
		                 ETOT = ex[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] * ex[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + ey[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*ey[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
                    // ETOT = hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] * hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
                  //  ETOT = hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] * hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + hy[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*hy[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + hz[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*hz[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
              //    ETOT = ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
                    ETOT = pow(ETOT,0.5);
		    fprintf(Snap,"%e ",ETOT);
                    #endif
            }
            fprintf(Snap,"\n");

          //  }
      }
    }
    if(YZ == 1){
      sprintf(filename,"videofiles/YZ_%s.%d.txt",snapshot_name,snapshot_count);
      Snap=fopen(filename,"w");
      for(j=0;j<NCELLY;j+=y_skip){
         for(k=0;k<NCELLZ;k+=z_skip){
           ETOT = ex[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] * ex[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)]+  ey[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)]*ey[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] + ez[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)]*ez[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)];
           ETOT = pow(ETOT,0.5);
           fprintf(Snap,"%e ",ETOT);

         }
         fprintf(Snap,"\n");

       }
    }
        fclose(Snap);
        snapshot_count++;

    return;
}



void SNAPSHOT_2D_N(){
    int i,j,k;
    char filename[100];
    real ETOT;
    //Is it the appropriate time step?
    FILE *Snap;

    //fprintf(Snap,"%d\t%d\n",NCELLX/x_skip,NCELLY/y_skip);
    int XZ = 0;
    int XY = 0;
    int YZ = 1;

    if(XZ == 1){
      sprintf(filename,"videofiles/XZ_%s.%d_N.txt",snapshot_name,snapshot_count);
      Snap=fopen(filename,"w");
    for(i=0;i<NCELLX;i+=x_skip){
       for(k=0;k<NCELLZ;k+=z_skip){
           // for(k=0;k<NCELLZ;k+=z_skip){

                    #ifdef DOUBLECOMPLEX
                    fprintf(Snap,"%f,",creal(ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]));
                    #endif
                    #ifndef DOUBLECOMPLEX
		                 // ETOT = ex[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] * ex[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + ey[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*ey[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
                    // ETOT = hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] * hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
                  //  ETOT = hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] * hx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + hy[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*hy[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + hz[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*hz[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
              //    ETOT = ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)]*ez[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
                    // ETOT = pow(ETOT,0.5);
                    ETOT = NDx[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + NDy[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)] + NDz[ThreeDMap(i,(int)(NCELLY/2),k,NCELLZ,NCELLY)];
		    fprintf(Snap,"%e ",ETOT/3.0);
                    #endif
            }
            fprintf(Snap,"\n");

          //  }
      }
    }
    if(YZ == 1){
      sprintf(filename,"videofiles/YZ_%s.%d_N.txt",snapshot_name,snapshot_count);
      Snap=fopen(filename,"w");
      for(j=0;j<NCELLY;j+=y_skip){
         for(k=0;k<NCELLZ;k+=z_skip){
           // ETOT = ex[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] * ex[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)]+  ey[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)]*ey[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] + ez[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)]*ez[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)];
           // ETOT = pow(ETOT,0.5);
           // printf("%e\t%d\t%d\n", NDx[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] ,j,k);
           ETOT = NDx[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] + NDy[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)] + NDz[ThreeDMap((int)(NCELLX/2),j,k,NCELLZ,NCELLY)];
           fprintf(Snap,"%e ",ETOT/3.0);

         }
         fprintf(Snap,"\n");

       }
    }
        fclose(Snap);
        snapshot_count++;

    return;
}


void SNAPSHOT_1D(){
    int i;
    char filename[100];
    FILE *Snap;
    sprintf(filename,"%s.%d.txt",aux_snapshot_name,snapshot_count);
    Snap=fopen(filename,"w");
    for(i=0;i<NCELLZ;i+=z_skip){
      #ifdef DOUBLECOMPLEX
        fprintf(Snap,"%e\n",creal(hx[ThreeDMap(1,1,i,NCELLZ,NCELLY)]));
      #endif

      #ifndef DOUBLECOMPLEX
      //fprintf(Snap,"%e\n",Py_d[2][2][i][0]);
    fprintf(Snap,"%e\n",ey[ThreeDMap(2,2,i,NCELLZ,NCELLY)]);
      #endif
    }
    fclose(Snap);




    sprintf(filename,"%s.%d_2.txt",aux_snapshot_name,snapshot_count);
    Snap=fopen(filename,"w");
    for(i=0;i<inc_Length;i+=z_skip){
      #ifdef DOUBLECOMPLEX
        fprintf(Snap,"%e\n",creal(hx[ThreeDMap(1,1,i,NCELLZ,NCELLY)]));
      #endif

      #ifndef DOUBLECOMPLEX
      fprintf(Snap,"%e\n",e_inc[i]);
    // fprintf(Snap,"%e\n",ey[2][2][i]);
      #endif
    }
    fclose(Snap);
    snapshot_count++;


    return;
}
