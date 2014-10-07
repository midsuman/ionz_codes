#ifndef ION_H

#ifdef PARALLEL
#include <mpi.h>
#include "ionz_mpi.h"
#endif

/* in allotarrays.c */
extern float **allocate_float_2d(long N1,int N2);
extern fftw_real  ***allocate_fftw_real_3d(int N1,int N2,int N3);

/* in ionz_funcs.c */
extern void Setting_Up_Memory_For_ionz(int Nnion);
extern void smooth(fftw_real ***ro_dum,float Radii);

/* in read_param.c */
struct params 
{
  int Nnion;
  int *nion;
  float a_expansion;
  float z;
  float Hubble_h;
  float omegam;
  float omegalam;
  float omegab; 
  int N1,N2,N3;
  float boxsize;
  float gridsize;
};
extern struct params input_param;
extern void read_param(char filename[2048]);

#define ION_H
#endif










