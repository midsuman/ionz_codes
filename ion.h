/**
 * @file   ion.h
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 20:49:08 2014
 * 
 * @brief  
 * Define all functions
 * 
 */

#ifndef ION_H_

#include "srfftw.h"
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef PARALLEL
#include <mpi.h>
#include "ionz_mpi.h"
#endif

#include "global.h"

/* in ionz_misc.c */
extern int make_radii_list(float *radii_p, float r_min, float r_max);
extern double Get_Current_time();

/* in ionz_io.c */


/* in allotarrays.c */
extern float **allocate_float_2d(long N1,int N2);
extern fftw_real  ***allocate_fftw_real_3d(int N1,int N2,int N3);

/* in ionz_funcs.c */
extern void Setting_Up_Memory_For_ionz(int Nnion);
extern void smooth(fftw_real ***ro_dum,float Radii);
extern void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3);

/* in read_param.c */
struct params 
{
  int Nnion;
  float *nion;
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

#define ION_H_
#endif










