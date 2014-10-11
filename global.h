/**
 * @file   global.h
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 20:39:04 2014
 * 
 * @brief  
 * Global variables
 * 
 */

#ifndef GLOBAL_H_
#include "srfftw.h"
#include <math.h>

struct GLOBALVARS {
  /// Hubble parameter
  float  vhh;
  /// Omega_matter
  float vomegam;
  /// Omega_lambda
  float vomegalam; // Cosmological Constant 
  /// Omega_baryon
  float vomegab;
  /// grid spacing in Mpc
  float LL; 
  /// pi constant
  float  pi=M_PI;

  /// arrays for storing data
  fftw_real ***ro; /// for density/potential
  rfftwnd_plan p_ro; /// for FFT
  rfftwnd_plan q_ro; /// for FFT

  ///end of declaration of global variables for output binary file
  fftw_real ***nh,***nhs,***ngamma,***ngammas,****nxion;
} globals;

#define GLOBAL_H_
#endif
