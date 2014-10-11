/**
 * @file   ionz_misc.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 23:55:03 2014
 * 
 * @brief  Misc functions
 * 
 * 
 */

#include "ion.h"
/** 
 * Get current time
 * 
 * 
 * @return current time
 */
double Get_Current_time() {
#ifdef PARALLEL
  double time;
  time = MPI_Wtime();
  return time;
#else
  time_t timer;
  timer = time(NULL);
  return (double) timer;
#endif
}

/** 
 * Generate a list of radii to compute semi-numerical reionization 
 * 
 * @param radii_p Pointer to radii list (pre-allocated)
 * @param r_min Minimum radius
 * @param r_max Maximum radius
 * @param dr_inc Increase in radius (0.0 < dr_inc <= 1.0; 0.1 by default)
 * @param max_dr Maximum increase in radius (2.0 by default)
 * @return total number of radii 
 */
int make_radii_list(float *radii_p, float r_min, float r_max, float dr_inc = constant.dr_inc, float max_dr = constants.max_dr) {
  float r,dr;
  int i;

  /// Test if nRadii <= max_Nradii
  i = 0;
  r = r_min;
  while (r < r_max) {
    dr = min(r*dr_inc,max_dr);
    r += dr;
    i++;
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(mympi.ThisTask == 0) 
    if(i > constants.max_Nradii) {
      printf("[%s][%s] ***Error: File:%s Line:%d\n",__DATE__,__TIME__,__FILE__,__LINE__);
      printf("The number of smoothing radii is greater than %d\nTerminate\n",constants.max_Nradii);
      exit(1);
    }
  r = r_min;
  i = 0;
  /// Store radii list
  while (r < r_max) {
    radii_p[i] = r;
    dr = min(r*dr_inc,max_dr);
    r += dr;
    i++;
  }  
  radii_p = realloc(radii_p,sizeof(float)*i);
  return i;
}
