#include "ion.h"

/* Chaichalit Boyd Srisawat 11 Oct 2014 */

/* double Get_Current_time(); */
/* >>> Return current time (for timing) */

/* int make_radii_list(float *radii_p, float r_min, float r_max, float dr)  */
/* >>> Make a list of radii for computing semi-numerical process */
/* >>> radii_p must be pre-allocated */
/* >>> r_min, r_max :: minimum, maximum radii */
/* >>> dr :: step in radius in fraction 0.0 <= dr <= 1.0 */
/* >>> Return the total number of radii used */

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
int make_radii_list(float *radii_p, float r_min, float r_max, float dr_inc) {
  float r,dr;
  int i;

  /// Test if nRadii <= max_Nradii
  i = 0;
  r = r_min;
  while (r < r_max) {
    dr = (r*dr_inc) < 2.0 ? (r*dr_inc):2.0;
    r += dr;
    i++;
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(ThisTask == 0) 
    if(i > constants.max_Nradii) {
      printf("[%s][%s] ***Error: File:%s Line:%d\n",__DATE__,__TIME__,__FILE__,__LINE__);
      printf("The number of smoothing radii is greater than %d\nTerminate\n",constants.max_Nradii);
      exit(1);
    }
  r = r_min;
  i = 0;
  while (r < r_max) {
    radii_p[i] = r;
    dr = (r*dr_inc) < 2.0 ? (r*dr_inc:2.0;
    r += dr;
    i++;
  }  
  radii_p = realloc(radii_p,sizeof(float)*i);
  return i;
}
