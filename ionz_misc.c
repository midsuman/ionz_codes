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
int make_radii_list(float *radii_p, float r_min, float r_max, float dt) {
  float r,dr;
  int i = 0;
  r = r_min;
  while (r < r_max) {
    radii_p[i] = r;
    dr = (r*.1) < 2.0 ? (r*.1):2.0;
    r += dr;
    i++;
  }
  radii_p = realloc(radii_p,sizeof(float)*i);
  return i;
}
