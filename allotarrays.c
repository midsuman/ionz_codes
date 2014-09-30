#include<stdlib.h>
#include<sfftw.h>
#include<srfftw.h>

fftw_real  ***allocate_fftw_real_3d(int N1,int N2,int N3)
{
  int ii,jj;
  fftw_real ***phia, *phi;

  phia=(fftw_real ***)malloc(N1 *  sizeof(fftw_real **));


  for(ii=0;ii<N1;++ii)
      phia[ii]=(fftw_real **)malloc(N2 *  sizeof(fftw_real*));

  phi = (fftw_real *) calloc(N1 * N2 * N3,sizeof(fftw_real));

  for(ii=0;ii<N1;++ii)
    for(jj=0;jj<N2;++jj)
      phia[ii][jj]=phi+ (N2*N3)*ii + N3*jj;

  return(phia);
}

float **allocate_float_2d(long N1,int N2)
/* Typically N1 is large and N2 is small. */ 
/* N1 is number of particles and N2 is number of dimensions namely 3 */
{
  float **xxa, *xx;
  long ii;
  xx = (float *) calloc((size_t)(N1 * N2),sizeof(float));

  xxa=(float**)malloc(N1 *  sizeof(float*));

  for(ii=0;ii<N1;++ii) xxa[ii]=xx + N2*ii ;

return(xxa);
}



