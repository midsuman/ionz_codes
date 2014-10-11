#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<srfftw.h>
#include"ion.h"



// arrays for storing data
void Setting_Up_Memory_For_ionz(int Nnion) {
  int jk;
  globals.nh=allocate_fftw_real_3d(globals.N1,globals.N2,globals.N3+2);
  globals.ngamma=allocate_fftw_real_3d(globals.N1,globals.N2,globals.N3+2);
  global.nxion=(fftw_real****)malloc(sizeof(fftw_real***)*globals.Nnion);
    
  for(jk=0;jk<globals.Nnion;++jk) {
      global.nxion[jk]=allocate_fftw_real_3d(globals.N1,globals.N2,globals.N3+2);
    }
  // allocate area for storing densities  DONE
    
  /* The last dimension gets padded because of using REAL FFT */
    

  // done making plan  
  
}


void smooth(fftw_real ***ro_dum,float Radii) {
  int i,j,k,index,x1,y1,z1,x2,y2,z2,a,b,c;
  float m,tempre,tempim,tot;
  fftw_real ***rosp;
  fftw_complex *A;
  fftw_complex *B;

  rfftwnd_plan p_ro; // for FFT
  rfftwnd_plan q_ro; // for FFT
  
  rosp = allocate_fftw_real_3d(globals.N1,globals.N2,globals.N3+2);
  /****************************************************/	
  /* Creating the plans for forward and reverse FFT's */
  
  p_ro = rfftw3d_create_plan(globals.N1,globals.N2,globals.N3,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);  
  q_ro = rfftw3d_create_plan(globals.N1,globals.N2,globals.N3,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE | FFTW_IN_PLACE);

  //generating the filtering function
  for(i=0;i<globals.N1;i++)
    for(j=0;j<globals.N2;j++)
      for(k=0;k<globals.N3;k++)
  	rosp[i][j][k]=0.0;
 
  //Radii is radius of the sphere in grid unit
  //generating a sphere at the centre of the box
   
  tot=0.;
  for(i=0;i<globals.N1;i++)
    for(j=0;j<globals.N2;j++)
      for(k=0;k<globals.N3;k++) {
	if((float)((globals.N1/2-i)*(globals.N1/2-i)+(globals.N2/2-j)*(globals.N2/2-j)+(globals.N3/2-k)*(globals.N3/2-k))<=Radii*Radii)
	  rosp[i][j][k]=1.0;//centre N1/2,N2/2,N3/2
	tot += rosp[i][j][k];
      }
  //Sphere generation complete 
  //Doing Fourier Transform of the sphere
  rfftwnd_one_real_to_complex(p_ro,&rosp[0][0][0], NULL);
  B=(fftw_complex*)&(rosp[0][0][0]);

  //We will multiply the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to one corner of the box from box centre while applying boundary condition below
  //----------------------------------------------------------------------

  //Doing Fourier Transform of the density field
  rfftwnd_one_real_to_complex(p_ro,&ro_dum[0][0][0], NULL);
  A=(fftw_complex*)&(ro_dum[0][0][0]);
  

  for(i=0;i<globals.N1;i++)
    for(j=0;j<globals.N2;j++)
      for(k=0;k<=globals.N3/2;k++)    { 
	index = i*globals.N2*(globals.N3/2 +1) + j*(globals.N3/2 +1) + k;
	tempre=(A[index].re*B[index].re-A[index].im*B[index].im)*powf((-1.),1.*(i+j+k))/tot;
	tempim=(A[index].im*B[index].re+A[index].re*B[index].im)*powf((-1.),1.*(i+j+k))/tot;
	//multiplying the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to box centre from one corner of the box after complex to real FT
	A[index].re=tempre;
	A[index].im=tempim;
      }

  rfftwnd_one_complex_to_real(q_ro,(fftw_complex *) &ro_dum[0][0][0], NULL);
  for(i=0;i<globals.N1;i++)
    for(j=0;j<globals.N2;j++)
      for(k=0;k<=globals.N3;k++)
  	ro_dum[i][j][k]=ro_dum[i][j][k]/(globals.N1*globals.N2*globals.N3);

  fftw_free(rosp);
  fftw_destroy_plan(p_ro);
  fftw_destroy_plan(q_ro);
  /* A and B are aliases so there is no need to free them... Boyd */
  // fftw_free(A);
  // fftw_free(B);
}

/**************************************************************************
                              FUNCTIONS
**************************************************************************/
