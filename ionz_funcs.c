#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<srfftw.h>
#include"ion.h"

/*  GLOBAL VARIABLES  */

// cosmological parameters read from input file  "input.ionz"

extern float  vhh, // Hubble parameter 
  vomegam, // Omega_matter) 
  vomegalam, // Cosmological Constant 
  vomegab, //Omega_baryon */
  tcmb, //  CMB temperature*/
  sigma_8_present ,//  Last updated value of sigma_8 (Presently WMAP)
  vnn; // Spectral index of primordial Power spectrum

extern int N1,N2,N3,// box dimension (grid) 
  MM, // Number of particles
  NF, // Fill every NF grid point 
  Nbin; // Number of pins to calculate final P(k) (output)
extern float   LL; // grid spacing in Mpc

extern float pi;


extern fftw_real ***ro;
extern  rfftwnd_plan p_ro; // for FFT
extern   rfftwnd_plan q_ro; // for FFT

extern fftw_real ***nh,***nhs,***ngamma,***ngammas,***rosp,****nxion;

//end of declaration of global variables for output binary file


// arrays for storing data

void Setting_Up_Memory_For_ionz(int Nnion)
{
  int jk;
  nh=allocate_fftw_real_3d(N1,N2,N3+2);
  nhs=allocate_fftw_real_3d(N1,N2,N3+2);
  ngamma=allocate_fftw_real_3d(N1,N2,N3+2);
  ngammas=allocate_fftw_real_3d(N1,N2,N3+2);
  rosp = allocate_fftw_real_3d(N1,N2,N3+2);
  nxion=(fftw_real****)malloc(sizeof(fftw_real***)*Nnion);
    
  for(jk=0;jk<Nnion;++jk)
    {
      nxion[jk]=allocate_fftw_real_3d(N1,N2,N3+2);
    }
  // allocate area for storing densities  DONE
    
  /* The last dimension gets padded because of using REAL FFT */
    
  /****************************************************/	
  /* Creating the plans for forward and reverse FFT's */
  
  p_ro = rfftw3d_create_plan(N1,N2,N3,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);  
  q_ro = rfftw3d_create_plan(N1,N2,N3,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE | FFTW_IN_PLACE);
    
  // done making plan  
  
}


void smooth(fftw_real ***ro_dum,float Radii)
{
  int i,j,k,index,x1,y1,z1,x2,y2,z2,a,b,c;
  float m,tempre,tempim,tot;
  fftw_complex *A;
  fftw_complex *B;
  
  //generating the filtering function
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
  	rosp[i][j][k]=0.0;
 
  //Radii is radius of the sphere in grid unit
  //generating a sphere at the centre of the box
   
  tot=0.;
  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<N3;k++)
  	{
	  if((float)((N1/2-i)*(N1/2-i)+(N2/2-j)*(N2/2-j)+(N3/2-k)*(N3/2-k))<=Radii*Radii)
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
  

  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<=N3/2;k++)
	{ 
	  /* a=i; */
	  /* if(i>N1/2) a=N1-i; */
	  
	  /* b=j; */
	  /* if(j>N2/2) b=N2-j; */
	  
	  /* c=k; */
	  /* if(k>N3/2) c=N3-k; */
	  
	  index = i*N2*(N3/2 +1) + j*(N3/2 +1) + k;
	  
	  
	  tempre=(A[index].re*B[index].re-A[index].im*B[index].im)*powf((-1.),1.*(i+j+k))/tot;
	  tempim=(A[index].im*B[index].re+A[index].re*B[index].im)*powf((-1.),1.*(i+j+k))/tot;
	  //multiplying the factor powf((-1.),(i+j+k)) with FT of the sphere to shift it to box centre from one corner of the box after complex to real FT
	  
	  A[index].re=tempre;
	  A[index].im=tempim;

	  
	 
	}

  rfftwnd_one_complex_to_real(q_ro,(fftw_complex *) &ro_dum[0][0][0], NULL);

  for(i=0;i<N1;i++)
    for(j=0;j<N2;j++)
      for(k=0;k<=N3;k++)
  	ro_dum[i][j][k]=ro_dum[i][j][k]/(N1*N2*N3);
  //free(rosp);
}

/**************************************************************************
                              FUNCTIONS
**************************************************************************/
