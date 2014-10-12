/**
 * @file   ionz_funcs.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sun Oct 12 00:59:24 2014
 * 
 * @brief  Several core functions
 * 
 * 
 */

#include "ion.h"

/** 
 * Set up memory for nh, ngamma and nxion
 * 
 * @param Nnion Number of nion
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 * @param nh density array
 * @param ngamma source array
 * @param nxion xfrac array
 */
void Setting_Up_Memory_For_ionz(int Nnion, int N1, int N2, int N3, fftw_real ***nh, fftw_real ***ngamma, fftw_real ****nxion) {
  int jk;
  free(nh);
  free(ngamma);
  free(nxion);
  printf("ok checked");
  nh=allocate_fftw_real_3d(N1,N2,N3+2);
  ngamma=allocate_fftw_real_3d(N1,N2,N3+2);
  nxion=(fftw_real****)malloc(sizeof(fftw_real***)*Nnion);
    
  for(jk=0;jk<Nnion;++jk) {
    nxion[jk]=allocate_fftw_real_3d(N1,N2,N3+2);
  }
  // allocate area for storing densities  DONE
  /* The last dimension gets padded because of using REAL FFT */
}

/** 
 * Smooth radius with FFT
 * 
 * @param ro_dum density field
 * @param Radii Smooth radius 
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void smooth(fftw_real ***ro_dum,float Radii,int N1,int N2, int N3) {
  int i,j,k,index;
  float tempre,tempim,tot;
  fftw_real ***rosp;
  fftw_complex *A;
  fftw_complex *B;

  rfftwnd_plan p_ro; // for FFT
  rfftwnd_plan q_ro; // for FFT
  
  rosp = allocate_fftw_real_3d(N1,N2,N3+2);
  /****************************************************/	
  /* Creating the plans for forward and reverse FFT's */
  
  p_ro = rfftw3d_create_plan(N1,N2,N3,FFTW_REAL_TO_COMPLEX,FFTW_ESTIMATE | FFTW_IN_PLACE);  
  q_ro = rfftw3d_create_plan(N1,N2,N3,FFTW_COMPLEX_TO_REAL,FFTW_ESTIMATE | FFTW_IN_PLACE);

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
      for(k=0;k<N3;k++) {
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
      for(k=0;k<=N3/2;k++)    { 
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

  fftw_free(rosp);
  rfftwnd_destroy_plan(p_ro);
  rfftwnd_destroy_plan(q_ro);
  /* A and B are aliases so there is no need to free them... Boyd */
  // fftw_free(A);
  // fftw_free(B);
}
void subgrid_reionization(fftw_real ***nh, fftw_real ***ngamma, fftw_real ****nxion, double robar, float *nion, int Nnion, int N1, int N2, int N3) {
  int ii,jj,kk,jk;
  double vion[Nnion],roion[Nnion];
  
  for(jk=0;jk<Nnion;jk++) {
    //calculating avg. ionization frction
    vion[jk]=0.0;
    roion[jk]=0.0;

    for(kk=0;kk<N3;kk++)
      for(jj=0;jj<N2;jj++)
	for(ii=0;ii<N1;ii++) {
	  if(nh[ii][jj][kk]>nion[jk]*ngamma[ii][jj][kk]) {
	    nxion[jk][ii][jj][kk]=nion[jk]*ngamma[ii][jj][kk]/nh[ii][jj][kk];
	  }    
	  else {
	    nxion[jk][ii][jj][kk]=1.;
	  }
	  vion[jk]+=nxion[jk][ii][jj][kk];
	  roion[jk]+=nxion[jk][ii][jj][kk]*nh[ii][jj][kk];	  
	}
    vion[jk]/=(1.*N1*N2*N3);
    roion[jk]/=(float)(robar*N1*N2*N3);
    if(mympi.ThisTask == 0)
      printf("Subgrid: obtained vol. avg. x_ion=%e mass avg. x_ion=%e\n",vion[jk],roion[jk]);
 }
}


/** 
 * Core function to do semi-numnerical reionization estimation
 * 
 * @param Radii Radius 
 * @param nh_p Pointer to baryon density field
 * @param ngamma_p Pointer to photon density field
 * @param nxion_p Pointer to Xfrac field
 * @param nion_p nion arrays
 * @param Nnion Total element of nion
 * @param N1 1st dimension grid
 * @param N2 2nd dimension grid
 * @param N3 3rd dimension grid
 */
void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3) {
  fftw_real ***nhs,***ngammas;
  int ii,jj,kk,jk;

  nhs=allocate_fftw_real_3d(N1,N2,N3+2);
  ngammas=allocate_fftw_real_3d(N1,N2,N3+2);
  for(ii=0;ii<N1;ii++)
    for(jj=0;jj<N2;jj++)
      for(kk=0;kk<N3;kk++) {
	//Filling smoothing arrays with the dark matter and source density data
	nhs[ii][jj][kk]=nh_p[ii][jj][kk];
	ngammas[ii][jj][kk]=ngamma_p[ii][jj][kk];	     
      }
      
  // printf("starting smoothing for radius of size %e (in units of grid size)\n",Radii);

  //Smoothing with real space spherical filter
  
  smooth(nhs,Radii,N1,N2,N3);
  smooth(ngammas,Radii,N1,N2,N3); 

  for(jk=0;jk<Nnion;++jk) {	 
    for(ii=0;ii<N1;ii++)
      for(jj=0;jj<N2;jj++)
	for(kk=0;kk<N3;kk++) {
	  //Checking the ionization condition
	  if(nhs[ii][jj][kk]<nion_p[jk]*ngammas[ii][jj][kk]) {
	    nxion_p[jk][ii][jj][kk]=1.;
	  }
	}
  }
  fftw_free(nhs);
  fftw_free(ngammas);
}



