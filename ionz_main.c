#include"srfftw.h"
#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include"ion.h"

#ifdef PARALLEL
#include <mpi.h>
#include "ionz_mpi.h"
#endif

#include<string.h>

#define  min(x,y)  ((x)<(y) ?(x):(y))
#define  max(x,y)  ((x)>(y) ?(x):(y))
 /*  GLOBAL VARIABLES  */



// cosmological parameters read from input file  "input.ionz"

float  vhh, // Hubble parameter 
  vomegam, // Omega_matter) 
  vomegalam, // Cosmological Constant 
  vomegab, //Omega_baryon */
  sigma_8_present ,//  Last updated value of sigma_8 (Presently WMAP)
  vnn; // Spectral index of primordial Power spectrum

int N1,N2,N3;// box dimension (grid) 

float   LL; // grid spacing in Mpc

float  pi;

// arrays for storing data

fftw_real ***ro; // for density/potential
rfftwnd_plan p_ro; // for FFT
rfftwnd_plan q_ro; // for FFT

//end of declaration of global variables for output binary file

fftw_real ***nh,***nhs,***ngamma,***ngammas,***rosp,****nxion, ****nxion_buffer;


  //Reading the Nbody dark matter density field 
  //Density are in c^2ray simulation unit (raw)
  //Density files are in binary format with three integers written in
  //the begining of the file, which represents the grid size N1xN2xN3
  //After that the density array of size N1xN2xN3 is written in row-major order (C order), each value is of size float (float32 in python) 

 
  //To avoid confusion we take input the density and the source data in same units


void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int n1, int n2, int n3 )
{
  int ii,jj,kk,jk;
     
  for(ii=0;ii<n1;ii++)
    for(jj=0;jj<n2;jj++)
      for(kk=0;kk<n3;kk++)
	{
	  //Filling smoothing arrays with the dark matter and source density data
	  nhs[ii][jj][kk]=nh_p[ii][jj][kk];
	  ngammas[ii][jj][kk]=ngamma_p[ii][jj][kk];
	     
	}
      
  // printf("starting smoothing for radius of size %e (in units of grid size)\n",Radii);

  //Smoothing with real space spherical filter

  smooth(nhs,Radii);
  smooth(ngammas,Radii); 

  for(jk=0;jk<Nnion;++jk)
    {
	 
      for(ii=0;ii<n1;ii++)
	for(jj=0;jj<n2;jj++)
	  for(kk=0;kk<n3;kk++)
	    {
	      //Checking the ionization condition
	      if(nhs[ii][jj][kk]<nion_p[jk]*ngammas[ii][jj][kk])
		{
		  nxion_p[jk][ii][jj][kk]=1.;
		}
	    }
    }

}
int make_radii_list(float *radii_p, float r_min, float r_max)
{
  float r,dr;
  int i = 0;
  r = r_min;
  while (r < r_max)
    {
      radii_p[i] = r;
      dr = (r*.1) < 2.0 ? (r*.1):2.0;
      r += dr;
      i++;
    }
  radii_p = realloc(radii_p,sizeof(float)*i);
  return i;
}
void pack_4d_array_mpi_transfer(fftw_real ****input, float *output, int n_nion, int n1, int n2, int n3)
{ 
  int ii,jj,kk,jk;
  for(jk=0;jk<n_nion;jk++)
    for(ii=0;ii<n1;ii++)
      for(jj=0;jj<n2;jj++)
	for(kk=0;kk<n3;kk++)
	  output[jk*n1*n2*n3 + ii*n2*n3 + jj*n3 + kk] = input[jk][ii][jj][kk];
}
void unpack_4d_array_mpi_transfer(float *input, fftw_real ****output, int n_nion,int n1, int n2, int n3)
{
  int ii,jj,kk,jk;
  for(jk=0;jk<n_nion;jk++)
    for(ii=0;ii<n1;ii++)
      for(jj=0;jj<n2;jj++)
	for(kk=0;kk<n3;kk++)
	  output[jk][ii][jj][kk]=input[jk*n1*n2*n3 + ii*n2*n3 + jj*n3 + kk];
}
void read_density(char filename[2048],int *N1_p, int *N2_p, int *N3_p, fftw_real ***nh_p, double *robar_p)
{  
  int ii,jj,kk;
  FILE *inp;
  // printf("start read_density\n");
  inp=fopen(filename,"r");
  *robar_p=0.;
  fread(N1_p,sizeof(int),1,inp);
  fread(N2_p,sizeof(int),1,inp);
  fread(N3_p,sizeof(int),1,inp);
  printf("N1=%d\n",*N1_p);
  for(ii=0;ii<*N1_p;ii++)
    for(jj=0;jj<*N2_p;jj++)
      for(kk=0;kk<*N3_p;kk++)
	{
	  // printf("reading %d %d %d\n", ii,jj,kk);
	  fread(&nh_p[ii][jj][kk],sizeof(float),1,inp);
	  
	  *robar_p += nh_p[ii][jj][kk];
	}
  fclose(inp);
  *robar_p /= (1.*(*N1_p)*(*N2_p)*(*N3_p));
  // printf("ok with density\n");
}


void read_sources(char *filename, int N1, int N2, int N3, fftw_real ***ngamma_p, double *robarhalo_p)
{
  FILE *inp;
  int ii,jj,kk,ll;
  int nhalo;
  float mass1,mass2;
  mass1 = 0.;
  mass2 = 0.;
  /* Clear out the array before reading source density ******/
  for (ii=0;ii<N1;ii++)
    for (jj=0;jj<N2;jj++)
      for (kk=0;kk<N3;kk++)
	ngamma_p[ii][jj][kk] = 0.0;
  

  /********************************/
  // Reading the source density data
  // This is also in the same unit as the dark matter density field (in C^2-Ray simulation unit)
  // This file is in ascii
  // Total no. of filled grid points with sources are written in the beginning of the file
  // Total 5 or 6 columns of data is written after that
  // First 3 columns are three FORTRAN indices of the array (need to subtract 1 from each to convert them into C indices)
  // Next 2 columns are low mass and high mass source contribution to that grid
  // There could be another column in the file which we don't need
  
  *robarhalo_p=0.;
  inp=fopen(filename,"r");
  
  fscanf(inp,"%d",&nhalo);
  //Total number of filled grid points with sources
  printf("nhalo =%d\n",nhalo);
  
  
  for(ll=0;ll<nhalo;ll++)
    {
      //If there are 6 columns in the file 

      //fscanf(inp,"%d%d%d%f%f%f",&ii,&jj,&kk,&mass1,&mass2,&dump);   
      //If there are 4 columns in the file 
      fscanf(inp,"%d%d%d%f",&ii,&jj,&kk,&mass1);
   
      //You can treat both mass ranges similarly
      //Or one can use different weights or functional response for each mass range
      //We show the simplest case here, treating both of them similarly
      ngamma_p[ii-1][jj-1][kk-1] = mass1+mass2;     
      *robarhalo_p += ngamma_p[ii-1][jj-1][kk-1]; 
    }
    
  fclose(inp);
  //Avg. source density
  *robarhalo_p/=(1.*N1*N2*N3);
}
main(int argc, char **argv)
{
  FILE  *inp,*outpp;
  int ii, jj, kk,ll,jk,mm,sfac;
  float  vaa,epsilon,box;
  int  Nnion,temp,flag_sup;
  float dr,r_min,r_max;
  char file1[300],file2[300],num[50];
  float *nion,xh1,dump,mass1,mass2;  
  int nhalo,output_flag,in_flag;
  float Radii,his_z,zval;
  double robar,robarhalo,vfac,*vion,*roion;
  float **rra,**vva,**halo,**data,*dummy,junk1=1.0,junk2=0.0;
  float *Radii_list;
  int n_radii;
  int *NjobsperTask;
  int *JobsTask;
  float *buffer, *buffer_final;
#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
  NTask = 1;
  ThisTask = 0;
#endif //PARALLEL

  pi=4.0*atan(1.0);
  
  system("date");
  //Reading the input simulation parameter file
  inp=fopen("input.ionz","r");
  /* get parameters for reion simulation */
  fscanf(inp,"%d",&Nnion);
  //Nnion is the number of nion values for which we will run the simulation
  //You will get a x_HI map for each value of nion
  //Allocating memory for some variables
  nion=(float*)calloc(Nnion,sizeof(float));
  vion=(double*)calloc(Nnion,sizeof(double));
  roion=(double*)calloc(Nnion,sizeof(double));
  //Reading nion values from the input file
  //nion is the efficiency paramter for ionization (a combination of f_esc,f_star etc)
  // We try with variuos values of nion
  for(ii=0;ii<Nnion;ii++)
    {
      fscanf(inp,"%f",&nion[ii]);
    }
  
  fscanf(inp,"%f%f%f%f%f",&vaa,&vomegam,&vomegalam,&vhh,&vomegab);
  //Reading cosmological parameters values
  //vaa initially scans the redshift of the Nbody simulation
  zval = vaa; //value of z stored
  vaa = 1/(1+vaa);// converts it into scale factor
  fscanf(inp,"%d%d%d%f",&N1,&N2,&N3,&box);
    
  fclose(inp);
 
  LL = box/(vhh*N1); // grid size in Mpc

  printf("scale factor= %f\n",vaa);

  // Allocating memory to different arrays
  Setting_Up_Memory_For_ionz(Nnion);
  
  read_density("/research/prace/sph_smooth_cubepm_130315_6_1728_47Mpc_ext2/nc306/7.859n_all.dat",&N1,&N2,&N3,nh,&robar);
  read_sources("/research/prace/47Mpc_RT/47Mpc_f2_gs_306/sources/7.859-coarsest_sources_used_wfgamma.dat",N1,N2,N3,ngamma,&robarhalo);
  

  //calculating max and min radius for smoothing in units of grid size
  r_min=1.;
  r_max=pow((1.*N1*N2*N3),(1./3.))/2.;
  Radii_list = malloc(sizeof(float)*1000); // max is 1000
  NjobsperTask = malloc(sizeof(float)*NTask);
  n_radii = make_radii_list(Radii_list,r_min,r_max);
  for(jj=0;jj<NTask;jj++)
    {
      NjobsperTask[jj] = n_radii/NTask;
      if(jj < n_radii%NTask)
	NjobsperTask[jj]++;
      if(jj == ThisTask)
	{
	  JobsTask = malloc(sizeof(int)*NjobsperTask[jj]);
	  for(ii=0;ii<NjobsperTask[jj];ii++)
	    {
	      JobsTask[ii] = ii*NTask+ThisTask;
	    }
	}
    }
  buffer = malloc(sizeof(float)*Nnion*N1*N2*N3);
  if(ThisTask == 0)
    buffer_final = malloc(sizeof(float)*Nnion*N1*N2*N3);
  // The max smoothing radius here is set as half of the diagonal of the box
  // This can be changed, one can choose a redshift dependent function instead
  // Or one can choose a model for redshift evolution of the mean free path of the UV photons
  // We are showing the most simple case here
  
  // printf("robar=%e  robarhalo=%e ratio= %e\n",robar,robarhalo,robar/robarhalo);

  //subgrid re-ionization

  for(jk=0;jk<Nnion;++jk)
    {
      
      for(ii=0;ii<N1;ii++)
  	for(jj=0;jj<N2;jj++)
  	  for(kk=0;kk<N3;kk++)
  	    {
  	      if(nh[ii][jj][kk]>nion[jk]*ngamma[ii][jj][kk])
  		{
  		  nxion[jk][ii][jj][kk]=nion[jk]*ngamma[ii][jj][kk]/nh[ii][jj][kk];
  		}
	      
  	      else
  		{
		  nxion[jk][ii][jj][kk]=1.;
  		}
  	    }
    }

  for(jk=0;jk<Nnion;++jk)
    {
      
      //calculating avg. ionization frction
      vion[jk]=0.0;
      roion[jk]=0.0;
      
      for(ii=0;ii<N1;ii++)
  	for(jj=0;jj<N2;jj++)
  	  for(kk=0;kk<N3;kk++)
  	    {
  	      vion[jk]+=nxion[jk][ii][jj][kk];
  	      roion[jk]+=nxion[jk][ii][jj][kk]*nh[ii][jj][kk];
  	    }
      vion[jk]/=(1.*N1*N2*N3);
      roion[jk]/=(float)(robar*N1*N2*N3);
      // printf("obtained vol. avg. x_ion=%e mass avg. x_ion=%e\n",vion[jk],roion[jk]);
    }


  //subgrid re-ionization  done

 
  //smoothing loop 
  //done for a range of length scales from r_min to r_max (in units of grid size)
  
  // system("date");
  /* smoothing */
  printf("Task: %d Njobs %d\n",ThisTask,NjobsperTask[ThisTask]);
  for(ii=0;ii<NjobsperTask[ThisTask];ii++)
    {
      printf("Task: %d do the job %d\n",ThisTask,JobsTask[ii]);
      reionization(Radii_list[JobsTask[ii]], nh, ngamma, nxion, nion, Nnion, N1, N2, N3 );    
    }  // system("date");
  

  pack_4d_array_mpi_transfer(nxion,buffer,Nnion, N1, N2, N3);
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Reduce(buffer,buffer_final,Nnion*N1*N2*N3,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);
  if(ThisTask == 0)
    printf("finish finding max\n");
  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      unpack_4d_array_mpi_transfer(buffer_final, nxion, Nnion, N1, N2, N3);
      for(jk=0;jk<Nnion;++jk)
	{
	  //calculating avg. ionization frction
	  vion[jk]=0.0;
	  roion[jk]=0.0;
      
	  // Defining the ionization map output file name
	  // This is based on the value of nion assigned to it

	  strcpy(file2,"xHI_map_");
	  sprintf(num,"%4.2f",nion[jk]);
	  strcat(file2,num);

	  // Writing the x_HI map in binary
	  // In the begining 3 integers are written which defines the size
	  // of the x_HI array

	  inp=fopen(file2,"w");
	  fwrite(&N1,sizeof(int),1,inp);
	  fwrite(&N2,sizeof(int),1,inp);
	  fwrite(&N3,sizeof(int),1,inp);
	  for(ii=0;ii<N1;ii++)
	    for(jj=0;jj<N2;jj++)
	      for(kk=0;kk<N3;kk++)
		{
		  xh1=(1.-nxion[jk][ii][jj][kk]);
		  xh1=(xh1 >0.0)? xh1: 0.0;
		  vion[jk]+=xh1;
		  nxion[jk][ii][jj][kk]=xh1; // store x_HI instead of x_ion
		  nhs[ii][jj][kk]=xh1*nh[ii][jj][kk]; // ro_HI on grid
		  roion[jk]+=nhs[ii][jj][kk];
		  fwrite(&nxion[jk][ii][jj][kk],sizeof(float),1,inp);
	      
		}
	    
	  fclose(inp);
	  roion[jk]/=(N1*N2*N3); // mean HI density in grid units
      
	  vion[jk]/=(1.*N1*N2*N3); // volume avg xHI
	  roion[jk]/=(float)robar; // divide by H density to get mass avg. xHI
	  printf("nion = %f obtained vol. avg. x_HI=%e mass avg. x_HI=%e\n",nion[jk],vion[jk],roion[jk]);
	  system("date");
	}
    }
 system("date");            
} /* END OF MAIN */

