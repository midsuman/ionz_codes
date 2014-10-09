#include "ion.h"

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
int Nnion;

float   LL; // grid spacing in Mpc

float  pi;

// arrays for storing data

fftw_real ***ro; // for density/potential
rfftwnd_plan p_ro; // for FFT
rfftwnd_plan q_ro; // for FFT

//end of declaration of global variables for output binary file

fftw_real ***nh,***nhs,***ngamma,***ngammas,***rosp,****nxion;


//Reading the Nbody dark matter density field 
//Density are in c^2ray simulation unit (raw)
//Density files are in binary format with three integers written in
//the begining of the file, which represents the grid size N1xN2xN3
//After that the density array of size N1xN2xN3 is written in row-major order (C order), each value is of size float (float32 in python) 

 
//To avoid confusion we take input the density and the source data in same units


void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int n1, int n2, int n3 )
{
  int ii,jj,kk,jk;
  double t_start,t_stop;   
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
void pack_3d_array_mpi_transfer(fftw_real ***input, float *output, int n1, int n2, int n3)
{ 
  int ii,jj,kk;
  for(ii=0;ii<n1;ii++)
    for(jj=0;jj<n2;jj++)
      for(kk=0;kk<n3;kk++)
	output[ii*n2*n3 + jj*n3 + kk] = input[ii][jj][kk];
}

void unpack_3d_array_mpi_transfer(float *input, fftw_real ***output, int n1, int n2, int n3)
{
  int ii,jj,kk;
  for(ii=0;ii<n1;ii++)
    for(jj=0;jj<n2;jj++)
      for(kk=0;kk<n3;kk++)
	output[ii][jj][kk]=input[ii*n2*n3 + jj*n3 + kk];
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

/* Read density in cubep3m format (Fortran binary) */
/* The output will be the mass of baryons in cubep3m grid mass */
void read_density(char filename[2048],int *N1_p, int *N2_p, int *N3_p, fftw_real ***nh_p, double *robar_p)
{  
  int ii,jj,kk;
  FILE *inp;
  // printf("start read_density\n");
  inp=fopen(filename,"rb");
  *robar_p=0.;
  fread(N1_p,sizeof(int),1,inp);
  fread(N2_p,sizeof(int),1,inp);
  fread(N3_p,sizeof(int),1,inp);
  if(ThisTask == 0)
    printf("N1=%d\n",*N1_p);
  for(kk=0;kk<*N3_p;kk++)
    for(jj=0;jj<*N2_p;jj++)
      for(ii=0;ii<*N1_p;ii++)
	{
	  fread(&nh_p[ii][jj][kk],sizeof(float),1,inp);
	  nh_p[ii][jj][kk] *= vomegab/vomegam;
	  *robar_p += nh_p[ii][jj][kk];
	}
  fclose(inp);
  *robar_p /= (1.*(*N1_p)*(*N2_p)*(*N3_p));
  // printf("ok with density\n");
}

/* Read source in C2Ray like */
/* The output will be the SFR in cubep3m grid mass/year */
void read_sources(char *filename, int N1, int N2, int N3, fftw_real ***ngamma_p, double *robarhalo_p)
{
  FILE *inp;
  int ii,jj,kk,ll;
  int n1,n2,n3;
  float dt = 11.6e6;
  /* Clear out the array before reading source density ******/
  /* for (ii=0;ii<N1;ii++) */
  /*   for (jj=0;jj<N2;jj++) */
  /*     for (kk=0;kk<N3;kk++) */
  /* 	ngamma_p[ii][jj][kk] = 0.0; */
  

  /********************************/
  // Reading the source density data
  // This is also in the same unit as the dark matter density field (in C^2-Ray simulation unit)
  // This file is in ascii
  // Total no. of filled grid points with sources are written in the beginning of the file
  // Total 5 or 6 columns of data is written after that
  // First 3 columns are three FORTRAN indices of the array (need to subtract 1 from each to convert them into C indices)
  // Next 2 columns are low mass and high mass source contribution to that grid
  // There could be another column in the file which we don't need
  
  inp=fopen(filename,"rb");
  *robarhalo_p=0.;
  fread(&n1,sizeof(int),1,inp);
  fread(&n2,sizeof(int),1,inp);
  fread(&n3,sizeof(int),1,inp);
  for(kk=0;kk<n3;kk++)
    for(jj=0;jj<n2;jj++)
      for(ii=0;ii<n3;ii++)
	{
	  fread(&ngamma_p[ii][jj][kk],sizeof(float),1,inp);
	  ngamma_p[ii][jj][kk] *= dt;
	  *robarhalo_p += ngamma_p[ii][jj][kk];
	}
  fclose(inp);
  *robarhalo_p /= (1.*(n1)*(n2)*(n3));
}


main(int argc, char **argv)
{
  FILE  *inp,*outpp;
  int ii, jj, kk,ll,jk,mm,sfac;
  float  Nnion,vaa,epsilon,box;
  int  temp,flag_sup;
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
  double t_start, t_stop;
  float *buffer, *buffer_final;
  int mpi_buffer=1000000;
  int cur_len;
  char densfilename[2000], sourcefilename[2000],z_out[1000];
#ifdef PARALLEL
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &ThisTask);
  MPI_Comm_size(MPI_COMM_WORLD, &NTask);
#else
  NTask = 1;
  ThisTask = 0;
#endif //PARALLEL
  if(argc != 4) {
    printf("Need 3 inputs\n");
    exit(0);
  }
  sprintf(densfilename,"%s",argv[1]);
  sprintf(sourcefilename,"%s",argv[2]);
  sprintf(z_out,"%s",argv[3]);
  if(ThisTask == 0)
    {
      system("date");
      printf("Start semi-numerical reionization process\n");
    }
  pi=4.0*atan(1.0);
  read_params("input.ionz");
  Nnion = input_param.Nnion;
  nion=(float*)calloc(Nnion,sizeof(float));
  for(ii=0;ii<Nnion;ii++)
    {
      nion[ii] = input_param.nion[ii];
    }
  vion=(double*)calloc(Nnion,sizeof(double));
  roion=(double*)calloc(Nnion,sizeof(double));
  N1 = input_param.N1;
  N2 = input_param.N2;
  N3 = input_param.N3;
  vomegam = input_param.omegam;
  vomegalam = input_param.omegalam;
  vomegab = input_param.omegab;
  
  

  // Allocating memory to different arrays
  Setting_Up_Memory_For_ionz(Nnion);
  t_start = MPI_Wtime();
  if(ThisTask == 0)
    {
      read_density(densfilename,&N1,&N2,&N3,nh,&robar);
      buffer = malloc(sizeof(float)*N1*N2*N3);
      pack_3d_array_mpi_transfer(nh,buffer,N1,N2,N3);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N1, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N2, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N3, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(ThisTask > 0)
    {
      buffer = malloc(sizeof(float)*N1*N2*N3);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if(ThisTask > 0)
    {
      unpack_3d_array_mpi_transfer(buffer,nh,N1,N2,N3);
    }
  free(buffer);
  MPI_Barrier(MPI_COMM_WORLD);
  if(ThisTask == 0)
    {
      read_sources(sourcefilename,N1,N2,N3,ngamma,&robarhalo);  
      buffer = malloc(sizeof(float)*N1*N2*N3);
      pack_3d_array_mpi_transfer(ngamma,buffer,N1,N2,N3);
    }
  if(ThisTask > 0)
    {
      buffer = malloc(sizeof(float)*N1*N2*N3);
    }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robarhalo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if(ThisTask > 0)
    {
      unpack_3d_array_mpi_transfer(buffer,ngamma,N1,N2,N3);
    }

  free(buffer);
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = MPI_Wtime();
  /* Sanity MPI check */
  if(ThisTask == 1)
    printf("N1=%d N2=%d N3=%d\n",N1,N2,N3);
  if(ThisTask == 0)
    printf("reading in data %lf s\n",t_stop-t_start);
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

  if(ThisTask == 0)
    printf("Start semi-numerical reionization\n");
  //subgrid re-ionization  done

 
  //smoothing loop 
  //done for a range of length scales from r_min to r_max (in units of grid size)
  
  MPI_Barrier(MPI_COMM_WORLD);
  t_start = MPI_Wtime();
  for(ii=0;ii<NjobsperTask[ThisTask];ii++)
    {
      reionization(Radii_list[JobsTask[ii]], nh, ngamma, nxion, nion, Nnion, N1, N2, N3 );    
    }  
  MPI_Barrier(MPI_COMM_WORLD);
  buffer = malloc(sizeof(float)*Nnion*N1*N2*N3);
  t_stop = MPI_Wtime();
  if(ThisTask == 0)
    printf("Finish reionizing process %lf s\n",t_stop-t_start);
  t_start = MPI_Wtime();
  pack_4d_array_mpi_transfer(nxion,buffer,Nnion, N1, N2, N3);
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = MPI_Wtime();
  if(ThisTask == 0)
    {
      printf("Finish packing data %lf s\n",t_stop-t_start);
      /* allocate buffer_final for final outputs */
      buffer_final = malloc(sizeof(float)*Nnion*N1*N2*N3);
    }
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
  t_start = MPI_Wtime();
  ii = 0;
  while (ii*mpi_buffer < Nnion*N1*N2*N3)
    {
      cur_len = min(mpi_buffer, Nnion*N1*N2*N3-ii*mpi_buffer);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(&buffer[ii*mpi_buffer],&buffer_final[ii*mpi_buffer],cur_len,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
      ii++;
    }

  MPI_Barrier(MPI_COMM_WORLD);

  t_stop = MPI_Wtime();
  if(ThisTask == 0)
    printf("Finish finding max:split %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#else
  t_start = MPI_Wtime();
  MPI_Reduce(buffer, buffer_final, Nnion*N1*N2*N3, MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = MPI_Wtime();
  if(ThisTask == 0)
    printf("Finish finding max:whole %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(ThisTask == 0) {
    t_start = MPI_Wtime();
    unpack_4d_array_mpi_transfer(buffer_final, nxion, Nnion, N1, N2, N3);
    t_stop = MPI_Wtime();
    printf("Finish unpacking data %lf s\n",t_stop-t_start);
    for(jk=0;jk<Nnion;++jk)	{
      //calculating avg. ionization frction
      vion[jk]=0.0;
      roion[jk]=0.0;
      
      // Defining the ionization map output file name
      // This is based on the value of nion assigned to it


      sprintf(file2,"xHI_map_%s_%10.2f",z_out,nion[jk]);
      printf("Saving %s\n",file2);
      // Writing the x_HI map in binary
      // In the begining 3 integers are written which defines the size
      // of the x_HI array

      inp=fopen(file2,"w");
      fwrite(&N1,sizeof(int),1,inp);
      fwrite(&N2,sizeof(int),1,inp);
      fwrite(&N3,sizeof(int),1,inp);
      for(kk=0;kk<N3;kk++)
	for(jj=0;jj<N2;jj++)
	  for(ii=0;ii<N1;ii++) {
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
  MPI_Finalize();

  return 0;
} /* END OF MAIN */

