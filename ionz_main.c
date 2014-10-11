#include "ion.h"

#define  min(x,y)  ((x)<(y) ?(x):(y))
#define  max(x,y)  ((x)>(y) ?(x):(y))



void reionization(float Radii,fftw_real ***nh_p, fftw_real ***ngamma_p, fftw_real ****nxion_p, float *nion_p, int Nnion, int N1, int N2, int N3) {
  fftw_real ***nhs,***ngammas;
  int ii,jj,kk,jk;
  double t_start,t_stop;   
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
  
  smooth(nhs,Radii);
  smooth(ngammas,Radii); 

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

void pack_3d_array_mpi_transfer(fftw_real ***input, float *output, int n1, int n2, int n3) { 
  int ii,jj,kk;
  for(kk=0;kk<n3;kk++)
    for(jj=0;jj<n2;jj++)
      for(ii=0;ii<n1;ii++)
	output[kk*n2*n1 + jj*n1 + ii] = input[ii][jj][kk];
}

void unpack_3d_array_mpi_transfer(float *input, fftw_real ***output, int n1, int n2, int n3) {
  int ii,jj,kk;
  for(kk=0;kk<n3;kk++)
    for(jj=0;jj<n2;jj++)
      for(ii=0;ii<n1;ii++)
	output[ii][jj][kk]=input[kk*n2*n1 + jj*n1 + ii];
}
void pack_4d_array_mpi_transfer(fftw_real ****input, float *output, int n_nion, int n1, int n2, int n3) { 
  int ii,jj,kk,jk;
  for(jk=0;jk<n_nion;jk++)
    for(kk=0;kk<n3;kk++)
      for(jj=0;jj<n2;jj++)
	for(ii=0;ii<n1;ii++)
	  output[jk*n1*n2*n3 + ii*n1*n2 + jj*n1 + ii] = input[jk][ii][jj][kk];
}
void unpack_4d_array_mpi_transfer(float *input, fftw_real ****output, int n_nion,int n1, int n2, int n3) {
  int ii,jj,kk,jk;
  for(jk=0;jk<n_nion;jk++)
    for(kk=0;kk<n3;kk++)
      for(jj=0;jj<n2;jj++)
	for(ii=0;ii<n1;ii++)
	  output[jk][ii][jj][kk]=input[jk*n1*n2*n3 + kk*n1*n2 + jj*n1 + ii];
}

/* Read density in cubep3m format (Fortran binary) */
/* The output will be the mass of baryons in cubep3m grid mass */


void read_density(char *filename, float *buffer_3d, double *robar_p, int N1, int N2, int N3, float vomegam, float vomegab) {  
  int ii;
  int n1,n2,n2;
  FILE *inp;
  // printf("start read_density\n");
  inp=fopen(filename,"rb");
  *robar_p=0.;
  fread(n1,sizeof(int),1,inp);
  fread(n2,sizeof(int),1,inp);
  fread(n3,sizeof(int),1,inp);
  if(n1 != N1 || n2 != N2 || n3 != N3) {
    printf("Grid dimensions in %s are not the same as in config file\n",filename);
    printf("Density file: %d:%d:%d   Config file: %d:%d:%d\n",n1,n2,n3,N1,N2,N3);
    printf("Terminate\n");
    exit(1);
  }
  fread(buffer_3d,sizeof(float),n1*n2*n3,inp);
  fclose(inp);
  
  for(ii=0;ii<n1*n2*n3;ii++) {
    buffer_3d[ii] *= vomegab/vomegam;
    *robar_p += buffer_3d[ii];
  }
  *robar_p /= (1.*(n1)*(n2)*(n3));
}

#ifdef XH_HISTORY
void read_xfrac(char *filename, buffer *buffer_4d, int Nnion, int N1, int N2, int N3) {
  FILE *inp;
  int ii,jk;
  int n1,n2,n3;
  inp=fopen(filename,"rb");
  fread(&n1,sizeof(int),1,inp);
  fread(&n2,sizeof(int),1,inp);
  fread(&n3,sizeof(int),1,inp);
  if(n1 != N1 || n2 != N2 || n3 != N3) {
    printf("Grid dimensions in %s are not the same as in config file\n",filename);
    printf("Xfrac file: %d:%d:%d   Config file: %d:%d:%d\n",n1,n2,n3,N1,N2,N3);
    printf("Terminate\n");
    exit(1);
  }
  fread(buffer_4d,sizeof(float),n1*n2*n3,inp);
  fclose(inp);
}
#endif

/* Read source in C2Ray like */
/* The output will be the SFR in cubep3m grid mass/year */
void read_sources(char *filename, float *buffer_3d, double *robarhalo_p, int N1, int N2, int N3) {
  FILE *inp;
  int ii,jj,kk,ll;
  int n1,n2,n3;
  float dt = 11.6e6;
  
  inp=fopen(filename,"rb");
  *robarhalo_p=0.;
  fread(&n1,sizeof(int),1,inp);
  fread(&n2,sizeof(int),1,inp);
  fread(&n3,sizeof(int),1,inp);
  if(n1 != N1 || n2 != N2 || n3 != N3) {
    printf("Grid dimensions in %s are not the same as in config file\n",filename);
    printf("Source file: %d:%d:%d   Config file: %d:%d:%d\n",n1,n2,n3,N1,N2,N3);
    printf("Terminate\n");
    exit(1);
  }
  fread(buffer_3d,sizeof(float),n1*n2*n3,inp);
  fclose(inp);
  for(ii=0;ii<n1*n2*n3;ii++) {
    buffer_3d[ii] *= dt;
    *robarhalo_p += buffer_3d;
  }
  *robarhalo_p /= (1.*(n1)*(n2)*(n3));
}


main(int argc, char **argv) {
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
  if(ThisTask == 0) {
    system("date");
    printf("Start semi-numerical reionization process\n");
  }
  pi=4.0*atan(1.0);
  read_params("input.ionz");
  Nnion = input_param.Nnion;
  nion=(float*)calloc(Nnion,sizeof(float));
  for(ii=0;ii<Nnion;ii++) {
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
  t_start =Get_Current_time();
  if(ThisTask == 0) {
    read_density(densfilename,&N1,&N2,&N3,nh,&robar);
    buffer = malloc(sizeof(float)*N1*N2*N3);
    pack_3d_array_mpi_transfer(nh,buffer,N1,N2,N3);
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robar, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N1, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N2, 1, MPI_INT, 0, MPI_COMM_WORLD);
  MPI_Bcast(&N3, 1, MPI_INT, 0, MPI_COMM_WORLD);
  if(ThisTask > 0) {
    buffer = malloc(sizeof(float)*N1*N2*N3);
  }
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
  if(ThisTask > 0) {
#endif
    unpack_3d_array_mpi_transfer(buffer,nh,N1,N2,N3);
#ifdef PARALLEL
  }
#endif
  free(buffer);
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  if(ThisTask == 0)  {
    read_sources(sourcefilename,N1,N2,N3,ngamma,&robarhalo);  
    buffer = malloc(sizeof(float)*N1*N2*N3);
    pack_3d_array_mpi_transfer(ngamma,buffer,N1,N2,N3);
  }
  if(ThisTask > 0) {
    buffer = malloc(sizeof(float)*N1*N2*N3);
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&robarhalo, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(buffer, N1*N2*N3, MPI_FLOAT, 0, MPI_COMM_WORLD);
#endif
  if(ThisTask > 0)    {
    unpack_3d_array_mpi_transfer(buffer,ngamma,N1,N2,N3);
  }

  free(buffer);
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
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
  for(jj=0;jj<NTask;jj++) {
    NjobsperTask[jj] = n_radii/NTask;
    if(jj < n_radii%NTask)
      NjobsperTask[jj]++;
    if(jj == ThisTask) {
      JobsTask = malloc(sizeof(int)*NjobsperTask[jj]);
      for(ii=0;ii<NjobsperTask[jj];ii++) {
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
    if(ThisTask == 0)
      printf("Subgrid: obtained vol. avg. x_ion=%e mass avg. x_ion=%e\n",vion[jk],roion[jk]);
 }

  if(ThisTask == 0)
    printf("Start semi-numerical reionization process\n");



#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t_start = Get_Current_time();
  for(ii=0;ii<NjobsperTask[ThisTask];ii++) {
    reionization(Radii_list[JobsTask[ii]], nh, ngamma, nxion, nion, Nnion, N1, N2, N3 );    
  }  
  fftw_free(ngamma);
  fftw_free(ngammas);
  fftw_free(nhs);
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  buffer = malloc(sizeof(float)*Nnion*N1*N2*N3);
  t_stop = Get_Current_time();
  if(ThisTask == 0)
    printf("Finish reionizing process %lf s\n",t_stop-t_start);
  t_start = Get_Current_time();

#ifdef PARALLEL
  pack_4d_array_mpi_transfer(nxion,buffer,Nnion, N1, N2, N3);
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  t_stop = Get_Current_time();
  if(ThisTask == 0) {
    printf("Finish packing data %lf s\n",t_stop-t_start);
#ifdef PARALLEL
    buffer_final = malloc(sizeof(float)*Nnion*N1*N2*N3);
#endif
  }
#ifdef PARALLEL
  MPI_Barrier(MPI_COMM_WORLD);
#ifdef CHUNKTRANSFER
  t_start = Get_Current_time();
  ii = 0;
  while (ii*mpi_buffer < Nnion*N1*N2*N3) {
    cur_len = min(mpi_buffer, Nnion*N1*N2*N3-ii*mpi_buffer);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(&buffer[ii*mpi_buffer],&buffer_final[ii*mpi_buffer],cur_len,MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
    ii++;
  }

  MPI_Barrier(MPI_COMM_WORLD);

  t_stop = Get_Current_time();
  if(ThisTask == 0)
    printf("Finish finding max:split %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#else
  t_start = Get_Current_time();
  MPI_Reduce(buffer, buffer_final, Nnion*N1*N2*N3, MPI_FLOAT,MPI_MAX,0,MPI_COMM_WORLD);      
  MPI_Barrier(MPI_COMM_WORLD);
  t_stop = Get_Current_time();
  if(ThisTask == 0)
    printf("Finish finding max:whole %lf s\n",t_stop-t_start); 
  MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif // PARALLEL

  if(ThisTask == 0) {
    for(jk=0;jk<Nnion;jk++) {
      //calculating avg. ionization frction
      vion[jk]=0.0;
      roion[jk]=0.0;
      // Defining the ionization map output file name
      // This is based on the value of nion assigned to it
      sprintf(file2,"xHI_map_%s_%4.2f",z_out,nion[jk]);
      printf("Saving %s\n",file2);
      // Writing the x_HI map in binary
      // In the begining 3 integers are written which defines the size
      // of the x_HI array
#ifdef PARALLEL
      
#endif

      inp=fopen(file2,"w");
      fwrite(&N1,sizeof(int),1,inp);
      fwrite(&N2,sizeof(int),1,inp);
      fwrite(&N3,sizeof(int),1,inp);
      
      
#ifdef PARALLEL
      fwrite(&buffer_final[jk*N1*N2*N3],sizeof(float),N1*N2*N3,inp);
#else
      fwrite(&buffer[jk*N1*N2*N3],sizeof(float),N1*N2*N3,inp);	
#endif
      for(ii=0;ii<N1*N2*N3;ii++) {
	    xh1=(1.-nxion[jk][ii][jj][kk]);
	    xh1=(xh1 >0.0)? xh1: 0.0;
	    vion[jk]+=xh1;
	    nxion[jk][ii][jj][kk]=xh1; // store x_HI instead of x_ion
	    nhs[ii][jj][kk]=xh1*nh[ii][jj][kk]; // ro_HI on grid
	    roion[jk]+=nhs[ii][jj][kk];	
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

