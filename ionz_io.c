/**
 * @file   ionz_io.c
 * @author Chaichalit Srisawat < boyd.srisawat@gmail.com>
 * @date   Sat Oct 11 21:44:53 2014
 * 
 * @brief  IO operations
 * 
 * 
 */

#include "ion.h"

/** 
 * Convert 3D-fftw_real array to 1D-float array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 3D-fftw_real array input
 * @param output 1D-float array output
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void pack_3d_array_mpi_transfer(fftw_real ***input, float *output, int N1, int N2, int N3) { 
  int ii,jj,kk;
  for(kk=0;kk<N3;kk++)
    for(jj=0;jj<N2;jj++)
      for(ii=0;ii<N1;ii++)
	output[kk*N2*N1 + jj*N1 + ii] = input[ii][jj][kk];
}

/** 
 * Convert 1D-float array to 3D-fftw_real array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 1D-float array input
 * @param output 3D-fftw_real array output
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void unpack_3d_array_mpi_transfer(float *input, fftw_real ***output, int N1, int N2, int N3) {
  int ii,jj,kk;
  for(kk=0;kk<N3;kk++)
    for(jj=0;jj<N2;jj++)
      for(ii=0;ii<N1;ii++)
	output[ii][jj][kk]=input[kk*N2*N1 + jj*N1 + ii];
}

/** 
 * Convert 3D-fftw_real array to 1D-float array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 3D-fftw_real array output
 * @param output 1D-float array input
 * @param Nion Number of the choices of total reionizing photons/neutral baryon
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void pack_4d_array_mpi_transfer(fftw_real ****input, float *output, int Nnion, int N1, int N2, int N3) { 
  int ii,jj,kk,jk;
  for(jk=0;jk<Nnion;jk++)
    for(kk=0;kk<N3;kk++)
      for(jj=0;jj<N2;jj++)
	for(ii=0;ii<N1;ii++)
	  output[jk*N1*N2*N3 + ii*N1*N2 + jj*N1 + ii] = input[jk][ii][jj][kk];
}
/** 
 * Convert 1D-float array to 3D-fftw_real array. This is used
 * to simplify MPI transfer process. 
 *
 * @param input 1D-float array input
 * @param output 3D-fftw_real array output
 * @param Nion Number of the choices of total reionizing photons/neutral baryon
 * @param N1 1st dimension grid number
 * @param N2 2nd dimension grid number
 * @param N3 3rd dimension grid number
 */
void unpack_4d_array_mpi_transfer(float *input, fftw_real ****output, int Nnion,int N1, int N2, int N3) {
  int ii,jj,kk,jk;
  for(jk=0;jk<Nnion;jk++)
    for(kk=0;kk<N3;kk++)
      for(jj=0;jj<N2;jj++)
	for(ii=0;ii<N1;ii++)
	  output[jk][ii][jj][kk]=input[jk*N1*N2*N3 + kk*N1*N2 + jj*N1 + ii];
}


/** 
 * Read in density of matter in cubep3m format (Fortran binary)
 *  
 *
 * @param filename Filename (char* input) 
 * @param buffer_3d Output buffer[N1*N2*N3]
 * @param robar_p Average density of baryons (output)
 * @param N1 1st dimension grid (input)
 * @param N2 2nd dimension grid (input)
 * @param N3 3rd dimension grid (input)
 * @param vomegam Omega_matter (input)
 * @param vomegab Omega_baryon (input)
 */
void read_density(char *filename, float *buffer_3d, double *robar_p, int N1, int N2, int N3, float vomegam, float vomegab) {  
  int ii;
  int n1,n2,n3;
  FILE *inp;
  // printf("start read_density\n");
  inp=fopen(filename,"rb");
  *robar_p=0.;
  fread(&n1,sizeof(int),1,inp);
  fread(&n2,sizeof(int),1,inp);
  fread(&n3,sizeof(int),1,inp);
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

/** 
 * 
 * Read in the history of ionized fraction (previous snapshot).
 * 
 * @param dirname Directory of the xh (input) 
 * @param buffer_4d 4D float output
 * @param nion_list List of minimum N_gamma/N_b (input) 
 * @param Nnion Number of elements in nion_list
 * @param N1 1st dimension grid (input)
 * @param N2 2nd dimension grid (input)
 * @param N3 3rd dimension grid (input)
 */
void read_xfrac(char *dirname, float *buffer_4d, float *nion_list, int Nnion, int N1, int N2, int N3) {
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


/** 
 * Read sources list. The data will be the SFR in cubep3m_gridmass/year
 * 
 * @param filename Filename (char* input) 
 * @param buffer_3d Output buffer[N1*N2*N3]
 * @param robar_p Average density of photons (output)
 * @param N1 1st dimension grid (input)
 * @param N2 2nd dimension grid (input)
 * @param N3 3rd dimension grid (input)
 */
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
    *robarhalo_p += buffer_3d[ii];
  }
  *robarhalo_p /= (1.*(n1)*(n2)*(n3));
}

