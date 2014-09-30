/* in allotarrays.c */
float **allocate_float_2d(long N1,int N2);
fftw_real  ***allocate_fftw_real_3d(int N1,int N2,int N3);
void Setting_Up_Memory_For_ionz(int Nnion);

void smooth(fftw_real ***ro_dum,float Radii);
void smooth_k_sharp(fftw_real ***ro_dum,float Radii);              












