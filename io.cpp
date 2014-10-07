#include <iostream>
class ionz_io
{
  void read_dens(char* Filename, int Format, int *N1, int *N2, int *N3, float ***dens)
  {    
    // 1: CUBEP3M
    if(Format == 1)
      read_dens_CUBEP3M(Filename, dens);
  }
  void read_sources(char* File, int Format, int N1, int N2, int N3, float ***sources)
  {
  }
  void read_dens_CUBEP3M(char* Filename, int *N1, int *N2, int *N3, float ***dens)
  {
    FILE *fp;
    fp = fopen(Filename, "rb");
    fclose(fp);
  }
}
