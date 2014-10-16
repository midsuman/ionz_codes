#include "read_param.h"
class read_param
{
 public:
  int Nnion;
  int *nion;
  float a_expansion;
  float omegam;
  float omegalam;
  float omegab; 
}
void read_param::read_params(char filename[2048])
  {
    FILE *inp;
    int ii;
    //Reading the input simulation parameter file
    inp=fopen(filename,"r");
    /* get parameters for reion simulation */
    fscanf(inp,"%d",&Nnion);
    //Nnion is the number of nion values for which we will run the simulation
    //You will get a x_HI map for each value of nion
    //Allocating memory for some variables

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
  }
