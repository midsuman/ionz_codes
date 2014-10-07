#include "ion.h"

struct params input_param;
void read_params(char filename[2048])
  {
    FILE *inp;
    int ii;
    //Reading the input simulation parameter file
    inp=fopen(filename,"r");
    /* get parameters for reion simulation */
    fscanf(inp,"%d",&input_param.Nnion);
    input_param.nion = malloc(sizeof(int)*input_param.Nnion);
    //Nnion is the number of nion values for which we will run the simulation
    //You will get a x_HI map for each value of nion
    //Allocating memory for some variables

    //Reading nion values from the input file
    //nion is the efficiency paramter for ionization (a combination of f_esc,f_star etc)
    // We try with variuos values of nion
    for(ii=0;ii<input_param.Nnion;ii++)
      {
	fscanf(inp,"%f",&(input_param.nion[ii]));
      }
    fscanf(inp,"%f%f%f%f%f",&input_param.z,&input_param.omegam,&input_param.omegalam,&input_param.Hubble_h,&input_param.omegab);
    //Reading cosmological parameters values
    //vaa initially scans the redshift of the Nbody simulation
    input_param.a_expansion = 1/(1+input_param.z);// converts it into scale factor
    fscanf(inp,"%d%d%d%f",&input_param.N1,&input_param.N2,&input_param.N3,&box);
    
    fclose(inp);
 
    input_param.gridsize = input_param.boxsize/(input_param.N1); // grid size in Mpc/h
  }
