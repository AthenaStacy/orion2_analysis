//Reads in Gadget2 data, converts to an Orion2 unifrom grid and outputs file to be read in during Orion2 initialization

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#define MAXREF 20
#define PI 3.14159265359

#define wpot 0
#define use_tsc 0
#define convert_vpot 0
#define readB 1
#define writeV 0
#define cell_center 1

#define const_gamma 0
#define which_sim 3

//#define which_snap 679  //bin_HR10_ideal
//#define which_snap 5979    //bin_zoom10_new_cut (which_sim == 0)
//#define which_snap 6600  //bin_zoom10_new_cut_ref3_6600 max dens ~ 1.e6 cm^-3
//#define which_snap 7000  //bin_zoom10_new_cut_ref3_6600 max dens ~ 1.e8 cm^-3 (?)

#define which_snap 7055  //bin_zoom10_new_cut_ref2 t = 12 x 1e10 seconds after snapshot 7000
//#define which_snap 7130  //bin_zoom10_new_cut_ref3_7130 max dens ~ 5.e11 cm^-3 
//#define which_snap 7131  //bin_zoom10_new_cut_ref3_wpot_7130 max dens ~ 5.e11 cm^-3

//#define ref_lev 256 //number of grid cells along each x-y-z direction
//#define width  20.0 //full width of box in pc to be "cut out" of Gadget2 simulation box and mapped onto the Orion2 grid 

//#define ref_lev 256
//#define width 1.6

//#define ref_lev 256
//#define width 12.8

//#define ref_lev 16
//#define width 1.0

/*
#define ref_init 128
#define extra 10
#define width_init  1.0
*/

#define ref_init 128
#define extra 10
#define width_init  0.1

/*
#define ref_init 200
#define extra 10
#define width_init  1.0
*/

/*
#define ref_init 256
#define extra 10
//#define width_init  (1.3/2/2/2)
#define width_init  1.0
*/

#define hfac_out 0.7

#define hubble_param 0.7


double width;
int ref_lev, idnum=1;
int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize);
int write_snapshot(char *fname, int files, char *outname, char *outname0, char *outname1, char *outname2, char *outname3, char* outname4, char *outname8, char *outname9, char *outname10, double delx, double dely, double delz, double delvx, double delvy, double delvz, int NumGas, int proc);
int reordering(void);
int unit_conversion(void);
int allocate_memory(void);
int allocate_memory2(void);
int allocate_memory3(void);
double calc_kernel_spline(int n, double x, double y, double z, double hsm, double grid_size, double grid_size_half);
double calc_kernel_tsc(int n, double x, double y, double z, double hsm, double grid_size, double grid_size_half);

//functions for calculating bfield from vpot
double calc_kernel_spline_vpot(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double cosmo_fac);
double del_kernel(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp);
double del_rho(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp);
double curl(int dir, double x1, double y1, double z1, double x2, double y2, double z2);

struct io_header_1
{
  int      npart[MAXREF];
  double   mass[MAXREF];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[MAXREF];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam; 
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header1;


struct io_header_old
{
  int      npart[MAXREF];
  double   mass[MAXREF];
  double   time;
  double   redshift;
  int      flag_sfr;
  int      flag_feedback;
  int      npartTotal[MAXREF];
  int      flag_cooling;
  int      num_files;
  double   BoxSize;
  double   Omega0;
  double   OmegaLambda;
  double   HubbleParam;
  char     fill[256- 6*4- 6*8- 2*8- 2*4- 6*4- 2*4 - 4*8];  /* fills to 256 Bytes */
} header_old;


int     NumPart, Ngas;

struct plist_data
{
  int Id;
  double disx, disy, disz;
}*plist;


struct particle_data 
{
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  int    Id;

  double disx, disy, disz;
  double  Rho, U, Pres, nh, Density, hsm, hsm_phys;
  //double Temp, sink;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII gam;
  double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double nh_test;
  double Bfieldx, Bfieldy, Bfieldz;
  double Afieldx, Afieldy, Afieldz, fcor;
  double mass_mapped;
  double mass_shift_mapped;
  double pot;
  double dummy;
  int convert;
} *P;


double  Time, zred;



/* Here we load a snapshot file. It can be distributed
 * onto several files (for files>1).
 * The particles are brought back into the order
 * implied by their ID's.
 * A unit conversion routine is called to do unit
 * conversion, and to evaluate the gas temperature.
 */
int main(int argc, char **argv)
{
  char path_in[200], path_out[200], input_fname[200], input_fname2[200], output_fname[200], basename[200], basename2[200], basenameout[200];
  char output_dens[200], output_vx[200], output_vy[200], output_vz[200], output_egy[200];
  char output_bx[200], output_by[200], output_bz[200];
  int  i, j, n, type, snapshot_number, files, Ngas, Ngas2, random, num_err_HI;
  double x,y,z,x1,y1, error, nh, nhmax, ref_lev_doub, grid_size, grid_size_half;
  double delx=0, dely=0, delz=0, delvx, delvy, delvz, boxsize, dis, disx, disy, disz, disAU;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;
  FILE *outfile, *outfile0, *outfile1, *outfile2, *outfile3, *outfile4, *outfile8, *outfile9, *outfile10, *infile;
  int npes, myrank, ierr;

  //For NON-cosmological runs
  //Time = 1.0;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  if(which_sim == 0)
    {
    sprintf(path_in, "/work/00863/minerva/bin_zoom_cut");
    sprintf(path_out, "/work/00863/minerva/orion");
    if(writeV == 1)
      sprintf(path_out, "/work/00863/minerva/orion/bfield_comp_vpot_fwards");
    sprintf(basename, "bin_zoom10_new_cut");
    sprintf(basename2, "bin_zoom10");
    }

   if(which_sim == 1)
    {
    sprintf(path_in, "/work/00863/minerva/");
    sprintf(path_out, "/work/00863/minerva/orion");
    if(writeV == 1)
      sprintf(path_out, "/work/00863/minerva/orion/bfield_comp_vpot_fwards");
    sprintf(basename, "bin_HR10");
    sprintf(basename2, "bin_zoom10");
    }

  if(which_sim == 3)
    {
    sprintf(path_in,  "/scratch/00863/minerva/");

    //sprintf(path_in,  "/work/00863/minerva/bin_zoom");
    sprintf(path_out, "/work/00863/minerva/orion/bin_zoom");

    sprintf(basename, "bin_zoom10_new_cut_ref3");
    //sprintf(basename, "bin_zoom10_new_cut_ref3_wpot");

    sprintf(basename2, "bin_zoom10");
    }
 
   if(which_sim == 4)
    {
    sprintf(path_in, "/work/00863/minerva/");
    sprintf(path_out, "/work/00863/minerva/orion");
    if(writeV == 1) 
      sprintf(path_out, "/work/00863/minerva/orion/bfield_comp_vpot_fwards");
    sprintf(basename, "bin_HR10_ideal");
    sprintf(basename2, "bin_HR10_ideal");
    }
 
/*
////////////////////////////////////////////////////////////////////////////////
//Use this section if we want central position to be taken from an earlier snapshot
  snapshot_number = snapcenter;
  files=1;
  boxsize = 140.0;

  sprintf(input_fname, "%s/%s_%04d", path_in, basename, snapshot_number);
  if(which_sim == 4)  sprintf(input_fname, "%s/%s_%03d", path_in, basename, snapshot_number);
  if(snapshot_number > 9999) sprintf(input_fname, "%s/%s_%05d", path_in, basename, snapshot_number);

  Ngas = read_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  unit_conversion();

  nhmax = 0;
  for(n=0;n<Ngas;n++)
     {
     nh = P[n].nh;
     if(nh > nhmax && P[n].sink > -1)
        {
        nhmax = nh;
        delx = P[n].Pos[0];  
        dely = P[n].Pos[1];
        delz = P[n].Pos[2];
        }
     }
   if(myrank == 0) printf("nhmax = %lg, delx = %lg, dely = %lg, delz = %lg\n", nhmax, delx, dely, delz);
/////////////////////////////////////////////////////////////////////////////////
*/
  for(j=7;j<=7;j=j++)
  {

  //snapshot_number = j;
  snapshot_number = which_snap;  //bin_HR10_ideal
  files=1;     

  double extra_doub, ref_init_doub;

  extra_doub = (double)extra;

  ref_lev = (ref_init+extra);
  ref_init_doub = (double)ref_init;
  width = width_init + (extra_doub/ref_init_doub)*width_init;


  boxsize = 140.0;

  ref_lev_doub = (double)ref_lev;
  grid_size = width/ref_lev_doub;
  grid_size_half = grid_size/2.0;

  if(myrank == 0) printf("ref_lev = %d, width = %lg, grid_size = %lg\n", ref_lev, width, grid_size);

  sprintf(input_fname, "%s/%s_%04d", path_in, basename, snapshot_number);
  if(which_sim == 4)  sprintf(input_fname, "%s/%s_%03d", path_in, basename, snapshot_number);  
  if(snapshot_number > 9999) sprintf(input_fname, "%s/%s_%05d", path_in, basename, snapshot_number);

  ncount_doub = xCOM = yCOM = zCOM = delvx = delvy = delvz = vxCOM = vyCOM = vzCOM = 0.;

  //read in Gadget2 file
  Ngas = read_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  //read in density and Bfield values from analytic calculation
  if(readB == 1)
    {
    sprintf(input_fname2, "%s/%s_bfield_%04d", path_out, basename2, snapshot_number);
    if(which_sim == 3)
      {
      sprintf(input_fname2, "%s/%s_bfield_ref3_%04d", path_out, basename2, snapshot_number);
      if(snapshot_number > 9999) sprintf(input_fname2, "%s/%s_bfield_ref3_%05d", path_out, basename2, snapshot_number);
      }
    if(!(infile=fopen(input_fname2,"r")))
      {
        printf("can't open file `%s`\n",input_fname2);
        exit(0);
      }

    if(convert_vpot != 1)
    {
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldx, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldy, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldz, sizeof(double), 1, infile);
    }

    if(convert_vpot == 1)
    {
    for(n=0;n<Ngas;n++)
        P[n].Bfieldx = P[n].Bfieldy = P[n].Bfieldz = 0;
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afieldx, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afieldy, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Afieldz, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++) P[n].fcor = 1;
    }

    for(n=0;n<Ngas;n++)
        fread(&P[n].nh_test, sizeof(double), 1, infile); 
    fclose(infile);
    }

  unit_conversion();

  // calculate central density maximum and velocity CoM for "re-centering" of Orion2 data 

  nhmax = 0;
  for(n=0;n<Ngas;n++) 
     {
     nh = P[n].nh;
     if(nh > nhmax && P[n].sink > -1)
        {
        nhmax = nh;
        delx = P[n].Pos[0];  //del's based upon location of maximum density
        dely = P[n].Pos[1];
        delz = P[n].Pos[2];
        }
      if(myrank == 0 && n % 1000 == 0)
        printf("n = %d, rho = %lg, rho_test = %lg, Bfieldx = %lg\n", n, P[n].Rho, P[n].nh_test, P[n].Bfieldx); 
     }

  num_err_HI = 0;
  double nh_min = 1.e50;
  for(n=0;n < Ngas; n++)
     {
     dis = pow(((P[n].Pos[0]-delx)*(P[n].Pos[0]-delx) + (P[n].Pos[1]-dely)*(P[n].Pos[1]-dely) + (P[n].Pos[2]-delz)*(P[n].Pos[2]-delz)), 0.5);
     dis=dis*1.e3*Time/(hubble_param);
     disAU=dis*206264.806;

     disx = ((P[n].Pos[0]-delx))*1.e3*Time/(hubble_param);
     disy = ((P[n].Pos[1]-dely))*1.e3*Time/(hubble_param);
     disz = ((P[n].Pos[2]-delz))*1.e3*Time/(hubble_param);

     P[n].Density = P[n].convert = 0;

     //if(dis < width)
     if(fabs(disx) < width/2.0 && fabs(disy) < width/2.0 && fabs(disz) < width/2.0)
     //if(P[n].nh > nhmax/10.0) 
       {
       vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
       vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
       vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
       ncount_doub = ncount_doub + P[n].Mass;
       }
     if(fabs(disx) < width && fabs(disy) < width && fabs(disz) < width)
       {
       P[n].convert = 1;
       idnum++;
       if(P[n].nh < nh_min) nh_min = P[n].nh;       
       }
      if(myrank == 0 && n % 1000 == 0)
        printf("n = %d, rho = %lg, idnum = %d\n", n, P[n].Rho, idnum);
     }   

    printf("final idnum = %d, nh_min = %lg\n", idnum, nh_min);

    if(!(plist=(struct plist_data *) malloc(idnum*sizeof(struct plist_data))))
       {
       fprintf(stderr,"failed to allocate memory.\n");
       exit(0);
       }
    printf("plist memory allocated\n");

    n=0;
    for(i = 0; i < Ngas; i++)
      {
      if(P[i].convert == 1)
        {
        plist[n].Id = P[i].Id;
        n++;
        }
      if(myrank == 0 && i % 1000 == 0 && n > 1)
        printf("i = %d, rho = %lg, Id = %d, id = %d\n", i, P[i].Rho, P[i].Id, plist[n-1].Id);
      }

   delvx = vxCOM/ncount_doub;
   delvy = vyCOM/ncount_doub;
   delvz = vzCOM/ncount_doub;
   //delvx = delvy = delvz = 0.0;

  printf("nhmax = %lg, delvx = %lg, delvy = %lg, delvz = %lg\n", nhmax, delvx, delvy, delvz);

  sprintf(output_fname, "%s/gadget2orion", path_out);
  if(which_sim == 1)
    sprintf(output_fname, "%s/gadget2orion_bfield_%04d", path_out, snapshot_number);
  if(which_sim == 3)
    sprintf(output_fname, "%s/gadget2orion_bfield_ref3_%04d", path_out, snapshot_number);

  sprintf(output_dens, "%s/gadget2dens_ideal_%04d", path_out, snapshot_number);
  sprintf(output_vx, "%s/gadget2vx_ideal_%04d", path_out, snapshot_number);
  sprintf(output_vy, "%s/gadget2vy_ideal_%04d", path_out, snapshot_number);
  sprintf(output_vz, "%s/gadget2vz_ideal_%04d", path_out, snapshot_number);
  sprintf(output_egy, "%s/gadget2egy_ideal_%04d", path_out, snapshot_number);
  sprintf(output_bx, "%s/gadget2bx_ideal_%04d", path_out, snapshot_number);
  sprintf(output_by, "%s/gadget2by_ideal_%04d", path_out, snapshot_number);
  sprintf(output_bz, "%s/gadget2bz_ideal_%04d", path_out, snapshot_number);


  if(myrank == 0)
    {
    outfile=fopen(output_fname,"w");
    fclose(outfile);

    outfile0=fopen(output_dens,"w");
    fclose(outfile0);

    outfile1=fopen(output_vx,"w");
    fclose(outfile1);

    outfile2=fopen(output_vy,"w");
    fclose(outfile2);

    outfile3=fopen(output_vz,"w");
    fclose(outfile3);

    outfile4=fopen(output_egy,"w");
    fclose(outfile4);

    outfile8=fopen(output_bx,"w");
    fclose(outfile1);

    outfile9=fopen(output_by,"w");
    fclose(outfile2);

    outfile10=fopen(output_bz,"w");
    fclose(outfile3);
    }


  Ngas2 = write_snapshot(input_fname, files, output_fname, output_dens, output_vx, output_vy, output_vz, output_egy, output_bx, output_by, output_bz, delx, dely, delz, delvx, delvy, delvz, Ngas, myrank);

  free(P);
  }

  ierr=MPI_Finalize();

}


/* this template shows how one may convert from Gadget's units
 * to cgs units.
 * In this example, the temperate of the gas is computed.
 * (assuming that the electron density in units of the hydrogen density
 * was computed by the code. This is done if cooling is enabled.)
 */
int unit_conversion(void)
{
  double GRAVITY, BOLTZMANN, PROTONMASS;
  double UnitLength_in_cm, UnitMass_in_g, UnitVelocity_in_cm_per_s;
  double UnitTime_in_s, UnitDensity_in_cgs, UnitPressure_in_cgs, UnitEnergy_in_cgs;  
  double G, Xh, HubbleParam;

  int i;
  double MeanWeight, u, gamma;
  double h2frac, muh2, muh2in, temp;

  /* physical constants in cgs units */
  GRAVITY   = 6.672e-8;
  BOLTZMANN = 1.3806e-16;
  PROTONMASS = 1.6726e-24;

  /* internal unit system of the code */
  UnitLength_in_cm= 3.085678e21;   /*  code length unit in cm/h */
  UnitMass_in_g= 1.9891e43;         /*  code mass unit in g/h */
  UnitVelocity_in_cm_per_s= 1.0e5;

  UnitTime_in_s= UnitLength_in_cm / UnitVelocity_in_cm_per_s;
  UnitDensity_in_cgs= UnitMass_in_g/ pow(UnitLength_in_cm,3);
  UnitPressure_in_cgs= UnitMass_in_g/ UnitLength_in_cm/ pow(UnitTime_in_s,2);
  UnitEnergy_in_cgs= UnitMass_in_g * pow(UnitLength_in_cm,2) / pow(UnitTime_in_s,2);

  G=GRAVITY/ pow(UnitLength_in_cm,3) * UnitMass_in_g * pow(UnitTime_in_s,2);


  Xh= 0.76e0;  /* mass fraction of hydrogen */
  HubbleParam= 0.7e0;

  //for NON-comoving sims
  //HubbleParam = 1.0;

  for(i=0; i<=NumPart; i++)
    {
      if(P[i].Type==0)  /* gas particle */
	{
/*	  MeanWeight= 4.0/(3*Xh+1+4*Xh*P[i].elec) * PROTONMASS;  */

       MeanWeight=1.2195;
       h2frac=2.0*P[i].H2I;

       muh2in=(0.24/4.0) + ((1.0-h2frac)*0.76) + (h2frac*.76/2.0);
       muh2=pow(muh2in, -1.0);

       if(muh2 >= 1.22)
         {
          MeanWeight=muh2;
         }

          MeanWeight=MeanWeight*PROTONMASS;

	  //MeanWeight= 1.22e0 * PROTONMASS;

	  /* convert internal energy to cgs units */

	  u  = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

#if(const_gamma)
	  gamma= 5.0/3.0;
#else
          gamma = P[i].gam;	 
#endif
	  /* get temperature in Kelvin */

	  temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Density * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
	  P[i].nh= P[i].Rho / MeanWeight;
          P[i].Pres = (gamma-1)*u*P[i].Rho;
          P[i].U = P[i].U * UnitEnergy_in_cgs/ UnitMass_in_g;

          //if(readB == 1)
          //  P[i].nh_test = P[i].nh_test* UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);

	  /*  printf("zred = %g", zred);*/
	}
    }
return(0);
}





/* this routine loads particle data from Gadget's default
 * binary file format. (A snapshot may be distributed
 * into multiple files.
 */
int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize)
{
  FILE *fd;
  FILE *outfile;
  char   buf[200];
  int    i,j,k,l,dummy,ntot_withmasses; 
  int    t,n,off,pc,pc_new=0,pc_sph;
  int NumPart_new = 0, Ngas_new = 0;
  int Idnew, nnew=1;
  double *pos, massnew, hsmnew, x, y, z, nnew_doub=1.0;
  double randomx, randomy, randomz;

  pos= (double*)malloc(sizeof(double) * 3);
	
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=0; i<files; i++, pc=pc_new)
    {
      if(files>1)
	sprintf(buf,"%s.%d",fname,i);
      else
	sprintf(buf,"%s",fname);

      if(!(fd=fopen(buf,"r")))
	{
	  printf("can't open file `%s`\n",buf);
	  exit(0);
	}

      printf("reading `%s' ...\n",buf); fflush(stdout);

      fread(&dummy, sizeof(dummy), 1, fd);
      fread(&header1, sizeof(header1), 1, fd);
      fread(&dummy, sizeof(dummy), 1, fd);

      if(files==1)
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
            {
	    NumPart+= header1.npart[k];
            printf("NumPart[%d] = %d\n", k, header1.npart[k]);
            }
	  Ngas= header1.npart[0];
	}
      else
	{
	  for(k=0, NumPart=0, ntot_withmasses=0; k<5; k++)
	    NumPart+= header1.npartTotal[k];
	  Ngas= header1.npartTotal[0];
	}

      for(k=0, ntot_withmasses=0; k<5; k++)
	{
	  if(header1.mass[k]==0)
	    ntot_withmasses+= header1.npart[k];
	}

      if(i==0)
	allocate_memory();

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
	{
	  for(n=0;n<header1.npart[k];n++)
	    {
	      fread(&P[pc_new].Pos[0], sizeof(double), 3, fd);
	      pc_new++;
	    }
	}
      SKIP;

      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              fread(&P[pc_new].Vel[0], sizeof(double), 3, fd);
              pc_new++;
            }
        }
      SKIP;


      SKIP;
      for(k=0,pc_new=pc;k<6;k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              fread(&P[pc_new].Id, sizeof(int), 1, fd);
              pc_new++;
            }
        }
      SKIP;

      if(ntot_withmasses>0)
      {
        SKIP;
     }

      for(k=0, pc_new=pc; k<6; k++)
        {
          for(n=0;n<header1.npart[k];n++)
            {
              P[pc_new].Type=k;
              if(header1.mass[k]==0)
                {
                fread(&P[pc_new].Mass, sizeof(double), 1, fd);
                }
              else
                P[pc_new].Mass= header1.mass[k];
              pc_new++;
            }
        }
      if(ntot_withmasses>0)
      {
        SKIP;
      }


		if(header1.npart[0]>0)
		{
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].U, sizeof(double), 1, fd);
			    pc_sph++;
			}
			SKIP;
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
			{
				fread(&P[pc_sph].Density, sizeof(double), 1, fd);
				pc_sph++;
			}
			SKIP;
			
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].hsm, sizeof(double), 1, fd);
				pc_sph++;
            }
			SKIP;

if(wpot == 1)
{		
        SKIP;
         for(k=0,pc_new=pc;k<6;k++)
           {
             for(n=0;n<header1.npart[k];n++)
               {
                 fread(&P[pc_new].pot, sizeof(double), 1, fd);
                 pc_new++;
               }
           }
         SKIP;
}

if(const_gamma != 1)
{	
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].H2I, sizeof(double), 1, fd);
				fread(&P[pc_sph].HII, sizeof(double), 1, fd);
				fread(&P[pc_sph].DII, sizeof(double), 1, fd);
				fread(&P[pc_sph].HDI, sizeof(double), 1, fd);
				fread(&P[pc_sph].HeII, sizeof(double), 1, fd);
				fread(&P[pc_sph].HeIII, sizeof(double), 1, fd);
				pc_sph++;
            }
			SKIP;
			
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].gam, sizeof(double), 1, fd);
				pc_sph++;
            }
			SKIP;
		
            printf("gam = %lg\n",P[100].gam);
	
			SKIP;
			for(n=0, pc_sph=pc; n<header1.npart[0];n++)
            {
				fread(&P[pc_sph].sink, sizeof(double), 1, fd); 
				pc_sph++;
            }
			SKIP;
} 
			
        }		
      fclose(fd);
  }

  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.7e \n",Time);

  //For NON-cosmological runs
  //Time = 1.0;
  printf("Time= %12.7e \n",Time);

  printf("L= %6.2f \n",header1.BoxSize);
  fflush(stdout);
  return(Ngas);
}


int write_snapshot(char *fname, int files, char *outname, char *outname0, char *outname1, char *outname2, char *outname3, char *outname4, char *outname8, char *outname9, char *outname10, double delx, double dely, double delz, double delvx, double delvy, double delvz, int NumGas, int myrank)
{
  FILE *outfile, *outfile0, *outfile1, *outfile2, *outfile3, *outfile4, *outfile8, *outfile9, *outfile10;
  char   buf[200];
  int i, nmin, nmax, j, k, l, n, send_int, send_size, tot_proc=1, tot_proc_sum;
  int imin, imax, jmin, jmax, kmin, kmax;
  int NumGas_per_proc, vartype, varnum=8;
  double ref_lev_doub, i_doub; 
  double hsm_conv, vel_conv, dens_conv, mass_conv, mass, rho, rho_shift, hsm;
  double kernel_conv, kernel_phys, kernel_shift, kernel_shift_phys; 
  double fac, fac_shift, kernel, rad, ratio, xgrid, ygrid, zgrid, grid_size, grid_size_half, grid_vol;
  double grid_arr[ref_lev], ntot=0, masstot=0, xmomtot=0, ymomtot=0, zmomtot=0, etot=0;
  double masstot_alt=0, xmomtot_alt=0, ymomtot_alt=0, zmomtot_alt=0, etot_alt=0, hfac = 1.5;
  double *Gdum, *Gdumx, *Gdumy, *Gdumz, *Grho, *Grho_shift, *Grho_tot, *Grho_shift_tot; 
  double *Gvelx, *Gvely, *Gvelz, *Gpres, *GH2I, *GHDI, *GHII; 
  double *Gtot, *Gtotx, *Gtoty, *Gtotz, *Gnpart, *Gnpart_tot;
  MPI_Status status;

  if(readB == 1)
    varnum=12;

  //For NON-cosmological runs
  //Time = 1.0;

  ref_lev_doub = (double)ref_lev;
  vel_conv = 1.e5*pow(Time,0.5);
  dens_conv = P[100].Rho/P[100].Density;
  mass_conv = (1.e10/hubble_param)*1.9891e33;
  kernel_conv =  pow(3.08567758e18*1.e3*Time/(hubble_param),-3);

  grid_size = width/ref_lev_doub;
  grid_size_half = grid_size/2.0;
  grid_vol = pow(grid_size*3.08567758e18, 3); 
 
  MPI_Allreduce(&tot_proc, &tot_proc_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  printf("myrank = %d, tot_proc_sum = %d\n", myrank, tot_proc_sum);

//subtract out COM velocities and shift particle positions to center of Orion2 box
//let densest point sit at (-gsize_half, -gsize_half, -gsize_half) instead of (0,0,0)
      for(n=0;n<NumGas;n++)
        {
        P[n].hsm_phys = P[n].hsm*1.e3*Time/hubble_param; 
        P[n].disx = 1.e3*Time/(hubble_param)*(P[n].Pos[0] - delx) - grid_size_half;  
        P[n].disy = 1.e3*Time/(hubble_param)*(P[n].Pos[1] - dely) - grid_size_half;
        P[n].disz = 1.e3*Time/(hubble_param)*(P[n].Pos[2] - delz) - grid_size_half;
        P[n].Vel[0] = P[n].Vel[0] - delvx;
        P[n].Vel[1] = P[n].Vel[1] - delvy;
        P[n].Vel[2] = P[n].Vel[2] - delvz;
        P[n].mass_mapped = P[n].mass_shift_mapped = 0;
        if(fabs(fabs(P[n].disx) - 0.0*P[n].hsm_phys) < width/2.0 && fabs(fabs(P[n].disy) - 0.0*P[n].hsm_phys) < width/2.0 && fabs(fabs(P[n].disz) - 0.0*P[n].hsm_phys) < width/2.0)
          masstot = masstot + P[n].Mass;
          ntot++;
          xmomtot = xmomtot + P[n].Mass*P[n].Vel[0];
          ymomtot = ymomtot + P[n].Mass*P[n].Vel[1];
          zmomtot = zmomtot + P[n].Mass*P[n].Vel[2];
          etot = etot + P[n].Mass*P[n].U;
        }

//////////write out some header info///////////////////
    if(myrank==0)
    {
        outfile=fopen(outname,"w");
        fwrite(&myrank, sizeof(int), 1, outfile);
        fwrite(&ref_lev, sizeof(int), 1, outfile);
        fwrite(&width, sizeof(double), 1, outfile);
        fclose(outfile);
    }
///////////////////////////////////////////////////////////

      for(i=0;i<ref_lev;i++)
        {
        i_doub = (double)i;
        grid_arr[i] = (width * i_doub/ref_lev_doub) - width/2.0;
        }	

      NumGas_per_proc = (int) (NumGas/tot_proc_sum); 
      nmin = NumGas_per_proc * (myrank);
      nmax = NumGas_per_proc * (myrank+1) -1;

      if(myrank == tot_proc_sum - 1)
        nmax = NumGas-1;
      if(nmin >= NumGas-1)
        nmin = nmax = NumGas-1;
      if(nmax >= NumGas-1)
        nmax = NumGas-1;


 printf("myrank = %d, nmin = %d, nmax = %d, ref_lev = %d, ref_lev_doub = %lg, masstot = %lg, xmomtot = %lg, ymomtot = %lg, zmomtot = %lg, etot = %lg\n", myrank, nmin, nmax, ref_lev, ref_lev_doub, masstot*1.e10/hubble_param, xmomtot*vel_conv*mass_conv, ymomtot*vel_conv*mass_conv, zmomtot*vel_conv*mass_conv, etot*mass_conv);


//////////////////////////////////////////////////////////////////////////////////////////////
/////conversion of vector potential to B-field!
#if (convert_vpot)

  Gdumx = (double *) malloc(Ngas * sizeof(double));
  Gdumy = (double *) malloc(Ngas * sizeof(double));
  Gdumz = (double *) malloc(Ngas * sizeof(double));
  Gtotx = (double *) malloc(Ngas * sizeof(double));
  Gtoty = (double *) malloc(Ngas * sizeof(double));
  Gtotz = (double *) malloc(Ngas * sizeof(double));

  double dwdx, dwdy, dwdz, ax, ay, az, rcheck;
  double cosmo_fac = Time/hubble_param;

  for(n=0;n<Ngas;n++) 
    {
    P[n].Pos[0] = P[n].Pos[0]*cosmo_fac;
    P[n].Pos[1] = P[n].Pos[1]*cosmo_fac;
    P[n].Pos[2] = P[n].Pos[2]*cosmo_fac;
    P[n].hsm = P[n].hsm*cosmo_fac;
    }

  double nh_min = 1.e50;
  for(n=0;n<Ngas;n++)
    if(P[n].convert == 1)
      {
      if(P[n].nh < nh_min) nh_min = P[n].nh;
      }
  if(myrank == 0) printf ("vec. pot convert, nh_min = %lg\n", nh_min);

  struct particle_data psave, psource;
  printf("start reordering\n");
  for(n=0;n<Ngas;n++)
      if(P[n].nh > nh_min)
      for(k=0;k<idnum;k++)
          if(P[n].Id == plist[k].Id)
             {
             psource = P[n];
             psave   = P[k];
             P[n] = psave;
             P[k] = psource;
             }
  printf("finished reordering\n");

  for(n=0;n<Ngas;n++)
    {
    Gdumx[n] = 0;
    Gtotx[n] = 0;
    } 
  for(n=nmin;n<=nmax;n++)
     {
     if(P[n].convert == 1)
       {
       for(i=0;i<idnum;i++)
         {
         if(i == n) continue;
         rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
         if(rcheck > P[n].hsm) continue;
         Gdumx[n] = Gdumx[n] + P[i].Mass * calc_kernel_spline_vpot(n, P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
          }
        }
     }

  send_size = Ngas;
  MPI_Allreduce(&Gdumx[0], &Gtotx[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(n=0;n<Ngas;n++)
    P[n].Density = Gtotx[n];
 
  if(myrank == 0) printf("finished density calculation, Ngas = %d, idnum = %d\n", Ngas, idnum);
  for(n=0;n<Ngas;n++)
    if(myrank == 0 && P[n].convert == 1)
      printf("n = %d, rho = %lg,  density = %lg\n", n, P[n].Rho, P[n].Density);

  for(n=0;n<Ngas;n++)
    {
    Gdumx[n] = 0;
    Gtotx[n] = 0;
    }
  for(n=nmin;n<=nmax;n++)
     {
     if(P[n].convert == 1)
       {
        for(i=0;i<idnum;i++)
           {
           if(i == n) continue;
           rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
           if(rcheck > P[n].hsm) continue;
           Gdumx[n] = Gdumx[n] + (P[n].hsm /  3. / P[n].Density) * P[i].Mass * del_rho(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0],  P[i].Pos[1], P[i].Pos[2], P[n].hsm, 0);
           }
         }   //end if(convert == 1) statement
      }  //end for(n=0;n < Ngas; n++) loop

  send_size = Ngas;
  MPI_Allreduce(&Gdumx[0], &Gtotx[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  for(n=0;n<Ngas;n++)
    P[n].fcor = Gtotx[n] + 1.0;

  printf("finished fcor calculation\n");
  for(n=0;n<Ngas;n++)
    if(myrank == 0 && P[n].convert == 1)      
      printf("n = %d, rho = %lg,  fcor = %lg\n", n, P[n].Rho, P[n].fcor);

  for(n=0;n<Ngas;n++)
    {
    Gdumx[n] = 0;
    Gdumy[n] = 0;
    Gdumz[n] = 0;
    Gtotx[n] = 0;
    Gtoty[n] = 0;
    Gtotz[n] = 0;
    }

  for(n=nmin;n<=nmax;n++)
     {
     if(P[n].convert == 1)
       {
        for(i=0;i<idnum;i++)
           {
           if(i == n) continue;
           rcheck = pow((P[n].Pos[0] - P[i].Pos[0])*(P[n].Pos[0] - P[i].Pos[0]) + (P[n].Pos[1] - P[i].Pos[1])*(P[n].Pos[1] - P[i].Pos[1]) + (P[n].Pos[2] - P[i].Pos[2])*(P[n].Pos[2] - P[i].Pos[2]), 0.5);
           if(rcheck > P[n].hsm) continue;
           dwdx = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[0] - P[i].Pos[0]);
           dwdy = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[1] - P[i].Pos[1]);
           dwdz = del_kernel(P[n].Pos[0], P[n].Pos[1], P[n].Pos[2], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2], P[n].hsm, P[n].Pos[2] - P[i].Pos[2]);
           ax = P[n].Afieldx - P[i].Afieldx;
           ay = P[n].Afieldy - P[i].Afieldy;
           az = P[n].Afieldz - P[i].Afieldz;
           Gdumx[n] = Gdumx[n] + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(0, ax, ay, az, dwdx, dwdy, dwdz);
           Gdumy[n] = Gdumy[n] + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(1, ax, ay, az, dwdx, dwdy, dwdz);
           Gdumz[n] = Gdumz[n] + (1./P[n].Density/P[n].fcor)*P[i].Mass*curl(2, ax, ay, az, dwdx, dwdy, dwdz);
           if(myrank == 0 && n%2000 == 0) 
             printf("n = %d, dwdx = %lg, dwdy = %lg, dwdz = %lg, ax = %lg, ay = %lg, az = %lg\n", n, dwdx, dwdy, dwdz, ax, ay, az);
           }
         }   //end if(convert == 1) statement
      }  //end for(n=0;n < Ngas; n++) loop

  printf("finished vector potential conversion\n");

      send_size = Ngas;
      MPI_Allreduce(&Gdumx[0], &Gtotx[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&Gdumy[0], &Gtoty[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      MPI_Allreduce(&Gdumz[0], &Gtotz[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

   for(n=0;n<Ngas;n++)
    {
    P[n].Bfieldx = Gtotx[n];
    P[n].Bfieldy = Gtoty[n];
    P[n].Bfieldz = Gtotz[n];
    if(myrank == 0 && P[n].convert == 1)
      printf("n = %d, rho = %lg, rho_test = %lg, density = %lg, fcor = %lg, Bfieldx = %lg, Bfield y = %lg, Bfieldz = %lg\n", n, P[n].Rho, P[n].nh_test, P[n].Density, P[n].fcor, P[n].Bfieldx, P[n].Bfieldy, P[n].Bfieldz);
    }

  for(n=0;n<Ngas;n++) 
    {
    P[n].Pos[0] = P[n].Pos[0]/cosmo_fac;
    P[n].Pos[1] = P[n].Pos[1]/cosmo_fac;
    P[n].Pos[2] = P[n].Pos[2]/cosmo_fac;
    P[n].hsm = P[n].hsm/cosmo_fac;
    }

free(Gdumx);
free(Gdumy);
free(Gdumz);
free(Gtotx);
free(Gtoty);
free(Gtotz);
#endif
/////////////////////////////////////////////////////////////////////////////////


  Gdum =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
  Grho =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
  Gtot =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
  Gnpart =(double *) malloc(pow(ref_lev,3) * sizeof(double));
  Gnpart_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));
  Grho_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));
  Grho_shift =(double *) malloc(pow(ref_lev,3) * sizeof(double));
  Grho_shift_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));

  printf("myrank = %d, Memory allocation done\n", myrank);

// For each processor, find particles whose smoothing lengths overlap with grid cell centers, and add their weighted values to the grid
      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             Grho[l] = 0.0;
             Grho_shift[l] = 0.0;
             Gtot[l] = 0.0;
             Gnpart[l] = 0.0;
             Gnpart_tot[l] = 0.0;
             Grho_tot[l] = 0.0;
             Grho_shift_tot[l] = 0.0;
             l++;
             }

    for(vartype=-1;vartype<varnum;vartype++)
     {

      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             Gdum[l] = 0.0;
             l++;
             }

      for(n=nmin;n<=nmax;n++)
        {
        hsm = P[n].hsm;
        mass = P[n].Mass*mass_conv;

        if(P[n].disx + hfac*P[n].hsm_phys < -width/2.0) continue;
        if(P[n].disx - hfac*P[n].hsm_phys >  width/2.0) continue;
        if(P[n].disy + hfac*P[n].hsm_phys < -width/2.0) continue;
        if(P[n].disy - hfac*P[n].hsm_phys >  width/2.0) continue;
        if(P[n].disz + hfac*P[n].hsm_phys < -width/2.0) continue;
        if(P[n].disz - hfac*P[n].hsm_phys >  width/2.0) continue;

        imin = (int)((P[n].disx - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        imax = (int)((P[n].disx + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        jmin = (int)((P[n].disy - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        jmax = (int)((P[n].disy + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        kmin = (int)((P[n].disz - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        kmax = (int)((P[n].disz + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);

        if(vartype == -1)
          for(i=imin;i<=imax;i++)
            for(j=jmin;j<=jmax;j++)
              for(k=kmin;k<=kmax;k++)
                 {
                 xgrid = -width/2. + ((double) i) * grid_size + grid_size_half;
                 ygrid = -width/2. + ((double) j) * grid_size + grid_size_half;
                 zgrid = -width/2. + ((double) k) * grid_size + grid_size_half;

                 kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                 if(P[n].hsm_phys < grid_size_half && use_tsc ==1) 
                    kernel  = calc_kernel_tsc(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;
                 P[n].mass_mapped = P[n].mass_mapped + mass*kernel_phys*grid_vol;

                 kernel = calc_kernel_spline(n, xgrid-grid_size_half, ygrid-grid_size_half, zgrid-grid_size_half, hfac_out*hsm, grid_size, grid_size_half);
                 if(P[n].hsm_phys < grid_size_half && use_tsc ==1)
                    kernel  = calc_kernel_tsc(n, xgrid-grid_size_half, ygrid-grid_size_half, zgrid-grid_size_half, hfac_out*hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;
                 P[n].mass_shift_mapped = P[n].mass_shift_mapped + mass*kernel_phys*grid_vol;
                 }

        if(imin < 0)           imin = 0;
        if(imax > ref_lev-1)   imax = ref_lev-1;
        if(jmin < 0)           jmin = 0;
        if(jmax > ref_lev-1)   jmax = ref_lev-1;
        if(kmin < 0)           kmin = 0;
        if(kmax > ref_lev-1)   kmax = ref_lev-1;

        if(imin > ref_lev-1 || jmin > ref_lev-1 || kmin > ref_lev-1 || imax < 0 || jmax < 0 || kmax < 0)
          continue;

        for(i=imin;i<=imax;i++)
          for(j=jmin;j<=jmax;j++)
            for(k=kmin;k<=kmax;k++)
               {
               l = k + (j)*(ref_lev) + (i)*(ref_lev)*(ref_lev);
 
               if(vartype == 0)
                 {
                 xgrid = grid_arr[i]+grid_size_half;
                 ygrid = grid_arr[j]+grid_size_half;
                 zgrid = grid_arr[k]+grid_size_half;

                 kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);                            if(P[n].hsm_phys < grid_size_half && use_tsc ==1) 
                    kernel  = calc_kernel_tsc(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;
                 rho = mass*kernel_phys;
                 fac = mass / P[n].mass_mapped;
                 if(rho > 0)  Grho[l]  = Grho[l] + rho*fac;               

                 kernel = calc_kernel_spline(n, xgrid-grid_size_half, ygrid-grid_size_half, zgrid-grid_size_half, hfac_out*hsm, grid_size, grid_size_half);
                 if(P[n].hsm_phys < grid_size_half && use_tsc ==1)
                    kernel  = calc_kernel_tsc(n, xgrid-grid_size_half, ygrid-grid_size_half, zgrid-grid_size_half, hfac_out*hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;
                 rho_shift = mass*kernel_phys;
                 fac_shift = mass / P[n].mass_shift_mapped;
                 if(rho_shift > 0)  Grho_shift[l]  = Grho_shift[l] + rho_shift*fac_shift;
                  }


               if(vartype > 0 && (Grho[l] > 0 || Grho_shift[l] > 0))
                 {
                 xgrid = grid_arr[i]+grid_size_half;
                 ygrid = grid_arr[j]+grid_size_half;
                 zgrid = grid_arr[k]+grid_size_half;

                 kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                 if(P[n].hsm_phys < grid_size_half && use_tsc ==1) 
                    kernel  = calc_kernel_tsc(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;
                 rho = mass*kernel_phys;
                 fac = mass / P[n].mass_mapped;

                 kernel = calc_kernel_spline(n, xgrid-grid_size_half, ygrid-grid_size_half, zgrid-grid_size_half, hfac_out*hsm, grid_size, grid_size_half);                          
                 if(P[n].hsm_phys < grid_size_half && use_tsc ==1)
                    kernel  = calc_kernel_tsc(n, xgrid-grid_size_half, ygrid-grid_size_half, zgrid-grid_size_half, hfac_out*hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;
                 rho_shift = mass*kernel_phys;
                 fac_shift = mass / P[n].mass_shift_mapped;


                  if(vartype == 1 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].Vel[0]*vel_conv*rho*fac; //Go ahead and convert velocity to cgs units 
                  if(vartype == 2 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].Vel[1]*vel_conv*rho*fac; 
                  if(vartype == 3 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].Vel[2]*vel_conv*rho*fac; 
                  if(vartype == 4 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].U*rho*fac;
                  if(vartype == 5 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].H2I*rho*fac;
                  if(vartype == 6 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].HDI*rho*fac;
                  if(vartype == 7 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].HII*rho*fac; 
#ifndef cell_center
if(myrank == 0 && l == 1) printf("Using cell edges!\n");
                  if(vartype == 8 && rho_shift > 0)
                    Gdum[l] = Gdum[l] + P[n].Bfieldx*rho_shift*fac_shift;
                  if(vartype == 9 && rho_shift > 0)
                    Gdum[l] = Gdum[l] + P[n].Bfieldy*rho_shift*fac_shift;
                  if(vartype == 10 && rho_shift > 0)
                    Gdum[l] = Gdum[l] + P[n].Bfieldz*rho_shift*fac_shift; 
#endif
#ifdef cell_center
if(myrank == 0 && l == 1) printf("Using cell centers!\n");
                  if(vartype == 8 && rho > 0)
                    //Gdum[l] = Gdum[l] + P[n].Bfieldx*rho*fac;
                    Gdum[l] = Gdum[l] + P[n].Bfieldx*rho;
                  if(vartype == 9 && rho > 0)
                    //Gdum[l] = Gdum[l] + P[n].Bfieldy*rho*fac;
                    Gdum[l] = Gdum[l] + P[n].Bfieldy*rho;
                  if(vartype == 10 && rho > 0)
                    //Gdum[l] = Gdum[l] + P[n].Bfieldz*rho*fac;
                    Gdum[l] = Gdum[l] + P[n].Bfieldz*rho;
#endif
                  if(vartype == 11 && rho > 0)
                    Gdum[l] = Gdum[l] + P[n].nh_test*rho*fac; 
                  }
               }

        }       

      printf("myrank = %d, main loop finished\n", myrank);

      MPI_Barrier(MPI_COMM_WORLD);

      send_size = pow(ref_lev,3);
      if(vartype == 0)
        {
        MPI_Allreduce(&Grho[0], &Grho_tot[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Grho_shift[0], &Grho_shift_tot[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

//////////write mapped data to file for Orion2!!!!!!!!!!!

     //write out density for each grid cell
     if(vartype == 0)
       MPI_Reduce(&Grho[0], &Gtot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

     if(vartype > 0) 
        MPI_Reduce(&Gdum[0], &Gtot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      printf("Rank 0 has the data now, myrank=%d\n", myrank);  fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
      if(myrank == 0 && vartype >= 0)
      {
      outfile=fopen(outname,"a");
      outfile0=fopen(outname0,"a");
      outfile1=fopen(outname1,"a");
      outfile2=fopen(outname2,"a");
      outfile3=fopen(outname3,"a");
      outfile4=fopen(outname4,"a");
      outfile8=fopen(outname8,"a");
      outfile9=fopen(outname9,"a");
      outfile10=fopen(outname10,"a");
      masstot = 0;
      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             if(vartype == 0)  masstot_alt = masstot_alt + Gtot[l]*grid_vol;  //total mass        
             if(vartype == 1)  xmomtot_alt = xmomtot_alt + Gtot[l]*grid_vol; //total momentum
             if(vartype == 2)  ymomtot_alt = ymomtot_alt + Gtot[l]*grid_vol; //total momentum
             if(vartype == 3)  zmomtot_alt = zmomtot_alt + Gtot[l]*grid_vol; //total momentum
             if(vartype == 4)  etot_alt = etot_alt + Gtot[l]*grid_vol;  //total internal energy
             if(vartype > 0 && vartype < 8)  Gtot[l] = (Gtot[l]/Grho_tot[l]);  //write out stuff that is NOT density!
#ifndef cell_center
             if(vartype == 8 || vartype == 9 || vartype == 10) Gtot[l] = (Gtot[l]/Grho_shift_tot[l]);
#endif
#ifdef cell_center
             if(vartype == 8 || vartype == 9 || vartype == 10) Gtot[l] = (Gtot[l]/Grho_tot[l]);
#endif
             if(i >= extra/2 && i < ref_lev-extra/2 && j >= extra/2 && j < ref_lev-extra/2 
                && k >= extra/2 && k < ref_lev-extra/2)  
               {
               if(vartype == 0)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile0);
               if(vartype == 1)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile1);
               if(vartype == 2)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile2);
               if(vartype == 3)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile3);
               if(vartype == 4)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile4);
               if(vartype == 8)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile8);
               if(vartype == 9)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile9);
               if(vartype == 10)
                 fwrite(&Gtot[l], sizeof(double), 1, outfile10);
               }
             fwrite(&Gtot[l], sizeof(double), 1, outfile);
             if(l % 1000 == 0)
               printf("vartype = %d, Gtot[l] = %lg\n", vartype, Gtot[l]);
             l++;
             }
      if(vartype == 0) printf("masstot_alt = %lg\n", masstot_alt/1.9891e33);
      if(vartype == 1) printf("xmomtot_alt = %lg\n", xmomtot_alt);
      if(vartype == 2) printf("ymomtot_alt = %lg\n", ymomtot_alt);
      if(vartype == 3) printf("zmomtot_alt = %lg\n", zmomtot_alt);
      if(vartype == 4) printf("etot_alt = %lg\n", etot_alt);
      fclose(outfile);
      fclose(outfile0);
      fclose(outfile1);
      fclose(outfile2);
      fclose(outfile3);
      fclose(outfile4);
      fclose(outfile8);
      fclose(outfile9);
      fclose(outfile10);
      }

  }

  free(Gdum);
  free(Grho);
  free(Gtot);
  free(Gnpart);
  free(Gnpart_tot);
  free(Grho_tot);

  return(Ngas);
}





/* this routine allocates the memory for the 
 * particle data.
 */
int allocate_memory(void)
{
  printf("allocating memory...\n");

  if(!(P=(struct particle_data *) malloc(NumPart*sizeof(struct particle_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
//  P--;   /* start with offset 1 */

 /* 
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  */
//  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
  return(0);
}

/*
int allocate_memory3(void)
{
  printf("allocating memory...\n");

  if(!(Gtot=(struct grid_data_tot *) malloc(pow(ref_lev,3)*sizeof(struct grid_data_tot))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  printf("allocating memory3...done\n");
  return(0);
}
*/

/* This routine brings the particles back into
 * the order of their ID's.
 * NOTE: The routine only works if the ID's cover
 * the range from 1 to NumPart !
 * In other cases, one has to use more general
 * sorting routines.
 */
/*
int reordering(void)
{
  int i,j;
  int idsource, idsave, dest;
  struct particle_data psave, psource;


  printf("reordering....\n");

  for(i=1; i<=NumPart; i++)
    {
      if(Id[i] != i)
	{
	  psource= P[i];
	  idsource=Id[i];
	  dest=Id[i];

	  do
	    {
	      psave= P[dest];
	      idsave=Id[dest];

	      P[dest]= psource;
	      Id[dest]= idsource;
	      
	      if(dest == i) 
		break;

	      psource= psave;
	      idsource=idsave;

	      dest=idsource;
	    }
	  while(1);
	}
    }

  printf("done.\n");

  Id++;   
  free(Id);

  printf("space for particle ID freed\n");
  return(0);
}
*/

double calc_kernel_tsc(int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size, double grid_size_half)
{
  double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
  double radx, rady, radz;

  rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
  rad = pow(rad,0.5);

  ratio = rad/P[n].hsm_phys;

  radx = fabs(P[n].disx - xgrid);
  rady = fabs(P[n].disy - ygrid);
  radz = fabs(P[n].disz - zgrid);

//triangular-shaped cloud kernel
  grid_size_simu = grid_size/1.e3/Time*(hubble_param); 

  xratio = fabs((P[n].disx-xgrid)/grid_size);
  yratio = fabs((P[n].disy-ygrid)/grid_size);
  zratio = fabs((P[n].disz-zgrid)/grid_size);

  Wx = Wy = Wz = 0;

  if(ratio >= 1. && rad < grid_size_half)
    //xratio = yratio = zratio = rad/(grid_size_half);
    xratio = yratio = zratio = 0;

  if(xratio < 0.5)
    Wx = 0.75 - pow(xratio, 2);
  if(xratio < 1.5 && xratio >= 0.5)
    Wx = 0.5*pow(1.5 - xratio,2);

  if(yratio < 0.5)
    Wy = 0.75 - pow(yratio, 2);
  if(yratio < 1.5 && yratio >= 0.5)
    Wy = 0.5*pow(1.5 - yratio,2);

  if(zratio < 0.5)
    Wz = 0.75 - pow(zratio, 2);
  if(zratio < 1.5 && zratio >= 0.5)
    Wz = 0.5*pow(1.5 - zratio,2);

  kernel = pow(1./grid_size_simu,3)*Wx*Wy*Wz; 

  //if(radx > grid_size_half || rady > grid_size_half || radz > grid_size_half)
  //  kernel = 0;

  return(kernel);
}


double calc_kernel_spline(int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size, double grid_size_half)
{
  double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
  double grid_size_max, hfac=2.0;
  double radx, rady, radz;

  grid_size_max = pow(grid_size_half*grid_size_half + grid_size_half*grid_size_half + grid_size_half*grid_size_half,0.5);

  rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
  rad = pow(rad,0.5);

  radx = fabs(P[n].disx - xgrid);
  rady = fabs(P[n].disy - ygrid);
  radz = fabs(P[n].disz - zgrid);

  ratio = rad/P[n].hsm_phys;


  if(ratio >= hfac && radx < hfac*grid_size_half && rady < hfac*grid_size_half && radz < hfac*grid_size_half)
     {
     ratio = rad/(grid_size_half);
     if(hsm < grid_size_half/1.e3/Time*(hubble_param))
       hsm = grid_size_half/1.e3/Time*(hubble_param);
     printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, P[n].hsm_phys);
     }


//Gadget kernel

/*
  if(ratio <= 0.5)
    kernel = fac*(8./3.14159/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
  if(ratio > 0.5 && ratio <= 1.)
    kernel = (8./3.14159/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
  if(ratio > 1.)
    kernel = 0.;
*/


  if(ratio < 1)
    kernel = (1./3.14159/pow(hsm,3)) * (1. - 1.5*pow(ratio,2) + 0.75*pow(ratio,3));
  if(ratio >= 1 && ratio <= 2)
    kernel = (1./3.14159/pow(hsm,3)) * (0.25*pow(2. - ratio, 3));
  if(ratio > 2)
    kernel = 0;


 //if(radx > grid_size_half || rady > grid_size_half || radz > grid_size_half)
 //   kernel = 0;

// if(fabs(P[n].disy) > width/2.0 || fabs(P[n].disx) > width/2.0 || fabs(P[n].disz) > width/2.0)
//    kernel = 0;

  return(kernel);
}

double calc_kernel_spline_vpot(int n, double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double grid_size_half)
{
  double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
  double hsm_phys, fac=1.0;
  double radx, rady, radz;

  rad = (x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z);
  rad = pow(rad,0.5);

  radx = fabs(x_part - x);
  rady = fabs(y_part - y);
  radz = fabs(z_part - z);

  ratio = rad/hsm;

  if(ratio >= 1. && radx < grid_size_half && rady < grid_size_half && radz < grid_size_half)
     {
     ratio = rad/(grid_size_half);
     if(hsm < grid_size_half)
       hsm = grid_size_half;
     printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, P[n].hsm);
     }

  if(ratio <= 0.5)
    kernel = (8./PI/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
  if(ratio > 0.5 && ratio <= 1.)
    kernel = (8./PI/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
  if(ratio > 1.)
    kernel = 0.;

  return(kernel);
}

double del_kernel(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp)
{
  double dkernel;
  double rad, del_kernel, ratio, drdu;

  rad = (x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z);
  rad = pow(rad,0.5);

  ratio = rad/hsm;
  drdu = 0.5 * pow((x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z), -0.5 ) * 2. * rad_comp;

  dkernel = 0;

  if(ratio <= 0.5)
    dkernel = (8./3.14159/pow(hsm,3)) * ( -12*pow(ratio,1)*(1./hsm)*drdu + 18*pow(ratio,2)*(1./hsm)*drdu );
  if(ratio > 0.5 && ratio <= 1.)
    dkernel = (8./3.14159/pow(hsm,3)) * 6.*pow(1. - ratio, 2)*(-1./hsm)*drdu;
  if(ratio > 1.)
    dkernel = 0.;

  return(dkernel);
}

double del_rho(double x_part, double y_part, double z_part, double x, double y, double z, double hsm, double rad_comp)
{
  double dkernel;
  double rad, del_kernel, ratio, drdu;

  rad = (x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z);
  rad = pow(rad,0.5);

  ratio = rad/hsm;
  drdu = 0.5 * pow((x_part - x)*(x_part - x) + (y_part - y)*(y_part - y)+ (z_part - z)*(z_part - z), -0.5 ) * 2. * rad_comp;

  dkernel = 0;

  if(ratio <= 0.5)
    dkernel = -3.*(8./3.14159/pow(hsm,4)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3))
             + (8./3.14159/pow(hsm,3)) * ( - 12.*pow(ratio,1)*(-rad/hsm/hsm) +  18*pow(ratio,2)*(-rad/hsm/hsm) );
  if(ratio > 0.5 && ratio <= 1.)
    dkernel = -3.*(8./3.14159/pow(hsm,4)) * 2.*pow(1. - ratio, 3) + (8./3.14159/pow(hsm,3)) * 6.*pow(1. - ratio, 2)*(rad/hsm/hsm);
  if(ratio > 1.)
    dkernel = 0.;

  return(dkernel);
}

double curl(int dir, double x1, double y1, double z1, double x2, double y2, double z2)
{

double xnew, ynew, znew;

xnew = y1*z2 - z1*y2;
ynew = z1*x2 - x1*z2;
znew = x1*y2 - y1*x2;

if(dir == 0) return(xnew);
else if(dir == 1) return(ynew);
else if(dir == 2) return(znew);
else return(0);
}

