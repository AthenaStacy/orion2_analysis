//Reads in Gadget2 data, converts to an Orion2 unifrom grid and outputs file to be read in during Orion2 initialization

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#define MAXREF 20
#define readB 1

#define pathname_in "/work/00863/minerva/bin_map_smallt"

#define pathname_out "/work/00863/minerva/orion"

#define ref_lev 128 //number of grid cells along each x-y-z direction
#define width  0.2 //size of box in pc to be "cut out" of Gadget2 simulation box and mapped onto the Orion2 grid 

//#define ref_lev 64
//#define width  10.0

#define hubble_param 0.7

int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize);
int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double delvx, double delvy, double delvz, int NumGas, int proc);
int reordering(void);
int unit_conversion(void);
int allocate_memory(void);
int allocate_memory2(void);
int allocate_memory3(void);
double calc_kernel_spline(int n, double x, double y, double z, double hsm, double grid_size, double grid_size_half);
double calc_kernel_tsc(int n, double x, double y, double z, double hsm, double grid_size, double grid_size_half);

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

struct particle_data 
{
  double  Pos[3];
  double  Vel[3];
  double  Mass;
  int    Type;
  double disx, disy, disz;
  double  Rho, U, Pres, nh, Density, hsm, hsm_phys;
  //double Temp, sink;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII gam;
  double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double nh_test, Bfieldx, Bfieldy, Bfieldz;
  double dummy;
} *P;

int *Id;


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
  int  j, n, type, snapshot_number, files, Ngas, Ngas2, random;
  double x,y,z,x1,y1, error, nh, nhmax, ref_lev_doub, grid_size, grid_size_half;
  double delx, dely, delz, delvx, delvy, delvz, boxsize, dis, disx, disy, disz, disAU;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;
  FILE *outfile, *infile;
  int npes, myrank, ierr;

  //For NON-cosmological runs
  //Time = 1.0;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  sprintf(path_in, pathname_in);
  sprintf(path_out, pathname_out);

  //sprintf(basename, "bin_HR10");
  sprintf(basename, "bin_HR10_map");
  //sprintf(basename, "iso_map");
  
  sprintf(basename2, "bin_HR10_map");

  //snapshot_number = 0;
  snapshot_number = 514;
  //snapshot_number = 519;
  files=1;     

  boxsize = 140.0;
  delx = dely = delz = 0;

  ref_lev_doub = (double)ref_lev;
  grid_size = width/ref_lev_doub;
  grid_size_half = grid_size/2.0;

  sprintf(input_fname, "%s/%s_%04d", path_in, basename, snapshot_number);

  ncount_doub = xCOM = yCOM = zCOM = delvx = delvy = delvz = vxCOM = vyCOM = vzCOM = 0.;

  //read in Gadget2 file
  Ngas = read_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  unit_conversion();

  //read in density and Bfield values from analytic calculation
  if(readB == 1)
    {
    sprintf(input_fname2, "%s/%s_bfield_%04d", path_out, basename2, snapshot_number);
    if(!(infile=fopen(input_fname2,"r")))
      {
        printf("can't open file `%s`\n",input_fname2);
        exit(0);
      }

    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldx, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldy, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].Bfieldz, sizeof(double), 1, infile);
    for(n=0;n<Ngas;n++)
        fread(&P[n].nh_test, sizeof(double), 1, infile); 
    fclose(infile);
    for(n=0;n<Ngas;n++)
      {
      error = 2. * fabs(P[n].Rho - P[n].nh_test) / (P[n].Rho + P[n].nh_test);
      if(myrank == 0 && error > 0.5)
        printf("ID = %d, rho = %lg, rho_test = %lg, error = %lg \n", Id[n], P[n].Rho, P[n].nh_test, error);
      }
    }


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
        printf("rho = %lg, rho_test = %lg, Bfieldx = %lg\n", P[n].Rho, P[n].nh_test, P[n].Bfieldx); 
     }

  for(n=0;n < Ngas; n++)
     {
     dis = pow(((P[n].Pos[0]-delx)*(P[n].Pos[0]-delx) + (P[n].Pos[1]-dely)*(P[n].Pos[1]-dely) + (P[n].Pos[2]-delz)*(P[n].Pos[2]-delz)), 0.5);
     dis=dis*1.e3*Time/(hubble_param);
     disAU=dis*206264.806;

     disx = fabs((P[n].Pos[0]-delx))*1.e3*Time/(hubble_param);
     disy = fabs((P[n].Pos[1]-dely))*1.e3*Time/(hubble_param);
     disz = fabs((P[n].Pos[2]-delz))*1.e3*Time/(hubble_param);
     //if(dis < width)
     if(disx < width && disy < width && disz < width)
     //if(P[n].nh > nhmax/10.0) 
       {
       vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
       vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
       vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
       ncount_doub = ncount_doub + P[n].Mass;
       }
     }


   delvx = vxCOM/ncount_doub;
   delvy = vyCOM/ncount_doub;
   delvz = vzCOM/ncount_doub;
   //delvx = delvy = delvz = 0.0;

  printf("delvx = %lg, delvy = %lg, delvz = %lg\n", delvx, delvy, delvz);

  sprintf(output_fname, "%s/gadget2orion", path_out);

  if(myrank == 0)
    {
    outfile=fopen(output_fname,"w");
    fclose(outfile);
    }


  Ngas2 = write_snapshot(input_fname, files, output_fname, delx, dely, delz, delvx, delvy, delvz, Ngas, myrank);

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
  UnitMass_in_g= 1.989e43;         /*  code mass unit in g/h */
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
       //h2frac=2.0*P[i].H2I;
       h2frac = 2.0*1.e-3;

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

	  gamma= 5.0/3.0;
          //gamma = P[i].gam;	 

	  /* get temperature in Kelvin */

	  temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Density * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
	  P[i].nh= P[i].Rho / MeanWeight;
          P[i].Pres = (gamma-1)*u*P[i].Rho;

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
              fread(&Id[pc_new], sizeof(int), 1, fd);
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
		
/*	
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
*/			
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


int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double delvx, double delvy, double delvz, int NumGas, int myrank)
{
  FILE *outfile;
  char   buf[200];
  int i, nmin, nmax, j, k, l, n, send_int, send_size, tot_proc=1, tot_proc_sum;
  int imin, imax, jmin, jmax, kmin, kmax;
  int NumGas_per_proc, vartype, varnum=8;
  double ref_lev_doub, i_doub, vel_conv, dens_conv, mass_conv, kernel_conv, kernel_phys, mass, rho, hsm; 
  double fac, kernel, rad, ratio, xgrid, ygrid, zgrid, grid_size, grid_size_half;
  double grid_arr[ref_lev], masstot, hfac = 2.0;
  double *Gdum, *Grho, *Grho_tot, *Gvelx, *Gvely, *Gvelz, *Gpres, *GH2I, *GHDI, *GHII, *Gtot, *Gnpart, *Gnpart_tot;
  MPI_Status status;

  if(readB == 1)
    varnum=12;

  //For NON-cosmological runs
  //Time = 1.0;

  Gdum =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
  Grho =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
  Gtot =  (double *) malloc(pow(ref_lev,3) * sizeof(double));
  Gnpart =(double *) malloc(pow(ref_lev,3) * sizeof(double));
  Gnpart_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));
  Grho_tot =(double *) malloc(pow(ref_lev,3) * sizeof(double));

  printf("myrank = %d, Memory allocation done\n", myrank);

  ref_lev_doub = (double)ref_lev;
  send_size = pow(ref_lev,3);
  vel_conv = 1.e5*pow(Time,0.5);
  dens_conv = P[100].Rho/P[100].Density;
  mass_conv = (1.e10/hubble_param)*1.989e33;
  kernel_conv =  pow(3.0857e18*1.e3*Time/(hubble_param),-3); 
  grid_size = width/ref_lev_doub;
  grid_size_half = grid_size/2.0;
  
  MPI_Allreduce(&tot_proc, &tot_proc_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  printf("myrank = %d, tot_proc_sum = %d\n", myrank, tot_proc_sum);

//subtract out COM velocities and shift particle positions to center of Orion2 box
      for(n=0;n<NumGas;n++)
        {
        P[n].hsm_phys = P[n].hsm*1.e3*Time/hubble_param; 
        P[n].disx = 1.e3*Time/(hubble_param)*(P[n].Pos[0] - delx) - grid_size_half;
        P[n].disy = 1.e3*Time/(hubble_param)*(P[n].Pos[1] - dely) - grid_size_half;
        P[n].disz = 1.e3*Time/(hubble_param)*(P[n].Pos[2] - delz) - grid_size_half;
        P[n].Vel[0] = P[n].Vel[0] - delvx;
        P[n].Vel[1] = P[n].Vel[1] - delvy;
        P[n].Vel[2] = P[n].Vel[2] - delvz;
        if(fabs(fabs(P[n].disx) - 0.0*P[n].hsm_phys) < width/2.0 && fabs(fabs(P[n].disy) - 0.0*P[n].hsm_phys) < width/2.0 && fabs(fabs(P[n].disz) - 0.0*P[n].hsm_phys) < width/2.0)
          masstot = masstot + P[n].Mass;
        }

      for(i=0;i<ref_lev;i++)
        {
        i_doub = (double)i;
        grid_arr[i] = (width * i_doub/ref_lev_doub) - width/2.0;
        }	

      NumGas_per_proc = (int) (NumGas/tot_proc_sum); 
      nmin = NumGas_per_proc * (myrank);
      nmax = NumGas_per_proc * (myrank+1) -1;

/*
      if(myrank == tot_proc_sum - 1)
        nmax = NumGas-1;
*/

  printf("myrank = %d, nmin = %d, nmax = %d, ref_lev = %d, ref_lev_doub = %lg\n", myrank, nmin, nmax, ref_lev, ref_lev_doub);

// For each processor, find particles whose smoothing lengths overlap with grid cell centers, and add their weighted values to the grid
      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             Grho[l] = 0.0;
             Gtot[l] = 0.0;
             Gnpart[l] = 0.0;
             Gnpart_tot[l] = 0.0;
             Grho_tot[l] = 0.0;
             l++;
             }

    for(vartype=0;vartype<varnum;vartype++)
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

        imin = (int)((P[n].disx - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        imax = (int)((P[n].disx + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        jmin = (int)((P[n].disy - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        jmax = (int)((P[n].disy + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        kmin = (int)((P[n].disz - hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);
        kmax = (int)((P[n].disz + hfac*P[n].hsm_phys + width/2.0)/width*ref_lev_doub);

        if(imin < 0)
          imin = 0;
        if(imax > ref_lev-1)
          imax = ref_lev-1;
        if(jmin < 0)
          jmin = 0;
        if(jmax > ref_lev-1)
          jmax = ref_lev-1;
        if(kmin < 0)
          kmin = 0;
        if(kmax > ref_lev-1)
          kmax = ref_lev-1;

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

                 kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);               
                 //kernel = calc_kernel_tsc(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half); 
                 kernel_phys = kernel*kernel_conv;

                 rho = mass*kernel_phys;
                 fac = P[n].Mass/P[n].Density*kernel;
                 //fac = 1./200.;

                 if(rho > 0)
                    {
                    Grho[l]  = Grho[l] + rho;
                    //Grho[l]  = Grho[l] + P[n].Rho*fac;               
                    Gnpart[l] = Gnpart[l] + 1.0;
                    }
                  }


               if(vartype > 0 && Gnpart[l] > 0)
                 {
                 xgrid = grid_arr[i]+grid_size_half;
                 ygrid = grid_arr[j]+grid_size_half;
                 zgrid = grid_arr[k]+grid_size_half;

                 kernel = calc_kernel_spline(n, xgrid, ygrid, zgrid, hsm, grid_size, grid_size_half);
                 kernel_phys = kernel*kernel_conv;

                 rho = mass*kernel_phys;
                 fac = P[n].Mass/P[n].Density*kernel;

                  if(vartype == 1)
                    Gdum[l] = Gdum[l] + P[n].Vel[0]*rho; //do not yet convert velocity to cgs units 
                  if(vartype == 2)
                    Gdum[l] = Gdum[l] + P[n].Vel[1]*rho; 
                  if(vartype == 3)
                    Gdum[l] = Gdum[l] + P[n].Vel[2]*rho; 
                  if(vartype == 4)
                    Gdum[l] = Gdum[l] + P[n].Pres*fac;
                  if(vartype == 5)
                    Gdum[l] = Gdum[l] + P[n].H2I*fac;
                  if(vartype == 6)
                    Gdum[l] = Gdum[l] + P[n].HDI*fac;
                  if(vartype == 7)
                    Gdum[l] = Gdum[l] + P[n].HII*fac; 
                  if(vartype == 8)
                    Gdum[l] = Gdum[l] + P[n].Bfieldx*fac;
                  if(vartype == 9)
                    Gdum[l] = Gdum[l] + P[n].Bfieldy*fac;
                  if(vartype == 10)
                    Gdum[l] = Gdum[l] + P[n].Bfieldz*fac; 
                  if(vartype == 11)
                    Gdum[l] = Gdum[l] + P[n].nh_test*fac; 
                  }
               }

        }       

      printf("myrank = %d, main loop finished\n", myrank);

      MPI_Barrier(MPI_COMM_WORLD);

      if(vartype == 0)
        {
        MPI_Allreduce(&Gnpart[0], &Gnpart_tot[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Allreduce(&Grho[0], &Grho_tot[0], send_size, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        }

//////////write mapped data to file for Orion2!!!!!!!!!!!

     //write out density for each grid cell
     if(vartype == 0)
       MPI_Reduce(&Grho[0], &Gtot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

     if(vartype > 0) 
        MPI_Reduce(&Gdum[0], &Gtot[0], send_size, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

      printf("Rank 0 has the data now, myrank=%d\n", myrank);  fflush(stdout);
      MPI_Barrier(MPI_COMM_WORLD);
      if(myrank == 0)
      {
      outfile=fopen(outname,"a");
      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             if(vartype == 1 || vartype == 2 || vartype == 3)   //write out velocitiesi
               Gtot[l] = (Gtot[l]/Grho_tot[l])*vel_conv;
             fwrite(&Gtot[l], sizeof(double), 1, outfile);
             if(l % 1000 == 0)
               printf("vartype = %d, Gtot[l] = %lg\n", vartype, Gtot[l]);
             l++;
             }
      fclose(outfile);
      }

  }

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

  
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
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


double calc_kernel_tsc(int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size, double grid_size_half)
{
  double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;

  rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
  rad = pow(rad,0.5);

  ratio = rad/P[n].hsm_phys;

//triangular-shaped cloud kernel
  grid_size_simu = grid_size/1.e3/Time*(hubble_param); 

  xratio = fabs((P[n].disx-xgrid)/grid_size);
  yratio = fabs((P[n].disy-ygrid)/grid_size);
  zratio = fabs((P[n].disz-zgrid)/grid_size);

  Wx = Wy = Wz = 0;

  if(ratio >= 1. && rad < grid_size_half)
    xratio = yratio = zratio = rad/(grid_size_half);

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

  return(kernel);
}


double calc_kernel_spline(int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size, double grid_size_half)
{
  double rad, kernel, ratio, xratio, yratio, zratio, Wx, Wy, Wz, grid_size_simu;
  double grid_size_max, fac=1.0;

  grid_size_max = pow(grid_size_half*grid_size_half + grid_size_half*grid_size_half + grid_size_half*grid_size_half,0.5);

  rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
  rad = pow(rad,0.5);

  ratio = rad/P[n].hsm_phys;


  if(ratio >= 1. && rad < grid_size_half)
     {
     ratio = rad/(grid_size_half);
     if(hsm < grid_size_half/1.e3/Time*(hubble_param))
       hsm = grid_size_half/1.e3/Time*(hubble_param);
     printf("Dense particle1! nh = %lg, hsm = %lg\n", P[n].nh, P[n].hsm_phys);
     }


//Gadget kernel

  if(ratio <= 0.5)
    kernel = fac*(8./3.14159/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
  if(ratio > 0.5 && ratio <= 1.)
    kernel = (8./3.14159/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
  if(ratio > 1.)
    kernel = 0.;

/*
  if(ratio < 1)
    kernel = (1./3.14159/pow(hsm,3)) * (1. - 1.5*pow(ratio,2) + 0.75*pow(ratio,3));
  if(ratio >= 1 && ratio <= 2)
    kernel = (1./3.14159/pow(hsm,3)) * (0.25*pow(2. - ratio, 3));
  if(ratio > 2)
    kernel = 0;
*/
  return(kernel);
}

