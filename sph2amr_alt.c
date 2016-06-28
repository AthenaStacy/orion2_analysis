//Adds gamma values of the particles to a new snapshot file

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"
#define MAXREF 20

#define ref_lev 128
#define dmax  1.0

int read_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double  boxsize);
int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double delvx, double delvy, double delvz, int NumGas, int proc);
int reordering(void);
int unit_conversion(void);
int allocate_memory(void);
int allocate_memory2(void);
double calc_kernel(int i, double x, double y, double z, double grid_size, double grid_size_half);

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
  int flag;
  double disx, disy, disz;
  double  Rho, U, Temp, Pres, nh, Density, hsm, hsm_phys;
  //double  elec, HI, HII, HeI, HeII,  HeIII, H2I, H2II, HM, hsm, DI, DII, HDI, DM, HDII, FosHII, sink, gam;
  double H2I, HII, DII, HDI, HeII, HeIII, gam, sink;
  double dummy;
} *P;

int *Id;

struct grid_data
{
  double rho; 
  double pres;
  double velx, vely, velz;
} *G;


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
  char path[200], input_fname[200], output_fname[200], basename[200], basenameout[200];
  int  j, n, type, snapshot_number, files, Ngas, Ngas2, random, ncount, ncounthalo1, ncount2;
  float x,y,z,x1,y1, nh, nhmax;
  double delx, dely, delz, delvx, delvy, delvz, boxsize;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;
  FILE *outfile;
  int npes, myrank, ierr;

  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  sprintf(path, "/work/00863/minerva");
  //sprintf(basename, "midres_small");

  //snapshot_number=300;                   /* number of snapshot */
  //snapshot_number=500;
  //snapshot_number=800;  
  //snapshot_number=910;

  //sprintf(basename, "midres");
  //snapshot_number=250;                   /* number of snapshot */

  sprintf(basename, "bin_HR10");
  snapshot_number = 1;

  files=1;     

  boxsize = 140.0;
  delx = dely = delz = 0;

  sprintf(input_fname, "%s/%s_%03d", path, basename, snapshot_number);
  Ngas = read_snapshot(input_fname, files, output_fname, delx, dely, delz, boxsize);

  printf("hello line 140\n");

  unit_conversion();

  nhmax = 0;
  for(n=1;n<=Ngas;n++) 
     {
     nh = P[n].nh;
     if(nh > nhmax && P[n].sink > -1)
        {
        nhmax = nh;
        delx = P[n].Pos[0];  //del's based upon location of maximum density
        dely = P[n].Pos[1];
        delz = P[n].Pos[2];
        }
     }

  ncount_doub = xCOM = yCOM = zCOM = delvx = delvy = delvz = vxCOM = vyCOM = vzCOM = 0.;
  for(n=1;n <= Ngas; n++)
     {
     nh = P[n].nh;
     if(nh > nhmax/10.0)
       {
       vxCOM = vxCOM + P[n].Vel[0]*P[n].Mass;
       vyCOM = vyCOM + P[n].Vel[1]*P[n].Mass;
       vzCOM = vzCOM + P[n].Vel[2]*P[n].Mass;
       xCOM = xCOM + P[n].Pos[0]*P[n].Mass;
       yCOM = yCOM + P[n].Pos[1]*P[n].Mass;
       zCOM = zCOM + P[n].Pos[2]*P[n].Mass;
       ncount_doub = ncount_doub + P[n].Mass;
       }
     }

   delvx = vxCOM/ncount_doub;
   delvy = vyCOM/ncount_doub;
   delvz = vzCOM/ncount_doub;

  printf("hello line 154, nhmax = %lg, delx = %lg, dely = %lg, delz = %lg\n", nhmax, delx, dely, delz);

  sprintf(output_fname, "/work/00863/minerva/orion/gadget2orion");
  Ngas2 = write_snapshot(input_fname, files, output_fname, delx, dely, delz, delvx, delvy, delvz, Ngas, myrank);

  ncount = 0;
  ncount2 = 0;

  printf("ncount = %d.\n", ncount);
  printf("ncount2 = %d.\n", ncount2);
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
  double h2frac, muh2, muh2in;

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


  for(i=1; i<=NumPart; i++)
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

	  //gamma= 5.0/3.0;
          gamma = P[i].gam;	 

	  /* get temperature in Kelvin */

	  P[i].Temp= MeanWeight/BOLTZMANN * (gamma-1) * u;
	  P[i].Rho= P[i].Density * UnitDensity_in_cgs * HubbleParam * HubbleParam * pow(1.e0+zred,3.e0);
	  P[i].nh= P[i].Rho / MeanWeight;
          P[i].Pres = (gamma-1)*u*P[i].Rho/gamma;

          if(i%100000 == 0)
            printf("Pres = %lg, pres = %lg\n ",P[i].Pres, P[i].nh*BOLTZMANN*P[i].Temp);

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
  int    t,n,off,pc,pc_new,pc_sph;
  int NumPart_new = 0, Ngas_new = 0;
  int Idnew, nnew=1;
  double *pos, massnew, hsmnew, x, y, z, nnew_doub=1.0;
  double randomx, randomy, randomz;

  pos= (double*)malloc(sizeof(double) * 3);
	
#define SKIP fread(&dummy, sizeof(dummy), 1, fd);

  for(i=0, pc=1; i<files; i++, pc=pc_new)
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
      fclose(fd);
  }

  Time= header1.time;
  zred= header1.redshift;
  printf("z= %6.2f \n",zred);
  printf("Time= %12.7e \n",Time);
  printf("L= %6.2f \n",header1.BoxSize);
  fflush(stdout);
  return(Ngas);
}


int write_snapshot(char *fname, int files, char *outname, double delx, double dely, double delz, double delvx, double delvy, double delvz, int NumGas, int myrank)
{
  FILE *outfile;
  char   buf[200];
  int    i,imin,imax, j,k,l,n,send_int,send_size,tot_proc=1, tot_proc_sum;
  int num_phys=5;
  double ref_lev_doub, i_doub, vel_conv, dens_conv, mass, hsm;  //grid size in pc
  double fac, kernel, rad, ratio, xgrid, ygrid, zgrid, grid_size, grid_size_half;
  double grid_arr[ref_lev], masstot;
  MPI_Status status;

  allocate_memory2();

  ref_lev_doub = (double)ref_lev;
  vel_conv = 1.e5*pow(Time,0.5);
  dens_conv = P[100].Rho/P[100].Density; 
  grid_size = dmax/ref_lev_doub;
  grid_size_half = grid_size/2.0;

//create grid!!
  
  MPI_Allreduce(&tot_proc, &tot_proc_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  printf("hello line 508, myrank = %d, tot_proc_sum = %d\n", myrank, tot_proc_sum);  fflush(stdout);

      for(n=1;n<=NumGas;n++)
        {
        P[n].hsm_phys = P[n].hsm*1.e3*header1.time/0.7; 
        P[n].disx = 1.e3*header1.time/(0.7)*(P[n].Pos[0] - delx);
        P[n].disy = 1.e3*header1.time/(0.7)*(P[n].Pos[1] - dely);
        P[n].disz = 1.e3*header1.time/(0.7)*(P[n].Pos[2] - delz);
        P[n].Vel[0] = P[n].Vel[0] - delvx;
        P[n].Vel[1] = P[n].Vel[1] - delvy;
        P[n].Vel[2] = P[n].Vel[2] - delvz;
        if(fabs(P[n].disx) < dmax/2.0 && fabs(P[n].disy) < dmax/2.0 && fabs(P[n].disz) < dmax/2.0)
          masstot = masstot + P[n].Mass;
        }

      if(myrank == 0)
        printf("masstot = %lg", masstot*1.e10/0.7);

      for(i=0;i<ref_lev;i++)
        {
        i_doub = (double)i;
        grid_arr[i] = (dmax * i_doub/ref_lev_doub) - dmax/2.0;
        }	

      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             G[l].rho = 0.0;
             G[l].velx = 0.0;
             G[l].vely = 0.0;
             G[l].velz = 0.0;
             G[l].pres = 0.0;
             l++;
             }

  printf("hello line 527\n");  fflush(stdout);

      l=myrank*pow(ref_lev,2);
      imin = ref_lev / tot_proc_sum * myrank;
      imax = ref_lev / tot_proc_sum * (myrank+1);
      for(i=imin;i<imax;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             xgrid = grid_arr[i]+grid_size_half;
             ygrid = grid_arr[j]+grid_size_half;
             zgrid = grid_arr[k]+grid_size_half;

             //printf("xgrid = %lg, ygrid = %lg, zgrid = %lg, grid_size=%lg\n", xgrid, ygrid, zgrid, grid_size);             
             
             for(n=1;n<=NumGas;n++)
               {
               mass = P[n].Mass;

               kernel = calc_kernel(n, xgrid, ygrid, zgrid, grid_size, grid_size_half);               

               fac = mass/P[n].Density*kernel;
               G[l].rho  = G[l].rho + mass*kernel*dens_conv;
               G[l].velx = G[l].velx + P[n].Vel[0]*vel_conv*fac; //convert velocity to cgs units
               G[l].vely = G[l].vely + P[n].Vel[1]*vel_conv*fac;
               G[l].velz = G[l].velz + P[n].Vel[2]*vel_conv*fac;
               G[l].pres = G[l].pres + P[n].Pres*fac;
               }

               if(i%10==0 && j%10==0 && k%10==0)
                 {
                 printf("xgrid = %lg, ygrid = %lg, zgrid = %lg, grid_size=%lg\n", xgrid, ygrid, zgrid, grid_size);
                 printf("n = %d, i=%d, j=%d, k=%d, l=%d, myrank=%d, rho=%lg, velx=%lg \n", n, i, j, k, l, myrank, G[l].rho, G[l].velx);
                 fflush(stdout);
                 }
 
             l++;
             }       

  printf("hello line 550\n");

//////////write new file!!!!!!!!!!!

      MPI_Barrier(MPI_COMM_WORLD);

      send_size = pow(ref_lev,3)/tot_proc_sum*num_phys;
      for(i=1;i<tot_proc_sum;i++)
        {
        send_int = i*pow(ref_lev,2);
        if(myrank == i)
          {
          MPI_Send(&G[send_int], send_size, MPI_DOUBLE, 0, myrank, MPI_COMM_WORLD);
          printf("myrank = %d, send_int = %d, send_size = %d, G[send_int].rho = %lg \n", myrank, send_int, send_size, G[send_int].rho);
          fflush(stdout);    
          }
        if(myrank == 0)
          MPI_Recv(&G[send_int], send_size, MPI_DOUBLE, i, i, MPI_COMM_WORLD, &status);
        } 

      MPI_Barrier(MPI_COMM_WORLD);

      printf("Rank 0 has the data now, myrank=%d\n", myrank);  fflush(stdout);

      if(myrank == 0)
      {
      outfile=fopen(outname,"w");
      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             fwrite(&G[l].rho, sizeof(double), 1, outfile);
             l++;
             if(l%1000 == 0)
               printf("l = %d, G[l].rho = %lg, G[l].pres = %lg \n", l, G[l].rho, G[l].pres);
             }
 
      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             fwrite(&G[l].velx, sizeof(double), 1, outfile);
             fwrite(&G[l].vely, sizeof(double), 1, outfile);
             fwrite(&G[l].velz, sizeof(double), 1, outfile);
             l++;
             }

      l=0;
      for(i=0;i<ref_lev;i++)
        for(j=0;j<ref_lev;j++)
          for(k=0;k<ref_lev;k++)
             {
             fwrite(&G[l].pres, sizeof(double), 1, outfile);
             l++;
             }
 
      fclose(outfile);
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
  
  P--;   /* start with offset 1 */

  
  if(!(Id=(int *) malloc(NumPart*sizeof(int))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  
  Id--;   /* start with offset 1 */

  printf("allocating memory...done\n");
  return(0);
}



int allocate_memory2(void)
{
  printf("allocating memory...\n");

  if(!(G=(struct grid_data *) malloc(pow(ref_lev,3)*sizeof(struct grid_data))))
    {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
    }
  printf("allocating memory2...done\n");
  return(0);
}



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


double calc_kernel(int n, double xgrid, double ygrid, double zgrid, double grid_size, double grid_size_half)
{
  double rad, mass, hsm, kernel, ratio;

  rad = (P[n].disx - xgrid)*(P[n].disx - xgrid) + (P[n].disy - ygrid)*(P[n].disy - ygrid)+ (P[n].disz - zgrid)*(P[n].disz - zgrid);
  rad = pow(rad,0.5);

  hsm = P[n].hsm;
  ratio = rad/P[n].hsm_phys;

  if(ratio >= 1. && (rad - grid_size_half) < P[n].hsm_phys && (rad - grid_size_half) > 0)
     {
     ratio = rad/(P[n].hsm_phys + grid_size_half);
     mass = P[n].Mass * pow(P[n].hsm_phys - rad + grid_size_half,3)/pow(P[n].hsm_phys,3);
     }

  if(ratio >= 1. && rad < grid_size_half)
     ratio = rad/(grid_size_half);

  if(ratio <= 0.5)
    kernel = (8./3.14159/pow(hsm,3)) * (1. - 6.*pow(ratio,2) + 6.*pow(ratio,3));
  if(ratio > 0.5 && ratio <= 1.)
    kernel = (8./3.14159/pow(hsm,3)) * 2.*pow(1. - ratio, 3);
  if(ratio > 1.)
    kernel = 0.;

  return(kernel);
}


