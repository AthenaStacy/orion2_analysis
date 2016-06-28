//Reads in grid data for Orion2 initialization, mapped from Gadget2 sim (gadget2orion).  
//Finds vector potential for b-field data
//Takes curl of vector potential to determine new b-field values
//Outputes new "gadget2orion" file with new divergence-free b-field values

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mpi.h"

#define wpot 1
#define MAXREF 20
#define readB 1
#define which_sim 3
#define print_vec 0
#define interp_edge 0
#define interp_cell 1

#define pi 3.14159265359

//#define gnum (200 + 0)
#define gnum (128 + 10)
//#define gnum (64 + 10)
//#define gnum 32
//#define gnum (32 + 10)
//#define width 1.0

struct grid_data
{
double rho[gnum][gnum][gnum];
double velx[gnum][gnum][gnum];
double vely[gnum][gnum][gnum];
double velz[gnum][gnum][gnum];
double eint[gnum][gnum][gnum];

double H2I[gnum][gnum][gnum];
double HDI[gnum][gnum][gnum];
double HII[gnum][gnum][gnum];

double Bfieldx[gnum][gnum][gnum];
double Bfieldy[gnum][gnum][gnum];
double Bfieldz[gnum][gnum][gnum];
double nh_test[gnum][gnum][gnum];

double Ax[gnum][gnum][gnum];
double Ay[gnum][gnum][gnum];
double Az[gnum][gnum][gnum];

double DivB[gnum][gnum][gnum];
double DivB_alt[gnum][gnum][gnum];
double DivB_shift[gnum][gnum][gnum]; //Bx(i+1,j,k) - Bx(i,j,k) + By(i,j+1,k) 

double Bx_new[gnum][gnum][gnum];
double By_new[gnum][gnum][gnum];
double Bz_new[gnum][gnum][gnum];
double DivB_new[gnum][gnum][gnum];

double Bx_new_face[gnum][gnum][gnum];
double By_new_face[gnum][gnum][gnum];
double Bz_new_face[gnum][gnum][gnum];

double grad_phi[gnum][gnum][gnum];
double grad_phi_face[gnum][gnum][gnum]; //grad_phi at cell faces (instead of cell centers)

double PosX[gnum][gnum][gnum];
double PosY[gnum][gnum][gnum];
double PosZ[gnum][gnum][gnum];
} *grid;

void read_grid(char *fname);
void write_grid(char *fname_out, char *fout_bx, char *fout_by, char *fout_bz);
int get_hilo(void);
void get_gphi(int n);
void get_gphi_alt(int n);
void get_Bnew(int n);
void get_Bnew_alt(int n);
void set_Bedge(int n);
void skip_check();

int gnum_file;
double width, DeltaX;

int igrid, jgrid, kgrid, i_lo, i_hi, j_lo, j_hi, k_lo, k_hi;
double gphi_xmin, gphi_xplus, gphi_ymin, gphi_yplus, gphi_zmin, gphi_zplus, gphi, gphi_check;
double Bx_new, By_new, Bz_new;
double Bx_left, Bx_right, By_left, By_right, Bz_left, Bz_right;

int myrank=0;

int main(int argc, char **argv)
{
  char path_in[200], path_out[200], input_fname[200], input_fname2[200], output_fname[200], output_bx[200], output_by[200], output_bz[200], basename[200], basename2[200], basenameout[200];
  int  j, type, snapshot_number, files, l;
  double x,y,z,x1,y1, error, nh, nhmax, ref_lev_doub, grid_size, grid_size_half;
  double delx, dely, delz, delvx, delvy, delvz, boxsize, dis, disx, disy, disz, disAU;
  double xCOM, yCOM, zCOM, vxCOM, vyCOM, vzCOM, ncount_doub;
  double pc_to_cm = 3.08567758e18;
  int npes, ierr, tot_proc=1, tot_proc_sum=1, gmin, gmax;
  int n=0, num=1;


  ierr = MPI_Init(&argc, &argv);
  ierr = MPI_Comm_size(MPI_COMM_WORLD, &npes);
  ierr = MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  MPI_Allreduce(&tot_proc, &tot_proc_sum, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

  printf("myrank = %d, tot_proc_sum = %d\n", myrank, tot_proc_sum);

  if(myrank == 0) printf("gnum = %d, width = %lg\n", gnum, width);


  sprintf(path_in, "/work/00863/minerva/popiii_Bscope");
  sprintf(path_out, "/work/00863/minerva/popiii_Bscope");

  //sprintf(path_in, "/global/scratch/minerva/orion_Btest");
  //sprintf(path_out, "/global/scratch/minerva/orion_Btest");


  sprintf(input_fname, "%s/gadget2orion_test", path_in);
  sprintf(output_fname, "%s/gadget2orion_clean", path_out);

  snapshot_number=7000;
  sprintf(output_bx, "%s/gadget2bx_clean_%04d", path_out, snapshot_number);
  sprintf(output_by, "%s/gadget2by_clean_%04d", path_out, snapshot_number);
  sprintf(output_bz, "%s/gadget2bz_clean_%04d", path_out, snapshot_number);

  if(!(grid=(struct grid_data *) malloc( num * sizeof(struct grid_data))))
     {
      fprintf(stderr,"failed to allocate memory.\n");
      exit(0);
     }

if(myrank == 0) printf("Read in the file! %s\n", input_fname);

  read_grid(input_fname);

  DeltaX = width / (double) gnum;
  DeltaX = DeltaX * pc_to_cm; 
  double DeltaXh = DeltaX / 2.;

  if(myrank == 0) printf("gnum = %d, width = %lg, DeltaX = %lg\n", gnum, width, DeltaX);

////////////////////////////////////////////////////////////////////////////////////////////
//calculate vector potential
/////////////////////////////////////////////////////////////////////////////////////////////

if(myrank == 0) printf("Calculate DivB!\n");

double Ax_out, Ay_out, Az_out;

for(igrid=0;igrid<gnum;igrid++)
 for(jgrid=0;jgrid<gnum;jgrid++)
   for(kgrid=0;kgrid<gnum;kgrid++)
     {
     grid[n].PosX[igrid][jgrid][kgrid] = ((double) igrid) / ((double) gnum) * width * pc_to_cm - (width / 2.0 * pc_to_cm) + DeltaXh;
     grid[n].PosY[igrid][jgrid][kgrid] = ((double) jgrid) / ((double) gnum) * width * pc_to_cm - (width / 2.0 * pc_to_cm) + DeltaXh;
     grid[n].PosZ[igrid][jgrid][kgrid] = ((double) kgrid) / ((double) gnum) * width * pc_to_cm - (width / 2.0 * pc_to_cm) + DeltaXh;

     grid[n].DivB[igrid][jgrid][kgrid] = 0;

     grid[n].Ax[igrid][jgrid][kgrid] = grid[n].Ay[igrid][jgrid][kgrid] = grid[n].Ay[igrid][jgrid][kgrid] = 0;
     }

//original read-in-from-file B-fields are assumed to be at cell centers.  
double Bx, By, Bz;
double B_avg=0, divB_avg=0, divB_cur;
double B_new_avg=0, B2_avg, B_avg_prev=0;
int gogo;

for(igrid=0; igrid<gnum; igrid++)
 for(jgrid=0; jgrid<gnum; jgrid++)
   for(kgrid=0; kgrid<gnum; kgrid++)
     {
     gogo = get_hilo();

     if(gogo == 1) continue;

     Bx = grid[n].Bfieldx[igrid][jgrid][kgrid];
     By = grid[n].Bfieldy[igrid][jgrid][kgrid];
     Bz = grid[n].Bfieldz[igrid][jgrid][kgrid];

     Bx_left = (grid[n].Bfieldx[i_lo][jgrid][kgrid] + grid[n].Bfieldx[igrid][jgrid][kgrid])/2.;
     Bx_right = (grid[n].Bfieldx[i_hi][jgrid][kgrid] + grid[n].Bfieldx[igrid][jgrid][kgrid])/2.;
     By_left = (grid[n].Bfieldy[igrid][j_lo][kgrid] + grid[n].Bfieldy[igrid][jgrid][kgrid])/2.;
     By_right = (grid[n].Bfieldy[igrid][j_hi][kgrid] + grid[n].Bfieldy[igrid][jgrid][kgrid])/2.;
     Bz_left = (grid[n].Bfieldz[igrid][jgrid][k_lo] + grid[n].Bfieldz[igrid][jgrid][kgrid])/2.;
     Bz_right = (grid[n].Bfieldz[igrid][jgrid][k_hi] + grid[n].Bfieldz[igrid][jgrid][kgrid])/2.;

if(interp_edge == 1)
{
     if(igrid == 0 /*&& i_lo == 0*/) 
       Bx_left = Bx - (grid[n].Bfieldx[0][jgrid][kgrid] + grid[n].Bfieldx[1][jgrid][kgrid])/2;    
     if(igrid == gnum-1 /*&& i_hi == gnum-1*/)
       Bx_right = Bx + (grid[n].Bfieldx[gnum-1][jgrid][kgrid] + grid[n].Bfieldx[gnum-2][jgrid][kgrid])/2;

     if(jgrid == 0 /*&& j_lo == 0*/)
       By_left  = By_new - (grid[n].Bfieldy[igrid][0][kgrid] + grid[n].Bfieldy[igrid][1][kgrid])/2;
     if(jgrid == gnum-1 /*&& j_hi == gnum-1*/)
       By_right = By + (grid[n].Bfieldy[igrid][gnum-1][kgrid] + grid[n].Bfieldy[igrid][gnum-2][kgrid])/2;

     if(kgrid == 0 /*&& k_lo == 0*/)
       Bz_left  = Bz - (grid[n].Bfieldz[igrid][jgrid][0] + grid[n].Bfieldz[igrid][jgrid][1])/2;
     if(kgrid == gnum-1 /*&& k_hi == gnum-1*/)
       Bz_right = Bz + (grid[n].Bfieldz[igrid][jgrid][gnum-1] + grid[n].Bfieldz[igrid][jgrid][gnum-2])/2;
}


     grid[n].DivB[igrid][jgrid][kgrid] 
      =  ((Bx_right - Bx_left)  +  (By_right - By_left) + (Bz_right - Bz_left)) / DeltaX;          

     Bx_left = grid[n].Bfieldx[i_lo][jgrid][kgrid];
     Bx_right = grid[n].Bfieldx[i_hi][jgrid][kgrid];
     By_left = grid[n].Bfieldy[igrid][j_lo][kgrid];
     By_right = grid[n].Bfieldy[igrid][j_hi][kgrid];
     Bz_left = grid[n].Bfieldz[igrid][jgrid][k_lo];
     Bz_right = grid[n].Bfieldz[igrid][jgrid][k_hi];

     grid[n].DivB_alt[igrid][jgrid][kgrid] 
     =  ((Bx_right - Bx_left)  +  (By_right - By_left) + (Bz_right - Bz_left)) / DeltaX;    
     grid[n].DivB_alt[igrid][jgrid][kgrid] = grid[n].DivB_alt[igrid][jgrid][kgrid]/2.;

     grid[n].DivB_shift[igrid][jgrid][kgrid]
     = ((Bx_right - Bx)  +  (By_right - By) + (Bz_right - Bz)) / DeltaX; 

     divB_avg = divB_avg + fabs(grid[n].DivB[igrid][jgrid][kgrid]);
     B_avg = B_avg + pow(Bx*Bx + By*By + Bz*Bz, 0.5); 
     }

if(myrank == 0) printf("We have DivB! divB_avg = %lg B_avg = %lg\n", divB_avg/pow(gnum,3), B_avg/pow(gnum,3));

if(myrank == 0)
for(igrid=0;igrid<gnum;igrid++)
 for(jgrid=0;jgrid<gnum;jgrid++)
   for(kgrid=0;kgrid<gnum;kgrid++)
     {
     if(igrid == jgrid && igrid == kgrid && jgrid == kgrid && igrid % 8 == 0)
        printf("Bx = %lg, By = %lg, Bz = %lg, divB = %lg, \n", grid[n].Bfieldx[igrid][jgrid][kgrid], 
        grid[n].Bfieldy[igrid][jgrid][kgrid], grid[n].Bfieldz[igrid][jgrid][kgrid], 
        grid[n].DivB[igrid][jgrid][kgrid]);
     }


///////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////
//B = B_sol + B_comp
//B_comp = grad(phi)
//div(B) = del^2 phi = div(grad(phi))
//B_sol = B - grad(phi)

//make initial guess for grad_phi and B_sol...
for(igrid=0;igrid<gnum;igrid++)
 for(jgrid=0;jgrid<gnum;jgrid++)
   for(kgrid=0;kgrid<gnum;kgrid++)
     {
     grid[n].grad_phi[igrid][jgrid][kgrid] = grid[n].grad_phi[igrid][jgrid][kgrid] = 0;
     grid[n].Bx_new[igrid][jgrid][kgrid] = grid[n].Bfieldx[igrid][jgrid][kgrid];
     grid[n].By_new[igrid][jgrid][kgrid] = grid[n].Bfieldy[igrid][jgrid][kgrid];
     grid[n].Bz_new[igrid][jgrid][kgrid] = grid[n].Bfieldz[igrid][jgrid][kgrid];
     }

//div (B - grad_phi) = 0
//div (grad_phi) = divB

int i, k, ntot, n111, n011, n211, n101, n121, n110, n112; 
int n_iter=0, a[7][7], index[7], a_shift[9][9], index_shift[9];
double sum, dfac = 1.;

ntot = pow(gnum,3);

for(i=0; i<7; i++)
  for(j=0; j<7; j++)
    {
    if(i!=j) a[i][j] = 0;
    if(i==j) a[i][j] = -6;
    }
a[3][0] = 1; a[3][1] = 1; a[3][2] = 1; a[3][4] = 1; a[3][5] = 1; a[3][6] = 1;
a[0][3] = 1; a[1][3] = 1; a[2][3] = 1; a[4][3] = 1; a[5][3] = 1; a[6][3] = 1;

for(i=0; i<9; i++)
  for(j=0; j<9; j++)
    {
    a_shift[i][j] = 0;
    if(i==j) a_shift[i][j] = 1;
    }
a_shift[3][3] = -2; a_shift[4][4] = -2; a_shift[5][5] = -2; 

double *rhs, *rhs_shift, *xcalc, *xcalc_old, *xcalc_shift, *xcalc_shift_old;
rhs = (double *) malloc(ntot * sizeof(double));
rhs_shift = (double *) malloc(ntot * sizeof(double)); 
xcalc = (double *) malloc(ntot * sizeof(double));
xcalc_old = (double *) malloc(ntot * sizeof(double));
xcalc_shift = (double *) malloc(ntot * sizeof(double));
xcalc_shift_old = (double *) malloc(ntot * sizeof(double));

l=0;
for(igrid=0; igrid<gnum; igrid++)
 for(jgrid=0; jgrid<gnum; jgrid++)
   for(kgrid=0; kgrid<gnum; kgrid++)
     {
     //grid[n].DivB_alt[igrid][jgrid][kgrid] = 1.0/DeltaX;
     rhs[l] = grid[n].DivB[igrid][jgrid][kgrid]  * DeltaX ;
     rhs_shift[l] = grid[n].DivB_shift[igrid][jgrid][kgrid]  * DeltaX ;
     xcalc_old[l] = grid[n].grad_phi[igrid][jgrid][kgrid];
     l++;
     }

int num_per_proc;
num_per_proc = (int) (gnum/tot_proc_sum);
gmin = num_per_proc * myrank;
gmax = num_per_proc * (myrank+1);

if(myrank == 0) printf("num_per_proc = %d\n", num_per_proc);

if(myrank == tot_proc_sum - 1) gmax = gnum;
if(gmin > gnum) gmin = gmax = gnum;
if(gmax > gnum) gmax = gnum;

printf("myrank = %d, gmin = %d, gmax = %d \n", myrank, gmin, gmax);

//gmin=0; gmax=gnum;
//gmin = gmin-1; gmax = gmax+1;

do{
n_iter++;

l=0;
for(igrid=0; igrid<gnum; igrid=igrid+1)
 for(jgrid=0; jgrid<gnum; jgrid=jgrid+1)
   for(kgrid=0; kgrid<gnum; kgrid=kgrid+1)
      {
      xcalc[l] = 0;
      l++;
      }

for(igrid=gmin; igrid<gmax; igrid=igrid+1)
 for(jgrid=0; jgrid<gnum; jgrid=jgrid+1)
   for(kgrid=0; kgrid<gnum; kgrid=kgrid+1)
     {

     gogo = get_hilo();

     if(gogo == 1) continue;

     index[0] = n011 = kgrid + (jgrid)*(gnum) + (i_lo)*(gnum)*(gnum);
     index[1] = n211 = kgrid + (jgrid)*(gnum) + (i_hi)*(gnum)*(gnum);
     index[2] = n101 = kgrid + (j_lo)*(gnum) + (igrid)*(gnum)*(gnum);
     index[3] = n111 = kgrid + (jgrid)*(gnum) + (igrid)*(gnum)*(gnum);
     index[4] = n121 = kgrid + (j_hi)*(gnum) + (igrid)*(gnum)*(gnum);
     index[5] = n110 = k_lo + (jgrid)*(gnum) + (igrid)*(gnum)*(gnum); 
     index[6] = n112 = k_hi + (jgrid)*(gnum) + (igrid)*(gnum)*(gnum); 

       sum=0; i=3;
       for(j=0; j<7; j++)
         if(j != i) sum = sum + a[i][j] * xcalc_old[index[j]];

     xcalc[index[3]] = (1. / a[3][3]) * (rhs[index[3]] - sum); 

     if(n_iter%100==0 && igrid==0 && jgrid==0 && kgrid==0) 
        {
        gphi = xcalc[index[3]];
        gphi_check = xcalc[index[0]] + xcalc[index[1]] + xcalc[index[2]] + xcalc[index[4]] + xcalc[index[5]] + xcalc[index[6]] - grid[n].DivB[igrid][jgrid][kgrid] * DeltaX;
        gphi_check = gphi_check/6.0;
        printf("ind0=%d, ind1=%d, ind2=%d, n_iter = %d, gphi = %lg, gphi_check = %lg\n", index[0], index[1], index[2], n_iter, gphi, gphi_check);
        }

     }

 MPI_Allreduce(&xcalc[0], &xcalc_old[0], ntot, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

/*
l=0;
for(igrid=gmin; igrid<gmax; igrid=igrid+1)
 for(jgrid=gmin; jgrid<gmax; jgrid=jgrid+1)
   for(kgrid=gmin; kgrid<gmax; kgrid=kgrid+1)
      {
      index[0] = kgrid + (jgrid)*(gnum) + (igrid)*(gnum)*(gnum);
      xcalc_old[index[0]] = xcalc[index[0]];
      l++;
      }
*/
}while(n_iter < 50000);

l=0;
for(igrid=0; igrid<gnum; igrid=igrid+1)
 for(jgrid=0; jgrid<gnum; jgrid=jgrid+1)
   for(kgrid=0; kgrid<gnum; kgrid=kgrid+1)
      {
      grid[n].grad_phi[igrid][jgrid][kgrid] = xcalc_old[l];
      l++;
      }
///////////////////////////////////////////////////////////////////////////////////////////////////////////
//For an alternate stencil where cell-centered B-fields are used to calculate DivB, also calculate a set of face-centeres phi values
//

///////////////////////////////////////////////////////////////////////////////////////////////////////////

/*
for(igrid=0; igrid<gnum; igrid=igrid+1)
 for(jgrid=0; jgrid<gnum; jgrid=jgrid+1)
   for(kgrid=0; kgrid<gnum; kgrid=kgrid+1)
     {
     gogo = get_hilo();

     get_gphi_alt(n);

     gphi = grid[n].grad_phi[igrid][jgrid][kgrid];
     gphi_check = (gphi_xplus + gphi_xmin)/(1.0) + (gphi_yplus + gphi_ymin)/(1.0) + (gphi_zplus + gphi_zmin)/(1.0) - grid[n].DivB_alt[igrid][jgrid][kgrid] * DeltaX ;
     gphi_check = gphi_check/6.0;

     if(igrid%20==0 && jgrid%20==0 && kgrid%20==0) printf("igrid=%d, jgrid=%d, kgrid=%d, xcalc = %lg, gphi = %lg, gphi_check = %lg\n", igrid, jgrid, kgrid, xcalc[0], gphi, gphi_check);
     }
*/


double Bx_new_left, Bx_new_right, By_new_left, By_new_right, Bz_new_left, Bz_new_right;
for(igrid=0; igrid<gnum; igrid++)
 for(jgrid=0; jgrid<gnum; jgrid++)
   for(kgrid=0; kgrid<gnum; kgrid++)
     {
     gogo = get_hilo();

     get_gphi_alt(n);
     gphi = grid[n].grad_phi[igrid][jgrid][kgrid];

     Bx_right = (grid[n].Bfieldx[i_hi][jgrid][kgrid] + grid[n].Bfieldx[igrid][jgrid][kgrid])/2.;
     By_right = (grid[n].Bfieldy[igrid][j_hi][kgrid] + grid[n].Bfieldy[igrid][jgrid][kgrid])/2.;
     Bz_right = (grid[n].Bfieldz[igrid][jgrid][k_hi] + grid[n].Bfieldz[igrid][jgrid][kgrid])/2.;

     grid[n].Bx_new_face[igrid][jgrid][kgrid] = Bx_right - (gphi_xplus - gphi);
     grid[n].By_new_face[igrid][jgrid][kgrid] = By_right - (gphi_yplus - gphi);
     grid[n].Bz_new_face[igrid][jgrid][kgrid] = Bz_right - (gphi_zplus - gphi);

     grid[n].Bx_new[igrid][jgrid][kgrid] = grid[n].Bx_new_face[igrid][jgrid][kgrid];
     grid[n].By_new[igrid][jgrid][kgrid] = grid[n].By_new_face[igrid][jgrid][kgrid];
     grid[n].Bz_new[igrid][jgrid][kgrid] = grid[n].Bz_new_face[igrid][jgrid][kgrid];
     }


////////////////////////////////////////////////////////////////////////////////////////////////
if(myrank == 0 /*&& fabs(gphi_xplus) > 1*/) printf("divB = %lg, grid[n].grad_phi = %lg, gphi_xplus = %lg, gphi_xmin = %lg\n", grid[n].DivB[10][10][10], grid[n].grad_phi[10][10][10], gphi_xplus, gphi_xmin);
if(myrank == 0) printf("divB = %lg, grid[n].grad_phi = %lg, gphi_xplus = %lg, gphi_xmin = %lg\n", grid[n].DivB[gnum-1][gnum-1][gnum-1], grid[n].grad_phi[gnum-1][gnum-1][gnum-1], gphi_xplus, gphi_xmin);

if(myrank == 0) printf("A Bx_new = %lg, By_new = %lg, Bz_new = %lg\n", grid[n].Bx_new[10][10][10], grid[n].By_new[10][10][10], grid[n].Bz_new[10][10][10]);
if(myrank == 0) printf("Bx = %lg, By = %lg, Bz = %lg\n", grid[n].Bfieldx[10][10][10], grid[n].Bfieldy[10][10][10], grid[n].Bfieldz[10][10][10]);
if(myrank == 0) printf("A Bx_new = %lg, By_new = %lg, Bz_new = %lg\n", grid[n].Bx_new[gnum-1][gnum-1][gnum-1], grid[n].By_new[gnum-1][gnum-1][gnum-1], grid[n].Bz_new[gnum-1][gnum-1][gnum-1]);
if(myrank == 0) printf("Bx = %lg, By = %lg, Bz = %lg\n", grid[n].Bfieldx[gnum-1][gnum-1][gnum-1], grid[n].Bfieldy[gnum-1][gnum-1][gnum-1], grid[n].Bfieldz[gnum-1][gnum-1][gnum-1]);


//////////////////////////////////////////////////////////////////////////////////////////////////
//Check the divergence of the new field!
B_avg_prev = B_new_avg;
divB_avg = B2_avg = B_new_avg = 0;
for(igrid=0; igrid<gnum; igrid++)
 for(jgrid=0; jgrid<gnum; jgrid++)
   for(kgrid=0; kgrid<gnum; kgrid++)
     {
     gogo = get_hilo();     

     get_gphi_alt(n);

     if(interp_cell == 1) get_Bnew(n);  else get_Bnew_alt(n);

     Bx = grid[n].Bx_new[igrid][jgrid][kgrid] + (gphi_xplus - gphi_xmin)/dfac;
     By = grid[n].By_new[igrid][jgrid][kgrid] + (gphi_yplus - gphi_ymin)/dfac;
     Bz = grid[n].Bz_new[igrid][jgrid][kgrid] + (gphi_zplus - gphi_zmin)/dfac;

     divB_cur =  ((Bx_right - Bx_left) +  (By_right - By_left) + (Bz_right - Bz_left)) / DeltaX;

     divB_cur = ((grid[n].Bx_new_face[igrid][jgrid][kgrid] - grid[n].Bx_new_face[i_lo][jgrid][kgrid])
                +(grid[n].By_new_face[igrid][jgrid][kgrid] - grid[n].By_new_face[igrid][j_lo][kgrid])
                +(grid[n].Bz_new_face[igrid][jgrid][kgrid] - grid[n].Bz_new_face[igrid][jgrid][k_lo]))
                /DeltaX;   

/*
     divB_cur = ((grid[n].Bx_new[igrid][jgrid][kgrid] - grid[n].Bx_new[i_lo][jgrid][kgrid])
                +(grid[n].By_new[igrid][jgrid][kgrid] - grid[n].By_new[igrid][j_lo][kgrid])
                +(grid[n].Bz_new[igrid][jgrid][kgrid] - grid[n].Bz_new[igrid][jgrid][k_lo]))
                /DeltaX;    
*/

     if(interp_cell != 1) divB_cur = divB_cur/2.;
     grid[n].DivB_new[igrid][jgrid][kgrid] = divB_cur;

     //grid[n].Bx_new_face[igrid][jgrid][kgrid] = Bx_right;
     //grid[n].By_new_face[igrid][jgrid][kgrid] = By_right;
     //grid[n].Bz_new_face[igrid][jgrid][kgrid] = Bz_right;

     divB_avg = divB_avg + fabs(divB_cur);
     B_new_avg = B_new_avg + pow(Bx_new*Bx_new + By_new*By_new + Bz_new*Bz_new, 0.5);
     B2_avg = B2_avg + pow(Bx*Bx + By*By + Bz*Bz, 0.5);
     }

if(myrank == 0) printf("n_iter = %d New DivB! divB_avg = %lg, B2_avg = %lg, B_new_avg = %lg B_avg = %lg\n", n_iter, divB_avg/pow(gnum,3), B2_avg/pow(gnum,3), B_new_avg/pow(gnum,3), B_avg/pow(gnum,3));

//}while((fabs(B_avg_prev - B_new_avg) / B_new_avg > 1.e-10) && n_iter < 10000);
//while(fabs(B_avg - B_new_avg) / B_avg > 0.1);
//while(divB_avg > 1.e-40 && n_iter < 100);

if(myrank == 0)
for(igrid=0;igrid<gnum;igrid++)
 for(jgrid=0;jgrid<gnum;jgrid++)
   for(kgrid=0;kgrid<gnum;kgrid++)
     {
     if(igrid == jgrid && igrid == kgrid && jgrid == kgrid && igrid % 8 == 0)
        printf("Bx = %lg, Bx_new = %lg, By = %lg, By_new = %lg, Bz = %lg, Bz_new = %lg\n",
        grid[n].Bfieldx[igrid][jgrid][kgrid], grid[n].Bx_new[igrid][jgrid][kgrid], 
        grid[n].Bfieldy[igrid][jgrid][kgrid], grid[n].By_new[igrid][jgrid][kgrid],
        grid[n].Bfieldz[igrid][jgrid][kgrid], grid[n].Bz_new[igrid][jgrid][kgrid]);

     gogo = get_hilo();
     if(interp_cell == 1) get_Bnew(n);  else get_Bnew_alt(n);

    //if(grid[n].DivB_new[igrid][jgrid][kgrid] > 1.e-40)
     // printf("high DivB = %lg, i = %d, j = %d, k = %d, Bx_left = %lg, Bx_right = %lg, By_left = %lg, By_right = %lg, Bz_left = %lg, Bz_right = %lg, \n", grid[n].DivB_new[igrid][jgrid][kgrid], igrid, jgrid, kgrid, Bx_left, Bx_right, By_left, By_right, Bz_left, Bz_right);
     }

if(myrank == 0) printf("n_iter = %d New DivB! divB_avg = %lg, B2_avg = %lg, B_new_avg = %lg B_avg = %lg\n", n_iter, divB_avg/pow(gnum,3), B2_avg/pow(gnum,3), B_new_avg/pow(gnum,3), B_avg/pow(gnum,3));
///////////////////////////////////////////////////////////////////////////////////////////////

  if(myrank == 0) write_grid(output_fname, output_bx, output_by, output_bz);

}


void read_grid(char *fname)
{
  FILE *gadget_data; 
  char   buf[200];
  int i2, j2, k2;
  int num=1, n=0, nProc;

  sprintf(buf,"%s",fname);

  if(!(gadget_data=fopen(buf,"r")))
    {
      printf("can't open file `%s`\n",buf);
      exit(0);
    }


//read in header info

    fread(&nProc, sizeof(int),1,gadget_data);
    fread(&gnum_file, sizeof(int),1,gadget_data);
    fread(&width, sizeof(double), 1, gadget_data);

    if(gnum - gnum_file != 0) 
      printf("ARRAY SIZE MISMATCH!!, gnum = %d, gnum_file = %d, width = %lg \n", gnum, gnum_file, width); 

     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].rho[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("rho[0][0][0] = %lg\n", grid[n].rho[0][0][0]);

     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].velx[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("vex[0][0][0] = %lg\n", grid[n].velx[0][0][0]);

     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].vely[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("vely[0][0][0] = %lg\n", grid[n].vely[0][0][0]);

    for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].velz[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("velz[0][0][0] = %lg\n", grid[n].velz[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].eint[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("eint[0][0][0] = %lg\n", grid[n].eint[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].H2I[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("H2I[0][0][0] = %lg\n", grid[n].H2I[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].HDI[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("HDI[0][0][0] = %lg\n", grid[n].HDI[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].HII[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("HII[0][0][0] = %lg\n", grid[n].HII[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].Bfieldx[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("Bfieldx[0][0][0] = %lg\n", grid[n].Bfieldx[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].Bfieldy[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("Bfieldy[0][0][0] = %lg\n", grid[n].Bfieldy[0][0][0]);

     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].Bfieldz[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("Bfieldz[0][0][0] = %lg\n", grid[n].Bfieldz[0][0][0]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fread(&grid[n].nh_test[i2][j2][k2], sizeof(double), 1, gadget_data);
     printf("nh_test[0][0][0] = %lg\n", grid[n].nh_test[0][0][0]);
     printf("nh_test[%d][%d][%d] = %lg\n", i2, j2, k2,grid[n].nh_test[gnum-1][gnum-1][gnum-1]);

fclose(gadget_data);
}


void write_grid(char *fname_out, char *fout_bx, char *fout_by, char *fout_bz)
{
  FILE *gadget_out, *gadget_out_bx, *gadget_out_by, *gadget_out_bz;
  char   buf[200], buf1[200], buf2[200], buf3[200];
  int i2, j2, k2;
  int num=1, n=0;

  printf("Writing file!\n");

  sprintf(buf,"%s",fname_out);
  sprintf(buf1,"%s",fout_bx);
  sprintf(buf2,"%s",fout_by);
  sprintf(buf3,"%s",fout_bz);

  if(!(gadget_out=fopen(buf,"w")))
    {
      printf("can't open file `%s`\n",buf);
      exit(0);
    }

    gadget_out_bx=fopen(buf1,"w");
    gadget_out_by=fopen(buf2,"w");
    gadget_out_bz=fopen(buf3,"w");
 

    double width_val; 
    int gnum_val;
    width_val = width;
    gnum_val = gnum;

    fwrite(&myrank, sizeof(int), 1, gadget_out);
    fwrite(&gnum_val, sizeof(int), 1, gadget_out);
    fwrite(&width_val, sizeof(double), 1, gadget_out);

     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].rho[i2][j2][k2], sizeof(double), 1, gadget_out);


     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].velx[i2][j2][k2], sizeof(double), 1, gadget_out);


     for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].vely[i2][j2][k2], sizeof(double), 1, gadget_out);


    for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].velz[i2][j2][k2], sizeof(double), 1, gadget_out);


      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].eint[i2][j2][k2], sizeof(double), 1, gadget_out);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].H2I[i2][j2][k2], sizeof(double), 1, gadget_out);


      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].HDI[i2][j2][k2], sizeof(double), 1, gadget_out);


      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].HII[i2][j2][k2], sizeof(double), 1, gadget_out);

if(print_vec == 0)
{
     
      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            {
            fwrite(&grid[n].Bx_new[i2][j2][k2], sizeof(double), 1, gadget_out);
            fwrite(&grid[n].Bx_new[i2][j2][k2], sizeof(double), 1, gadget_out_bx);
            }

      printf("Bx_new[10][10][10] = %lg\n", grid[n].Bx_new[10][10][10]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            {
            fwrite(&grid[n].By_new[i2][j2][k2], sizeof(double), 1, gadget_out);
            fwrite(&grid[n].By_new[i2][j2][k2], sizeof(double), 1, gadget_out_by);
            }

      printf("By_new[10][10][10] = %lg\n", grid[n].By_new[10][10][10]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            {
            fwrite(&grid[n].Bz_new[i2][j2][k2], sizeof(double), 1, gadget_out);
            fwrite(&grid[n].Bz_new[i2][j2][k2], sizeof(double), 1, gadget_out_bz);
            }

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].Bx_new_face[i2][j2][k2], sizeof(double), 1, gadget_out);

      printf("Bx_new_face[10][10][10] = %lg\n", grid[n].Bx_new_face[10][10][10]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].By_new_face[i2][j2][k2], sizeof(double), 1, gadget_out);

      printf("By_new_face[10][10][10] = %lg\n", grid[n].By_new_face[10][10][10]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].Bz_new_face[i2][j2][k2], sizeof(double), 1, gadget_out);

      printf("Bz_new_face[10][10][10] = %lg\n", grid[n].Bz_new_face[10][10][10]);


}

if(print_vec == 1)
{
      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].Ax[i2][j2][k2], sizeof(double), 1, gadget_out);

     printf("Ax[10][10][10] = %lg\n", grid[n].Ax[10][10][10]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].Ay[i2][j2][k2], sizeof(double), 1, gadget_out);

     printf("Ay[10][10][10] = %lg\n", grid[n].Ay[10][10][10]);

      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].Az[i2][j2][k2], sizeof(double), 1, gadget_out);

    printf("Az[10][10][10] = %lg\n", grid[n].Az[10][10][10]);
}



      for(i2=0;i2<gnum;i2++)
       for(j2=0;j2<gnum;j2++)
         for(k2=0;k2<gnum;k2++)
            fwrite(&grid[n].DivB_new[i2][j2][k2], sizeof(double), 1, gadget_out);

fclose(gadget_out);
fclose(gadget_out_bx);
fclose(gadget_out_by);
fclose(gadget_out_bz);

kgrid=gnum/2-1;
//int shift=5;
int shift=0;
int pvalue = gnum/2-1, gogo;

     for(igrid=0;igrid<gnum;igrid++)
       for(jgrid=0;jgrid<gnum;jgrid++)
         for(kgrid=0;kgrid<gnum;kgrid++)
           {
           gogo = get_hilo();

           get_gphi_alt(n);

           gphi = grid[n].grad_phi[igrid][jgrid][kgrid]; 
           gphi_check = (gphi_xplus + gphi_xmin)/(1.0) + (gphi_yplus + gphi_ymin)/(1.0) + (gphi_zplus + gphi_zmin)/(1.0) - grid[n].DivB[igrid][jgrid][kgrid] * DeltaX ; 
           gphi_check = gphi_check/6.0;

           //if(((igrid-shift)%20 == 0 && (jgrid-shift)%20 == 0) || (igrid == pvalue && jgrid == pvalue))
           if(igrid%20==0 && jgrid%20==0 && kgrid%20==0)
             {
             printf("x = %lg, y = %lg, z = %lg, igrid = %d, jgrid = %d, kgrid = %d, i_hi = %d, j_hi = %d, k_hi = %d\n", grid[n].PosX[igrid][jgrid][kgrid], grid[n].PosY[igrid][jgrid][kgrid], grid[n].PosZ[igrid][jgrid][kgrid], igrid, jgrid, kgrid, i_hi, j_hi, k_hi); 
             printf("i2 = %d, j2 = %d, k2 = %d, DivB = %lg, DivB_new = %lg\n", igrid, jgrid, kgrid, grid[n].DivB_alt[igrid][jgrid][kgrid], grid[n].DivB_new[igrid][kgrid][jgrid]); 
             printf("Bx = %lg By = %lg Bz = %lg, gphi = %lg, gphi_check = %lg\n",  grid[n].Bfieldx[igrid][jgrid][kgrid], grid[n].Bfieldy[igrid][jgrid][kgrid], grid[n].Bfieldz[igrid][jgrid][kgrid], gphi, gphi_check );
             printf("face value: Bfieldx_face = %lg, Bfieldy_face = %lg, Bfieldz_face = %lg\n", grid[n].Bx_new_face[igrid][jgrid][kgrid], grid[n].By_new_face[igrid][jgrid][kgrid],  grid[n].Bz_new_face[igrid][jgrid][kgrid] );
             }
           }
}

int get_hilo(void)
{

     int go_on=0;


     i_lo  = igrid-1; i_hi = igrid+1;
     j_lo  = jgrid-1; j_hi = jgrid+1;
     k_lo  = kgrid-1; k_hi = kgrid+1;

     if(i_lo < 0) i_lo = gnum-1;  if(i_hi >= gnum) i_hi = 0;
     if(j_lo < 0) j_lo = gnum-1;  if(j_hi >= gnum) j_hi = 0;
     if(k_lo < 0) k_lo = gnum-1;  if(k_hi >= gnum) k_hi = 0;

/*
     if(i_lo < 0) i_lo = 0;  if(i_hi >= gnum) i_hi = gnum-1;
     if(j_lo < 0) j_lo = 0;  if(j_hi >= gnum) j_hi = gnum-1;
     if(k_lo < 0) k_lo = 0;  if(k_hi >= gnum) k_hi = gnum-1;
*/

     if(igrid == gnum-1 && i_hi == gnum-1) go_on = 1;
     if(jgrid == gnum-1 && j_hi == gnum-1) go_on = 1;
     if(kgrid == gnum-1 && k_hi == gnum-1) go_on = 1;
     if(igrid == 0 && i_lo == 0) go_on = 1;
     if(jgrid == 0 && j_lo == 0) go_on = 1;
     if(kgrid == 0 && k_lo == 0) go_on = 1;

     go_on = 0;

return(go_on);
}
                                                           

void get_gphi(int n)
{
     //if(igrid == 0 || jgrid == 0 || kgrid == 0) grid[n].grad_phi[igrid][jgrid][kgrid] = 0.;
     //if(igrid == gnum-1 || jgrid == gnum-1 || kgrid == gnum-1) grid[n].grad_phi[igrid][jgrid][kgrid] = 0.;

     gphi_xmin  = (grid[n].grad_phi[i_lo][jgrid][kgrid] + grid[n].grad_phi[igrid][jgrid][kgrid])/2.;
     gphi_xplus = (grid[n].grad_phi[i_hi][jgrid][kgrid] + grid[n].grad_phi[igrid][jgrid][kgrid])/2.;
     gphi_ymin  = (grid[n].grad_phi[igrid][j_lo][kgrid] + grid[n].grad_phi[igrid][jgrid][kgrid])/2.;
     gphi_yplus = (grid[n].grad_phi[igrid][j_hi][kgrid] + grid[n].grad_phi[igrid][jgrid][kgrid])/2.;
     gphi_zmin  = (grid[n].grad_phi[igrid][jgrid][k_lo] + grid[n].grad_phi[igrid][jgrid][kgrid])/2.;
     gphi_zplus = (grid[n].grad_phi[igrid][jgrid][k_hi] + grid[n].grad_phi[igrid][jgrid][kgrid])/2.;

     double gphi_here = grid[n].grad_phi[igrid][jgrid][kgrid];

if(interp_edge == 1)
{
    if(igrid == 0 /*&& i_lo == 0*/)
      gphi_xmin = gphi_here - (grid[n].grad_phi[0][jgrid][kgrid] + grid[n].grad_phi[1][jgrid][kgrid])/2;
    if(igrid == gnum-1 /*&& i_hi == gnum-1*/)
      gphi_xplus = gphi + (grid[n].grad_phi[gnum-1][jgrid][kgrid] + grid[n].grad_phi[gnum-2][jgrid][kgrid])/2;

    if(jgrid == 0 /*&& j_lo == 0*/)
      gphi_ymin  = gphi_here - (grid[n].grad_phi[igrid][0][kgrid] + grid[n].grad_phi[igrid][1][kgrid])/2;
    if(jgrid == gnum-1 /*&& j_hi == gnum-1*/)
      gphi_yplus = gphi + (grid[n].grad_phi[igrid][gnum-1][kgrid] + grid[n].grad_phi[igrid][gnum-2][kgrid])/2;

    if(kgrid == 0 /*&& k_lo == 0*/)
      gphi_zmin  = gphi_here - (grid[n].grad_phi[igrid][jgrid][0] + grid[n].grad_phi[igrid][jgrid][1])/2;
    if(kgrid == gnum-1 /*&& k_hi == gnum-1*/)
      gphi_zplus = gphi_here + (grid[n].grad_phi[igrid][jgrid][gnum-1] + grid[n].grad_phi[igrid][jgrid][gnum-2])/2;
}


//return(0);
}

void get_gphi_alt(int n)
{

     gphi_xmin  = (grid[n].grad_phi[i_lo][jgrid][kgrid]);
     gphi_xplus = (grid[n].grad_phi[i_hi][jgrid][kgrid]);
     gphi_ymin  = (grid[n].grad_phi[igrid][j_lo][kgrid]);
     gphi_yplus = (grid[n].grad_phi[igrid][j_hi][kgrid]);
     gphi_zmin  = (grid[n].grad_phi[igrid][jgrid][k_lo]);
     gphi_zplus = (grid[n].grad_phi[igrid][jgrid][k_hi]);

     double gphi_here = grid[n].grad_phi[igrid][jgrid][kgrid];


if(interp_edge == 1)
{
     if(igrid == 0 /*&& i_lo == 0*/)
       gphi_xmin = gphi_here - (grid[n].grad_phi[0][jgrid][kgrid] + grid[n].grad_phi[1][jgrid][kgrid])/2;
     if(igrid == gnum-1 /*&& i_hi == gnum-1*/)
       gphi_xplus = gphi + (grid[n].grad_phi[gnum-1][jgrid][kgrid] + grid[n].grad_phi[gnum-2][jgrid][kgrid])/2;

     if(jgrid == 0 /*&& j_lo == 0*/)
       gphi_ymin  = gphi_here - (grid[n].grad_phi[igrid][0][kgrid] + grid[n].grad_phi[igrid][1][kgrid])/2;
     if(jgrid == gnum-1 /*&& j_hi == gnum-1*/)
       gphi_yplus = gphi + (grid[n].grad_phi[igrid][gnum-1][kgrid] + grid[n].grad_phi[igrid][gnum-2][kgrid])/2;

     if(kgrid == 0 /*&& k_lo == 0*/)
       gphi_zmin  = gphi_here - (grid[n].grad_phi[igrid][jgrid][0] + grid[n].grad_phi[igrid][jgrid][1])/2;
     if(kgrid == gnum-1 /*&& k_hi == gnum-1*/)
       gphi_zplus = gphi_here + (grid[n].grad_phi[igrid][jgrid][gnum-1] + grid[n].grad_phi[igrid][jgrid][gnum-2])/2;
}

}

void get_Bnew(int n)
{

     Bx_new = grid[n].Bx_new[igrid][jgrid][kgrid];
     By_new = grid[n].By_new[igrid][jgrid][kgrid];
     Bz_new = grid[n].Bz_new[igrid][jgrid][kgrid];

     Bx_left = (grid[n].Bx_new[i_lo][jgrid][kgrid] + grid[n].Bx_new[igrid][jgrid][kgrid])/2.;
     Bx_right = (grid[n].Bx_new[i_hi][jgrid][kgrid] + grid[n].Bx_new[igrid][jgrid][kgrid])/2.;

     By_left = (grid[n].By_new[igrid][j_lo][kgrid] + grid[n].By_new[igrid][jgrid][kgrid])/2.;
     By_right = (grid[n].By_new[igrid][j_hi][kgrid] + grid[n].By_new[igrid][jgrid][kgrid])/2.;

     Bz_left = (grid[n].Bz_new[igrid][jgrid][k_lo] + grid[n].Bz_new[igrid][jgrid][kgrid])/2.;
     Bz_right = (grid[n].Bz_new[igrid][jgrid][k_hi] + grid[n].Bz_new[igrid][jgrid][kgrid])/2.;


    if(interp_edge == 1)
     {
     if(igrid == 0 /*&& i_lo == 0*/)
       Bx_left = Bx_new - (grid[n].Bx_new[0][jgrid][kgrid] + grid[n].Bx_new[1][jgrid][kgrid])/2;
     if(igrid == gnum-1 /*&& i_hi == gnum-1*/)
       Bx_right = Bx_new + (grid[n].Bx_new[gnum-1][jgrid][kgrid] + grid[n].Bx_new[gnum-2][jgrid][kgrid])/2;

     if(jgrid == 0 /*&& j_lo == 0*/)
       By_left  = By_new - (grid[n].By_new[igrid][0][kgrid] + grid[n].By_new[igrid][1][kgrid])/2;
     if(jgrid == gnum-1 /*&& j_hi == gnum-1*/)
       By_right = By_new + (grid[n].By_new[igrid][gnum-1][kgrid] + grid[n].By_new[igrid][gnum-2][kgrid])/2;

     if(kgrid == 0 /*&& k_lo == 0*/)
       Bz_left  = Bz_new - (grid[n].Bz_new[igrid][jgrid][0] + grid[n].Bz_new[igrid][jgrid][1])/2;
     if(kgrid == gnum-1 /*&& k_hi == gnum-1*/)
       Bz_right = Bz_new + (grid[n].Bz_new[igrid][jgrid][gnum-1] + grid[n].Bz_new[igrid][jgrid][gnum-2])/2;
     }
    if(interp_edge == 2)
     {
     if(igrid == 0 /*&& i_lo == 0*/)
       Bx_left = 1.e-16;
     if(igrid == gnum-1 /*&& i_hi == gnum-1*/)
       Bx_right = 1.e-16;

     if(jgrid == 0 /*&& j_lo == 0*/)
       By_left  = 1.e-16;
     if(jgrid == gnum-1 /*&& j_hi == gnum-1*/)
       By_right = 1.e-16;

     if(kgrid == 0 /*&& k_lo == 0*/)
       Bz_left  = 1.e-16;
     if(kgrid == gnum-1 /*&& k_hi == gnum-1*/)
       Bz_right = 1.e-16;

/*
     if(igrid == 0 && i_lo == 0)
       Bx_left = 1.e-17;
     if(igrid == gnum-1 && i_hi == gnum-1)
       Bx_right = 1.e-17;

     if(jgrid == 0 && j_lo == 0)
       By_left  = 1.e-17;
     if(jgrid == gnum-1 && j_hi == gnum-1)
       By_right = 1.e-17;

     if(kgrid == 0 && k_lo == 0)
       Bz_left  = 1.e-17;
     if(kgrid == gnum-1 && k_hi == gnum-1)
       Bz_right = 1.e-17;
*/
     }
 

}

void get_Bnew_alt(int n)
{ 
     Bx_new = grid[n].Bx_new[igrid][jgrid][kgrid];
     By_new = grid[n].By_new[igrid][jgrid][kgrid];
     Bz_new = grid[n].Bz_new[igrid][jgrid][kgrid];

     Bx_left  = grid[n].Bx_new[i_lo][jgrid][kgrid];
     Bx_right = grid[n].Bx_new[i_hi][jgrid][kgrid];

     By_left  = grid[n].By_new[igrid][j_lo][kgrid];
     By_right = grid[n].By_new[igrid][j_hi][kgrid];

     Bz_left  = grid[n].Bz_new[igrid][jgrid][k_lo];
     Bz_right = grid[n].Bz_new[igrid][jgrid][k_hi];

    if(interp_edge == 1)
     {
     if(igrid == 0 /*&& i_lo == 0*/) 
       Bx_left = Bx_new - (grid[n].Bx_new[0][jgrid][kgrid] + grid[n].Bx_new[1][jgrid][kgrid])/2;    
     if(igrid == gnum-1 /*&& i_hi == gnum-1*/)
       Bx_right = Bx_new + (grid[n].Bx_new[gnum-1][jgrid][kgrid] + grid[n].Bx_new[gnum-2][jgrid][kgrid])/2;

     if(jgrid == 0 /*&& j_lo == 0*/)
       By_left  = By_new - (grid[n].By_new[igrid][0][kgrid] + grid[n].By_new[igrid][1][kgrid])/2;
     if(jgrid == gnum-1 /*&& j_hi == gnum-1*/)
       By_right = By_new + (grid[n].By_new[igrid][gnum-1][kgrid] + grid[n].By_new[igrid][gnum-2][kgrid])/2;

     if(kgrid == 0 /*&& k_lo == 0*/)
       Bz_left  = Bz_new - (grid[n].Bz_new[igrid][jgrid][0] + grid[n].Bz_new[igrid][jgrid][1])/2;
     if(kgrid == gnum-1 /*&& k_hi == gnum-1*/)
       Bz_right = Bz_new + (grid[n].Bz_new[igrid][jgrid][gnum-1] + grid[n].Bz_new[igrid][jgrid][gnum-2])/2;
     }

     if(interp_cell == 2) { Bx_right = Bx_new; By_right = By_new; Bz_right = Bz_new; }
}

void set_Bedge(int n)
{
     grid[n].Bx_new[0][jgrid][kgrid] = grid[n].Bfieldx[0][jgrid][kgrid];
     grid[n].Bx_new[gnum-1][jgrid][kgrid] = grid[n].Bfieldx[gnum-1][jgrid][kgrid];

     grid[n].Bx_new[igrid][0][kgrid] = grid[n].Bfieldx[igrid][0][kgrid];
     grid[n].Bx_new[igrid][gnum-1][kgrid] = grid[n].Bfieldx[igrid][gnum-1][kgrid];

     grid[n].Bx_new[igrid][jgrid][0] = grid[n].Bfieldx[igrid][jgrid][0];
     grid[n].Bx_new[igrid][jgrid][gnum-1] = grid[n].Bfieldx[igrid][jgrid][gnum-1];

     grid[n].By_new[0][jgrid][kgrid] = grid[n].Bfieldy[0][jgrid][kgrid];
     grid[n].By_new[gnum-1][jgrid][kgrid] = grid[n].Bfieldy[gnum-1][jgrid][kgrid];

     grid[n].By_new[igrid][0][kgrid] = grid[n].Bfieldy[igrid][0][kgrid];
     grid[n].By_new[igrid][gnum-1][kgrid] = grid[n].Bfieldy[igrid][gnum-1][kgrid];

     grid[n].By_new[igrid][jgrid][0] = grid[n].Bfieldy[igrid][jgrid][0];
     grid[n].By_new[igrid][jgrid][gnum-1] = grid[n].Bfieldy[igrid][jgrid][gnum-1];

     grid[n].Bz_new[0][jgrid][kgrid] = grid[n].Bfieldz[0][jgrid][kgrid];
     grid[n].Bz_new[gnum-1][jgrid][kgrid] = grid[n].Bfieldz[gnum-1][jgrid][kgrid];

     grid[n].Bz_new[igrid][0][kgrid] = grid[n].Bfieldz[igrid][0][kgrid];
     grid[n].Bz_new[igrid][gnum-1][kgrid] = grid[n].Bfieldz[igrid][gnum-1][kgrid];

     grid[n].Bz_new[igrid][jgrid][0] = grid[n].Bfieldz[igrid][jgrid][0];
     grid[n].Bz_new[igrid][jgrid][gnum-1] = grid[n].Bfieldz[igrid][jgrid][gnum-1];
}

void skip_check(void)
{
}





 
