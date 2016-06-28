
#define pi 3.14159265359

double vec_pot(int j, int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size_cm, double grid_size_half_cm, double Bx0, double By0, double Bz0)
{


double Jx, Jy, Jz; 
double Bx_new, By_new, Bz_new;
double mu0 = 1.0;
double A_send;
double DeltaX; 
double Ax, Ay, Az; 
double Bx[2][2][2], By[2][2][2], Bz[2][2][2];


DeltaX = grid_size_cm;

// [line integral] B dot dl = mu0 [surface integral] J dot dS

//Jx = Bx / grid_size_cm / mu0 * 4. *pi;
//Jy = By / grid_size_cm / mu0 * 4. *pi;
//Jz = Bz / grid_size_cm / mu0 * 4. *pi;

//Ax = (Bz - By) * grid_size_cm ;
//Ay = (Bx - Bz) * grid_size_cm ;
//Az = (By - Bx) * grid_size_cm ;

//A(r) = (1 / 4pi) * [volume integral] (del X B(r) ) / r

Ax = ((Bz[1][1][1] - Bz[1][0][1]) - (By[1][1][1] - By[1][1][0])) * grid_size_cm ;
Ay = ((Bx[1][1][1] - Bx[1][1][0]) - (Bz[1][1][1] - Bz[0][1][1])) * grid_size_cm ;
Ax = ((By[1][1][1] - By[0][1][1]) - (Bx[1][1][1] - Bx[1][0][1])) * grid_size_cm ;

if(j == 0) A_send = Ax;
if(j == 1) A_send = Ay;
if(j == 2) A_send = Az;

return(A_send);

}
