
#define pi 3.14159265359

double curlA(int j, int n, double xgrid, double ygrid, double zgrid, double hsm, double grid_size_cm, double grid_size_half_cm, double Ax0, double Ay0, double Az0)
{


double Jx, Jy, Jz; 
double Bx_new, By_new, Bz_new;
double mu0 = 1.0;
double B_send;
double DeltaX; 
double Bx, By, Bz; 
double Ax[2][2][2], Ay[2][2][2], Az[2][2][2];


DeltaX = grid_size_cm;

// [line integral] B dot dl = mu0 [surface integral] J dot dS

//A(r) = (1 / 4pi) * [volume integral] (del X B(r) ) / r

Bx = ((Az[1][1][1] - Az[1][0][1]) - (Ay[1][1][1] - Ay[1][1][0])) / grid_size_cm ;
By = ((Ax[1][1][1] - Ax[1][1][0]) - (Az[1][1][1] - Az[0][1][1])) / grid_size_cm ;
Bx = ((Ay[1][1][1] - Ay[0][1][1]) - (Ax[1][1][1] - Ax[1][0][1])) / grid_size_cm ;

if(j == 0) B_send = Bx;
if(j == 1) B_send = By;
if(j == 2) B_send = Bz;

return(B_send);

}
