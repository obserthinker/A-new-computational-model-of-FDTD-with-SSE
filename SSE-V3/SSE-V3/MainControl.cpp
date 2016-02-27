#include "Input.h"
#include "Compute.h"
#include <iostream>
#include <fstream>

using namespace std;

typedef __declspec(align(16)) float ALN16;

int COUNT;
int Nt, Nx, Ny;
float coe_H, coe_Ez, coe_MUR, dt, dz;
float *E_nbd_up, *E_nbd_down, *E_nbd_right, *E_nbd_left;
float *E_bd_up, *E_bd_down, *E_bd_right, *E_bd_left;
ALN16 *Ez, *Hy, *Hx;
const float PI = 3.141592653589793;
const float mu = (4.0*PI)*1e-7;
const float epsilon = 8.85e-12;

void main()
{
	Input();
	compute();
}
