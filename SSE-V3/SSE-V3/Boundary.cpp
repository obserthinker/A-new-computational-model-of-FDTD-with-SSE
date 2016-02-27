#include "Boundary.h"
#include <iostream>
#include <fstream>

using namespace std;

extern float* Ez;
extern float *E_nbd_up, *E_nbd_down, *E_nbd_right, *E_nbd_left;
extern float *E_bd_up, *E_bd_down, *E_bd_right, *E_bd_left;
extern int Nx, Ny;
extern float dt, dz, coe_MUR;

const float C = 3e8f;

void Bd_proc_up();
void Bd_proc_down();
void Bd_proc_left_right();
void Bd_proc_point();

void Boundary_Init()
{
	int i;
	E_bd_left = (float *)malloc((Ny+1)*sizeof(float));
	memset(E_bd_left, 0, (Ny + 1)*sizeof(float));

	E_bd_right = (float *)malloc((Ny+1)*sizeof(float));
	memset(E_bd_right, 0, (Ny + 1)*sizeof(float));

	E_bd_up = (float *)malloc((Nx+1)*sizeof(float));
	memset(E_bd_up, 0, (Nx + 1)*sizeof(float));

	E_bd_down = (float *)malloc((Nx+1)*sizeof(float));
	memset(E_bd_down, 0, (Nx + 1)*sizeof(float));

	E_nbd_left = (float *)malloc((Ny+1)*sizeof(float));
	memset(E_nbd_left, 0, (Ny + 1)*sizeof(float));

	E_nbd_right = (float *)malloc((Ny+1)*sizeof(float));
	memset(E_nbd_right, 0, (Ny + 1)*sizeof(float));

	E_nbd_up = (float *)malloc((Nx+1)*sizeof(float));
	memset(E_nbd_up, 0, (Nx+1)*sizeof(float));

	E_nbd_down = (float *)malloc((Nx+1)*sizeof(float));
	memset(E_nbd_down, 0, (Nx + 1)*sizeof(float));

	coe_MUR = (C*dt - dz) / (C*dt + dz);
}

//has not been modified yet
void Boundary_PEC()
{
	int i, j;
	
	for (i = 0; i < Nx + 1;i++){
		if (i == 0 || i == Nx){
			for (j = 0; j < Nx + 1; j++){
				Ez[i*(Nx+1) + j] = 0.f;
			}
		}
		else{
			Ez[i*(Nx+1) + 0] = 0.f;
			Ez[i*(Nx+1) + Nx] = 0.f;
		}
	}
}

void Boundary_MUR1()
{
	Bd_proc_up();
	Bd_proc_down();
	Bd_proc_left_right();
	Bd_proc_point();
}

void Bd_proc_down()
{
	int j;
	for (j = 1; j < Nx; j++){
		Ez[j] = E_nbd_down[j] + coe_MUR*(Ez[(Nx + 1) + j] - E_bd_down[j]);
		E_nbd_down[j] = Ez[Nx + 1 + j];
		E_bd_down[j] = Ez[j];
	}
}

void Bd_proc_up()
{
	int j;
	for (j = 1; j < Nx; j++){
		Ez[Ny*(Nx + 1) + j] = E_nbd_up[j] + coe_MUR*(Ez[(Ny - 1)*(Nx + 1) + j] - E_bd_up[j]);
		E_nbd_up[j] = Ez[(Ny - 1)*(Nx + 1) + j];
		E_bd_up[j] = Ez[Ny*(Nx + 1) + j];
	}
}

void Bd_proc_left_right()
{
	int i;
	for (i = 1; i < Ny; i++){
		//left
		Ez[i*(Nx + 1)] = E_nbd_left[i] + coe_MUR*(Ez[i*(Nx + 1) + 1] - E_bd_left[i]);
		E_nbd_left[i] = Ez[i*(Nx + 1) + 1];
		E_bd_left[i] = Ez[i*(Nx + 1)];
		//right
		Ez[i*(Nx + 1) + Nx] = E_nbd_right[i] + coe_MUR*(Ez[i*(Nx + 1) + Nx - 1] - E_bd_right[i]);
		E_nbd_right[i] = Ez[i*(Nx + 1) + Nx - 1];
		E_bd_right[i] = Ez[i*(Nx + 1) + Nx];
	}
}

void Bd_proc_point()
{
	int i, j;
	Ez[0] = (Ez[1] + Ez[Nx + 1]) / 2;
	Ez[Nx] = (Ez[Nx - 1] + Ez[Nx + 1 + Nx]) / 2;
	Ez[Ny*(Nx + 1)] = (Ez[(Ny - 1)*(Nx + 1)] + Ez[Ny*(Nx + 1) + 1]) / 2;
	Ez[Ny*(Nx + 1) + Nx] = (Ez[(Ny - 1)*(Nx + 1) + Nx] + Ez[Ny*(Nx + 1) + Nx - 1]) / 2;
}