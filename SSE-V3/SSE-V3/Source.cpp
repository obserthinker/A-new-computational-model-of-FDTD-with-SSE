#include "Source.h"
#include <cmath>
#include <iostream>

using namespace std;

typedef __declspec(align(16)) float ALN16;

extern float dz, dt, omega;
extern const float PI;
extern ALN16 *Ez;
extern int Nx, Ny;
const float C = 3e8;

void Src_Init()
{
	/*float period, lambda;
	period = 2 * PI / omega;
	lambda = period * C;
	dz = lambda / 12;*/
	dz = 0.015;
	dt = dz / (2*C);
}

void Src_compute(int i)
{
	float time = i * dt;
	float T, T0, vt, val_src;
	int src_position;

	T = 5e-10;
	T0 = 3 * T;
	vt = (time - T0) / T;

	val_src = exp(-pow(vt, 2));

	src_position = (Ny / 2)*(Nx + 1) + (Nx / 2);
	Ez[src_position] = val_src;
}