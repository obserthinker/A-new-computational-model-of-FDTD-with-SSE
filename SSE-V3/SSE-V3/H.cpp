#include "H.h"
#include <xmmintrin.h>
#include <cmath>
#include <iostream>
#include <fstream>

typedef __declspec(align(16)) float ALN16;

#define SSE

using namespace std;

extern ALN16 *Ez,*Hy,*Hx;
extern int Nx, Ny;
extern int COUNT;
extern float dt, dz, coe_H;
//extern const float mu;
const float PI = 3.141592653589793;
const float mu = (4.0*PI)*1e-7;

void H_Init()
{
	int N_Hx, N_Hy, i;

	N_Hx = Nx*(Ny - 1);
	N_Hy = (Nx - 1)*Ny;
	Hx = (ALN16*)_aligned_malloc(N_Hx * sizeof(float), 16);
	for (i = 0; i < N_Hx; i++) {
		Hx[i] = 0;
	}
	Hy = (ALN16*)_aligned_malloc(N_Hy * sizeof(float), 16);
	for (i = 0; i < N_Hy; i++) {
		Hy[i] = 0.f;
	}
	
	coe_H = dt / (mu* dz);
}

void H_cmp()
{
	int i, j;
	
#ifndef SSE
	//Hy
	for (i = 0; i < Ny; i++){
		for (j = 0; j < Nx-1; j++){
			//cout << i << "\t" << j << endl;
			Hy[i*(Nx - 1) + j] += coe_H*(Ez[(i + 1)*(Nx + 1) + j + 1] - Ez[i*(Nx + 1) + j + 1]);
			//cout << Hy[i*(Nx - 1) + j] << endl;
		}
	}
	//getchar();
	//Hx
	for (i = 0; i < Ny-1; i++){
		for (j = 0; j < Nx; j++){
			//cout << i << "\t" << j << endl;
			Hx[i*Nx + j] += coe_H*(Ez[(i+1)*(Nx + 1) + j] - Ez[(i+1)*(Nx + 1) + j + 1]);
			//cout << Hx[i*Nx + j] << endl;
		}
	}
	//getchar();
#else
	__m128 v_coe_H, v_Ez_down, v_temp, v_Ez_right;
	__m128 *v_Hy, *v_Hx, *v_Ez;
	__m128 tail;
	int row, limit;

	v_Hy = (__m128 *)Hy;
	v_Hx = (__m128 *)Hx;
	v_Ez = (__m128 *)Ez;
	v_coe_H = _mm_load1_ps(&coe_H);

	//Hx
	__m128 v_Hx_temp, v_Ez_left;
	limit = (int)(Nx / 4) * 4;
	for (row = 0; row < Ny - 1; row++) {
		tail = _mm_loadu_ps(Hx + row*Nx + Nx - 4);
		for (i = 0; i < limit; i += 4) {
			COUNT++;
			v_Hx_temp = _mm_loadu_ps(Hx + row*Nx + i);
			v_Ez_right = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i + 1);
			v_Ez_left = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i);
			v_temp = _mm_sub_ps(v_Ez_left, v_Ez_right);
			v_temp = _mm_mul_ps(v_temp, v_coe_H);
			v_Hx_temp = _mm_add_ps(v_Hx_temp, v_temp);
			_mm_storeu_ps(Hx + row*Nx + i, v_Hx_temp);
		}
		if (limit < Nx) {
			for (i = limit; i < Nx; i++) {
				COUNT++;
				*(Hx + row*Nx + i) += coe_H*(*(Ez + (row + 1)*(Nx + 1) + i) - *(Ez + (row + 1)*(Nx + 1) + i + 1));
			}
		}
		/*_mm_storeu_ps(Hx + row*Nx + Nx - 4, tail);
		i = Nx - 4;
		v_Hx_temp = _mm_loadu_ps(Hx + row*Nx + i);
		v_Ez_right = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i + 1);
		v_Ez_left = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i);
		v_temp = _mm_sub_ps(v_Ez_left, v_Ez_right);
		v_temp = _mm_mul_ps(v_temp, v_coe_H);
		v_Hx_temp = _mm_add_ps(v_Hx_temp, v_temp);
		_mm_storeu_ps(Hx + row*Nx + i, v_Hx_temp);*/
	}

	//¼ÆËãHy
	__m128 v_Hy_temp, v_Ez_up;
	limit = (int)((Nx - 1) / 4) * 4;
	for (row = 0; row < Ny; row++) {
		tail = _mm_loadu_ps(Hy + row*(Nx - 1) + Nx - 5);
		for (i = 0; i < limit; i += 4) {
			COUNT++;
			v_Hy_temp = _mm_loadu_ps(Hy + row*(Nx - 1) + i);
			v_Ez_down = _mm_loadu_ps(Ez + row*(Nx + 1) + i + 1);
			v_Ez_up = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i + 1);
			v_temp = _mm_sub_ps(v_Ez_up, v_Ez_down);
			v_temp = _mm_mul_ps(v_temp, v_coe_H);
			v_Hy_temp = _mm_add_ps(v_Hy_temp, v_temp);
			_mm_storeu_ps(Hy + row*(Nx - 1) + i, v_Hy_temp);
		}	
		if (limit < Nx-1) {
			for (i = limit; i < Nx-1; i++) {
				COUNT++;
				*(Hy + row*(Nx - 1) + i) += coe_H*(*(Ez + (row + 1)*(Nx + 1) + i + 1) - *(Ez + row*(Nx + 1) + i + 1));
			}
		}
		/*_mm_storeu_ps(Hy + row*(Nx - 1) + Nx - 5, tail);
		i = Nx - 5;
		v_Hy_temp = _mm_loadu_ps(Hy + row*(Nx - 1) + i);
		v_Ez_down = _mm_loadu_ps(Ez + row*(Nx + 1) + i + 1);
		v_Ez_up = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i + 1);
		v_temp = _mm_sub_ps(v_Ez_up, v_Ez_down);
		v_temp = _mm_mul_ps(v_temp, v_coe_H);
		v_Hy_temp = _mm_add_ps(v_Hy_temp, v_temp);
		_mm_storeu_ps(Hy + row*(Nx - 1) + i, v_Hy_temp);*/
	}
#endif
}