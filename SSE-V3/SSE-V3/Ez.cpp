#include "Ez.h"
#include "Boundary.h"
#include <cmath>
#include <xmmintrin.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>

typedef __declspec(align(16)) float ALN16;

#define SSE

using namespace std;

extern int COUNT;
extern float coe_Ez, dt, dz;
extern int Nx, Ny;
extern ALN16 *Ez, *Hy, *Hx;

const float epsilon = 8.85e-12;
int head_sign, tail_sign;

void Ez_vector_body();
float H_dif(int i, int j);

void Ez_Init()
{
	int i, N;

	N = (Nx + 1)*(Ny + 1);
	Ez = (ALN16*)_aligned_malloc(N*sizeof(float), 16);
	cout << "分配的电场空间网格数 : " << N << endl;
	
	for (i = 0; i < N; i++){
			Ez[i] = 0.f;
	}

	coe_Ez = dt / (epsilon* dz);
}

void Ez_cmp()
{
	int i, j, N;
	int m, n;
#ifndef SSE
	for (i = 1; i < Ny; i++) {
		for (j = 1; j < Nx; j++) {
			Ez[i*(Nx + 1) + j] += coe_Ez*H_dif(i, j);
		}
	}

#else 
	__m128 *v_Ez, *v_Hx, *v_Hy;
	__m128 v_temp, v_Hx_left, v_Hy_up, v_coe_E, v_temp_Hy;
	__m128 v_Ez_temp, v_Hy_down, v_Hx_right;
	__m128 tail;
	int row, limit;

	v_Ez = (__m128 *)Ez;
	v_Hx = (__m128 *)Hx;
	v_Hy = (__m128 *)Hy;
	v_coe_E = _mm_load1_ps(&coe_Ez);

	//共有Nx+1个未知量，内存中0-Nx，去掉头尾，计算1——Nx-1，共Nx-1个需要计算
	limit = (int)((Nx-1) / 4) * 4;//limit为内存地址，退出循环时i=limit+1
	for (row = 1; row < Ny; row++) {//共有Ny+1行，只计算1——Ny-1行
		tail = _mm_loadu_ps(Ez + row*(Nx + 1) + (Nx - 3));
		for (i = 1; i < limit; i += 4) {
			COUNT++;
			v_Ez_temp = _mm_loadu_ps(Ez + row*(Nx + 1) + i);

			v_Hy_up = _mm_loadu_ps(Hy + row*(Nx - 1) + (i - 1));
			v_Hy_down = _mm_loadu_ps(Hy + (row - 1)*(Nx - 1) + (i - 1));
			v_Hx_left = _mm_loadu_ps(Hx + (row - 1)*Nx + (i - 1));
			v_Hx_right = _mm_loadu_ps(Hx + (row - 1)*Nx + i);

			v_temp_Hy = _mm_sub_ps(v_Hy_up, v_Hy_down);
			v_temp = _mm_sub_ps(v_Hx_left, v_Hx_right);
			v_temp = _mm_add_ps(v_temp, v_temp_Hy);
			v_temp = _mm_mul_ps(v_temp, v_coe_E);
			v_Ez_temp = _mm_add_ps(v_Ez_temp, v_temp);
			_mm_storeu_ps(Ez + row*(Nx + 1) + i, v_Ez_temp);
		}
			if (limit < Nx-1) {
				for (i = limit+1; i < Nx; i++) {
					COUNT++;
					*(Ez + row*(Nx + 1) + i) += coe_Ez * ((*(Hx + (row - 1)*Nx + (i - 1)) - *(Hx + (row - 1)*Nx + i)) + (*(Hy + row*(Nx - 1) + (i - 1)) - *(Hy + (row - 1)*(Nx - 1) + (i - 1))));
				}
			}
		/*_mm_storeu_ps(Ez + row*(Nx + 1) + (Nx - 3), tail);
		i = Nx - 3;
		v_Ez_temp = _mm_loadu_ps(Ez + row*(Nx + 1) + i);

		v_Hy_up = _mm_loadu_ps(Hy + row*(Nx - 1) + (i - 1));
		v_Hy_down = _mm_loadu_ps(Hy + (row - 1)*(Nx - 1) + (i - 1));
		v_Hx_left = _mm_loadu_ps(Hx + (row - 1)*Nx + (i - 1));
		v_Hx_right = _mm_loadu_ps(Hx + (row - 1)*Nx + i);

		v_temp_Hy = _mm_sub_ps(v_Hy_up, v_Hy_down);
		v_temp = _mm_sub_ps(v_Hx_left, v_Hx_right);
		v_temp = _mm_add_ps(v_temp, v_temp_Hy);
		v_temp = _mm_mul_ps(v_temp, v_coe_E);
		v_Ez_temp = _mm_add_ps(v_Ez_temp, v_temp);
		_mm_storeu_ps(Ez + row*(Nx + 1) + i, v_Ez_temp);*/
	}
#endif
}

float H_dif(int i, int j)
{
	float dif_Hy, dif_Hx;
	//Hy(i,j)	-	Hy(i-1,j)	
	dif_Hy = Hy[i*(Nx - 1) + j-1] - Hy[(i - 1)*(Nx - 1) + j-1];	
	//Hx(i,j-1)	-	Hx(i,j)
	dif_Hx = Hx[(i-1)*Nx + j - 1] - Hx[(i-1)*Nx + j];

	return (dif_Hx + dif_Hy);
}


