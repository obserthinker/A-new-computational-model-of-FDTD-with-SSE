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
//#define NEW

using namespace std;
extern int COUNT;
extern float coe_Ez, dt, dz;
extern int Nx;
extern ALN16 *Ez, *Hy, *Hx;

const float epsilon = 8.85e-12;
int head_sign, tail_sign;

void Ez_cmp_head();
void Ez_cmp_tail();
void Ez_vector_body();
float Hy_dif(int i, int j);

void Ez_Init()
{
	int i, N;

	N = (Nx + 1)*(Nx + 1);
	Ez = (ALN16*)_aligned_malloc(N*sizeof(float), 16);
	cout << "分配的电场空间网格数 : " << N << endl;
	
	for (i = 0; i < N; i++){
			Ez[i] = 0.f;
	}

	coe_Ez = dt / (epsilon* dz);
}

void Ez_cmp()
{
#ifdef SSE
#ifndef NEW
	Ez_cmp_head();
	Ez_cmp_tail();
#endif
#endif
	Ez_vector_body();
}


void Ez_vector_body()
{
	int i, j, N;
	int m, n;
#ifndef SSE
	for (i = 1; i < Nx; i++) {
		for (j = 1; j < Nx; j++) {
			Ez[i*(Nx + 1) + j] += coe_Ez*Hy_dif(i, j);
		}
	}
#else 
	__m128 *v_Ez, *v_Hx, *v_Hy;
	__m128 v_temp, v_Hx_left, v_Hy_up, v_coe_E, v_temp_Hy;
	//v_Hy_down=Hy[i]
	v_Ez = (__m128 *)Ez;
	v_Hx = (__m128 *)Hx;
	v_Hy = (__m128 *)Hy;
	v_coe_E = _mm_load1_ps(&coe_Ez);
#ifndef NEW
	for (i = head_sign + 1, j = (head_sign + 1) / 4; i <= tail_sign; i += 4, j++) {
		COUNT++;
		v_Hy_up = _mm_loadu_ps(Hy + i - (Nx + 1));
		v_Hx_left = _mm_loadu_ps(Hx + i - 1);
		v_temp_Hy = _mm_sub_ps(v_Hy[j], v_Hy_up);
		v_temp = _mm_sub_ps(v_Hx_left, v_Hx[j]);
		v_temp = _mm_add_ps(v_temp, v_temp_Hy);
		v_temp = _mm_mul_ps(v_temp, v_coe_E);
		v_Ez[j] = _mm_add_ps(v_Ez[j], v_temp);
	}
#else
	int row, limit;
	__m128 v_Ez_temp, v_Hy_down, v_Hx_right;
	//共有Nx+1个未知量，内存中0-Nx，去掉头尾，计算1——Nx-1，共Nx-1个需要计算
	limit = (int)((Nx-1) / 4) * 4;//limit为内存地址，退出循环时i=limit+1
	for (row = 1; row < Nx; row++) {//共有Nx+1行，只计算1——Nx-1行
		for (i = 1; i < limit; i += 4) {
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
					*(Ez + row*(Nx + 1) + i) += coe_Ez * ((*(Hx + (row - 1)*Nx + (i - 1)) - *(Hx + (row - 1)*Nx + i)) + (*(Hy + row*(Nx - 1) + (i - 1)) - *(Hy + (row - 1)*(Nx - 1) + (i - 1))));
				}
			}
		}
#endif
#endif
}

void Ez_cmp_head()
{
	int j, trunc;
	for (head_sign = Nx; head_sign <= Nx + 5; head_sign++){
		if ((head_sign+1) % 4 == 0){
			trunc = head_sign - (Nx + 1);
			//cout << "find head sign: " << head_sign<<endl;
			break;
		}
	}

	for (j = 0; j <= trunc; j++){
		COUNT++;
		Ez[Nx + 1 + j] += coe_Ez*Hy_dif(1, j);
	}
}

void Ez_cmp_tail()
{
	int  j, trunc;
	for (tail_sign = (Nx - 1)*(Nx + 1) + Nx; tail_sign >= (Nx - 1)*(Nx + 1) + Nx - 5; tail_sign--){
		if ((tail_sign + 1) % 4 == 0){
			trunc = tail_sign - (Nx - 1)*(Nx + 1);
			//cout << "find tail sign: " <<tail_sign <<endl;
			break;
		}
	}

	for (j = trunc + 1; j < Nx; j++){
		COUNT++;
		Ez[(Nx-1)*(Nx + 1) + j] += coe_Ez*Hy_dif(Nx-1, j);
	}
}

float Hy_dif(int i, int j)
{
	float dif_Hy, dif_Hx;
	//Hy(i,j)	-	Hy(i-1,j)	
	dif_Hy = Hy[i*(Nx + 1) + j] - Hy[(i - 1)*(Nx + 1) + j];
	//Hx(i,j-1)	-	Hx(i,j)
	dif_Hx = Hx[i*(Nx + 1) + j - 1] - Hx[i*(Nx + 1) + j];

	return (dif_Hx + dif_Hy);
}


