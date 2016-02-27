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
	Ez_cmp_head();
	Ez_cmp_tail();
#endif
	Ez_vector_body();
}


void Ez_vector_body()
{
	int i, j, N;
#ifndef SSE
	for (i = 1; i < Nx; i++){
		for (j = 1; j < Nx; j++){
			Ez[i*(Nx+1) + j] += coe_Ez*Hy_dif(i, j);
		}
	}
#else 
	__m128 *v_Ez, *v_Hx, *v_Hy;
	__m128 v_temp, v_Hx_left, v_Hy_up, v_coe_E,v_temp_Hy;
	//v_Hy_down=Hy[i]
	v_Ez = (__m128 *)Ez;
	v_Hx = (__m128 *)Hx;
	v_Hy = (__m128 *)Hy;
	v_coe_E = _mm_load1_ps(&coe_Ez);
	
	for (i = head_sign+1,j=(head_sign+1)/4;i<=tail_sign;i+=4,j++){
		v_Hy_up = _mm_loadu_ps(Hy +i- (Nx + 1));
		v_Hx_left = _mm_loadu_ps(Hx + i - 1);
		v_temp_Hy = _mm_sub_ps(v_Hy[j], v_Hy_up);
		v_temp = _mm_sub_ps(v_Hx_left, v_Hx[j]);
		v_temp = _mm_add_ps(v_temp, v_temp_Hy);
		v_temp = _mm_mul_ps(v_temp, v_coe_E);
		v_Ez[j] = _mm_add_ps(v_Ez[j], v_temp);
	}
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


