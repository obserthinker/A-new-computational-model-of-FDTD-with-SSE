#include "H.h"
#include <xmmintrin.h>
#include <cmath>
#include <iostream>
#include <fstream>

typedef __declspec(align(16)) float ALN16;

#define SSE
//#define NEW

using namespace std;
extern int COUNT;
extern ALN16 *Ez,*Hy,*Hx;
extern int Nx;
extern float dt, dz, coe_H;
//extern const float mu;
const float PI = 3.141592653589793;
const float mu = (4.0*PI)*1e-7;

void H_Init()
{
	int N,i;
#ifndef NEW
	N = Nx*(Nx + 1);
#else
	N = Nx*(Nx - 1);
#endif
	Hy = (ALN16*)_aligned_malloc(N*sizeof(float), 16);
	Hx = (ALN16*)_aligned_malloc(N*sizeof(float), 16);
	for (i = 0; i < N; i++){
		Hy[i] = 0.f;
		Hx[i] = 0.f;
	}
	coe_H = dt / (mu* dz);
}

void H_compute()
{
	int i, j, N;
	
#ifndef SSE
	//Hy
	for (i = 0; i < Nx; i++){
		for (j = 1; j < Nx; j++){
			//cout << "before= "<<"Hy("<<i<<" , "<<j<<")"<<Hy[i*(Nx + 1) + j] << "\t";
			Hy[i*(Nx + 1) + j] += coe_H*(Ez[(i + 1)*(Nx + 1) + j] - Ez[i*(Nx + 1) + j]);
			//cout << coe_Hy*(Ez[(i + 1)*(Nx + 1) + j + 1] - Ez[i*(Nx + 1) + j + 1]) << "\t" << Hy[i*(Nx + 1) + j] << endl;
		}
	}
	//Hx
	for (i = 1; i < Nx; i++){
		for (j = 0; j < Nx; j++){
			//cout << "before= " << "Hx(" << i << " , " << j << ")" << Hy[i*Nx + j] << "\t";
			Hx[i*(Nx + 1) + j] += coe_H*(Ez[i*(Nx + 1) + j] - Ez[i*(Nx + 1) + j + 1]);
			//cout << coe_Hy*(Ez[(i + 1)*(Nx + 1) + j] - Ez[(i + 1)*(Nx + 1) + j + 1]) << "\t" << Hy[i*Nx + j] << endl;
		}
	}
#else
	__m128 v_coe_H, v_Ez_down, v_temp, v_Ez_right;
	__m128 *v_Hy, *v_Hx, *v_Ez;

	v_Hy = (__m128 *)Hy;
	v_Hx = (__m128 *)Hx;
	v_Ez = (__m128 *)Ez;
	v_coe_H = _mm_load1_ps(&coe_H);
	N = Nx*(Nx + 1);

	//计算Hy
#ifndef NEW
	for (i = 0, j = 0; i < N; i += 4, j++){
		COUNT++;
		v_Ez_down = _mm_loadu_ps(Ez + i + Nx + 1);
		v_temp = _mm_sub_ps(v_Ez_down, v_Ez[j]);
		v_temp = _mm_mul_ps(v_temp, v_coe_H);
		v_Hy[j] = _mm_add_ps(v_Hy[j], v_temp);
	}
#else
	int row, limit;
	__m128 v_Hy_temp, v_Ez_up;
	limit = (int)((Nx - 1) / 4) * 4;
	for (row = 0; row < Nx; row++) {
		for (i = 0; i < limit; i += 4) {
			v_Hy_temp = _mm_loadu_ps(Hy + row*(Nx - 1) + i);
			v_Ez_down = _mm_loadu_ps(Ez + row*(Nx + 1) + i + 1);//此处的v_Ez_down与上面相反
			v_Ez_up = _mm_loadu_ps(Ez + (row + 1)*(Nx + 1) + i + 1);
			v_temp = _mm_sub_ps(v_Ez_up, v_Ez_down);
			v_temp = _mm_mul_ps(v_temp, v_coe_H);
			v_Hy_temp = _mm_add_ps(v_Hy_temp, v_temp);
			_mm_storeu_ps(Hy + row*(Nx - 1) + i, v_Hy_temp);
		}	
		if (limit < Nx-1) {
			for (i = limit; i < Nx-1; i++) {
				*(Hy + row*(Nx - 1) + i) += coe_H*(*(Ez + (row + 1)*(Nx + 1) + i + 1) - *(Ez + row*(Nx + 1) + i + 1));
			}
		}
	}
#endif

#ifndef NEW
	for (i = 0, j = 0; i < N; i += 4, j++){
		COUNT++;
		v_Ez_right = _mm_loadu_ps(Ez + i + 1);
		v_temp = _mm_sub_ps(v_Ez[j], v_Ez_right);
		v_temp = _mm_mul_ps(v_temp, v_coe_H);
		v_Hx[j] = _mm_add_ps(v_Hx[j], v_temp);
	}
#else
	__m128 v_Hx_temp, v_Ez_left;
	limit = (int)(Nx / 4) * 4;
	for (row = 0; row < Nx - 1; row++) {//row指Hx的行数
		for (i = 0; i < limit; i += 4) {
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
				*(Hx + row*Nx + i) += coe_H*(*(Ez + (row + 1)*(Nx + 1) + i) - *(Ez + (row + 1)*(Nx + 1) + i + 1));
			}
		}
	}
#endif
#endif
}