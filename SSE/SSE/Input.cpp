#include "Input.h"
#include "Boundary.h"
#include "Ez.h"
#include "H.h"
#include "Source.h"
#include "Save2File.h"
#include <iostream>

#define test

using namespace std;

float space_length, time_length, omega;

extern float coe_MUR, dt, dz, coe_Ez, coe_H;
const float PI = 3.141592653589793;
extern const float mu;
extern const float epsilon;

extern int Nt, Nx;

void Get_Coe()
{
#ifndef test
#else
	space_length = 0.3;
	time_length = 3;
#endif
}

void Init_check()
{
	cout << "dz = " << dz << endl;
	cout << "dt = " << dt << endl;
	cout << "space grid number : " << Nx << endl;
	cout << "time grid number : " << Nt << endl;
	//cout << "space length: " << space_length << endl;
	//cout << "time length: " << time_length << endl;
	//cout << "omega: " << omega << endl;
	cout << "Coe_Ez: " << coe_Ez << endl;
	cout << "Coe_Hy: " << coe_H << endl;
	cout << "Coe_MUR: " << coe_MUR << endl;
}

void Init()
{
	Src_Init();
#ifndef test
	Nx = (int)ceilf(space_length / dz);
	Nt = (int)ceilf(time_length / dt);
#else
	Nx = 1000;
	Nt = 3000;

	//1000 1000
	//1000 5000
	//1001 1000
	//1002 1000
	//1003 1000
	//1004 1000
	//1005 1000
	//1006 1000
	//1007 1000
#endif
	Ez_Init();
	H_Init();
	Boundary_Init();
	File_Init();
	Init_check();
}


void Input()
{
	//Get_Coe();
	Init();
	cout << "press enter to continue" << endl;
	//getchar();
	cout << "computing..." << endl;
}

