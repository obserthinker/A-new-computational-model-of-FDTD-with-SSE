#include "Compute.h"
#include "Ez.h"
#include "H.h"
#include "Boundary.h"
#include "Source.h"
#include "Save2File.h"
#include <iostream>

using namespace std;

extern int Nt,Nx;
extern float* Ez;

void compute()
{
	int i;

	for (i = 0; i < Nt; i++){
		H_compute();
		Ez_cmp();
		Boundary_MUR1();
		//Boundary_cmp_PEC();
		Src_compute(i);
		
		//Save2File(i);
	}
}
