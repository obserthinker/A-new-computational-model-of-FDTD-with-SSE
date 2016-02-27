#include "Compute.h"
#include "Ez.h"
#include "H.h"
#include "Boundary.h"
#include "Source.h"
#include "Save2File.h"
#include <iostream>

using namespace std;

extern int Nt;
extern int COUNT;

void compute()
{
	int i;
	COUNT = 0;
	for (i = 0; i < Nt; i++){
		H_cmp();
		Ez_cmp();
		Boundary_MUR1();
		//Boundary_PEC();
		Src_compute(i);
		cout << "count: "<<COUNT << endl;
		//Save2File();
	}
}
