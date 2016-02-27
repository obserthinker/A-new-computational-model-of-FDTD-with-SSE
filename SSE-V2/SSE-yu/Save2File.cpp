#include "Save2File.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

#define NEW

typedef __declspec(align(16)) float ALN16;

using namespace std;

extern int Nx;
extern ALN16 *Ez, *Hy, *Hx;

void File_Init()
{
	fstream outEz, outHy;

	outEz.open("Ez_SSE.txt", ios::out);
	outEz.close();
	outHy.open("Hy.txt", ios::out);
	outHy.close();
	outHy.open("Hx.txt", ios::out);
	outHy.close();
}

void Save_Ez(int i);
void Save_Hy();
void Save_Hx();

void Save2File(int i)
{
	Save_Ez(i);
	Save_Hy();
	//Save_Hx();
}

void Save_Ez(int c)
{
	fstream outEz;

	outEz.open("Ez_SSE.txt", ios::app);
	int i=0, j=0;
	for (i = 0; i < Nx + 1; i++){
		for (j = 0; j < Nx + 1; j++){
			outEz << Ez[i*(Nx+1)+j] << "\t";
		}
		outEz << endl;
	}
	//outEz << "\n";
	outEz.close();
}
void Save_Hy()
{
	fstream outHy;
	int i,j;
	outHy.open("Hy.txt", ios::app);
#ifdef NEW
	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Nx-1; j++) {
			outHy << Hy[i*(Nx-1) + j] << "\t";
		}
#else
	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Nx - 1; j++) {
			outHy << Hy[i*(Nx + 1) + j] << "\t";
		}
#endif
		outHy << endl;
	}
	outHy << "\n";
	outHy.close();
}

void Save_Hx()
{
	fstream outHx;
	int i, j;
	outHx.open("Hx.txt", ios::app);
	for (i = 0; i < Nx-1; i++){
		for (j = 0; j < Nx; j++){
			outHx << Hx[i*Nx + j] << "\t";
		}
		outHx << endl;
	}
	outHx.close();
}