#include "Save2File.h"
#include <iostream>
#include <fstream>
#include <stdlib.h>

typedef __declspec(align(16)) float ALN16;

using namespace std;

extern int Nx, Ny;
extern ALN16 *Ez, *Hy, *Hx;

void File_Init()
{
	fstream OUT;

	OUT.open("Ez.txt", ios::out);
	OUT.close();
	OUT.open("Hx.txt", ios::out);

	OUT.close();
}

void Save_Ez();
void Save_Hy();
void Save_Hx();

void Save2File()
{
	Save_Ez();
	//Save_Hy();
	Save_Hx();
}

void Save_Ez()
{
	fstream outEz;

	outEz.open("Ez.txt", ios::app);
	int i = 0, j = 0;
	for (i = 0; i < Ny + 1; i++) {
		for (j = 0; j < Nx + 1; j++) {
			outEz << Ez[i*(Nx + 1) + j] << "\t";
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
	for (i = 0; i < Nx; i++) {
		for (j = 0; j < Nx-1; j++) {
			outHy << Hy[i*(Nx-1) + j] << "\t";
		}
		outHy << endl;
	}
	outHy << "\n";
	outHy.close();
}

void Save_Hx()
{
	ofstream outHx;
	int i, j;
	outHx.open("Hx.txt", ios::app);
	for (i = 0; i < Ny-1; i++){
		for (j = 0; j < Nx; j++){
			outHx << Hx[i*Nx + j] << "\t";
		}
		outHx << endl;
	}
	outHx.close();
}