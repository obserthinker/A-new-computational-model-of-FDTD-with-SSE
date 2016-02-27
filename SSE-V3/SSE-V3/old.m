function v=old(Nx,Ny)

%v = (Ny-1)*(Nx-1)+(Ny-1)*Nx+Ny*(Nx-1);

E=ceil((Ny-1)*(Nx+1)/4);

Hx=ceil(Ny*(Nx+1)/4);
Hy=Hx;

v=E+Hx+Hy;


