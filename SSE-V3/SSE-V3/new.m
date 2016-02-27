function v=new(Nx,Ny)

%Ez
q=floor((Nx-1)/4);
r=(Nx-1)-q*4;
E = (Ny-1)*(q+r);

%Hx
q=floor(Nx/4);
r=Nx-4*q;
Hx=(Ny-1)*(q+r);

%Hy
q=floor((Nx-1)/4);
r=(Nx-1)-4*q;
Hy=Ny*(q+r);

v=E+Hx+Hy;
