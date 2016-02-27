% clear;
clc;
close all;

dz=0.015;%设置空间步长
% dz=1;

Nx = 31;
Ny = 21;
ezx= (0:dz:(Nx-1)*dz);%设置空间采样点位置
ezy= (0:dz:(Ny-1)*dz);%设置空间采样点位置

for i=1:100%i的取值范围要求在Ex和Hy的矩阵行数范围之内，代表时间取样点
    Z=Ez(1+(i-1)*Ny:i*Ny,1:Nx);
%     Z=Hx(1+(i-1)*Ny:i*Ny,1:Nx);
    
    mesh(ezx,ezy,Z)
    xlabel('X');
    ylabel('Y');
    zlabel('Z');
    view(-64,28)%view E
%     view(0,0)
    %axis([0 (N-1)*dz 0 (N-1)*dz 0 0.1])%设置坐标轴范围
    axis equal;
    axis([0 (Nx-1)*dz 0 (Ny-1)*dz 0 1])%设置坐标轴范围
    
    pause(0.01)
end